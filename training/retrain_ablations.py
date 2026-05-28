"""Retrain-based ablation study for the classifier and regressor.

Trains fresh models from scratch under each critical ablation condition,
evaluates on the held-out test set, and prints a summary table.

Ablation conditions:
  baseline_retrain   - full model retrained from scratch (sanity check)
  no_primer_enc      - zero sequence channels 5:10 (primer one-hot) during training
  no_target_enc      - zero sequence channels 0:5  (target one-hot) during training
  no_fwd_seq         - zero sequence positions 0:28 (forward primer region)
  no_rev_seq         - zero sequence positions 28:56 (reverse primer region)
  no_pgc2            - Janus without the pgc2 layer (largest AUROC drop in probe study)
  no_s4d             - Janus without the s4d layer  (largest R² drop in probe study)
  mlp_only           - train without SSM branch (features only)
  ssm_only           - train without MLP branch (sequences only)

Usage (from repo root):
    python training/retrain_ablations.py
"""

import numpy as np
import pandas as pd
import torch
import torch.nn as nn
from torch.utils.data import DataLoader, TensorDataset
from sklearn.metrics import roc_auc_score, r2_score
from sklearn.preprocessing import StandardScaler
from tqdm import tqdm

from qprimer_designer.models import load_models
from qprimer_designer.models.architectures import (
    Janus, MLP, CombinedModel, CombinedModelClassifier,
)
from qprimer_designer.utils import encode_batch_parallel

# ── Column mapping ───────────────────────────────────────────────────────────
FEAT_RENAME = {
    "f_length": "len_f", "f_Tm": "Tm_f", "f_GC": "GC_f",
    "f_indel": "indel_f", "f_mm": "mm_f",
    "r_length": "len_r", "r_Tm": "Tm_r", "r_GC": "GC_r",
    "r_indel": "indel_r", "r_mm": "mm_r",
    "prod_length": "prod_len",
}
SEQ_RENAME = {
    "f_penc": "pseq_f", "f_tenc": "tseq_f",
    "r_penc": "pseq_r", "r_tenc": "tseq_r",
}
FEATURE_COLUMNS = [
    "len_f", "Tm_f", "GC_f", "indel_f", "mm_f",
    "len_r", "Tm_r", "GC_r", "indel_r", "mm_r",
    "prod_len", "prod_Tm",
]
SEQUENCE_COLUMNS = ["pseq_f", "tseq_f", "pseq_r", "tseq_r"]

NUM_EPOCHS = 50
BATCH_SIZE = 512
LR = 1e-3
WEIGHT_DECAY = 0.01
PATIENCE = 8


# ── Architectural variants ───────────────────────────────────────────────────

class JanusNoPGC2(Janus):
    """Janus with pgc2 skipped (pass-through)."""
    def forward(self, u):
        x = self.encoder(u)
        x = self.pgc1(x)
        z = self.norm(x)
        x = self.dropout(self.s4d(z)) + x
        x = x.mean(dim=1)
        return self.decoder(x)


class JanusNoS4D(Janus):
    """Janus with s4d skipped (no state-space convolution)."""
    def forward(self, u):
        x = self.encoder(u)
        x = self.pgc1(x)
        x = self.pgc2(x)
        x = x.mean(dim=1)
        return self.decoder(x)


# ── Branch-only model wrappers ───────────────────────────────────────────────

class MLPOnlyClassifier(nn.Module):
    def __init__(self, mlp_in, mlp_hid, mlp_out, com_hid):
        super().__init__()
        self.mlp = MLP(mlp_in, mlp_out, mlp_hid)
        self.head = nn.Sequential(
            nn.Linear(mlp_out, com_hid), nn.ReLU(),
            nn.Linear(com_hid, 1), nn.Sigmoid(),
        )

    def forward(self, feat, seq):
        return self.head(self.mlp(feat))


class SSMOnlyClassifier(nn.Module):
    def __init__(self, ssm_in, ssm_out, ssm_mod, com_hid):
        super().__init__()
        self.ssm = Janus(ssm_in, ssm_out, ssm_mod)
        self.head = nn.Sequential(
            nn.Linear(ssm_out, com_hid), nn.ReLU(),
            nn.Linear(com_hid, 1), nn.Sigmoid(),
        )

    def forward(self, feat, seq):
        return self.head(self.ssm(seq))


class MLPOnlyRegressor(nn.Module):
    def __init__(self, mlp_in, mlp_hid, mlp_out, com_hid):
        super().__init__()
        self.mlp = MLP(mlp_in, mlp_out, mlp_hid)
        self.head = nn.Sequential(
            nn.Linear(mlp_out, com_hid), nn.ReLU(),
            nn.Linear(com_hid, 1),
        )

    def forward(self, feat, seq):
        return self.head(self.mlp(feat))


class SSMOnlyRegressor(nn.Module):
    def __init__(self, ssm_in, ssm_out, ssm_mod, com_hid):
        super().__init__()
        self.ssm = Janus(ssm_in, ssm_out, ssm_mod)
        self.head = nn.Sequential(
            nn.Linear(ssm_out, com_hid), nn.ReLU(),
            nn.Linear(com_hid, 1),
        )

    def forward(self, feat, seq):
        return self.head(self.ssm(seq))


# ── Model construction from baseline dims ────────────────────────────────────

def _cls_dims(model):
    mlp = model.mlp
    ssm = model.ssm
    return {
        "mlp_in":  mlp.model[0].in_features,
        "mlp_hid": mlp.model[0].out_features,
        "mlp_out": mlp.model[2].out_features,
        "ssm_in":  ssm.encoder.in_features,
        "ssm_mod": ssm.encoder.out_features,
        "ssm_out": ssm.decoder.out_features,
        "com_hid": model.combiner[0].out_features,
    }


def _reg_dims(model):
    mlp = model.mlp
    dl = model.dl
    return {
        "mlp_in":  mlp.model[0].in_features,
        "mlp_hid": mlp.model[0].out_features,
        "mlp_out": mlp.model[2].out_features,
        "dl_in":   dl.encoder.in_features,
        "dl_mod":  dl.encoder.out_features,
        "dl_out":  dl.decoder.out_features,
        "com_hid": model.combiner[0].out_features,
    }


def build_classifier(cd, ssm_class=Janus):
    m = CombinedModelClassifier(
        mlp_dims=(cd["mlp_in"], cd["mlp_out"], cd["mlp_hid"]),
        ssm_dims=(cd["ssm_in"], cd["ssm_out"], cd["ssm_mod"]),
        combined_hidden=cd["com_hid"],
        final_output=1,
    )
    if ssm_class is not Janus:
        m.ssm = ssm_class(cd["ssm_in"], cd["ssm_out"], cd["ssm_mod"])
    return m


def build_regressor(rd, ssm_class=Janus):
    m = CombinedModel(
        DL=Janus,
        mlp_dims=(rd["mlp_in"], rd["mlp_out"], rd["mlp_hid"]),
        dl_dims=(rd["dl_in"], rd["dl_out"], rd["dl_mod"]),
        combined_hidden=rd["com_hid"],
        final_output=1,
    )
    if ssm_class is not Janus:
        m.dl = ssm_class(rd["dl_in"], rd["dl_out"], rd["dl_mod"])
    return m


# ── Sequence masking ─────────────────────────────────────────────────────────

def _apply_mask(seq, pos_slice=None, chan_slice=None):
    out = seq.clone()
    ps = pos_slice if pos_slice is not None else slice(None)
    cs = chan_slice if chan_slice is not None else slice(None)
    out[:, ps, cs] = 0.0
    return out


# ── Training ─────────────────────────────────────────────────────────────────

def _train(model, train_loader, val_loader, device, criterion, maximize_metric=False, desc=""):
    model.to(device)
    opt = torch.optim.AdamW(model.parameters(), lr=LR, weight_decay=WEIGHT_DECAY)
    best_val = float("-inf") if maximize_metric else float("inf")
    best_state = None
    patience_left = PATIENCE

    pbar = tqdm(range(NUM_EPOCHS), desc=desc, ncols=90, leave=False)
    for epoch in pbar:
        model.train()
        for seq_b, feat_b, y_b in train_loader:
            seq_b, feat_b, y_b = seq_b.to(device), feat_b.to(device), y_b.to(device).float()
            loss = criterion(model(feat_b, seq_b).squeeze(), y_b)
            opt.zero_grad()
            loss.backward()
            opt.step()

        model.eval()
        val_losses, val_preds, val_true = [], [], []
        with torch.no_grad():
            for seq_b, feat_b, y_b in val_loader:
                seq_b, feat_b, y_b = seq_b.to(device), feat_b.to(device), y_b.to(device).float()
                preds = model(feat_b, seq_b).squeeze()
                val_losses.append(criterion(preds, y_b).item())
                val_preds.append(preds.cpu().numpy())
                val_true.append(y_b.cpu().numpy())

        val_preds = np.concatenate(val_preds)
        val_true = np.concatenate(val_true)
        if maximize_metric:
            metric = r2_score(val_true, val_preds)
            improved = metric > best_val
            pbar.set_postfix(val_r2=f"{metric:.4f}", best=f"{best_val:.4f}")
        else:
            metric = sum(val_losses)
            improved = metric < best_val
            pbar.set_postfix(val_loss=f"{metric:.4f}", best=f"{best_val:.4f}")

        if improved:
            best_val = metric
            best_state = {k: v.cpu().clone() for k, v in model.state_dict().items()}
            patience_left = PATIENCE
        else:
            patience_left -= 1
            if patience_left == 0:
                pbar.set_description(f"{desc} (early stop @{epoch + 1})")
                break

    pbar.close()
    model.load_state_dict(best_state)
    return model.to(device)


def train_classifier(model, train_loader, val_loader, device, desc="classifier"):
    return _train(model, train_loader, val_loader, device, nn.BCELoss(), maximize_metric=False, desc=desc)


def train_regressor(model, train_loader, val_loader, device, desc="regressor"):
    return _train(model, train_loader, val_loader, device, nn.MSELoss(), maximize_metric=True, desc=desc)


# ── Evaluation ────────────────────────────────────────────────────────────────

def evaluate(cls_model, reg_model, loader, device):
    cls_model.eval()
    reg_model.eval()
    cls_preds, reg_preds, ys = [], [], []
    with torch.no_grad():
        for seq_b, feat_b, y_b in loader:
            seq_b, feat_b = seq_b.to(device), feat_b.to(device)
            cls_preds.append(cls_model(feat_b, seq_b).squeeze().cpu().numpy())
            reg_preds.append(reg_model(feat_b, seq_b).squeeze().cpu().numpy())
            ys.append(y_b.numpy())
    cls_preds = np.concatenate(cls_preds)
    reg_preds = np.concatenate(reg_preds)
    y = np.concatenate(ys)
    return roc_auc_score((y > 0).astype(int), cls_preds), r2_score(y, reg_preds)


def make_loaders(feat_tr, seq_tr, y_cls_tr, y_reg_tr,
                 feat_val, seq_val, y_cls_val, y_reg_val,
                 feat_te, seq_te, y_cls_te, y_reg_te, seq_mask=None):
    if seq_mask is not None:
        seq_tr  = seq_mask(seq_tr)
        seq_val = seq_mask(seq_val)
        seq_te  = seq_mask(seq_te)

    def _loader(seq, feat, y, shuffle):
        return DataLoader(TensorDataset(seq, feat, y), batch_size=BATCH_SIZE, shuffle=shuffle)

    cls_loaders = (
        _loader(seq_tr,  feat_tr,  y_cls_tr,  shuffle=True),
        _loader(seq_val, feat_val, y_cls_val, shuffle=False),
        _loader(seq_te,  feat_te,  y_cls_te,  shuffle=False),
    )
    reg_loaders = (
        _loader(seq_tr,  feat_tr,  y_reg_tr,  shuffle=True),
        _loader(seq_val, feat_val, y_reg_val, shuffle=False),
        _loader(seq_te,  feat_te,  y_reg_te,  shuffle=False),
    )
    return cls_loaders, reg_loaders


# ── Main ─────────────────────────────────────────────────────────────────────

def main():
    rename = {**FEAT_RENAME, **SEQ_RENAME}

    train_df = pd.read_csv("training/0716_dataset_train.csv").rename(columns=rename)
    val_df   = pd.read_csv("training/0716_dataset_valid.csv").rename(columns=rename)
    test_df  = pd.read_csv("training/0716_dataset_test.csv").rename(columns=rename)
    print(f"Train: {len(train_df)}, Val: {len(val_df)}, Test: {len(test_df)}")

    scaler = StandardScaler()
    feat_tr  = torch.tensor(scaler.fit_transform(train_df[FEATURE_COLUMNS]), dtype=torch.float32)
    feat_val = torch.tensor(scaler.transform(val_df[FEATURE_COLUMNS]),       dtype=torch.float32)
    feat_te  = torch.tensor(scaler.transform(test_df[FEATURE_COLUMNS]),      dtype=torch.float32)

    seq_tr  = encode_batch_parallel(train_df[SEQUENCE_COLUMNS]).float()
    seq_val = encode_batch_parallel(val_df[SEQUENCE_COLUMNS]).float()
    seq_te  = encode_batch_parallel(test_df[SEQUENCE_COLUMNS]).float()

    y_tr  = torch.tensor(train_df["score"].values, dtype=torch.float32)
    y_val = torch.tensor(val_df["score"].values,   dtype=torch.float32)
    y_te  = torch.tensor(test_df["score"].values,  dtype=torch.float32)

    y_tr_bin  = (y_tr  > 0).float()
    y_val_bin = (y_val > 0).float()
    y_te_bin  = (y_te  > 0).float()

    _, baseline_cls, baseline_reg, device = load_models(device="cuda")
    print(f"Device: {device}\n")

    cd = _cls_dims(baseline_cls)
    rd = _reg_dims(baseline_reg)

    # Evaluate pre-trained baseline for reference
    _, base_reg_loaders = make_loaders(
        feat_tr, seq_tr, y_tr_bin, y_tr,
        feat_val, seq_val, y_val_bin, y_val,
        feat_te, seq_te, y_te_bin, y_te,
    )
    base_auroc, base_r2 = evaluate(baseline_cls, baseline_reg, base_reg_loaders[2], device)
    print(f"Pre-trained baseline: AUROC={base_auroc:.4f}, R²={base_r2:.4f}\n")

    RUNS = {
        "baseline_retrain": {
            "cls": lambda: build_classifier(cd),
            "reg": lambda: build_regressor(rd),
        },
        "no_pgc2": {
            "cls": lambda: build_classifier(cd, ssm_class=JanusNoPGC2),
            "reg": lambda: build_regressor(rd, ssm_class=JanusNoPGC2),
        },
        "no_s4d": {
            "cls": lambda: build_classifier(cd, ssm_class=JanusNoS4D),
            "reg": lambda: build_regressor(rd, ssm_class=JanusNoS4D),
        },
        "no_primer_enc": {
            "seq_mask": lambda s: _apply_mask(s, chan_slice=slice(5, 10)),
            "cls": lambda: build_classifier(cd),
            "reg": lambda: build_regressor(rd),
        },
        "no_target_enc": {
            "seq_mask": lambda s: _apply_mask(s, chan_slice=slice(0, 5)),
            "cls": lambda: build_classifier(cd),
            "reg": lambda: build_regressor(rd),
        },
        "no_fwd_seq": {
            "seq_mask": lambda s: _apply_mask(s, pos_slice=slice(0, 28)),
            "cls": lambda: build_classifier(cd),
            "reg": lambda: build_regressor(rd),
        },
        "no_rev_seq": {
            "seq_mask": lambda s: _apply_mask(s, pos_slice=slice(28, 56)),
            "cls": lambda: build_classifier(cd),
            "reg": lambda: build_regressor(rd),
        },
        "mlp_only": {
            "cls": lambda: MLPOnlyClassifier(cd["mlp_in"], cd["mlp_hid"], cd["mlp_out"], cd["com_hid"]),
            "reg": lambda: MLPOnlyRegressor(rd["mlp_in"], rd["mlp_hid"], rd["mlp_out"], rd["com_hid"]),
        },
        "ssm_only": {
            "cls": lambda: SSMOnlyClassifier(cd["ssm_in"], cd["ssm_out"], cd["ssm_mod"], cd["com_hid"]),
            "reg": lambda: SSMOnlyRegressor(rd["dl_in"], rd["dl_out"], rd["dl_mod"], rd["com_hid"]),
        },
    }

    results = {}
    for name, cfg in RUNS.items():
        seq_mask = cfg.get("seq_mask", None)
        cls_loaders, reg_loaders = make_loaders(
            feat_tr, seq_tr, y_tr_bin, y_tr,
            feat_val, seq_val, y_val_bin, y_val,
            feat_te, seq_te, y_te_bin, y_te,
            seq_mask=seq_mask,
        )

        cls_model = cfg["cls"]()
        reg_model = cfg["reg"]()

        cls_model = train_classifier(cls_model, cls_loaders[0], cls_loaders[1], device, desc=f"{name} cls")
        reg_model = train_regressor(reg_model, reg_loaders[0], reg_loaders[1], device, desc=f"{name} reg")

        auroc, r2 = evaluate(cls_model, reg_model, reg_loaders[2], device)
        results[name] = (auroc, r2)
        print(f"  AUROC={auroc:.4f} ({auroc - base_auroc:+.4f})  R²={r2:.4f} ({r2 - base_r2:+.4f})\n")

    # ── Summary table ────────────────────────────────────────────────────────
    print(f"\n{'=' * 65}")
    print(f"  RETRAIN ABLATION SUMMARY  (reference: pre-trained baseline)")
    print(f"{'=' * 65}")
    print(f"  {'Condition':<22s} {'AUROC':>7s} {'ΔAUROC':>8s} {'R²':>7s} {'ΔR²':>8s}")
    print(f"  {'-' * 55}")
    print(f"  {'pretrained_baseline':<22s} {base_auroc:>7.4f} {'':>8s} {base_r2:>7.4f} {'':>8s}")
    for name, (auroc, r2) in results.items():
        da = auroc - base_auroc
        dr = r2 - base_r2
        print(f"  {name:<22s} {auroc:>7.4f} {da:>+8.4f} {r2:>7.4f} {dr:>+8.4f}")

    rows = [{"condition": "pretrained_baseline", "auroc": base_auroc, "r2": base_r2,
             "delta_auroc": 0.0, "delta_r2": 0.0}]
    for name, (auroc, r2) in results.items():
        rows.append({"condition": name, "auroc": auroc, "r2": r2,
                     "delta_auroc": auroc - base_auroc, "delta_r2": r2 - base_r2})
    pd.DataFrame(rows).to_csv("training/ablation_retrain_results.csv", index=False)
    print("\nResults saved to training/ablation_retrain_results.csv")


if __name__ == "__main__":
    main()
