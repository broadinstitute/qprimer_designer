"""Inference-time ablation study for the classifier and regressor.

Runs four ablation studies and prints tabular results:
1. Branch ablation (MLP-only vs SSM-only vs full)
2. Layer ablation (skip individual SSM layers)
3. Feature ablation (zero numerical feature groups)
4. Input ablation (zero sequence channels/positions)
"""

import numpy as np
import pandas as pd
import torch
from torch.utils.data import DataLoader, TensorDataset
from sklearn.metrics import roc_auc_score, r2_score

from qprimer_designer.models import load_models, PcrDataset
from qprimer_designer.utils import encode_batch_parallel

# ── Column mapping from training CSVs to model expectations ─
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


def evaluate(cls_model, reg_model, feat, seq, y_bin, y_raw, batch_size=512):
    loader = DataLoader(TensorDataset(seq, feat), batch_size=batch_size, shuffle=False)
    cls_preds, reg_preds = [], []
    with torch.no_grad():
        for seq_b, feat_b in loader:
            cls_preds.append(cls_model(feat_b, seq_b).squeeze().cpu().numpy())
            reg_preds.append(reg_model(feat_b, seq_b).squeeze().cpu().numpy())
    cls_preds = np.concatenate(cls_preds)
    reg_preds = np.concatenate(reg_preds)
    return roc_auc_score(y_bin, cls_preds), r2_score(y_raw, reg_preds)


def print_table(title, results, base_auroc, base_r2):
    print(f"\n{'=' * 65}")
    print(f"  {title}")
    print(f"{'=' * 65}")
    print(f"  {'Condition':<22s} {'AUROC':>7s} {'ΔAUROC':>8s} {'R²':>7s} {'ΔR²':>8s}")
    print(f"  {'-' * 52}")
    print(f"  {'Baseline':<22s} {base_auroc:>7.4f} {'':>8s} {base_r2:>7.4f} {'':>8s}")
    for name, (auroc, r2) in results.items():
        da = auroc - base_auroc
        dr = r2 - base_r2
        print(f"  {name:<22s} {auroc:>7.4f} {da:>+8.4f} {r2:>7.4f} {dr:>+8.4f}")


def main():
    rename = {**FEAT_RENAME, **SEQ_RENAME}
    test_df = pd.read_csv("training/0716_dataset_test.csv").rename(columns=rename)
    scaler, classifier, regressor, device = load_models(device="cpu")
    classifier.eval()
    regressor.eval()

    features = scaler.transform(test_df[FEATURE_COLUMNS])
    sequences = encode_batch_parallel(test_df[SEQUENCE_COLUMNS])
    scores = test_df["score"].values
    y_bin = (scores > 0).astype(int)

    feat_tensor = torch.tensor(features, dtype=torch.float32)
    seq_tensor = sequences.float()

    base_auroc, base_r2 = evaluate(classifier, regressor, feat_tensor, seq_tensor, y_bin, scores)

    # ── 1. Branch ablation ──────────────────────────────────
    branch_results = {}

    def zero_hook(module, inp, out):
        return torch.zeros_like(out)

    # MLP-only: zero the SSM branch output
    hooks = [
        classifier.ssm.decoder.register_forward_hook(zero_hook),
        regressor.dl.decoder.register_forward_hook(zero_hook),
    ]
    branch_results["MLP-only (zero SSM)"] = evaluate(
        classifier, regressor, feat_tensor, seq_tensor, y_bin, scores
    )
    for h in hooks:
        h.remove()

    # SSM-only: zero the MLP branch output
    hooks = [
        classifier.mlp.model[2].register_forward_hook(zero_hook),
        regressor.mlp.model[2].register_forward_hook(zero_hook),
    ]
    branch_results["SSM-only (zero MLP)"] = evaluate(
        classifier, regressor, feat_tensor, seq_tensor, y_bin, scores
    )
    for h in hooks:
        h.remove()

    print_table("1. BRANCH ABLATION", branch_results, base_auroc, base_r2)

    # ── 2. Layer ablation (SSM internals) ───────────────────
    layer_results = {}

    LAYER_ABLATIONS = {
        "No encoder": ("encoder", "zero"),
        "No pgc1":    ("pgc1",    "skip"),
        "No pgc2":    ("pgc2",    "skip"),
        "No norm":    ("norm",    "skip"),
        "No s4d":     ("s4d",     "skip"),
    }

    def make_zero_hook():
        def hook(module, inp, out):
            return torch.zeros_like(out)
        return hook

    def make_skip_hook():
        def hook(module, inp, out):
            return inp[0]
        return hook

    for label, (layer_name, mode) in LAYER_ABLATIONS.items():
        hook_fn = make_zero_hook() if mode == "zero" else make_skip_hook()
        hooks = [
            getattr(classifier.ssm, layer_name).register_forward_hook(hook_fn),
            getattr(regressor.dl, layer_name).register_forward_hook(hook_fn),
        ]
        layer_results[label] = evaluate(
            classifier, regressor, feat_tensor, seq_tensor, y_bin, scores
        )
        for h in hooks:
            h.remove()

    print_table("2. LAYER ABLATION (SSM internals)", layer_results, base_auroc, base_r2)

    # ── 3. Feature ablation ─────────────────────────────────
    feat_idx = {name: i for i, name in enumerate(FEATURE_COLUMNS)}

    FEATURE_GROUPS = {
        "No fwd primer feats": ["len_f", "Tm_f", "GC_f", "indel_f", "mm_f"],
        "No rev primer feats": ["len_r", "Tm_r", "GC_r", "indel_r", "mm_r"],
        "No product feats":    ["prod_len", "prod_Tm"],
        "No Tm (all)":         ["Tm_f", "Tm_r", "prod_Tm"],
        "No GC (all)":         ["GC_f", "GC_r"],
        "No mismatches":       ["indel_f", "mm_f", "indel_r", "mm_r"],
    }

    feat_results = {}
    for label, cols in FEATURE_GROUPS.items():
        idxs = [feat_idx[c] for c in cols]
        ablated_feat = feat_tensor.clone()
        ablated_feat[:, idxs] = 0.0
        feat_results[label] = evaluate(
            classifier, regressor, ablated_feat, seq_tensor, y_bin, scores
        )

    print_table("3. FEATURE ABLATION (numerical features)", feat_results, base_auroc, base_r2)

    # ── 4. Input ablation (sequences) ───────────────────────
    # Shape: (N, 56, 10)
    # Positions: 0-27 = forward, 28-55 = reverse
    # Channels: 0-4 = target encoding, 5-9 = primer encoding
    INPUT_ABLATIONS = {
        "No fwd sequences":   (slice(0, 28),  slice(None)),
        "No rev sequences":   (slice(28, 56), slice(None)),
        "No primer encoding": (slice(None),   slice(5, 10)),
        "No target encoding": (slice(None),   slice(0, 5)),
    }

    input_results = {}
    for label, (pos_slice, chan_slice) in INPUT_ABLATIONS.items():
        ablated_seq = seq_tensor.clone()
        ablated_seq[:, pos_slice, chan_slice] = 0.0
        input_results[label] = evaluate(
            classifier, regressor, feat_tensor, ablated_seq, y_bin, scores
        )

    print_table("4. INPUT ABLATION (sequence channels)", input_results, base_auroc, base_r2)

    print(f"\n{'=' * 65}")
    print(f"  Baseline: AUROC={base_auroc:.4f}, R²={base_r2:.4f}")
    print(f"{'=' * 65}")


if __name__ == "__main__":
    main()
