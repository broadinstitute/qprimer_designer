#!/usr/bin/env python
# coding: utf-8

"""
evaluate_primers.py

Purpose
-------
Evaluate primer–target pairs using trained ML models (classifier + regressor).
Produces:

1) Raw per-pair scores
2) Aggregated per-primer-pair scores
3) Final ranked table with coverage, activity, and score

Key assumptions
---------------
- Input CSV is produced by prepare_input.py
- Sequence length for one-hot encoding is fixed (default: 28)
- Classifier output > 0.5 is considered a "hit"
- ON-target and OFF-target modes differ in aggregation logic
"""

# ------------------------------------------------------------
# Imports
# ------------------------------------------------------------

import argparse
import time
import os
import sys
import ast
import warnings
from collections import defaultdict
from multiprocessing import Pool

import numpy as np
import pandas as pd
import joblib
import torch
import torch.nn as nn
from torch.utils.data import Dataset, DataLoader
from einops import rearrange, repeat
from Bio import SeqIO
from sklearn.preprocessing import MultiLabelBinarizer

torch.manual_seed(42)


# ------------------------------------------------------------
# Encoding utilities
# ------------------------------------------------------------

def one_hot_encode(seq: str, length: int = 28) -> np.ndarray:
    """
    One-hot encode a nucleotide sequence with gap support.

    Encoding order: A, T, C, G, gap
    Sequence is left-padded/truncated to fixed length.
    """
    mapping = {
        "A": [1, 0, 0, 0, 0],
        "T": [0, 1, 0, 0, 0],
        "C": [0, 0, 1, 0, 0],
        "G": [0, 0, 0, 1, 0],
        "N": [0, 0, 0, 0, 0],
        "-": [0, 0, 0, 0, 1],
    }
    seq = seq.ljust(length, "N")[-length:]
    return np.array([mapping[c.upper()] for c in seq])


def encode_row(row: dict) -> np.ndarray:
    """
    Encode one primer–target pair into a (56, 10) tensor.
    """
    fenc = one_hot_encode(row["pseq_f"])
    ftenc = one_hot_encode(row["tseq_f"])
    renc = one_hot_encode(row["pseq_r"])
    rtenc = one_hot_encode(row["tseq_r"])

    prienc = np.append(fenc, renc, axis=0)
    tarenc = np.append(ftenc, rtenc, axis=0)
    return np.append(tarenc, prienc, axis=1)


def one_hot_encode_pbs_gap_parallel(df_seqs: pd.DataFrame, threads: int) -> torch.Tensor:
    """
    Parallel one-hot encoding of sequence columns using multiprocessing.
    """
    rows = df_seqs.to_dict("records")
    with Pool(processes=threads) as pool:
        encoded = pool.map(encode_row, rows)
    return torch.tensor(np.array(encoded), dtype=torch.float32)


# ------------------------------------------------------------
# Model components (unchanged, documented)
# ------------------------------------------------------------

class PGC(nn.Module):
    """Position-wise gated convolution block."""
    def __init__(self, model_dim, expansion_factor=1.0, dropout=0.0):
        super().__init__()
        self.conv = nn.Conv1d(
            int(model_dim * expansion_factor),
            int(model_dim * expansion_factor),
            kernel_size=3,
            padding=1,
            groups=int(model_dim * expansion_factor),
        )
        self.in_proj = nn.Linear(model_dim, int(model_dim * expansion_factor * 2))
        self.in_norm = nn.RMSNorm(int(model_dim * expansion_factor * 2), eps=1e-8)
        self.out_norm = nn.RMSNorm(model_dim, eps=1e-8)
        self.out_proj = nn.Linear(int(model_dim * expansion_factor), model_dim)
        self.dropout = nn.Dropout(dropout)

    def forward(self, u):
        xv = self.in_norm(self.in_proj(u))
        x, v = xv.chunk(2, dim=-1)
        x = self.conv(x.transpose(-1, -2)).transpose(-1, -2)
        return self.out_norm(self.out_proj(v * x))


class DropoutNd(nn.Module):
    def __init__(self, p: float = 0.5, tie=True, transposed=True):
        """
        tie: tie dropout mask across sequence lengths (Dropout1d/2d/3d)
        """
        super().__init__()
        if p < 0 or p >= 1:
            raise ValueError("dropout probability has to be in [0, 1), " "but got {}".format(p))
        self.p = p
        self.tie = tie
        self.transposed = transposed
        self.binomial = torch.distributions.binomial.Binomial(probs=1-self.p)

    def forward(self, X):
        """X: (batch, dim, lengths...)."""
        if self.training:
            if not self.transposed: X = rearrange(X, 'b ... d -> b d ...')
            # binomial = torch.distributions.binomial.Binomial(probs=1-self.p)
            # This is incredibly slow because of CPU -> GPU copying
            mask_shape = X.shape[:2] + (1,)*(X.ndim-2) if self.tie else X.shape
            # mask = self.binomial.sample(mask_shape)
            mask = torch.rand(*mask_shape, device=X.device) < 1.-self.p
            X = X * mask * (1.0/(1-self.p))
            if not self.transposed: X = rearrange(X, 'b d ... -> b ... d')
            return X
        return X


class S4DKernel(nn.Module):
    """Generate convolution kernel from diagonal SSM parameters."""

    def __init__(self, model_dim, N=64, dt_min=0.001, dt_max=0.1, lr=None):
        super().__init__()
        # Generate dt
        H = model_dim
        log_dt = torch.rand(H) * (
            math.log(dt_max) - math.log(dt_min)
        ) + math.log(dt_min)

        C = torch.randn(H, N // 2, dtype=torch.cfloat)
        self.C = nn.Parameter(torch.view_as_real(C))
        self.register("log_dt", log_dt, lr)

        log_A_real = torch.log(0.5 * torch.ones(H, N//2))
        A_imag = math.pi * repeat(torch.arange(N//2), 'n -> h n', h=H)
        self.register("log_A_real", log_A_real, lr)
        self.register("A_imag", A_imag, lr)

    def forward(self, L):
        """
        returns: (..., c, L) where c is number of channels (default 1)
        """

        # Materialize parameters
        dt = torch.exp(self.log_dt) # (H)
        C = torch.view_as_complex(self.C) # (H N)
        A = -torch.exp(self.log_A_real) + 1j * self.A_imag # (H N)

        # Vandermonde multiplication
        dtA = A * dt.unsqueeze(-1)  # (H N)
        K = dtA.unsqueeze(-1) * torch.arange(L, device=A.device) # (H N L)
        C = C * (torch.exp(dtA)-1.) / A
        K = 2 * torch.einsum('hn, hnl -> hl', C, torch.exp(K)).real

        return K

    def register(self, name, tensor, lr=None):
        """Register a tensor with a configurable learning rate and 0 weight decay"""

        if lr == 0.0:
            self.register_buffer(name, tensor)
        else:
            self.register_parameter(name, nn.Parameter(tensor))

            optim = {"weight_decay": 0.0}
            if lr is not None: optim["lr"] = lr
            setattr(getattr(self, name), "_optim", optim)


class S4D(nn.Module):
    def __init__(self, model_dim, state_dim=64, dropout=0.0, transposed=True, **kernel_args):
        super().__init__()

        self.h = model_dim
        self.n = state_dim
        self.output_dim = self.h
        self.transposed = transposed

        self.D = nn.Parameter(torch.randn(self.h))
        # SSM Kernel
        self.kernel = S4DKernel(self.h, N=self.n, **kernel_args)
        # Pointwise
        self.activation = nn.GELU()
        dropout_fn = DropoutNd
        self.dropout = dropout_fn(dropout) if dropout > 0.0 else nn.Identity()

        # position-wise output transform to mix features
        self.output_linear = nn.Sequential(
            nn.Conv1d(self.h, 2*self.h, kernel_size=1),
            nn.GLU(dim=-2),
        )

    def forward(self, u, **kwargs): # absorbs return_output and transformer src mask
        """ Input and output shape (B, H, L) """
        if not self.transposed: u = u.transpose(-1, -2)
        L = u.size(-1)
        # Compute SSM Kernel
        k = self.kernel(L=L) # (H L)

        # Convolution
        k_f = torch.fft.rfft(k, n=2*L)  # (H L)
        u_f = torch.fft.rfft(u, n=2*L) # (B H L)
        y = torch.fft.irfft(u_f*k_f, n=2*L)[..., :L] # (B H L)

        # Compute D term in state space equation - essentially a skip connection
        y = y + u * self.D.unsqueeze(-1)

        y = self.dropout(self.activation(y))
        y = self.output_linear(y)
        if not self.transposed: y = y.transpose(-1, -2)
        return y


class Janus(nn.Module):
    def __init__(self, input_dim, output_dim, model_dim, state_dim=64, dropout=0.2, transposed=False, **kernel_args):
        super().__init__()
        self.encoder = nn.Linear(input_dim, model_dim)
        self.pgc1 = PGC(model_dim, expansion_factor=0.25, dropout=dropout)
        self.pgc2 = PGC(model_dim, expansion_factor=2, dropout=dropout)
        self.s4d = S4D(model_dim, state_dim=state_dim, dropout=dropout, transposed=transposed, **kernel_args)
        self.norm = nn.RMSNorm(model_dim)
        self.decoder = nn.Linear(model_dim, output_dim)
        self.dropout = nn.Dropout(dropout)

    def forward(self, u):
        x = self.encoder(u)
        x = self.pgc1(x)
        x = self.pgc2(x)
        z = x
        z = self.norm(z)
        x = self.dropout(self.s4d(z)) + x
        x = x.mean(dim=1)
        #x = self.dropout(x)
        x = self.decoder(x)
        return x


class MLP(nn.Module):
    def __init__(self, input_dim, output_dim, hidden_dim):
        super(MLP, self).__init__()
        self.model = nn.Sequential(
            nn.Linear(input_dim, hidden_dim),
            nn.ReLU(),
            nn.Linear(hidden_dim, output_dim)
        )
    def forward(self, x):
        return self.model(x)
   

class CombinedModel(nn.Module):
    def __init__(self, DL, mlp_dims, dl_dims, combined_hidden, final_output):
        super(CombinedModel, self).__init__()

        # Individual models
        self.mlp = MLP(*mlp_dims)
        self.dl = DL(*dl_dims)

        # Combining ml
        combined_input_dim = mlp_dims[1] + dl_dims[1]
        self.combiner = nn.Sequential(
            nn.Linear(combined_input_dim, combined_hidden),
            nn.ReLU(),
            nn.Linear(combined_hidden, final_output)
        )

    def forward(self, mlp_input, dl_input):
        mlp_out = self.mlp(mlp_input.float())  # Output from ml
        dl_out = self.dl(dl_input)  # Output from dl

        # Concatenate outputs
        combined = torch.cat((mlp_out, dl_out), dim=1)

        # Final prediction
        final_output = self.combiner(combined)
        return final_output
    

class CombinedModelClassifier(nn.Module):
    def __init__(self, mlp_dims, ssm_dims, combined_hidden, final_output):
        super(CombinedModelClassifier, self).__init__()

        # Individual models
        self.mlp = MLP(*mlp_dims)
        self.ssm = Janus(*ssm_dims)

        # Combining MLP
        combined_input_dim = mlp_dims[1] + ssm_dims[1]
        self.combiner = nn.Sequential(
            nn.Linear(combined_input_dim, combined_hidden),
            nn.ReLU(),
            nn.Linear(combined_hidden, final_output),
            nn.Sigmoid()
        )

    def forward(self, mlp_input, ssm_input):
        mlp_out = self.mlp(mlp_input.float())  # Output from MLP
        ssm_out = self.ssm(ssm_input)  # Output from SSM

        # Concatenate outputs
        combined = torch.cat((mlp_out, ssm_out), dim=1)

        # Final prediction
        final_output = self.combiner(combined)
        return final_output
    

class PcrDataset(Dataset):
    """Simple dataset wrapper for ML inference."""
    def __init__(self, encoded_input, custom_features, scores):
        self.encoded_input = encoded_input
        self.custom_features = custom_features
        self.scores = scores

    def __len__(self):
        return len(self.encoded_input)

    def __getitem__(self, idx):
        return (
            self.encoded_input[idx],
            self.custom_features[idx],
            self.scores[idx],
        )

# ------------------------------------------------------------
# Main
# ------------------------------------------------------------

def main():
    parser = argparse.ArgumentParser(
        prog="python -u evaluate_primers.py",
        description="Evaluate primer–target matches using ML",
    )
    parser.add_argument("--in", dest="input", required=True)
    parser.add_argument("--out", dest="output", required=True)
    parser.add_argument("--ref", dest="reference", required=True)
    parser.add_argument("--reftype", dest="reftype", required=True, choices=["on", "off"])
    parser.add_argument("--ml", dest="ml_dir", required=True)
    parser.add_argument("--threads", dest="threads", default=1, type=int)

    args = parser.parse_args()

    tnames = [s.id for s in SeqIO.parse(args.reference, "fasta")]

    # Load models and scaler
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        scaler = joblib.load(f"{args.ml_dir}/standard_scaler.joblib")

    device = torch.device("cuda" if torch.cuda.is_available() else "cpu")

    classifier = torch.load(f"{args.ml_dir}/combined_classifier.pth", map_location=device, weights_only=False)
    regressor = torch.load(f"{args.ml_dir}/combined_regressor.pth", map_location=device, weights_only=False)
    classifier.to(device).eval()
    regressor.to(device).eval()

    print(f"Evaluating {args.input} with {device}...")
    startTime = time.time()

    if os.path.getsize(args.input) == 0:
        open(args.output, "w").close()
        sys.exit()

    feats = [
        "len_f","Tm_f","GC_f","indel_f","mm_f",
        "len_r","Tm_r","GC_r","indel_r","mm_r",
        "prod_len","prod_Tm",
    ]

    header_flag = True
    mode = 'w'
    clstbl, regtbl = [], []

    for i, chunk in enumerate(pd.read_csv(args.input, chunksize=20000)):
        chunk['targets'] = chunk['targets'].apply(ast.literal_eval)

        inps_fe = chunk[feats]
        inps_fe = scaler.transform(inps_fe)
        inps_se = chunk[['pseq_f','tseq_f','pseq_r','tseq_r']]
        inps_se = one_hot_encode_pbs_gap_parallel(inps_se, int(args.threads))
        dataset = PcrDataset(inps_se, inps_fe, np.array([0]*len(chunk)))
        loader = DataLoader(dataset, batch_size=64, shuffle=False)

        predict_cls, predict_reg = [], []
        with torch.no_grad():
            for seq_in, fea_in, _ in loader:
                seq_in, fea_in = seq_in.to(device).float(), fea_in.to(device).float()
                out_cls = classifier(fea_in, seq_in)
                out_reg = regressor(fea_in, seq_in)
                if len(seq_in)==1:
                    predict_cls.append(np.array([out_cls.squeeze().detach().cpu().numpy()]))
                    predict_reg.append(np.array([out_reg.squeeze().detach().cpu().numpy()]))
                else:
                    predict_cls.append(out_cls.squeeze().detach().cpu().numpy())
                    predict_reg.append(out_reg.squeeze().detach().cpu().numpy())
            predict_cls = np.concatenate(predict_cls)
            predict_reg = np.round(np.concatenate(predict_reg), decimals=3)
        chunk.loc[:,'classifier'] = predict_cls
        chunk.loc[:,'regressor'] = predict_reg
        chunk.to_csv(f'{args.output}.full', mode=mode, header=header_flag, index=False)

        mlb = MultiLabelBinarizer()
        onehot = mlb.fit_transform(chunk['targets'])
        target_cols = list(mlb.classes_)
        
        for label, l in zip(['classifier','regressor'], [clstbl,regtbl]):
            targets_df = pd.DataFrame(onehot, columns=target_cols, index=chunk.index)
            targets_df = targets_df.mul(chunk[label], axis=0)
            evaltbl = pd.concat([chunk[['pname_f','pname_r']], targets_df], axis=1)
            agg_dict = {c: "max" for c in evaltbl.columns[2:]}
            evaltbl = evaltbl.groupby(['pname_f','pname_r']).agg(agg_dict).reset_index()
            l.append(evaltbl)
        header_flag = False
        mode = 'a'
    
    clstbl = pd.concat(clstbl, ignore_index=True)
    regtbl = pd.concat(regtbl, ignore_index=True).reindex(columns=clstbl.columns)
    agg_dict = {c: "max" for c in regtbl.columns[2:]}
    regtbl = regtbl.reset_index().groupby(['pname_f','pname_r']).agg(agg_dict)
    clstbl = clstbl.reset_index().groupby(['pname_f','pname_r']).agg(agg_dict) 
    clstbl.fillna(0).round(3).to_csv(f'{args.output}.cl')
    regtbl.fillna(0).round(3).to_csv(f'{args.output}.re')

    coverage = (clstbl>.5).sum(axis=1).reset_index(name='coverage')
    if args.reftype=='on':
        coverage['coverage'] = coverage['coverage'] / len(tnames)
        activity = (regtbl * (clstbl>.5)).mean(axis=1).reset_index(name='activity')
    else:
        activity = (regtbl * (clstbl>.5)).max(axis=1).reset_index(name='activity')
    res = coverage.merge(activity, on=['pname_f','pname_r'])
    res['score'] = (res['coverage'] * res['activity'])
    res = res.dropna().sort_values('score', ascending=False)
    res.round(4).to_csv(args.output, index=False)

    runtime = time.time() - startTime
    print(f"Wrote {len(res)} lines to {args.output} (%.1f sec)" % runtime)


if __name__ == "__main__":
    main()

