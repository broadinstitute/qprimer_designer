"""ML model architectures for primer evaluation.

These models are trained to predict primer-target binding affinity
using a combination of numerical features and sequence encodings.
"""

import math

import torch
import torch.nn as nn
from einops import rearrange, repeat
from torch.utils.data import Dataset


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
    """N-dimensional dropout with optional mask tying."""

    def __init__(self, p: float = 0.5, tie=True, transposed=True):
        super().__init__()
        if p < 0 or p >= 1:
            raise ValueError(f"dropout probability must be in [0, 1), got {p}")
        self.p = p
        self.tie = tie
        self.transposed = transposed

    def forward(self, X):
        if self.training:
            if not self.transposed:
                X = rearrange(X, 'b ... d -> b d ...')
            mask_shape = X.shape[:2] + (1,) * (X.ndim - 2) if self.tie else X.shape
            mask = torch.rand(*mask_shape, device=X.device) < 1.0 - self.p
            X = X * mask * (1.0 / (1 - self.p))
            if not self.transposed:
                X = rearrange(X, 'b d ... -> b ... d')
        return X


class S4DKernel(nn.Module):
    """Generate convolution kernel from diagonal SSM parameters."""

    def __init__(self, model_dim, N=64, dt_min=0.001, dt_max=0.1, lr=None):
        super().__init__()
        H = model_dim
        log_dt = torch.rand(H) * (math.log(dt_max) - math.log(dt_min)) + math.log(dt_min)

        C = torch.randn(H, N // 2, dtype=torch.cfloat)
        self.C = nn.Parameter(torch.view_as_real(C))
        self.register("log_dt", log_dt, lr)

        log_A_real = torch.log(0.5 * torch.ones(H, N // 2))
        A_imag = math.pi * repeat(torch.arange(N // 2), 'n -> h n', h=H)
        self.register("log_A_real", log_A_real, lr)
        self.register("A_imag", A_imag, lr)

    def forward(self, L):
        dt = torch.exp(self.log_dt)
        C = torch.view_as_complex(self.C)
        A = -torch.exp(self.log_A_real) + 1j * self.A_imag

        dtA = A * dt.unsqueeze(-1)
        K = dtA.unsqueeze(-1) * torch.arange(L, device=A.device)
        C = C * (torch.exp(dtA) - 1.0) / A
        K = 2 * torch.einsum('hn, hnl -> hl', C, torch.exp(K)).real

        return K

    def register(self, name, tensor, lr=None):
        if lr == 0.0:
            self.register_buffer(name, tensor)
        else:
            self.register_parameter(name, nn.Parameter(tensor))
            optim = {"weight_decay": 0.0}
            if lr is not None:
                optim["lr"] = lr
            setattr(getattr(self, name), "_optim", optim)


class S4D(nn.Module):
    """Diagonal State Space Model (S4D) layer."""

    def __init__(self, model_dim, state_dim=64, dropout=0.0, transposed=True, **kernel_args):
        super().__init__()
        self.h = model_dim
        self.n = state_dim
        self.output_dim = self.h
        self.transposed = transposed

        self.D = nn.Parameter(torch.randn(self.h))
        self.kernel = S4DKernel(self.h, N=self.n, **kernel_args)
        self.activation = nn.GELU()
        self.dropout = DropoutNd(dropout) if dropout > 0.0 else nn.Identity()

        self.output_linear = nn.Sequential(
            nn.Conv1d(self.h, 2 * self.h, kernel_size=1),
            nn.GLU(dim=-2),
        )

    def forward(self, u, **kwargs):
        if not self.transposed:
            u = u.transpose(-1, -2)
        L = u.size(-1)

        k = self.kernel(L=L)
        k_f = torch.fft.rfft(k, n=2 * L)
        u_f = torch.fft.rfft(u, n=2 * L)
        y = torch.fft.irfft(u_f * k_f, n=2 * L)[..., :L]

        y = y + u * self.D.unsqueeze(-1)
        y = self.dropout(self.activation(y))
        y = self.output_linear(y)

        if not self.transposed:
            y = y.transpose(-1, -2)
        return y


class Janus(nn.Module):
    """Janus architecture combining PGC and S4D for sequence modeling."""

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
        z = self.norm(x)
        x = self.dropout(self.s4d(z)) + x
        x = x.mean(dim=1)
        x = self.decoder(x)
        return x


class MLP(nn.Module):
    """Simple multi-layer perceptron."""

    def __init__(self, input_dim, output_dim, hidden_dim):
        super().__init__()
        self.model = nn.Sequential(
            nn.Linear(input_dim, hidden_dim),
            nn.ReLU(),
            nn.Linear(hidden_dim, output_dim),
        )

    def forward(self, x):
        return self.model(x)


class CombinedModel(nn.Module):
    """Combined MLP + Janus model for regression."""

    def __init__(self, DL, mlp_dims, dl_dims, combined_hidden, final_output):
        super().__init__()
        self.mlp = MLP(*mlp_dims)
        self.dl = DL(*dl_dims)

        combined_input_dim = mlp_dims[1] + dl_dims[1]
        self.combiner = nn.Sequential(
            nn.Linear(combined_input_dim, combined_hidden),
            nn.ReLU(),
            nn.Linear(combined_hidden, final_output),
        )

    def forward(self, mlp_input, dl_input):
        mlp_out = self.mlp(mlp_input.float())
        dl_out = self.dl(dl_input)
        combined = torch.cat((mlp_out, dl_out), dim=1)
        return self.combiner(combined)


class CombinedModelClassifier(nn.Module):
    """Combined MLP + Janus model for binary classification."""

    def __init__(self, mlp_dims, ssm_dims, combined_hidden, final_output):
        super().__init__()
        self.mlp = MLP(*mlp_dims)
        self.ssm = Janus(*ssm_dims)

        combined_input_dim = mlp_dims[1] + ssm_dims[1]
        self.combiner = nn.Sequential(
            nn.Linear(combined_input_dim, combined_hidden),
            nn.ReLU(),
            nn.Linear(combined_hidden, final_output),
            nn.Sigmoid(),
        )

    def forward(self, mlp_input, ssm_input):
        mlp_out = self.mlp(mlp_input.float())
        ssm_out = self.ssm(ssm_input)
        combined = torch.cat((mlp_out, ssm_out), dim=1)
        return self.combiner(combined)


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
