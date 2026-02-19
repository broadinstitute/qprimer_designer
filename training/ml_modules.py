#!/usr/bin/env python
# coding: utf-8

# In[1]:


import numpy as np
import torch
import torch.nn as nn
import torch.nn.functional as F
from torch.utils.data import Dataset, DataLoader
from einops import rearrange, repeat
from tqdm.auto import tqdm
torch.manual_seed(42)
device = torch.device('cuda' if torch.cuda.is_available() else "cpu")
device


# In[ ]:


class PGC(nn.Module):
    def __init__(self,model_dim,expansion_factor = 1.0,dropout = 0.0):
        super().__init__()
        self.model_dim = model_dim
        self.expansion_factor = expansion_factor
        self.dropout = dropout
        self.conv = nn.Conv1d(int(model_dim * expansion_factor), int(model_dim * expansion_factor),
                              kernel_size=3, padding=1, groups=int(model_dim * expansion_factor))
        self.in_proj = nn.Linear(model_dim, int(model_dim * expansion_factor * 2))
        self.out_norm = nn.RMSNorm(int(model_dim), eps=1e-8)
        self.in_norm = nn.RMSNorm(int(model_dim * expansion_factor * 2), eps=1e-8)
        self.out_proj = nn.Linear(int(model_dim * expansion_factor), model_dim)
        self.dropout = nn.Dropout(dropout)
    def forward(self, u):
        xv = self.in_norm(self.in_proj(u))
        x,v = xv.chunk(2,dim=-1)
        x_conv = self.conv(x.transpose(-1,-2)).transpose(-1,-2)
        gate =  v * x_conv
        x = self.out_norm(self.out_proj(gate))
        return x
    
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


# In[ ]:


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


# In[ ]:


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
        mlp_out = self.mlp(mlp_input)  # Output from ml
        dl_out = self.dl(dl_input)  # Output from dl

        # Concatenate outputs
        combined = torch.cat((mlp_out, dl_out), dim=1)

        # Final prediction
        final_output = self.combiner(combined)
        return final_output


# In[ ]:


class PcrDataset(Dataset):
    def __init__(self, encoded_input, custom_features, scores):
        self.encoded_input = encoded_input
        self.custom_features = custom_features
        self.scores = scores
    def __len__(self):
        return len(self.encoded_input)
    def __getitem__(self, idx):
        return self.encoded_input[idx], self.custom_features[idx], self.scores[idx]


# In[ ]:


class LSTMModel(nn.Module):
    def __init__(self, input_size, hidden_size, num_layers, output_size, bidirectional=False):
        super(LSTMModel, self).__init__()
        self.hidden_size = hidden_size
        self.num_layers = num_layers
        self.bidirectional = bidirectional
        self.num_directions = 2 if bidirectional else 1

        self.lstm = nn.LSTM(
            input_size,
            hidden_size,
            num_layers,
            batch_first=True,
            bidirectional=bidirectional
        )

        # Fully connected layer to change shape into output size
        self.fc = nn.Linear(hidden_size * self.num_directions, output_size)
        self.sigout = nn.Sigmoid()

    def forward(self, x):
        device = x.device

        h0 = torch.zeros(self.num_layers * self.num_directions, x.size(0), self.hidden_size).to(device)
        c0 = torch.zeros(self.num_layers * self.num_directions, x.size(0), self.hidden_size).to(device)

        out, _ = self.lstm(x, (h0, c0))  # out: (batch_size, seq_len, hidden_size * num_directions)

        out = self.fc(out[:, -1, :])  # Last time step
        out = self.sigout(out)
        return out


# In[ ]:


def rev_com_enc(enc):
    revcoms = {'A':'T',
               'T':'A',
               'C':'G',
               'G':'C',
               'a':'t',
               't':'a',
               'c':'g',
               'g':'c',
               '-':'-',
               'N':'N',
               'n':'n'
              }
    return ''.join([revcoms[char] for char in enc[::-1]])

def one_hot_encode(seq, length=28):
    mapping = { 'A':[1, 0, 0, 0, 0],
                'T':[0, 1, 0, 0, 0],
                'C':[0, 0, 1, 0, 0],
                'G':[0, 0, 0, 1, 0],
                'N':[0, 0, 0, 0, 0],
                '-':[0, 0, 0, 0, 1] }
    seq = seq.ljust(length, 'N') # (6, ATCG) -> NNATCG
    return np.array([mapping[char.upper()] for char in seq])


# In[ ]:


# def one_hot_encode_full_gap(df_seqs, maxl=1421):
#     primer_encoded = []
#     target_encoded = []
#     for (tname,pname),row in df_seqs.iterrows():
#         fseq, fst, rseq, rst, tseq = row[['f_seq','f_start','r_seq','r_start','target_seq']]
#         fenc, ftenc, renc, rtenc = row[['f_penc','f_tenc','r_penc','r_tenc']]
#         pseq = 'N'*fst + fenc + 'N'*(rst-(fst+len(fseq))) + renc + 'N'*(len(tseq)-(rst+len(rseq)))
#         tseq = tseq[:fst] + ftenc + tseq[fst+len(fseq):rst] + rtenc + tseq[rst+len(rseq):]
#         primer_encoded.append(one_hot_encode(pseq, maxl))
#         target_encoded.append(one_hot_encode(tseq, maxl))
#     final_encoded = np.append(np.array(target_encoded), np.array(primer_encoded), axis=2)
#     print(final_encoded.shape)
#     return torch.tensor(final_encoded, dtype=torch.float32)

# def one_hot_encode_pbs_gap(df_seqs):
#     primer_encoded = []
#     target_encoded = []
#     for (tname,pname),row in df_seqs.iterrows():
#         fenc, ftenc, renc, rtenc = row[['f_penc','f_tenc','r_penc','r_tenc']].apply(one_hot_encode)
#         prienc = np.append(fenc,renc,axis=0)
#         tarenc = np.append(ftenc,rtenc,axis=0)
#         primer_encoded.append(prienc)
#         target_encoded.append(tarenc)
#     primer_encoded = np.array(primer_encoded)
#     target_encoded = np.array(target_encoded)
#     final_encoded = np.append(target_encoded, primer_encoded, axis=2)
#     print(final_encoded.shape)
#     return torch.tensor(final_encoded, dtype=torch.float32)

