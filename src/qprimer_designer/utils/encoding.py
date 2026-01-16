"""Sequence encoding utilities for ML models."""

import numpy as np
from multiprocessing import Pool


def one_hot_encode(seq: str, length: int = 28) -> np.ndarray:
    """
    One-hot encode a nucleotide sequence with gap support.

    Encoding order: A, T, C, G, gap
    Sequence is left-padded/truncated to fixed length.

    Args:
        seq: DNA/RNA sequence
        length: Fixed output length (default: 28)

    Returns:
        numpy array of shape (length, 5)
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


def encode_primer_pair(row: dict) -> np.ndarray:
    """
    Encode one primer-target pair into a (56, 10) tensor.

    Args:
        row: Dictionary with keys: pseq_f, tseq_f, pseq_r, tseq_r

    Returns:
        numpy array of shape (56, 10)
    """
    fenc = one_hot_encode(row["pseq_f"])
    ftenc = one_hot_encode(row["tseq_f"])
    renc = one_hot_encode(row["pseq_r"])
    rtenc = one_hot_encode(row["tseq_r"])

    prienc = np.append(fenc, renc, axis=0)
    tarenc = np.append(ftenc, rtenc, axis=0)
    return np.append(tarenc, prienc, axis=1)


def encode_batch_parallel(df_seqs, threads: int = 1):
    """
    Parallel one-hot encoding of sequence columns using multiprocessing.

    Args:
        df_seqs: DataFrame with columns pseq_f, tseq_f, pseq_r, tseq_r
        threads: Number of parallel workers

    Returns:
        torch.Tensor of shape (n_samples, 56, 10)
    """
    import torch

    rows = df_seqs.to_dict("records")
    with Pool(processes=threads) as pool:
        encoded = pool.map(encode_primer_pair, rows)
    return torch.tensor(np.array(encoded), dtype=torch.float32)
