"""Sequence manipulation utilities."""

import re

from Bio.SeqUtils import MeltingTemp, gc_fraction

# IUPAC ambiguity codes (everything except A, T, C, G, N)
_IUPAC_AMBIGUITY = re.compile(r'[RYWSMKBDHVryswmkbdhv]')

# Complement translation table
COMPLEMENT_TABLE = str.maketrans({
    'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C',
    'a': 't', 't': 'a', 'c': 'g', 'g': 'c',
    '-': '-', 'N': 'N', 'n': 'n'
})


def reverse_complement_dna(seq: str) -> str:
    """Return reverse-complement of a DNA sequence using fast translation table."""
    return seq.translate(COMPLEMENT_TABLE)[::-1]


def complement_dna(seq: str) -> str:
    """Return complement (not reverse) of a DNA sequence."""
    return seq.translate(COMPLEMENT_TABLE)


def get_tm(seq: str) -> float:
    """Calculate primer melting temperature (Primer3-like ionic conditions)."""
    return MeltingTemp.Tm_NN(seq, Na=50, Mg=1.5, dNTPs=0.6)


def get_gc_fraction(seq: str) -> float:
    """Calculate GC content as a fraction (0-1)."""
    return gc_fraction(seq)


def sanitize_iupac(seq: str) -> str:
    """Replace IUPAC ambiguity codes with N.

    RNAduplex and primer generation only support A, T, C, G, and N.
    Ambiguity codes (R, Y, W, S, M, K, B, D, H, V) are replaced with N
    so that primer/probe candidates overlapping those positions are
    filtered out by the existing 'N' checks.
    """
    return _IUPAC_AMBIGUITY.sub('N', seq)


def has_homopolymer(seq: str, max_len: int) -> bool:
    """Check if sequence has homopolymer run > max_len.

    Args:
        seq: DNA sequence to check
        max_len: Maximum allowed homopolymer length

    Returns:
        True if sequence contains homopolymer run exceeding max_len
    """
    count = 1
    for i in range(1, len(seq)):
        if seq[i] == seq[i-1]:
            count += 1
            if count > max_len:
                return True
        else:
            count = 1
    return False

