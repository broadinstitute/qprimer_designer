"""Sequence manipulation utilities."""

from Bio.Seq import Seq
from Bio.SeqUtils import MeltingTemp, gc_fraction


def reverse_complement_dna(seq: str) -> str:
    """Return reverse-complement of a DNA sequence."""
    return str(Seq(seq).reverse_complement())


def get_tm(seq: str) -> float:
    """Calculate primer melting temperature (Primer3-like ionic conditions)."""
    return MeltingTemp.Tm_NN(seq, Na=50, Mg=1.5, dNTPs=0.6)


def get_gc_fraction(seq: str) -> float:
    """Calculate GC content as a fraction (0-1)."""
    return gc_fraction(seq)


# Complement translation table
COMPLEMENT_TABLE = str.maketrans({
    'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C',
    'a': 't', 't': 'a', 'c': 'g', 'g': 'c',
    '-': '-', 'N': 'N', 'n': 'n'
})


def complement_dna(seq: str) -> str:
    """Return complement (not reverse) of a DNA sequence."""
    return seq.translate(COMPLEMENT_TABLE)


def fast_reverse_complement(seq: str) -> str:
    """Fast reverse complement using translation table."""
    return seq.translate(COMPLEMENT_TABLE)[::-1]
