"""Tests for sequence utility functions."""

import pytest
from qprimer_designer.utils.sequences import (
    reverse_complement_dna,
    get_gc_fraction,
)


class TestReverseComplement:
    """Tests for reverse complement functions."""

    def test_reverse_complement_simple(self):
        """Test reverse complement of a simple sequence."""
        assert reverse_complement_dna("ATCG") == "CGAT"

    def test_reverse_complement_empty(self):
        """Test reverse complement of empty string."""
        assert reverse_complement_dna("") == ""


class TestGCFraction:
    """Tests for GC content calculation."""

    def test_gc_fraction_all_gc(self):
        """Test GC fraction of all G/C sequence."""
        assert get_gc_fraction("GCGCGC") == 1.0

    def test_gc_fraction_no_gc(self):
        """Test GC fraction of all A/T sequence."""
        assert get_gc_fraction("ATATAT") == 0.0

    def test_gc_fraction_half(self):
        """Test GC fraction of 50% GC sequence."""
        assert get_gc_fraction("ATCG") == 0.5

    def test_gc_fraction_empty(self):
        """Test GC fraction of empty string."""
        assert get_gc_fraction("") == 0.0
