"""Extended tests for sequence utility functions."""

import pytest

from qprimer_designer.utils.sequences import (
    reverse_complement_dna,
    complement_dna,
    get_tm,
    get_gc_fraction,
    sanitize_iupac,
    has_homopolymer,
    COMPLEMENT_TABLE,
)


class TestComplementDna:
    """Tests for complement_dna (not reverse)."""

    def test_simple(self):
        assert complement_dna("ATCG") == "TAGC"

    def test_empty(self):
        assert complement_dna("") == ""

    def test_lowercase(self):
        assert complement_dna("atcg") == "tagc"

    def test_with_gap(self):
        assert complement_dna("A-T") == "T-A"

    def test_with_n(self):
        assert complement_dna("ANC") == "TNG"

    def test_all_same(self):
        assert complement_dna("AAAA") == "TTTT"

    def test_palindromic(self):
        """A palindromic sequence's complement reversed should equal itself."""
        seq = "ATAT"
        assert complement_dna(seq) == "TATA"


class TestReverseComplementExtended:
    """Extended tests for reverse_complement_dna."""

    def test_palindrome(self):
        """ACGT is its own reverse complement."""
        assert reverse_complement_dna("ACGT") == "ACGT"

    def test_lowercase(self):
        assert reverse_complement_dna("atcg") == "cgat"

    def test_with_gap(self):
        assert reverse_complement_dna("A-CG") == "CG-T"

    def test_with_n(self):
        assert reverse_complement_dna("ANG") == "CNT"

    def test_long_sequence(self):
        seq = "ATCGATCGATCGATCGATCG"
        rc = reverse_complement_dna(seq)
        # Reverse complement of reverse complement should be original
        assert reverse_complement_dna(rc) == seq

    def test_single_base(self):
        assert reverse_complement_dna("A") == "T"
        assert reverse_complement_dna("T") == "A"
        assert reverse_complement_dna("C") == "G"
        assert reverse_complement_dna("G") == "C"


class TestGetTm:
    """Tests for melting temperature calculation."""

    def test_returns_float(self):
        tm = get_tm("ATCGATCGATCGATCGATCG")
        assert isinstance(tm, float)

    def test_gc_rich_higher_tm(self):
        """GC-rich sequences should have higher Tm than AT-rich."""
        tm_gc = get_tm("GCGCGCGCGCGCGCGCGCGC")
        tm_at = get_tm("ATATATATATATATATATAT")
        assert tm_gc > tm_at

    def test_longer_higher_tm(self):
        """Longer sequences generally have higher Tm."""
        tm_short = get_tm("ATCGATCG")
        tm_long = get_tm("ATCGATCGATCGATCGATCG")
        assert tm_long > tm_short

    def test_reasonable_range(self):
        """A 20-mer primer should have Tm in a biologically reasonable range."""
        tm = get_tm("ATCGATCGATCGATCGATCG")
        assert 40.0 < tm < 80.0


class TestSanitizeIupac:
    """Tests for IUPAC ambiguity code sanitization."""

    def test_no_ambiguity(self):
        assert sanitize_iupac("ATCGATCG") == "ATCGATCG"

    def test_replace_r(self):
        assert sanitize_iupac("ARCG") == "ANCG"

    def test_replace_y(self):
        assert sanitize_iupac("AYCG") == "ANCG"

    def test_replace_w(self):
        assert sanitize_iupac("AWCG") == "ANCG"

    def test_replace_s(self):
        assert sanitize_iupac("ASCG") == "ANCG"

    def test_replace_multiple(self):
        assert sanitize_iupac("RYSW") == "NNNN"

    def test_lowercase_ambiguity(self):
        assert sanitize_iupac("arcg") == "aNcg"

    def test_preserves_n(self):
        assert sanitize_iupac("ANCG") == "ANCG"

    def test_empty(self):
        assert sanitize_iupac("") == ""

    def test_all_iupac_codes(self):
        """All IUPAC ambiguity codes should be replaced."""
        result = sanitize_iupac("RYSWMKBDHV")
        assert result == "N" * 10


class TestHasHomopolymer:
    """Tests for homopolymer detection."""

    def test_no_homopolymer(self):
        assert has_homopolymer("ATCGATCG", max_len=3) is False

    def test_has_homopolymer(self):
        assert has_homopolymer("ATCGAAAAT", max_len=3) is True

    def test_exact_max_len(self):
        """Run exactly at max_len should NOT trigger."""
        assert has_homopolymer("AAATCG", max_len=3) is False

    def test_one_over_max_len(self):
        """Run one over max_len SHOULD trigger."""
        assert has_homopolymer("AAAATCG", max_len=3) is True

    def test_max_len_1(self):
        """With max_len=1, any repeated base triggers."""
        assert has_homopolymer("AATCG", max_len=1) is True
        assert has_homopolymer("ATCG", max_len=1) is False

    def test_single_base(self):
        assert has_homopolymer("A", max_len=3) is False

    def test_empty(self):
        assert has_homopolymer("", max_len=3) is False

    def test_all_same(self):
        assert has_homopolymer("AAAAA", max_len=4) is True
        assert has_homopolymer("AAAAA", max_len=5) is False

    def test_end_of_sequence(self):
        """Homopolymer at end of sequence."""
        assert has_homopolymer("ATCGGGGG", max_len=4) is True

    def test_different_bases(self):
        """Should detect runs of any nucleotide."""
        assert has_homopolymer("TTTTTT", max_len=5) is True
        assert has_homopolymer("CCCCCC", max_len=5) is True
        assert has_homopolymer("GGGGGG", max_len=5) is True


class TestComplementTable:
    """Tests for the complement translation table."""

    def test_standard_bases(self):
        assert "A".translate(COMPLEMENT_TABLE) == "T"
        assert "T".translate(COMPLEMENT_TABLE) == "A"
        assert "C".translate(COMPLEMENT_TABLE) == "G"
        assert "G".translate(COMPLEMENT_TABLE) == "C"

    def test_lowercase_bases(self):
        assert "a".translate(COMPLEMENT_TABLE) == "t"
        assert "t".translate(COMPLEMENT_TABLE) == "a"
        assert "c".translate(COMPLEMENT_TABLE) == "g"
        assert "g".translate(COMPLEMENT_TABLE) == "c"

    def test_gap(self):
        assert "-".translate(COMPLEMENT_TABLE) == "-"

    def test_n(self):
        assert "N".translate(COMPLEMENT_TABLE) == "N"
        assert "n".translate(COMPLEMENT_TABLE) == "n"
