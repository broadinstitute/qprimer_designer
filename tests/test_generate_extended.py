"""Extended tests for generate command logic."""

from unittest.mock import patch

import pytest

from qprimer_designer.commands.generate import (
    generate_primers_single,
    generate_primers_multi,
    count_primer_pairs,
)


class TestGeneratePrimersMulti:
    """Tests for generate_primers_multi (with mocked dG)."""

    @patch("qprimer_designer.commands.generate.compute_self_dimer_dg", return_value=0.0)
    def test_basic_multi(self, mock_dg):
        """Test multi-target primer generation."""
        targets = ["ATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCG"]
        for_filt, rev_filt, features = generate_primers_multi(
            targets, step=1, min_pri_len=20, max_pri_len=20,
            min_amp_len=30, max_amp_len=200,
            max_tm=100, min_tm=0, max_gc=100, min_dg=-100,
        )
        assert len(for_filt) > 0
        assert len(rev_filt) > 0
        # Each filtered primer should have features
        for seq in for_filt:
            assert "Tm" in features[seq]
            assert "GC" in features[seq]
            assert "len" in features[seq]
            assert "dG" in features[seq]

    @patch("qprimer_designer.commands.generate.compute_self_dimer_dg", return_value=-100.0)
    def test_dg_filter_removes_all(self, mock_dg):
        """All primers filtered when dG is below threshold."""
        targets = ["ATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCG"]
        for_filt, rev_filt, features = generate_primers_multi(
            targets, step=1, min_pri_len=20, max_pri_len=20,
            min_amp_len=30, max_amp_len=200,
            max_tm=100, min_tm=0, max_gc=100, min_dg=-6,
        )
        assert len(for_filt) == 0
        assert len(rev_filt) == 0

    @patch("qprimer_designer.commands.generate.compute_self_dimer_dg", return_value=0.0)
    def test_tm_filter(self, mock_dg):
        """Very narrow Tm range should filter out most primers."""
        targets = ["ATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCG"]
        for_filt, rev_filt, features = generate_primers_multi(
            targets, step=1, min_pri_len=20, max_pri_len=20,
            min_amp_len=30, max_amp_len=200,
            max_tm=0, min_tm=0,  # Impossible range
            max_gc=100, min_dg=-100,
        )
        # With Tm range [0, 0], almost no 20-mers will pass
        assert len(for_filt) == 0

    @patch("qprimer_designer.commands.generate.compute_self_dimer_dg", return_value=0.0)
    def test_multiple_targets(self, mock_dg):
        """Test with multiple target sequences."""
        targets = [
            "ATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCG",
            "GCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTA",
        ]
        for_filt, rev_filt, features = generate_primers_multi(
            targets, step=1, min_pri_len=20, max_pri_len=20,
            min_amp_len=30, max_amp_len=200,
            max_tm=100, min_tm=0, max_gc=100, min_dg=-100,
        )
        assert len(for_filt) > 0


class TestCountPrimerPairsExtended:
    """Extended tests for count_primer_pairs."""

    def test_multiple_valid_pairs(self):
        fors = {"A": 0, "B": 10, "C": 20}
        revs = {"X": 100, "Y": 150}
        count = count_primer_pairs(fors, revs, min_amp_len=60, max_amp_len=200)
        assert count > 0

    def test_max_amp_len_constraint(self):
        """Amplicons exceeding max_amp_len should not be counted."""
        fors = {"A": 0}
        revs = {"B": 500}
        count = count_primer_pairs(fors, revs, min_amp_len=60, max_amp_len=200)
        assert count == 0
