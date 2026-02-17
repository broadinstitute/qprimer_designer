"""Tests for generate_probe command logic."""

import pytest
from unittest.mock import patch

from qprimer_designer.commands.generate_probe import generate_probes


class TestGenerateProbes:
    """Tests for generate_probes (with mocked dG computation)."""

    @patch("qprimer_designer.commands.generate_probe.compute_self_dimer_dg", return_value=0.0)
    def test_basic_generation(self, mock_dg):
        """Test basic probe generation from a target sequence."""
        target_seqs = ["ATCGATCGATCGATCGATCGATCGATCGATCGATCGATCG"]
        filtered, features = generate_probes(
            target_seqs,
            len_min=8, len_max=10,
            max_tm=100, min_tm=0,  # Wide range to not filter by Tm
            max_gc=1.0, min_dg=-100,
            homopolymer_max=10,
            avoid_5prime_G=False,
            max_num=10000,
        )
        assert len(filtered) > 0
        for seq in filtered:
            assert 8 <= len(seq) <= 10
            assert "N" not in seq

    @patch("qprimer_designer.commands.generate_probe.compute_self_dimer_dg", return_value=0.0)
    def test_avoids_n_probes(self, mock_dg):
        """Probes with N should be excluded."""
        target_seqs = ["ATCGNNNATCGATCGATCGATCGATCGATCG"]
        filtered, features = generate_probes(
            target_seqs,
            len_min=8, len_max=10,
            max_tm=100, min_tm=0,
            max_gc=1.0, min_dg=-100,
            homopolymer_max=10,
            avoid_5prime_G=False,
            max_num=10000,
        )
        for seq in filtered:
            assert "N" not in seq

    @patch("qprimer_designer.commands.generate_probe.compute_self_dimer_dg", return_value=0.0)
    def test_avoid_5prime_g(self, mock_dg):
        """Probes starting with G should be excluded when avoid_5prime_G=True."""
        target_seqs = ["ATCGATCGATCGATCGATCGATCGATCGATCGATCGATCG"]
        filtered, features = generate_probes(
            target_seqs,
            len_min=8, len_max=10,
            max_tm=100, min_tm=0,
            max_gc=1.0, min_dg=-100,
            homopolymer_max=10,
            avoid_5prime_G=True,
            max_num=10000,
        )
        for seq in filtered:
            assert seq[0] != 'G'

    @patch("qprimer_designer.commands.generate_probe.compute_self_dimer_dg", return_value=0.0)
    def test_homopolymer_filter(self, mock_dg):
        """Probes with long homopolymer runs should be excluded."""
        # Create a target with a long run of A's
        target_seqs = ["AAAAAAAAAATCGATCGATCGATCGATCGATCGATCGATCG"]
        filtered, features = generate_probes(
            target_seqs,
            len_min=8, len_max=10,
            max_tm=100, min_tm=0,
            max_gc=1.0, min_dg=-100,
            homopolymer_max=3,
            avoid_5prime_G=False,
            max_num=10000,
        )
        # Any probe containing AAAA (4+ A's) should be excluded
        from qprimer_designer.utils.sequences import has_homopolymer
        for seq in filtered:
            assert not has_homopolymer(seq, 3)

    @patch("qprimer_designer.commands.generate_probe.compute_self_dimer_dg", return_value=0.0)
    def test_max_num_subsampling(self, mock_dg):
        """Output should be limited to max_num probes."""
        target_seqs = ["ATCGATCGATCGATCGATCGATCGATCGATCGATCGATCG"]
        filtered, features = generate_probes(
            target_seqs,
            len_min=8, len_max=10,
            max_tm=100, min_tm=0,
            max_gc=1.0, min_dg=-100,
            homopolymer_max=10,
            avoid_5prime_G=False,
            max_num=5,
        )
        assert len(filtered) <= 5

    @patch("qprimer_designer.commands.generate_probe.compute_self_dimer_dg", return_value=-100.0)
    def test_dg_filter(self, mock_dg):
        """Probes with dG below threshold should be excluded."""
        target_seqs = ["ATCGATCGATCGATCGATCGATCGATCGATCGATCGATCG"]
        filtered, _ = generate_probes(
            target_seqs,
            len_min=8, len_max=10,
            max_tm=100, min_tm=0,
            max_gc=1.0, min_dg=-6,  # All probes return dG=-100, which is < -6
            homopolymer_max=10,
            avoid_5prime_G=False,
            max_num=10000,
        )
        assert len(filtered) == 0

    @patch("qprimer_designer.commands.generate_probe.compute_self_dimer_dg", return_value=0.0)
    def test_features_populated(self, mock_dg):
        """Features dict should have Tm, GC, len, position, dG for each probe."""
        target_seqs = ["ATCGATCGATCGATCGATCGATCGATCGATCGATCGATCG"]
        filtered, features = generate_probes(
            target_seqs,
            len_min=8, len_max=8,
            max_tm=100, min_tm=0,
            max_gc=1.0, min_dg=-100,
            homopolymer_max=10,
            avoid_5prime_G=False,
            max_num=10000,
        )
        for seq in filtered:
            assert "Tm" in features[seq]
            assert "GC" in features[seq]
            assert "len" in features[seq]
            assert "position" in features[seq]
            assert "dG" in features[seq]
