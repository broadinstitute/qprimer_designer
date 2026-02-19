"""Tests for filter_primers command logic."""

import pytest

from qprimer_designer.commands.filter_primers import _find_position_valid_probes


class TestFindPositionValidProbes:
    """Tests for _find_position_valid_probes."""

    def test_basic_valid_probe(self):
        """Probe within amplicon bounds should pass."""
        pair_full_rows = [
            {
                "targets": "['target1']",
                "starts": "[100]",
                "prod_len": 200,
            }
        ]
        probes_by_target = {"target1": {"probe_a"}}
        probe_positions = {("target1", "probe_a"): 150}
        probe_seqs_dict = {"probe_a": "ATCGATCGATCGATCGATCGATCG"}  # len 24
        buffer = 20

        result = _find_position_valid_probes(
            pair_full_rows, probes_by_target, probe_positions,
            probe_seqs_dict, buffer
        )
        assert "probe_a" in result

    def test_probe_outside_buffer(self):
        """Probe too close to amplicon edge should fail."""
        pair_full_rows = [
            {
                "targets": "['target1']",
                "starts": "[100]",
                "prod_len": 200,
            }
        ]
        probes_by_target = {"target1": {"probe_a"}}
        # Probe starts at position 105, but buffer is 20, so amp_start = 120
        probe_positions = {("target1", "probe_a"): 105}
        probe_seqs_dict = {"probe_a": "ATCG"}
        buffer = 20

        result = _find_position_valid_probes(
            pair_full_rows, probes_by_target, probe_positions,
            probe_seqs_dict, buffer
        )
        assert "probe_a" not in result

    def test_no_candidate_probes(self):
        """No probes mapped to target should return empty."""
        pair_full_rows = [
            {
                "targets": "['target1']",
                "starts": "[100]",
                "prod_len": 200,
            }
        ]
        probes_by_target = {}
        probe_positions = {}
        probe_seqs_dict = {}
        buffer = 20

        result = _find_position_valid_probes(
            pair_full_rows, probes_by_target, probe_positions,
            probe_seqs_dict, buffer
        )
        assert result == []

    def test_probe_not_in_seqs_dict(self):
        """Probe without sequence should be skipped."""
        pair_full_rows = [
            {
                "targets": "['target1']",
                "starts": "[100]",
                "prod_len": 200,
            }
        ]
        probes_by_target = {"target1": {"probe_a"}}
        probe_positions = {("target1", "probe_a"): 150}
        probe_seqs_dict = {}  # probe_a not in dict
        buffer = 20

        result = _find_position_valid_probes(
            pair_full_rows, probes_by_target, probe_positions,
            probe_seqs_dict, buffer
        )
        assert result == []

    def test_probe_must_be_valid_for_all_targets(self):
        """Probe must pass position check for every target."""
        pair_full_rows = [
            {
                "targets": "['target1', 'target2']",
                "starts": "[100, 200]",
                "prod_len": 300,
            }
        ]
        probes_by_target = {
            "target1": {"probe_a"},
            "target2": {"probe_a"},
        }
        # Valid for target1 but no mapping for target2
        probe_positions = {("target1", "probe_a"): 150}
        probe_seqs_dict = {"probe_a": "ATCG"}
        buffer = 10

        result = _find_position_valid_probes(
            pair_full_rows, probes_by_target, probe_positions,
            probe_seqs_dict, buffer
        )
        # probe_a is missing from target2 positions, so it should fail
        assert "probe_a" not in result
