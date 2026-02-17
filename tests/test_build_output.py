"""Tests for build_output command logic."""

import pandas as pd
import pytest

from qprimer_designer.commands.build_output import find_valid_probes_offtarget


class TestFindValidProbesOfftarget:
    """Tests for find_valid_probes_offtarget."""

    def test_no_off_target_data(self):
        """All probes valid when off-target data is empty."""
        pair_key = ("primer_f", "primer_r")
        probes = ["probe1", "probe2"]
        offtarget_full = [pd.DataFrame()]
        offtarget_mapping = [pd.DataFrame(columns=["probe_name", "probe_seq", "target_id", "start_pos", "orientation"])]
        probe_seqs = {"probe1": "ATCGATCG", "probe2": "GCTAGCTA"}
        buffer = 10

        result = find_valid_probes_offtarget(
            pair_key, probes, offtarget_full, offtarget_mapping, probe_seqs, buffer
        )
        assert sorted(result) == ["probe1", "probe2"]

    def test_probe_not_in_seqs(self):
        """Probe without sequence should be skipped."""
        pair_key = ("primer_f", "primer_r")
        probes = ["probe_missing"]
        offtarget_full = [pd.DataFrame()]
        offtarget_mapping = [pd.DataFrame(columns=["probe_name", "probe_seq", "target_id", "start_pos", "orientation"])]
        probe_seqs = {}
        buffer = 10

        result = find_valid_probes_offtarget(
            pair_key, probes, offtarget_full, offtarget_mapping, probe_seqs, buffer
        )
        assert result == []

    def test_probe_in_off_target_amplicon_excluded(self):
        """Probe within off-target amplicon should be excluded."""
        pair_key = ("primer_f", "primer_r")
        probes = ["probe1"]
        probe_seqs = {"probe1": "ATCGATCG"}  # len 8

        offtarget_full = [pd.DataFrame({
            "pname_f": ["primer_f"],
            "pname_r": ["primer_r"],
            "targets": ["['off_target1']"],
            "starts": ["[100]"],
            "prod_len": [300],
        })]
        offtarget_mapping = [pd.DataFrame({
            "probe_name": ["probe1"],
            "probe_seq": ["ATCGATCG"],
            "target_id": ["off_target1"],
            "start_pos": [150],
            "orientation": ["+"],
        })]
        buffer = 10

        result = find_valid_probes_offtarget(
            pair_key, probes, offtarget_full, offtarget_mapping, probe_seqs, buffer
        )
        assert "probe1" not in result

    def test_probe_not_mapped_to_off_target(self):
        """Probe not mapped to any off-target should be valid."""
        pair_key = ("primer_f", "primer_r")
        probes = ["probe1"]
        probe_seqs = {"probe1": "ATCGATCG"}

        offtarget_full = [pd.DataFrame({
            "pname_f": ["primer_f"],
            "pname_r": ["primer_r"],
            "targets": ["['off_target1']"],
            "starts": ["[100]"],
            "prod_len": [300],
        })]
        # Probe not in the mapping
        offtarget_mapping = [pd.DataFrame({
            "probe_name": ["other_probe"],
            "probe_seq": ["GCTAGCTA"],
            "target_id": ["off_target1"],
            "start_pos": [150],
            "orientation": ["+"],
        })]
        buffer = 10

        result = find_valid_probes_offtarget(
            pair_key, probes, offtarget_full, offtarget_mapping, probe_seqs, buffer
        )
        assert "probe1" in result

    def test_empty_probes_list(self):
        """Empty input probes should return empty list."""
        result = find_valid_probes_offtarget(
            ("f", "r"), [], [pd.DataFrame()],
            [pd.DataFrame(columns=["probe_name", "probe_seq", "target_id", "start_pos", "orientation"])],
            {}, 10
        )
        assert result == []
