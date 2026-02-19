"""Tests for export_report command logic."""

import pandas as pd
import pytest

from qprimer_designer.commands.export_report import (
    filter_by_primer_id,
    build_alignment_string,
    build_sensitivity_table,
    build_specificity_table,
    get_target_name,
)


class TestBuildAlignmentString:
    """Tests for build_alignment_string."""

    def test_forward_alignment(self):
        df = pd.DataFrame({
            "pseq_f": ["ATCGATCG"],
            "len_f": [8],
            "match_f": ["||||||||"],
            "starts": [100],
            "tseq_f": ["ATCGATCG"],
        })
        result = build_alignment_string(df, "f")
        assert len(result) == 1
        assert "ATCGATCG" in result.iloc[0]

    def test_reverse_alignment(self):
        df = pd.DataFrame({
            "pseq_r": ["GCTAGCTA"],
            "len_r": [8],
            "match_r": ["||||||||"],
            "starts": [100],
            "tseq_r": ["GCTAGCTA"],
            "prod_len": [200],
        })
        result = build_alignment_string(df, "r")
        assert len(result) == 1
        assert "GCTAGCTA" in result.iloc[0]


class TestBuildSensitivityTable:
    """Tests for build_sensitivity_table."""

    def test_with_on_target_data(self):
        df = pd.DataFrame({
            "eval_type": ["on", "on", "on"],
            "classifier": [1, 1, 0],
            "regressor": [0.9, 0.8, 0.1],
            "target": ["HIV", "HIV", "HIV"],
        })
        result = build_sensitivity_table(df)
        assert not result.empty
        assert "Coverage" in result.columns
        assert result["Coverage"].iloc[0] == 2  # two classifier==1

    def test_no_on_target_data(self):
        df = pd.DataFrame({
            "eval_type": ["off"],
            "classifier": [1],
            "regressor": [0.5],
            "target": ["FLU"],
        })
        result = build_sensitivity_table(df)
        assert result.empty


class TestBuildSpecificityTable:
    """Tests for build_specificity_table."""

    def test_with_off_target_data(self):
        df = pd.DataFrame({
            "eval_type": ["off", "off"],
            "classifier": [1, 0],
            "regressor": [0.9, 0.1],
            "target": ["FLU", "FLU"],
        })
        result = build_specificity_table(df)
        assert not result.empty
        assert "Coverage" in result.columns

    def test_no_off_target_data(self):
        df = pd.DataFrame({
            "eval_type": ["on"],
            "classifier": [1],
            "regressor": [0.5],
            "target": ["HIV"],
        })
        result = build_specificity_table(df)
        assert result.empty
