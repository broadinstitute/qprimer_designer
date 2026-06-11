"""Tests for command modules (pure logic, no external tools)."""

import os
import tempfile
from unittest.mock import patch, MagicMock

import pandas as pd
import pytest

from qprimer_designer.commands.generate import (
    generate_primers_single,
    count_primer_pairs,
)
from qprimer_designer.commands.prepare_features import infer_primer_orientation
from qprimer_designer.commands.select_multiplex import (
    infer_target_name,
    find_off_score_columns,
    select_top_rows,
    _parse_sam_cross_hits,
)
from qprimer_designer.commands.prepare_input import _ensure_list
from qprimer_designer.commands.export_report import (
    filter_by_primer_id,
    get_target_name,
)


# --- Generate command tests ---


class TestGeneratePrimersSingle:
    """Tests for generate_primers_single."""

    def test_basic_generation(self):
        target = "ATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCG"
        fors, revs = generate_primers_single(target, step=1, min_pri_len=20, max_pri_len=20, min_amp_len=30)
        assert isinstance(fors, dict)
        assert isinstance(revs, dict)
        assert len(fors) > 0
        assert len(revs) > 0

    def test_no_n_primers(self):
        """Primers containing N should be excluded."""
        target = "ATCGATCNNNNCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCG"
        fors, revs = generate_primers_single(target, step=1, min_pri_len=20, max_pri_len=20, min_amp_len=20)
        for seq in fors:
            assert "N" not in seq
        for seq in revs:
            assert "N" not in seq

    def test_primer_length_range(self):
        target = "ATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCG"
        fors, revs = generate_primers_single(target, step=1, min_pri_len=18, max_pri_len=22, min_amp_len=20)
        for seq in fors:
            assert 18 <= len(seq) <= 22
        for seq in revs:
            assert 18 <= len(seq) <= 22

    def test_step_size(self):
        """Larger step size should produce fewer primers."""
        target = "ATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCG"
        fors_step1, _ = generate_primers_single(target, step=1, min_pri_len=20, max_pri_len=20, min_amp_len=20)
        fors_step5, _ = generate_primers_single(target, step=5, min_pri_len=20, max_pri_len=20, min_amp_len=20)
        assert len(fors_step1) >= len(fors_step5)

    def test_short_target(self):
        """Very short target should produce no primers."""
        target = "ATCG"
        fors, revs = generate_primers_single(target, step=1, min_pri_len=20, max_pri_len=20, min_amp_len=20)
        assert len(fors) == 0
        assert len(revs) == 0


class TestCountPrimerPairs:
    """Tests for count_primer_pairs."""

    def test_basic_count(self):
        fors = {"ATCG": 0, "GCTA": 10}
        revs = {"TTTT": 80, "CCCC": 100}
        count = count_primer_pairs(fors, revs, min_amp_len=60, max_amp_len=200)
        assert isinstance(count, int)
        assert count >= 0

    def test_no_pairs(self):
        """No valid amplicons when positions too close."""
        fors = {"ATCG": 0}
        revs = {"GCTA": 10}
        count = count_primer_pairs(fors, revs, min_amp_len=60, max_amp_len=200)
        assert count == 0

    def test_empty(self):
        assert count_primer_pairs({}, {}, min_amp_len=60, max_amp_len=200) == 0

    def test_all_valid(self):
        """All forward/reverse combos should be valid when positions are well-spaced."""
        fors = {"A": 0}
        revs = {"B": 100}
        count = count_primer_pairs(fors, revs, min_amp_len=50, max_amp_len=200)
        assert count == 1


# --- Prepare features tests ---


class TestInferPrimerOrientation:
    """Tests for infer_primer_orientation."""

    def test_forward_f(self):
        assert infer_primer_orientation("primer1_f") == "f"

    def test_forward_for(self):
        assert infer_primer_orientation("primer1_for") == "f"

    def test_reverse_r(self):
        assert infer_primer_orientation("primer1_r") == "r"

    def test_reverse_rev(self):
        assert infer_primer_orientation("primer1_rev") == "r"

    def test_case_insensitive(self):
        assert infer_primer_orientation("PRIMER1_F") == "f"
        assert infer_primer_orientation("PRIMER1_R") == "r"

    def test_unknown_raises(self):
        with pytest.raises(ValueError, match="Cannot infer primer orientation"):
            infer_primer_orientation("primer1")

    def test_unknown_suffix_raises(self):
        with pytest.raises(ValueError, match="Cannot infer primer orientation"):
            infer_primer_orientation("primer1_x")


# --- Pick representatives tests ---



# --- Select multiplex tests ---


class TestInferTargetName:
    """Tests for infer_target_name."""

    def test_csv(self):
        assert infer_target_name("final/HIV.csv") == "HIV"

    def test_tsv(self):
        assert infer_target_name("final/HIV.tsv") == "HIV"

    def test_csv_gz(self):
        assert infer_target_name("final/HIV.csv.gz") == "HIV"

    def test_tsv_gz(self):
        assert infer_target_name("final/HIV.tsv.gz") == "HIV"

    def test_no_extension(self):
        assert infer_target_name("final/HIV") == "HIV"

    def test_case_insensitive_extension(self):
        assert infer_target_name("final/HIV.CSV") == "HIV"

    def test_path_with_dirs(self):
        assert infer_target_name("/home/user/data/results/HIV.csv") == "HIV"


class TestFindOffScoreColumns:
    """Tests for find_off_score_columns."""

    def test_basic(self):
        df = pd.DataFrame(columns=["sco_target", "sco_off1", "sco_off2", "other"])
        result = find_off_score_columns(df)
        assert result == ["sco_off1", "sco_off2"]

    def test_no_off_cols(self):
        df = pd.DataFrame(columns=["sco_target", "other"])
        result = find_off_score_columns(df)
        assert result == []

    def test_no_sco_cols(self):
        df = pd.DataFrame(columns=["col1", "col2"])
        result = find_off_score_columns(df)
        assert result == []


class TestSelectTopRows:
    """Tests for select_top_rows."""

    def test_basic_selection(self):
        df = pd.DataFrame({
            "sco_target": [0.9, 0.8, 0.7, 0.6],
            "other": ["a", "b", "c", "d"],
        })
        result = select_top_rows(df, top_n=2)
        assert len(result) == 2
        assert result["sco_target"].iloc[0] == 0.9

    def test_missing_sco_target_raises(self):
        df = pd.DataFrame({"other": [1, 2, 3]})
        with pytest.raises(ValueError, match="missing required column"):
            select_top_rows(df, top_n=2)

    def test_with_off_target_columns(self):
        df = pd.DataFrame({
            "sco_target": [0.9, 0.9, 0.7],
            "sco_off1": [0.1, 0.5, 0.0],
            "sco_off2": [0.2, 0.3, 0.0],
        })
        result = select_top_rows(df, top_n=2)
        assert len(result) == 2
        # First should have lower off-target worst score
        assert result.iloc[0]["sco_off1"] == 0.1

    def test_top_n_larger_than_df(self):
        df = pd.DataFrame({"sco_target": [0.9, 0.8]})
        result = select_top_rows(df, top_n=10)
        assert len(result) == 2

    def test_all_nan_sco_target_raises(self):
        df = pd.DataFrame({"sco_target": [None, None]})
        with pytest.raises(ValueError, match="No valid rows"):
            select_top_rows(df, top_n=2)


class TestParseSamCrossHits:
    """Tests for _parse_sam_cross_hits."""

    def test_empty_file(self, tmp_path):
        sam = tmp_path / "test.sam"
        sam.write_text("")
        result = _parse_sam_cross_hits(str(sam), "HIV")
        assert result == []

    def test_nonexistent_file(self, tmp_path):
        result = _parse_sam_cross_hits(str(tmp_path / "nonexistent.sam"), "HIV")
        assert result == []

    def test_skips_unmapped(self, tmp_path):
        sam = tmp_path / "test.sam"
        sam.write_text("probe1\t4\t*\t0\t0\t*\t*\t0\t0\tATCG\t*\n")
        result = _parse_sam_cross_hits(str(sam), "HIV")
        assert result == []

    def test_skips_same_target(self, tmp_path):
        sam = tmp_path / "test.sam"
        sam.write_text("probe1\t0\tHIV_1\t100\t255\t20M\t*\t0\t0\tATCG\t*\n")
        result = _parse_sam_cross_hits(str(sam), "HIV")
        assert result == []

    def test_returns_cross_target_hit(self, tmp_path):
        sam = tmp_path / "test.sam"
        sam.write_text("probe1\t0\tFLU_1\t100\t255\t20M\t*\t0\t0\tATCG\t*\n")
        result = _parse_sam_cross_hits(str(sam), "HIV")
        assert len(result) == 1
        assert result[0] == ("probe1", "FLU", "FLU_1")

    def test_skips_header(self, tmp_path):
        sam = tmp_path / "test.sam"
        sam.write_text("@HD\tVN:1.0\n@SQ\tSN:ref\tLN:100\nprobe1\t0\tFLU_1\t100\t255\t20M\t*\t0\t0\tATCG\t*\n")
        result = _parse_sam_cross_hits(str(sam), "HIV")
        assert len(result) == 1


# --- Prepare input tests ---


class TestEnsureList:
    """Tests for _ensure_list."""

    def test_list_passthrough(self):
        assert _ensure_list([1, 2, 3]) == [1, 2, 3]

    def test_nan_returns_empty(self):
        import numpy as np
        assert _ensure_list(np.nan) == []
        assert _ensure_list(float('nan')) == []

    def test_scalar_wrapped(self):
        assert _ensure_list(42) == [42]
        assert _ensure_list("hello") == ["hello"]

    def test_none_is_nan(self):
        """pd.isna(None) is True."""
        assert _ensure_list(None) == []


# --- Export report tests ---


class TestFilterByPrimerId:
    """Tests for filter_by_primer_id."""

    def test_on_target(self):
        df = pd.DataFrame({
            "pname_f": ["pid1_for", "pid2_for", "pid1_for"],
            "pname_r": ["pid1_rev", "pid2_rev", "pid2_rev"],
            "score": [0.9, 0.8, 0.7],
        })
        result = filter_by_primer_id(df, "pid1", eval_type="on")
        assert len(result) == 1
        assert result.iloc[0]["score"] == 0.9

    def test_off_target(self):
        df = pd.DataFrame({
            "pname_f": ["pid1_for", "pid1_rev", "pid2_for"],
            "pname_r": ["pid1_rev", "pid1_for", "pid2_rev"],
            "score": [0.9, 0.8, 0.7],
        })
        result = filter_by_primer_id(df, "pid1", eval_type="off")
        assert len(result) == 2

    def test_invalid_eval_type(self):
        df = pd.DataFrame({"pname_f": ["a"], "pname_r": ["b"]})
        with pytest.raises(ValueError, match="Unknown eval_type"):
            filter_by_primer_id(df, "pid1", eval_type="invalid")

    def test_no_match(self):
        df = pd.DataFrame({
            "pname_f": ["other_for"],
            "pname_r": ["other_rev"],
        })
        result = filter_by_primer_id(df, "pid1", eval_type="on")
        assert len(result) == 0


class TestGetTargetName:
    """Tests for get_target_name."""

    def test_basic(self):
        assert get_target_name("eval.HIV.csv") == "HIV"

    def test_with_path(self):
        assert get_target_name("/data/eval.FLU.csv") == "FLU"
