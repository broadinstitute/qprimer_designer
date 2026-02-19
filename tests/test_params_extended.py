"""Extended tests for parameter parsing utilities."""

import os
import tempfile

import pytest

from qprimer_designer.utils.params import (
    parse_params,
    get_primer_params,
    get_probe_params,
    get_evaluation_params,
)


class TestParseParamsExtended:
    """Extended tests for parse_params."""

    def test_inline_comments(self):
        """Test that inline comments are stripped."""
        with tempfile.NamedTemporaryFile(mode='w', suffix='.txt', delete=False) as f:
            f.write("PARAM_A = 10  # This is a comment\n")
            path = f.name

        try:
            params = parse_params(path)
            assert params["PARAM_A"] == 10.0
        finally:
            os.unlink(path)

    def test_string_values(self):
        """Test parsing of non-numeric string values."""
        with tempfile.NamedTemporaryFile(mode='w', suffix='.txt', delete=False) as f:
            f.write("MODE = sensitive\n")
            path = f.name

        try:
            params = parse_params(path)
            assert params["MODE"] == "sensitive"
        finally:
            os.unlink(path)

    def test_float_values(self):
        """Test parsing of float values."""
        with tempfile.NamedTemporaryFile(mode='w', suffix='.txt', delete=False) as f:
            f.write("DG_MIN = -6.5\n")
            path = f.name

        try:
            params = parse_params(path)
            assert params["DG_MIN"] == -6.5
        finally:
            os.unlink(path)

    def test_comment_line_with_equals(self):
        """Comment lines containing = should not be parsed."""
        with tempfile.NamedTemporaryFile(mode='w', suffix='.txt', delete=False) as f:
            f.write("# MAX_VALUE = 100\n")
            f.write("PARAM_A = 5\n")
            path = f.name

        try:
            params = parse_params(path)
            assert "MAX_VALUE" not in params
            assert params["PARAM_A"] == 5.0
        finally:
            os.unlink(path)

    def test_whitespace_handling(self):
        """Test that extra whitespace is stripped."""
        with tempfile.NamedTemporaryFile(mode='w', suffix='.txt', delete=False) as f:
            f.write("  PARAM_A  =  10  \n")
            path = f.name

        try:
            params = parse_params(path)
            assert params["PARAM_A"] == 10.0
        finally:
            os.unlink(path)

    def test_negative_values(self):
        """Test parsing of negative values."""
        with tempfile.NamedTemporaryFile(mode='w', suffix='.txt', delete=False) as f:
            f.write("DG_MIN = -8\n")
            path = f.name

        try:
            params = parse_params(path)
            assert params["DG_MIN"] == -8.0
        finally:
            os.unlink(path)

    def test_multiple_equals_in_line(self):
        """Test parsing when value contains = sign."""
        with tempfile.NamedTemporaryFile(mode='w', suffix='.txt', delete=False) as f:
            f.write("FORMULA = a=b\n")
            path = f.name

        try:
            params = parse_params(path)
            assert params["FORMULA"] == "a=b"
        finally:
            os.unlink(path)


class TestGetPrimerParams:
    """Tests for get_primer_params."""

    def test_defaults(self):
        """Test default values when no params provided."""
        result = get_primer_params({})
        assert result["max_num"] == 10000
        assert result["step"] == 1
        assert result["min_pri_len"] == 20
        assert result["max_pri_len"] == 20
        assert result["min_amp_len"] == 60
        assert result["max_amp_len"] == 200
        assert result["max_tm"] == 60
        assert result["min_tm"] == 55
        assert result["max_gc"] == 60
        assert result["min_dg"] == -6

    def test_custom_values(self):
        """Test with custom parameter values."""
        params = {
            "MAX_PRIMER_CANDIDATES": 500,
            "TILING_STEP": 3,
            "PRIMER_LEN_MIN": 18,
            "PRIMER_LEN_MAX": 25,
            "AMPLEN_MIN": 100,
            "AMPLEN_MAX": 300,
            "TM_MAX": 65,
            "TM_MIN": 58,
            "GC_MAX": 70,
            "DG_MIN": -8,
        }
        result = get_primer_params(params)
        assert result["max_num"] == 500
        assert result["step"] == 3
        assert result["min_pri_len"] == 18
        assert result["max_pri_len"] == 25
        assert result["min_amp_len"] == 100
        assert result["max_amp_len"] == 300

    def test_types(self):
        """Test that returned types are correct."""
        result = get_primer_params({})
        assert isinstance(result["max_num"], int)
        assert isinstance(result["step"], int)
        assert isinstance(result["max_tm"], float)
        assert isinstance(result["min_dg"], float)


class TestGetProbeParams:
    """Tests for get_probe_params."""

    def test_defaults(self):
        """Test default probe parameters."""
        result = get_probe_params({})
        assert result["len_min"] == 24
        assert result["len_max"] == 28
        assert result["min_tm"] == 65
        assert result["max_tm"] == 70
        assert result["homopolymer_max"] == 3
        assert result["avoid_5prime_g"] is True
        assert result["max_gc"] == 60
        assert result["min_dg"] == -6
        assert result["max_num"] == 10000

    def test_avoid_5prime_g_true_string(self):
        """Test boolean parsing of avoid_5prime_g."""
        result = get_probe_params({"PROBE_AVOID_5PRIME_G": "True"})
        assert result["avoid_5prime_g"] is True

    def test_avoid_5prime_g_false_string(self):
        result = get_probe_params({"PROBE_AVOID_5PRIME_G": "False"})
        assert result["avoid_5prime_g"] is False

    def test_avoid_5prime_g_yes_string(self):
        result = get_probe_params({"PROBE_AVOID_5PRIME_G": "yes"})
        assert result["avoid_5prime_g"] is True

    def test_avoid_5prime_g_no_string(self):
        result = get_probe_params({"PROBE_AVOID_5PRIME_G": "no"})
        assert result["avoid_5prime_g"] is False

    def test_avoid_5prime_g_numeric(self):
        """When value is numeric (from parse_params), bool(1.0) should be True."""
        result = get_probe_params({"PROBE_AVOID_5PRIME_G": 1.0})
        assert result["avoid_5prime_g"] is True


class TestGetEvaluationParams:
    """Tests for get_evaluation_params."""

    def test_defaults(self):
        """Test default evaluation parameters."""
        result = get_evaluation_params({})
        assert result["primer_len"] == 20
        assert result["min_amp_len"] == 60
        assert result["max_amp_len"] == 200
        assert result["min_off_len"] == 60
        assert result["max_off_len"] == 2000
        assert result["num_select"] == 100
        assert result["min_dg"] == -8

    def test_custom_values(self):
        """Test with custom evaluation parameters."""
        params = {
            "PRIMER_LEN_MIN": 18,
            "AMPLEN_MIN": 80,
            "AMPLEN_MAX": 300,
            "OFFLEN_MIN": 100,
            "OFFLEN_MAX": 3000,
            "NUM_TOP_SENSITIVITY": 50,
            "DG_MIN": -10,
        }
        result = get_evaluation_params(params)
        assert result["primer_len"] == 18
        assert result["min_amp_len"] == 80
        assert result["max_amp_len"] == 300
        assert result["min_off_len"] == 100
        assert result["max_off_len"] == 3000
        assert result["num_select"] == 50
        assert result["min_dg"] == -10
