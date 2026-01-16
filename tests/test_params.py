"""Tests for parameter parsing utilities."""

import pytest
import tempfile
import os
from qprimer_designer.utils.params import parse_params


class TestParseParams:
    """Tests for parameter file parsing."""

    def test_parse_simple_params(self):
        """Test parsing a simple params file."""
        with tempfile.NamedTemporaryFile(mode='w', suffix='.txt', delete=False) as f:
            f.write("PARAM_A = 10\n")
            f.write("PARAM_B = 20\n")
            f.name

        try:
            params = parse_params(f.name)
            assert params["PARAM_A"] == 10
            assert params["PARAM_B"] == 20
        finally:
            os.unlink(f.name)

    def test_parse_params_with_comments(self):
        """Test parsing params file with comments."""
        with tempfile.NamedTemporaryFile(mode='w', suffix='.txt', delete=False) as f:
            f.write("## This is a comment\n")
            f.write("PARAM_A = 10\n")
            f.write("# Another comment\n")
            f.write("PARAM_B = 20\n")
            f.name

        try:
            params = parse_params(f.name)
            assert params["PARAM_A"] == 10
            assert params["PARAM_B"] == 20
            assert len(params) == 2
        finally:
            os.unlink(f.name)

    def test_parse_params_with_empty_lines(self):
        """Test parsing params file with empty lines."""
        with tempfile.NamedTemporaryFile(mode='w', suffix='.txt', delete=False) as f:
            f.write("PARAM_A = 10\n")
            f.write("\n")
            f.write("PARAM_B = 20\n")
            f.name

        try:
            params = parse_params(f.name)
            assert params["PARAM_A"] == 10
            assert params["PARAM_B"] == 20
        finally:
            os.unlink(f.name)
