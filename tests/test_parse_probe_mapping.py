"""Tests for parse_probe_mapping command."""

import os
from unittest.mock import MagicMock

import pandas as pd
import pytest

from qprimer_designer.commands.parse_probe_mapping import run


class TestParseProbeMapping:
    """Tests for the SAM parsing logic."""

    def test_basic_parsing(self, tmp_path):
        """Test parsing a simple SAM file."""
        sam = tmp_path / "test.sam"
        sam.write_text(
            "probe1\t0\ttarget1\t101\t255\t24M\t*\t0\t0\tATCGATCGATCGATCGATCGATCG\t*\n"
            "probe2\t16\ttarget2\t201\t255\t24M\t*\t0\t0\tGCTAGCTAGCTAGCTAGCTAGCTA\t*\n"
        )
        out = tmp_path / "out.csv"

        args = MagicMock()
        args.sam = str(sam)
        args.out = str(out)

        run(args)

        df = pd.read_csv(out)
        assert len(df) == 2
        assert df.iloc[0]["probe_name"] == "probe1"
        assert df.iloc[0]["start_pos"] == 100  # 0-based
        assert df.iloc[0]["orientation"] == "+"
        assert df.iloc[1]["probe_name"] == "probe2"
        assert df.iloc[1]["orientation"] == "-"

    def test_skips_header(self, tmp_path):
        """Test that SAM header lines are skipped."""
        sam = tmp_path / "test.sam"
        sam.write_text(
            "@HD\tVN:1.0\n"
            "@SQ\tSN:ref\tLN:1000\n"
            "probe1\t0\ttarget1\t101\t255\t24M\t*\t0\t0\tATCGATCGATCGATCGATCGATCG\t*\n"
        )
        out = tmp_path / "out.csv"

        args = MagicMock()
        args.sam = str(sam)
        args.out = str(out)

        run(args)

        df = pd.read_csv(out)
        assert len(df) == 1

    def test_empty_sam(self, tmp_path):
        """Test parsing an empty SAM file produces an empty output."""
        sam = tmp_path / "test.sam"
        sam.write_text("")
        out = tmp_path / "out.csv"

        args = MagicMock()
        args.sam = str(sam)
        args.out = str(out)

        run(args)

        # Empty SAM produces an empty/near-empty CSV
        assert out.exists()
        assert out.stat().st_size <= 10  # Just a newline or empty
