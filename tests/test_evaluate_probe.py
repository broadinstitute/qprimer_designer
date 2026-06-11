"""Tests for the evaluate-probe command."""

import argparse
import textwrap
from pathlib import Path

import pandas as pd
import pytest


@pytest.fixture
def eval_probe_data(tmp_path):
    """Create minimal input files for evaluate-probe."""
    # Probe FASTA
    probe_fa = tmp_path / "probes.fa"
    probe_fa.write_text(">probe1\nATCGATCGATCGATCGATCG\n")

    # Reference FASTA — target contains probe exactly at position 10
    ref_fa = tmp_path / "ref.fa"
    ref_fa.write_text(
        ">target1\n"
        "AAAAAAAAAA" "ATCGATCGATCGATCGATCG" "AAAAAAAAAA\n"
    )

    # eval.full CSV
    eval_full = tmp_path / "test.eval.full"
    eval_full.write_text(
        "pname_f,pname_r,targets,starts,prod_len\n"
        "fwd1,rev1,\"['target1']\",\"[0]\",40\n"
    )

    out_csv = tmp_path / "probe.csv"
    return {
        "probe_fa": probe_fa,
        "ref_fa": ref_fa,
        "eval_full": eval_full,
        "out": out_csv,
        "tmp_path": tmp_path,
    }


def _make_args(data, max_mismatches=3.0, max_indels=0):
    return argparse.Namespace(
        probe_fa=str(data["probe_fa"]),
        eval_full=str(data["eval_full"]),
        ref=str(data["ref_fa"]),
        out=str(data["out"]),
        max_mismatches=max_mismatches,
        max_indels=max_indels,
    )


class TestEvaluateProbeRun:
    """Tests for evaluate_probe.run()."""

    def test_exact_match_found(self, eval_probe_data):
        from qprimer_designer.commands.evaluate_probe import run
        args = _make_args(eval_probe_data)
        run(args)

        df = pd.read_csv(eval_probe_data["out"])
        assert len(df) > 0
        fwd_hits = df[df["orientation"] == "+"]
        assert any(fwd_hits["start_pos"] == 10)
        assert any(fwd_hits["mismatches"] == 0.0)

    def test_empty_probe_fa(self, eval_probe_data):
        from qprimer_designer.commands.evaluate_probe import run
        eval_probe_data["probe_fa"].write_text("")
        args = _make_args(eval_probe_data)
        run(args)

        df = pd.read_csv(eval_probe_data["out"])
        assert len(df) == 0

    def test_missing_eval_full(self, eval_probe_data):
        from qprimer_designer.commands.evaluate_probe import run
        eval_probe_data["eval_full"].unlink()
        args = _make_args(eval_probe_data)
        run(args)

        df = pd.read_csv(eval_probe_data["out"])
        assert len(df) == 0

    def test_tseq_and_match_columns_present(self, eval_probe_data):
        from qprimer_designer.commands.evaluate_probe import run
        args = _make_args(eval_probe_data)
        run(args)

        df = pd.read_csv(eval_probe_data["out"])
        assert "tseq" in df.columns
        assert "match" in df.columns

    def test_strict_threshold_filters(self, eval_probe_data):
        """With max_mismatches=0 and a probe that has mismatches, no hits."""
        from qprimer_designer.commands.evaluate_probe import run
        # Write a probe that doesn't exactly match
        eval_probe_data["probe_fa"].write_text(">probe1\nCCCCCCCCCCCCCCCCCCCC\n")
        args = _make_args(eval_probe_data, max_mismatches=0)
        run(args)

        out_text = eval_probe_data["out"].read_text().strip()
        assert out_text == "" or out_text.count("\n") == 0  # empty or header-only

    def test_target_not_in_ref(self, eval_probe_data):
        """eval.full references a target not in ref FASTA → gracefully skip."""
        from qprimer_designer.commands.evaluate_probe import run
        eval_probe_data["eval_full"].write_text(
            "pname_f,pname_r,targets,starts,prod_len\n"
            "fwd1,rev1,\"['nonexistent']\",\"[0]\",40\n"
        )
        args = _make_args(eval_probe_data)
        run(args)

        out_text = eval_probe_data["out"].read_text().strip()
        assert out_text == "" or out_text.count("\n") == 0  # empty or header-only

    def test_reverse_complement_hit(self, eval_probe_data):
        """Probe RC should also produce hits."""
        from qprimer_designer.commands.evaluate_probe import run
        from qprimer_designer.utils.sequences import reverse_complement_dna

        probe_seq = "ATCGATCGATCGATCGATCG"
        probe_rc = reverse_complement_dna(probe_seq)
        # Put the RC in the reference instead
        eval_probe_data["ref_fa"].write_text(
            f">target1\nAAAAAAAAAA{probe_rc}AAAAAAAAAA\n"
        )
        args = _make_args(eval_probe_data)
        run(args)

        df = pd.read_csv(eval_probe_data["out"])
        rc_hits = df[df["orientation"] == "-"]
        assert len(rc_hits) > 0


class TestEvaluateProbeRegister:
    """Tests for CLI registration."""

    def test_register_adds_subcommand(self):
        from qprimer_designer.commands.evaluate_probe import register
        parser = argparse.ArgumentParser()
        subparsers = parser.add_subparsers()
        register(subparsers)
        # Should parse without error
        args = parser.parse_args([
            "evaluate-probe",
            "--probe-fa", "p.fa",
            "--eval-full", "e.full",
            "--ref", "r.fa",
            "--out", "o.csv",
        ])
        assert args.probe_fa == "p.fa"
        assert args.max_mismatches == 3.0
