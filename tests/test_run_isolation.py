"""Tests for gui.run_isolation.prepare_run_dir (stdlib-only; no Streamlit)."""

import sys
from pathlib import Path

# gui/ lives at the repo root (not under src/), so make it importable.
sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from gui.run_isolation import prepare_run_dir  # noqa: E402


def _make_target_seqs(root: Path) -> Path:
    """Create a minimal shared target_seqs/{original,msa} tree."""
    target_seqs = root / "target_seqs"
    (target_seqs / "original").mkdir(parents=True)
    (target_seqs / "msa").mkdir(parents=True)
    (target_seqs / "original" / "virusA.fa").write_text(">virusA\nACGT\n")
    return target_seqs


def test_writes_control_files_and_symlinks(tmp_path):
    target_seqs = _make_target_seqs(tmp_path)
    runs_dir = tmp_path / "runs"

    scratch = prepare_run_dir(
        run_id="20260616_120000_abcd",
        snakefile_content="SNAKEFILE-BODY",
        params_content="PARAMS-BODY",
        target_seqs_dir=target_seqs,
        runs_dir=runs_dir,
    )

    assert scratch.is_dir()
    assert (scratch / "Snakefile").read_text() == "SNAKEFILE-BODY"
    assert (scratch / "params.txt").read_text() == "PARAMS-BODY"

    # Symlinks make the Snakefile's relative paths resolve.
    assert (scratch / "target_seqs").is_symlink()
    assert (scratch / "runs").is_symlink()
    assert (scratch / "target_seqs" / "original" / "virusA.fa").read_text().startswith(">virusA")
    assert (scratch / "target_seqs").resolve() == target_seqs.resolve()
    assert (scratch / "runs").resolve() == runs_dir.resolve()


def test_runs_dir_created_when_missing(tmp_path):
    target_seqs = _make_target_seqs(tmp_path)
    runs_dir = tmp_path / "nested" / "runs"
    assert not runs_dir.exists()

    prepare_run_dir(
        run_id="r1",
        snakefile_content="s",
        params_content="p",
        target_seqs_dir=target_seqs,
        runs_dir=runs_dir,
    )
    assert runs_dir.is_dir()


def test_each_run_gets_a_private_dir_but_shared_output(tmp_path):
    target_seqs = _make_target_seqs(tmp_path)
    runs_dir = tmp_path / "runs"

    s1 = prepare_run_dir("run-1", "s", "p", target_seqs, runs_dir)
    s2 = prepare_run_dir("run-2", "s", "p", target_seqs, runs_dir)

    # Distinct scratch dirs => isolated Snakefile / params.txt / .snakemake lock.
    assert s1 != s2

    # But outputs written via the per-run `runs` symlink land in the shared dir.
    (s1 / "runs" / "run-1").mkdir()
    (s2 / "runs" / "run-2").mkdir()
    assert (runs_dir / "run-1").is_dir()
    assert (runs_dir / "run-2").is_dir()


def test_run_id_with_unsafe_chars_is_sanitized(tmp_path):
    target_seqs = _make_target_seqs(tmp_path)
    runs_dir = tmp_path / "runs"
    # Slashes / spaces must not escape the temp dir or break creation.
    scratch = prepare_run_dir("a/b c..d", "s", "p", target_seqs, runs_dir)
    assert scratch.is_dir()
    assert "/" not in scratch.name
