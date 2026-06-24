"""Per-run working-directory isolation for the Streamlit GUI.

The Streamlit app can serve several users from one process/instance (Cloud Run
runs it with concurrency > 1). The Snakemake pipeline, however, was written to
run from a single working directory: it reads a ``Snakefile`` + ``params.txt``
and takes a ``.snakemake`` lock in the current directory, and most rules
reference inputs/outputs with paths *relative* to the cwd (``target_seqs/...``,
``runs/<RUN_ID>/...``). Two runs sharing one cwd would clobber each other's
control files and lock.

``prepare_run_dir`` gives every run its own scratch cwd containing just the
generated ``Snakefile`` + ``params.txt`` (so the ``.snakemake`` lock is private),
while the shared, read-mostly inputs and the persisted outputs stay where they
are via symlinks:

* ``<scratch>/target_seqs`` -> the shared ``target_seqs`` dir (baked reference
  sequences + user uploads), so the template's relative ``target_seqs/...``
  references resolve;
* ``<scratch>/runs`` -> the shared ``RUNS_DIR``, so outputs land in the shared
  (and, in future, persisted) location.

This module deliberately depends only on the standard library so it can be unit
tested without importing Streamlit / torch.
"""

from __future__ import annotations

import os
import re
import tempfile
from pathlib import Path


def _slug(value: str) -> str:
    """Filesystem-safe fragment of a run id for the scratch dir name."""
    return re.sub(r"[^A-Za-z0-9_.-]", "_", value)[:40] or "run"


def prepare_run_dir(
    run_id: str,
    snakefile_content: str,
    params_content: str,
    target_seqs_dir: os.PathLike | str,
    runs_dir: os.PathLike | str,
) -> Path:
    """Create an isolated scratch working directory for one pipeline run.

    Writes ``Snakefile`` and ``params.txt`` into a fresh temp dir and symlinks
    the shared ``target_seqs`` inputs and ``runs`` outputs into it. Returns the
    scratch dir, which the caller should use as the ``cwd`` for the Snakemake
    subprocess. The shared ``runs_dir`` is created if missing.
    """
    target_seqs_dir = Path(target_seqs_dir).resolve()
    runs_dir = Path(runs_dir)
    runs_dir.mkdir(parents=True, exist_ok=True)

    scratch = Path(tempfile.mkdtemp(prefix=f"qprimer_run_{_slug(run_id)}_"))
    (scratch / "Snakefile").write_text(snakefile_content)
    (scratch / "params.txt").write_text(params_content)

    # Symlink shared inputs/outputs so the Snakefile's relative paths resolve
    # while the .snakemake lock (created in cwd) stays private to this run.
    os.symlink(target_seqs_dir, scratch / "target_seqs", target_is_directory=True)
    os.symlink(runs_dir.resolve(), scratch / "runs", target_is_directory=True)

    return scratch
