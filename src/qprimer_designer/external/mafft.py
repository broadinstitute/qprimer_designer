"""MAFFT wrapper for multiple sequence alignment."""

import shutil
import subprocess
import tempfile
from pathlib import Path


def find_mafft() -> str:
    """
    Find mafft executable in PATH.

    Returns:
        Path to mafft executable

    Raises:
        FileNotFoundError: If mafft is not found in PATH
    """
    path = shutil.which("mafft")
    if path is None:
        raise FileNotFoundError(
            "mafft not found in PATH. "
            "Install: conda install -c bioconda mafft"
        )
    return path


def align_sequences(
    input_fasta: str | Path,
    output_fasta: str | Path,
    threads: int = 1,
    auto: bool = True,
    quiet: bool = True,
) -> None:
    """
    Perform multiple sequence alignment using MAFFT.

    Args:
        input_fasta: Path to input FASTA file
        output_fasta: Path to output aligned FASTA file
        threads: Number of threads to use
        auto: Use automatic algorithm selection (default: True)
        quiet: Suppress MAFFT output (default: True)

    Raises:
        FileNotFoundError: If mafft is not found or input file doesn't exist
        ValueError: If input file is empty or invalid
        subprocess.CalledProcessError: If alignment fails
    """
    find_mafft()  # Verify it exists

    # FIX: Add input file validation (was missing)
    input_fasta = Path(input_fasta)
    if not input_fasta.exists():
        raise FileNotFoundError(f"Input FASTA file not found: {input_fasta}")

    # FIX: Validate input file is not empty
    if input_fasta.stat().st_size == 0:
        raise ValueError(f"Input FASTA file is empty: {input_fasta}")

    # FIX: Basic FASTA format validation
    content = input_fasta.read_text()
    if not content.strip().startswith(">"):
        raise ValueError(f"Input file does not appear to be in FASTA format: {input_fasta}")

    cmd = ["mafft", "--thread", str(threads)]

    if auto:
        cmd.append("--auto")
    if quiet:
        cmd.append("--quiet")

    cmd.append(str(input_fasta))

    # ORIGINAL: Output file opened BEFORE subprocess runs
    # This leaves an empty file if mafft fails, breaking downstream processing
    # with open(output_fasta, "w") as outfile:
    #     subprocess.run(cmd, stdout=outfile, check=True)

    # FIX: Use temporary file + atomic move pattern
    # Only create final output file after successful completion
    output_fasta = Path(output_fasta)
    with tempfile.NamedTemporaryFile(
        mode="w",
        suffix=".fasta",
        dir=output_fasta.parent,
        delete=False
    ) as tmp_file:
        tmp_path = Path(tmp_file.name)
        try:
            # FIX: Add stderr capture and timeout
            result = subprocess.run(
                cmd,
                stdout=tmp_file,
                stderr=subprocess.PIPE,  # FIX: Capture stderr for error messages
                check=True,
                timeout=600,  # FIX: 10 minute timeout
            )
            tmp_file.flush()

            # FIX: Atomic move only after successful completion
            shutil.move(str(tmp_path), str(output_fasta))

        except (subprocess.CalledProcessError, subprocess.TimeoutExpired) as e:
            # FIX: Clean up temp file on failure
            if tmp_path.exists():
                tmp_path.unlink()
            # Include stderr in error message if available
            if hasattr(e, "stderr") and e.stderr:
                raise subprocess.CalledProcessError(
                    e.returncode if hasattr(e, "returncode") else 1,
                    cmd,
                    stderr=e.stderr
                ) from e
            raise
