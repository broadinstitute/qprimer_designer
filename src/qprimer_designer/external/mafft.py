"""MAFFT wrapper for multiple sequence alignment."""

import shutil
import subprocess
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
        FileNotFoundError: If mafft is not found
        subprocess.CalledProcessError: If alignment fails
    """
    find_mafft()  # Verify it exists

    cmd = ["mafft", "--thread", str(threads)]

    if auto:
        cmd.append("--auto")
    if quiet:
        cmd.append("--quiet")

    cmd.append(str(input_fasta))

    with open(output_fasta, "w") as outfile:
        subprocess.run(cmd, stdout=outfile, check=True)
