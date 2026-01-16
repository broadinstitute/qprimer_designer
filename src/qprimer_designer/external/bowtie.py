"""Bowtie2 wrapper for sequence alignment."""

import shutil
import subprocess
from pathlib import Path


def find_bowtie2() -> str:
    """
    Find bowtie2 executable in PATH.

    Returns:
        Path to bowtie2 executable

    Raises:
        FileNotFoundError: If bowtie2 is not found in PATH
    """
    path = shutil.which("bowtie2")
    if path is None:
        raise FileNotFoundError(
            "bowtie2 not found in PATH. "
            "Install: conda install -c bioconda bowtie2"
        )
    return path


def find_bowtie2_build() -> str:
    """
    Find bowtie2-build executable in PATH.

    Returns:
        Path to bowtie2-build executable

    Raises:
        FileNotFoundError: If bowtie2-build is not found in PATH
    """
    path = shutil.which("bowtie2-build")
    if path is None:
        raise FileNotFoundError(
            "bowtie2-build not found in PATH. "
            "Install: conda install -c bioconda bowtie2"
        )
    return path


def build_index(fasta_path: str | Path, index_prefix: str | Path, threads: int = 1) -> None:
    """
    Build a bowtie2 index from a FASTA file.

    Args:
        fasta_path: Path to input FASTA file
        index_prefix: Prefix for output index files
        threads: Number of threads to use

    Raises:
        FileNotFoundError: If bowtie2-build is not found
        subprocess.CalledProcessError: If index building fails
    """
    find_bowtie2_build()  # Verify it exists

    subprocess.run(
        ["bowtie2-build", "--threads", str(threads), str(fasta_path), str(index_prefix)],
        check=True,
        capture_output=True,
    )


def align_primers(
    index_prefix: str | Path,
    query_fasta: str | Path,
    output_sam: str | Path,
    threads: int = 1,
    local: bool = True,
    very_sensitive: bool = True,
) -> None:
    """
    Align primer sequences to a reference using bowtie2.

    Args:
        index_prefix: Bowtie2 index prefix
        query_fasta: FASTA file with primer sequences
        output_sam: Output SAM file path
        threads: Number of threads to use
        local: Use local alignment mode (default: True)
        very_sensitive: Use very-sensitive preset (default: True)

    Raises:
        FileNotFoundError: If bowtie2 is not found
        subprocess.CalledProcessError: If alignment fails
    """
    find_bowtie2()  # Verify it exists

    cmd = [
        "bowtie2",
        "-x", str(index_prefix),
        "-U", str(query_fasta),
        "-S", str(output_sam),
        "-p", str(threads),
        "--no-hd",  # No header
        "-a",  # Report all alignments
    ]

    if local:
        cmd.append("--local")
    if very_sensitive:
        cmd.append("--very-sensitive-local" if local else "--very-sensitive")

    subprocess.run(cmd, check=True, capture_output=True)
