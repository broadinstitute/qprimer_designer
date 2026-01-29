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
        FileNotFoundError: If bowtie2-build is not found or input file doesn't exist
        subprocess.CalledProcessError: If index building fails (includes stderr in message)
    """
    find_bowtie2_build()  # Verify it exists

    # FIX: Add input file validation (was missing)
    fasta_path = Path(fasta_path)
    if not fasta_path.exists():
        raise FileNotFoundError(f"Input FASTA file not found: {fasta_path}")

    # ORIGINAL: subprocess.run without timeout
    # subprocess.run(
    #     ["bowtie2-build", "--threads", str(threads), str(fasta_path), str(index_prefix)],
    #     check=True,
    #     capture_output=True,
    # )

    # TODO: Implement cleanup of partial .bt2 index files on failure.
    # If bowtie2-build fails partway through, partial index files may remain.
    # Could use try/finally with glob to find and remove partial files:
    #   index_files = list(Path(index_prefix).parent.glob(f"{Path(index_prefix).name}*.bt2*"))
    #   for f in index_files: f.unlink()

    # FIX: Add timeout parameter and expose stderr in exception
    cmd = ["bowtie2-build", "--threads", str(threads), str(fasta_path), str(index_prefix)]
    try:
        subprocess.run(
            cmd,
            check=True,
            capture_output=True,
            timeout=300,  # 5 minute timeout
        )
    except subprocess.CalledProcessError as e:
        # FIX: Include stderr in exception for better error messages
        stderr_msg = e.stderr.decode() if e.stderr else ""
        raise subprocess.CalledProcessError(
            e.returncode, cmd, output=e.stdout, stderr=e.stderr
        ) from None  # Re-raise with stderr accessible via e.stderr


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
        FileNotFoundError: If bowtie2 is not found or query file doesn't exist
        subprocess.CalledProcessError: If alignment fails
    """
    find_bowtie2()  # Verify it exists

    # FIX: Add input file validation (was missing)
    query_fasta = Path(query_fasta)
    if not query_fasta.exists():
        raise FileNotFoundError(f"Query FASTA file not found: {query_fasta}")

    # ORIGINAL: Missing -f flag for FASTA input - bowtie2 defaults to FASTQ
    # cmd = [
    #     "bowtie2",
    #     "-x", str(index_prefix),
    #     "-U", str(query_fasta),
    #     "-S", str(output_sam),
    #     "-p", str(threads),
    #     "--no-hd",  # No header
    #     "-a",  # Report all alignments
    # ]

    # FIX: Added "-f" flag to indicate FASTA input format (CRITICAL FIX)
    cmd = [
        "bowtie2",
        "-f",  # FIX: FASTA input format (was missing - caused "not FASTQ" error)
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

    # ORIGINAL: subprocess.run without timeout
    # subprocess.run(cmd, check=True, capture_output=True)

    # FIX: Add timeout parameter and expose stderr in exception
    try:
        subprocess.run(cmd, check=True, capture_output=True, timeout=300)
    except subprocess.CalledProcessError as e:
        # FIX: Include stderr in exception for better error messages
        raise subprocess.CalledProcessError(
            e.returncode, cmd, output=e.stdout, stderr=e.stderr
        ) from None
