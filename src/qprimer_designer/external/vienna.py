"""ViennaRNA RNAduplex wrapper for computing duplex free energy."""

import re
import shutil
import subprocess


def find_rnaduplex() -> str:
    """
    Find RNAduplex executable in PATH.

    Returns:
        Path to RNAduplex executable

    Raises:
        FileNotFoundError: If RNAduplex is not found in PATH
    """
    path = shutil.which("RNAduplex")
    if path is None:
        raise FileNotFoundError(
            "RNAduplex not found in PATH. "
            "Install ViennaRNA: conda install -c bioconda viennarna"
        )
    return path


def compute_dimer_dg(seq1: str, seq2: str) -> float:
    """
    Compute primer-primer dimer free energy (ΔG) using RNAduplex.

    Uses DNA parameter mode for accurate DNA duplex calculations.

    Args:
        seq1: First DNA sequence
        seq2: Second DNA sequence

    Returns:
        Free energy (ΔG) in kcal/mol

    Raises:
        FileNotFoundError: If RNAduplex is not found
        ValueError: If ΔG cannot be parsed from output
        subprocess.CalledProcessError: If RNAduplex fails
    """
    find_rnaduplex()  # Verify it exists

    inp = f"{seq1}\n{seq2}\n"
    result = subprocess.run(
        ["RNAduplex", "--noconv", "--paramFile=DNA"],
        input=inp,
        capture_output=True,
        text=True,
        check=True,
    )

    match = re.search(r"\(\s*([-+]?\d*\.?\d+)\s*\)", result.stdout)
    if not match:
        raise ValueError(
            f"Could not parse ΔG from RNAduplex output:\n{result.stdout}"
        )
    return float(match.group(1))


def compute_self_dimer_dg(seq: str) -> float:
    """
    Compute self-dimer free energy for a single sequence.

    Args:
        seq: DNA sequence

    Returns:
        Self-dimer ΔG in kcal/mol
    """
    return compute_dimer_dg(seq, seq)
