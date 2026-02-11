"""ViennaRNA RNAduplex wrapper for computing duplex free energy."""

import re
import shutil
import subprocess


# Valid DNA characters (standard bases only; RNAduplex does not support IUPAC ambiguity codes)
VALID_DNA_CHARS = set("ATCGatcgNn-")


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


def _validate_sequence(seq: str, name: str = "sequence") -> str:
    """
    Validate and clean a DNA sequence.

    Args:
        seq: DNA sequence to validate
        name: Name for error messages

    Returns:
        Cleaned sequence (whitespace stripped)

    Raises:
        TypeError: If sequence is not a string
        ValueError: If sequence is empty or contains invalid characters
    """
    if seq is None:
        raise TypeError(f"{name} cannot be None")
    if not isinstance(seq, str):
        raise TypeError(f"{name} must be a string, got {type(seq).__name__}")

    seq = seq.strip()
    if not seq:
        raise ValueError(f"{name} cannot be empty")

    if any(c.isspace() for c in seq):
        raise ValueError(f"{name} contains internal whitespace")

    invalid_chars = set(seq) - VALID_DNA_CHARS
    if invalid_chars:
        raise ValueError(
            f"{name} contains invalid characters: {sorted(invalid_chars)}. "
            f"Valid characters are: ATCGN (case-insensitive) and -"
        )

    return seq


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
        TypeError: If sequences are not strings
        ValueError: If sequences are empty or contain invalid characters
        ValueError: If ΔG cannot be parsed from output
        subprocess.CalledProcessError: If RNAduplex fails
    """
    find_rnaduplex()  # Verify it exists

    seq1 = _validate_sequence(seq1, "seq1")
    seq2 = _validate_sequence(seq2, "seq2")

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

    Raises:
        TypeError: If sequence is not a string
        ValueError: If sequence is empty or contains invalid characters
    """
    return compute_dimer_dg(seq, seq)
