"""ViennaRNA RNAduplex wrapper for computing duplex free energy."""

import functools
import re
import shutil
import subprocess
import warnings


# Valid DNA characters (standard bases only; RNAduplex does not support IUPAC ambiguity codes)
VALID_DNA_CHARS = set("ATCGatcgNn-")

# IUPAC ambiguity codes that can be safely replaced with N
_IUPAC_AMBIGUITY_CHARS = set("RYWSMKBDHVryswmkbdhv")


@functools.lru_cache(maxsize=1)
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
        # IUPAC ambiguity codes are replaced with N rather than raising
        iupac_chars = invalid_chars & _IUPAC_AMBIGUITY_CHARS
        truly_invalid = invalid_chars - _IUPAC_AMBIGUITY_CHARS

        if truly_invalid:
            raise ValueError(
                f"{name} contains invalid characters: {sorted(truly_invalid)}. "
                f"Valid characters are: ATCGN (case-insensitive) and -"
            )

        if iupac_chars:
            warnings.warn(
                f"{name} contains IUPAC ambiguity codes {sorted(iupac_chars)}, "
                f"replacing with N",
                stacklevel=3,
            )
            for ch in iupac_chars:
                seq = seq.replace(ch, 'N')

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


def compute_batch_dimer_dg(pairs: list) -> list:
    """
    Compute dimer ΔG for multiple sequence pairs in a single RNAduplex process.

    Args:
        pairs: List of (seq1, seq2) tuples

    Returns:
        List of ΔG values (float) in the same order as input pairs

    Raises:
        FileNotFoundError: If RNAduplex is not found
        ValueError: If output count doesn't match input count
    """
    if not pairs:
        return []

    find_rnaduplex()

    validated = []
    for i, (s1, s2) in enumerate(pairs):
        validated.append((
            _validate_sequence(s1, f"pair[{i}].seq1"),
            _validate_sequence(s2, f"pair[{i}].seq2"),
        ))

    inp_lines = []
    for s1, s2 in validated:
        inp_lines.append(s1)
        inp_lines.append(s2)
    inp = "\n".join(inp_lines) + "\n"

    result = subprocess.run(
        ["RNAduplex", "--noconv", "--paramFile=DNA"],
        input=inp,
        capture_output=True,
        text=True,
        check=True,
    )

    dg_pattern = re.compile(r"\(\s*([-+]?\d*\.?\d+)\s*\)")
    dg_values = []
    for line in result.stdout.split("\n"):
        m = dg_pattern.search(line)
        if m:
            dg_values.append(float(m.group(1)))

    if len(dg_values) != len(pairs):
        raise ValueError(
            f"Expected {len(pairs)} ΔG values from RNAduplex, got {len(dg_values)}"
        )

    return dg_values


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
