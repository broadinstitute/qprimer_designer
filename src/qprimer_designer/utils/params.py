"""Parameter file parsing utilities."""

from pathlib import Path
from typing import Any


def parse_params(param_file: str | Path) -> dict[str, Any]:
    """
    Parse params.txt file.

    Numeric values are parsed as float first, then cast where needed.

    Args:
        param_file: Path to the parameters file

    Returns:
        Dictionary of parameter name -> value
    """
    params = {}
    with open(param_file) as f:
        for line in f:
            if "=" in line:
                name, value = line.split("=", 1)
                name = name.strip()
                value = value.strip()
                # Remove inline comments
                if "#" in value:
                    value = value.split("#")[0].strip()
                try:
                    params[name] = float(value)
                except ValueError:
                    params[name] = value
    return params


def parse_primer_len(value: str | float | int) -> list[int]:
    """
    Parse PRIMER_LEN parameter value.

    Supports single value or comma-separated list.
    Examples:
        "20" -> [20]
        "19,20,21" -> [19, 20, 21]
        20.0 -> [20]

    Args:
        value: Parameter value (string, float, or int)

    Returns:
        List of primer lengths
    """
    # Convert to string if numeric
    if isinstance(value, (int, float)):
        return [int(value)]

    # Parse comma-separated values
    value_str = str(value).strip()
    if "," in value_str:
        lengths = []
        for item in value_str.split(","):
            item = item.strip()
            if item:
                try:
                    lengths.append(int(float(item)))
                except ValueError:
                    continue
        # Remove duplicates while preserving order
        seen = set()
        unique_lengths = []
        for length in lengths:
            if length not in seen:
                seen.add(length)
                unique_lengths.append(length)
        return unique_lengths if unique_lengths else [20]

    # Single value
    try:
        return [int(float(value_str))]
    except ValueError:
        return [20]  # Default fallback


def get_primer_params(params: dict) -> dict:
    """Extract primer generation parameters from parsed params dict."""
    return {
        "max_num": int(params.get("MAX_PRIMER_CANDIDATES", 10000)),
        "step": int(params.get("TILING_STEP", 1)),
        "primer_lens": parse_primer_len(params.get("PRIMER_LEN", 20)),
        "min_amp_len": int(params.get("AMPLEN_MIN", 60)),
        "max_amp_len": int(params.get("AMPLEN_MAX", 200)),
        "max_tm": float(params.get("TM_MAX", 60)),
        "min_tm": float(params.get("TM_MIN", 55)),
        "max_gc": float(params.get("GC_MAX", 60)),
        "min_dg": float(params.get("DG_MIN", -8)),
    }


def get_evaluation_params(params: dict) -> dict:
    """Extract evaluation parameters from parsed params dict."""
    return {
        "primer_len": int(params.get("PRIMER_LEN", 20)),
        "min_amp_len": int(params.get("AMPLEN_MIN", 60)),
        "max_amp_len": int(params.get("AMPLEN_MAX", 200)),
        "min_off_len": int(params.get("OFFLEN_MIN", 60)),
        "max_off_len": int(params.get("OFFLEN_MAX", 2000)),
        "num_select": int(params.get("NUM_TOP_SENSITIVITY", 100)),
        "min_dg": float(params.get("DG_MIN", -8)),
    }
