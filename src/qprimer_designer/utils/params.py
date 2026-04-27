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
            # BUG FIX: Skip comment lines before checking for '='
            # Original code only checked `if "=" in line:` which incorrectly
            # parsed comment lines like '# MAX_VALUE = 100' as parameters.
            # Original:
            # if "=" in line:
            if "=" in line and not line.lstrip().startswith("#"):
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


def parse_list_param(params: dict, key: str) -> list[str]:
    """Parse a comma-separated parameter value into a list of strings.

    Returns an empty list if the key is missing or the value is empty.
    """
    val = params.get(key, "")
    if isinstance(val, (int, float)):
        return [str(int(val))]
    if not val or not val.strip():
        return []
    return [item.strip() for item in val.split(",") if item.strip()]


def get_primer_params(params: dict) -> dict:
    """Extract primer generation parameters from parsed params dict."""
    return {
        "max_num": int(params.get("MAX_PRIMER_CANDIDATES", 10000)),
        "step": int(params.get("TILING_STEP", 1)),
        "min_pri_len": int(params.get("PRIMER_LEN_MIN", 20)),
        "max_pri_len": int(params.get("PRIMER_LEN_MAX", 20)),
        "min_amp_len": int(params.get("AMPLEN_MIN", 60)),
        "max_amp_len": int(params.get("AMPLEN_MAX", 200)),
        "max_tm": float(params.get("TM_MAX", 60)),
        "min_tm": float(params.get("TM_MIN", 55)),
        "max_gc": float(params.get("GC_MAX", 60)),
        "min_dg": float(params.get("DG_MIN", -6)),
    }


def get_probe_params(params: dict) -> dict:
    """Extract probe generation parameters from parsed params dict."""
    # Parse avoid_5prime_G as boolean
    avoid_5prime_g_str = params.get("PROBE_AVOID_5PRIME_G", "True")
    if isinstance(avoid_5prime_g_str, str):
        avoid_5prime_g = avoid_5prime_g_str.lower() in ("true", "1", "yes")
    else:
        avoid_5prime_g = bool(avoid_5prime_g_str)

    return {
        "len_min": int(params.get("PROBE_LEN_MIN", 24)),
        "len_max": int(params.get("PROBE_LEN_MAX", 28)),
        "min_tm": float(params.get("PROBE_TM_MIN", 65)),
        "max_tm": float(params.get("PROBE_TM_MAX", 70)),
        "homopolymer_max": int(params.get("PROBE_HOMOPOLYMER_MAX", 3)),
        "avoid_5prime_g": avoid_5prime_g,
        "max_gc": float(params.get("GC_MAX", 60)),
        "min_dg": float(params.get("DG_MIN", -6)),
        "max_num": int(params.get("MAX_PRIMER_CANDIDATES", 10000)),
    }


def get_evaluation_params(params: dict) -> dict:
    """Extract evaluation parameters from parsed params dict."""
    return {
        "primer_len": int(params.get("PRIMER_LEN_MIN", 20)),  # Use min as default single value
        "min_amp_len": int(params.get("AMPLEN_MIN", 60)),
        "max_amp_len": int(params.get("AMPLEN_MAX", 200)),
        "min_off_len": int(params.get("OFFLEN_MIN", 60)),
        "max_off_len": int(params.get("OFFLEN_MAX", 2000)),
        "num_select": int(params.get("NUM_TOP_SENSITIVITY", 100)),
        "min_dg": float(params.get("DG_MIN", -8)),
    }
