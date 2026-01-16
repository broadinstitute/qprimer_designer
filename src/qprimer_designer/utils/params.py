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
                try:
                    params[name] = float(value)
                except ValueError:
                    params[name] = value
    return params


def get_primer_params(params: dict) -> dict:
    """Extract primer generation parameters from parsed params dict."""
    return {
        "max_num": int(params.get("MAX_PRIMER_CANDIDATES", 10000)),
        "step": int(params.get("TILING_STEP", 1)),
        "primer_len": int(params.get("PRIMER_LEN", 20)),
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
