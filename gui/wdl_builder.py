"""Generate WDL inputs JSON and params.txt from GUI configuration."""

import json
from datetime import datetime
from pathlib import Path


def build_inputs_json(
    targets: list[str],
    cross: list[str],
    host: list[str],
    panel: list[str],
    target_dir: str = "target_seqs/original",
    params_file: str = "params.txt",
    multiplex: bool = False,
    evaluate: bool = False,
    probe_mode: bool = False,
    forward_seq: str | None = None,
    reverse_seq: str | None = None,
    pset_path: str | None = None,
    threads: int = 2,
    probe_mismatch_penalty: int = 10,
) -> str:
    """Build a WDL inputs JSON string from pipeline configuration.

    Returns the JSON string ready to write to a file.
    """
    inputs: dict = {
        "qprimer_pipeline.multiplex": multiplex,
        "qprimer_pipeline.evaluate": evaluate,
        "qprimer_pipeline.probe_mode": probe_mode,
        "qprimer_pipeline.targets": targets,
        "qprimer_pipeline.cross": cross,
        "qprimer_pipeline.host": host,
        "qprimer_pipeline.panel": panel,
        "qprimer_pipeline.target_dir": target_dir,
        "qprimer_pipeline.params_file": params_file,
        "qprimer_pipeline.threads": threads,
        "qprimer_pipeline.probe_mismatch_penalty": probe_mismatch_penalty,
    }

    if forward_seq:
        inputs["qprimer_pipeline.forward_seq"] = forward_seq
    if reverse_seq:
        inputs["qprimer_pipeline.reverse_seq"] = reverse_seq
    if pset_path:
        inputs["qprimer_pipeline.pset_file"] = pset_path

    return json.dumps(inputs, indent=4) + "\n"


def build_params_txt(params: dict) -> str:
    """Generate params.txt content from a dict.

    Parameters are grouped under section comment headers matching the
    template format.
    """
    sections = [
        (
            "## Parameters in picking representative sequences",
            ["DESIGN_WINDOW"],
        ),
        (
            "## Parameters in generating primers",
            [
                "MAX_PRIMER_CANDIDATES",
                "TM_MIN",
                "TM_MAX",
                "GC_MAX",
                "DG_MIN",
                "PRIMER_LEN_MIN",
                "PRIMER_LEN_MAX",
                "TILING_STEP",
            ],
        ),
        (
            "## Parameters for probe generation",
            [
                "PROBE_LEN_MIN",
                "PROBE_LEN_MAX",
                "PROBE_TM_MIN",
                "PROBE_TM_MAX",
                "PROBE_HOMOPOLYMER_MAX",
                "PROBE_AVOID_5PRIME_G",
            ],
        ),
        (
            "## Probe mode configuration (singleplex only)",
            [
                "PROBE_MAX_MISMATCHES",
                "PROBE_AMPLICON_BUFFER",
                "MIN_PROBES_PER_PAIR",
                "MAX_PROBES_PER_PAIR",
            ],
        ),
        (
            "## Parameters in input preparation",
            [
                "AMPLEN_MIN",
                "AMPLEN_MAX",
                "OFFLEN_MIN",
                "OFFLEN_MAX",
            ],
        ),
        (
            "## Parameters in evaluation",
            ["NUM_TOP_SENSITIVITY"],
        ),
    ]

    lines: list[str] = []
    for header, keys in sections:
        section_lines: list[str] = []
        for key in keys:
            if key in params:
                val = params[key]
                if isinstance(val, bool):
                    section_lines.append(f"{key} = {val}")
                elif isinstance(val, float) and val == int(val):
                    section_lines.append(f"{key} = {int(val)}")
                else:
                    section_lines.append(f"{key} = {val}")
        if section_lines:
            lines.append(header)
            lines.extend(section_lines)
            lines.append("")

    return "\n".join(lines) + "\n" if lines else ""
