"""Generate Snakefile and params.txt from GUI configuration."""

import re
from pathlib import Path


TEMPLATE_PATH = Path(__file__).resolve().parent.parent / "workflows" / "Snakefile.template"


def build_snakefile(
    targets: list[str],
    cross: list[str],
    host: list[str],
    panel: list[str],
    target_dir: str = "target_seqs/original",
    params_file: str = "params.txt",
) -> str:
    """Read Snakefile.template and substitute header variable assignments.

    Only the 6 variable assignments in the header block (TARGETS, CROSS,
    HOST, PANEL, TARGET_DIR, PARAMS) are replaced. Everything else passes
    through unchanged.
    """
    template = TEMPLATE_PATH.read_text()

    def _list_literal(items: list[str]) -> str:
        return "[" + ", ".join(f"'{i}'" for i in items) + "]"

    replacements = {
        r"^TARGETS\s*=\s*\[.*?\]": f"TARGETS = {_list_literal(targets)}",
        r"^CROSS\s*=\s*\[.*?\]": f"CROSS   = {_list_literal(cross)}",
        r"^HOST\s*=\s*\[.*?\]": f"HOST    = {_list_literal(host)}",
        r"^PANEL\s*=\s*\[.*?\]": f"PANEL = {_list_literal(panel)}",
        r'^TARGET_DIR\s*=\s*".*?"': f'TARGET_DIR = "{target_dir}"',
        r'^PARAMS\s*=\s*".*?"': f'PARAMS = "{params_file}"',
    }

    result = template
    for pattern, replacement in replacements.items():
        result = re.sub(pattern, replacement, result, count=1, flags=re.MULTILINE)

    return result


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
        # Only include sections that have at least one matching key
        section_lines: list[str] = []
        for key in keys:
            if key in params:
                val = params[key]
                # Format: booleans as True/False, ints without decimal,
                # floats with decimals
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
