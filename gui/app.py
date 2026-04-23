"""Streamlit GUI for the qPrimer Designer pipeline."""

import json
import os
import re
import select
import shutil
import signal
import subprocess
import sys
import time
from datetime import date, datetime
from pathlib import Path

import streamlit as st

# ---------------------------------------------------------------------------
# Resolve project root (parent of gui/)
# ---------------------------------------------------------------------------
PROJECT_ROOT = Path(__file__).resolve().parent.parent
FASTA_DIR = PROJECT_ROOT / "target_seqs" / "original"
FINAL_DIR = PROJECT_ROOT / "final"
EVALUATE_DIR = PROJECT_ROOT / "evaluate"
SCHEMATIC_PATH = PROJECT_ROOT / "schematic.png"

# Ensure the upload directory exists
FASTA_DIR.mkdir(parents=True, exist_ok=True)

# Add src/ to path so we can import qprimer_designer utilities
sys.path.insert(0, str(PROJECT_ROOT / "src"))
sys.path.insert(0, str(PROJECT_ROOT))

from gui.snakefile_builder import build_params_txt, build_snakefile
from qprimer_designer.utils.params import parse_params

# ---------------------------------------------------------------------------
# Page config
# ---------------------------------------------------------------------------
st.set_page_config(
    page_title="qPrimer Designer",
    layout="wide",
)

# ---------------------------------------------------------------------------
# Navigation helpers
# ---------------------------------------------------------------------------

def _navigate(page: str):
    """Set the current page in session state."""
    st.session_state.page = page


def _clear_pipeline_state():
    """Clear all pipeline run state for a fresh workflow."""
    for key in ("pipeline_return_code", "pipeline_log", "pipeline_running",
                "pipeline_completed_rules", "pipeline_should_start",
                "_run_page_fresh"):
        st.session_state.pop(key, None)


def _current_page() -> str:
    return st.session_state.get("page", "home")


# Workflow step definitions: (page_key, label)
_WORKFLOW_STEPS = [
    ("select_target", "Target"),
    ("select_offtarget", "Off-Target"),
    ("config", "Configuration"),
    ("run", "Run"),
    ("results", "Results"),
]


def _render_workflow_progress(current_step_key: str):
    """Render a horizontal progress stepper at the top of workflow pages."""
    parts = []
    found = False
    for key, label in _WORKFLOW_STEPS:
        if key == current_step_key:
            parts.append(f"<b>{label}</b>")
            found = True
        elif not found:
            parts.append(label)
        else:
            parts.append(f"<span style='color: #ccc;'>{label}</span>")

    st.markdown(
        "<div style='text-align: center; font-size: 1.05em; margin-bottom: 0.5em;'>"
        + " &nbsp;→&nbsp; ".join(parts)
        + "</div>",
        unsafe_allow_html=True,
    )
    st.divider()


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

# Common virus presets: (display name, TaxID)
VIRUS_PRESETS: list[tuple[str, int]] = [
    ("SARS-CoV-2", 2697049),
    ("SARS-CoV", 694009),
    ("MERS-CoV", 1335626),
    ("Influenza A", 11320),
    ("Influenza B", 11520),
    ("RSV-A (Respiratory Syncytial Virus A)", 208893),
    ("RSV-B (Respiratory Syncytial Virus B)", 208895),
    ("MPOX (Monkeypox virus)", 10244),
    ("Ebola virus (Zaire)", 186538),
    ("Marburg virus", 11269),
    ("Dengue virus 1", 11053),
    ("Dengue virus 2", 11060),
    ("Dengue virus 3", 11069),
    ("Dengue virus 4", 11070),
    ("Zika virus", 64320),
    ("Chikungunya virus", 37124),
    ("HIV-1", 11676),
    ("HIV-2", 11709),
    ("Hepatitis B virus", 10407),
    ("Hepatitis C virus", 11103),
    ("Human metapneumovirus", 162145),
    ("Measles virus", 11234),
    ("Mumps virus", 11161),
    ("Rubella virus", 11041),
    ("Norovirus GII", 142786),
    ("Rotavirus A", 28875),
    ("Adenovirus (Human mastadenovirus C)", 129951),
    ("Enterovirus D68", 42789),
    ("West Nile virus", 11082),
    ("Yellow Fever virus", 11089),
    ("Japanese Encephalitis virus", 11072),
    ("Lassa virus", 11620),
    ("Crimean-Congo hemorrhagic fever virus", 1980459),
    ("Nipah virus", 121227),
    ("Hendra virus", 63330),
    ("Rabies virus", 11292),
    ("Variola virus (Smallpox)", 10255),
    ("Human parainfluenza virus 1", 12730),
    ("Human parainfluenza virus 3", 11216),
    ("Human bocavirus 1", 329641),
    ("Human rhinovirus A", 147711),
]


def _available_fasta() -> list[str]:
    """Return sorted list of FASTA stem names in target_seqs/original/."""
    if not FASTA_DIR.exists():
        return []
    return sorted({
        p.stem for p in FASTA_DIR.iterdir()
        if p.suffix in (".fa", ".fasta", ".fna") and p.is_file()
    })


def _fasta_info(stem: str) -> dict:
    """Return sequence count and file size for a FASTA file."""
    path = FASTA_DIR / f"{stem}.fa"
    if not path.exists():
        # try other extensions
        for ext in (".fasta", ".fna"):
            alt = FASTA_DIR / f"{stem}{ext}"
            if alt.exists():
                path = alt
                break
    if not path.exists():
        return {"seqs": 0, "size": 0}
    with open(path) as f:
        n_seqs = sum(1 for line in f if line.startswith(">"))
    return {"seqs": n_seqs, "size": path.stat().st_size}


def _format_size(nbytes: int) -> str:
    for unit in ("B", "KB", "MB"):
        if nbytes < 1024:
            return f"{nbytes:.0f} {unit}"
        nbytes /= 1024
    return f"{nbytes:.1f} GB"


# Default parameter values (matches params.txt.template)
DEFAULT_PARAMS: dict = {
    "DESIGN_WINDOW": 500,
    "MAX_PRIMER_CANDIDATES": 10000,
    "TM_MIN": 55,
    "TM_MAX": 60,
    "GC_MAX": 60,
    "DG_MIN": -6,
    "PRIMER_LEN_MIN": 19,
    "PRIMER_LEN_MAX": 21,
    "TILING_STEP": 5,
    "PROBE_LEN_MIN": 24,
    "PROBE_LEN_MAX": 28,
    "PROBE_TM_MIN": 65,
    "PROBE_TM_MAX": 70,
    "PROBE_HOMOPOLYMER_MAX": 3,
    "PROBE_AVOID_5PRIME_G": True,
    "PROBE_MAX_MISMATCHES": 2,
    "PROBE_AMPLICON_BUFFER": 20,
    "MIN_PROBES_PER_PAIR": 1,
    "MAX_PROBES_PER_PAIR": 5,
    "AMPLEN_MIN": 60,
    "AMPLEN_MAX": 200,
    "OFFLEN_MIN": 50,
    "OFFLEN_MAX": 5000,
    "NUM_TOP_SENSITIVITY": 200,
}


def _init_params():
    """Initialize session-state params dict from defaults if not set."""
    if "params" not in st.session_state:
        st.session_state.params = dict(DEFAULT_PARAMS)


def _find_tool(name: str) -> str | None:
    """Find a tool in PATH or the current Python environment's bin directory."""
    found = shutil.which(name)
    if found:
        return found
    candidate = Path(sys.executable).parent / name
    if candidate.exists():
        return str(candidate)
    return None


def _check_tool(name: str) -> bool:
    return _find_tool(name) is not None


# ---------------------------------------------------------------------------
# Sidebar
# ---------------------------------------------------------------------------

def _render_sidebar():
    with st.sidebar:
        st.title("qPrimer Designer")
        st.caption(f"Working directory: `{PROJECT_ROOT}`")

        st.divider()

        page = _current_page()

        if page == "home":
            st.markdown("**Home**")
        else:
            if st.button("< Back to Home", use_container_width=True):
                _navigate("home")
                st.rerun()

            st.divider()

            mode = st.session_state.get("mode", "Singleplex")
            st.markdown(f"**Mode:** {mode}")

            probe = st.session_state.get("probe_enabled", False)
            st.markdown(f"**Probe mode:** {'On' if probe else 'Off'}")

            fastas = _available_fasta()
            st.markdown(f"**FASTA files:** {len(fastas)}")

            if mode == "Singleplex":
                targets = st.session_state.get("targets", [])
                st.markdown(f"**Targets:** {', '.join(targets) if targets else 'none'}")
            elif mode == "Multiplex":
                panel = st.session_state.get("panel", [])
                st.markdown(f"**Panel:** {', '.join(panel) if panel else 'none'}")



# ---------------------------------------------------------------------------
# Home page
# ---------------------------------------------------------------------------

def _page_home():
    st.markdown(
        "<h1 style='text-align: center;'>qPrimer Designer</h1>",
        unsafe_allow_html=True,
    )
    st.markdown(
        "<p style='text-align: center; font-size: 1.2em; color: gray;'>"
        "ML-guided PCR primer design and evaluation"
        "</p>",
        unsafe_allow_html=True,
    )

    st.divider()

    # Schematic
    if SCHEMATIC_PATH.exists():
        col_pad_l, col_img, col_pad_r = st.columns([1, 2, 1])
        with col_img:
            st.image(str(SCHEMATIC_PATH), use_container_width=True)

    st.markdown("")

    # Description
    st.markdown(
        """
qPrimer Designer is an AI-powered tool for PCR diagnostics primer design and
evaluation. It enables a **scalable, data-driven approach for adaptive
diagnostic design** against rapidly evolving pathogens — accessible anywhere
in the world, without requiring deep specialist expertise.

The tool is built on **machine learning models** that predict PCR activity for
a given primer pair and target sequence. These models were trained on
experimental data from over **50,000 unique primer pair–target sequence
combinations**, and significantly outperform predictions based on mismatch
counts or free energy alone
([Baek, Hsu et al., NeurIPS 2025 Workshop on AI for Science](https://openreview.net/forum?id=GJaNhMEZGM#discussion)).

Its core capabilities include **Design** — generating new primer sets for
given pathogen sequences — and **Evaluate** — assessing the performance of
existing primer sets. Pathogen sequences can be uploaded directly as FASTA
files or fetched from NCBI Virus by specifying a taxonomy ID and metadata
parameters. In addition to the target pathogen, users can specify
cross-reactivity sequences from other pathogens or host genomes to evaluate
and minimize non-specific amplification. Finally, the **Monitor** feature
tracks the validity of existing designs against newly observed pathogen
sequences by periodically running fetch and evaluate, and delivering results
via email reports.
"""
    )

    st.divider()

    # Action button illustrations (inline SVG)
    _ICON_DESIGN = """
    <svg viewBox="0 0 112 80" width="112" height="80" xmlns="http://www.w3.org/2000/svg">
      <!-- Target sequence -->
      <line x1="5" y1="20" x2="107" y2="20" stroke="#888" stroke-width="3" stroke-linecap="round"/>
      <!-- Unselected pair left -->
      <rect x="4" y="44" width="10" height="3" rx="1" fill="#4CAF50" opacity="0.25"/>
      <polygon points="15,45.5 13,43.5 13,47.5" fill="#4CAF50" opacity="0.25"/>
      <polygon points="18,45.5 20,43.5 20,47.5" fill="#E8451E" opacity="0.25"/>
      <rect x="19" y="44" width="10" height="3" rx="1" fill="#E8451E" opacity="0.25"/>
      <!-- Selected primer pair -->
      <rect x="38" y="41" width="14" height="5" rx="2" fill="#4CAF50"/>
      <polygon points="54,43.5 50,39.5 50,47.5" fill="#4CAF50"/>
      <polygon points="58,43.5 62,39.5 62,47.5" fill="#E8451E"/>
      <rect x="60" y="41" width="14" height="5" rx="2" fill="#E8451E"/>
      <!-- Selection box -->
      <rect x="33" y="36" width="46" height="16" rx="3" fill="none" stroke="#2E86AB" stroke-width="2"/>
      <!-- Unselected pair right -->
      <rect x="84" y="44" width="10" height="3" rx="1" fill="#4CAF50" opacity="0.25"/>
      <polygon points="95,45.5 93,43.5 93,47.5" fill="#4CAF50" opacity="0.25"/>
      <polygon points="98,45.5 100,43.5 100,47.5" fill="#E8451E" opacity="0.25"/>
      <rect x="99" y="44" width="10" height="3" rx="1" fill="#E8451E" opacity="0.25"/>
      <!-- Unselected pairs below -->
      <rect x="22" y="62" width="10" height="3" rx="1" fill="#4CAF50" opacity="0.25"/>
      <polygon points="33,63.5 31,61.5 31,65.5" fill="#4CAF50" opacity="0.25"/>
      <polygon points="36,63.5 38,61.5 38,65.5" fill="#E8451E" opacity="0.25"/>
      <rect x="37" y="62" width="10" height="3" rx="1" fill="#E8451E" opacity="0.25"/>
      <rect x="60" y="62" width="10" height="3" rx="1" fill="#4CAF50" opacity="0.25"/>
      <polygon points="71,63.5 69,61.5 69,65.5" fill="#4CAF50" opacity="0.25"/>
      <polygon points="74,63.5 76,61.5 76,65.5" fill="#E8451E" opacity="0.25"/>
      <rect x="75" y="62" width="10" height="3" rx="1" fill="#E8451E" opacity="0.25"/>
    </svg>
    """

    _ICON_EVALUATE = """
    <svg viewBox="0 0 100 80" width="100" height="80" xmlns="http://www.w3.org/2000/svg">
      <!-- Target sequence line -->
      <line x1="6" y1="20" x2="94" y2="20" stroke="#888" stroke-width="3" stroke-linecap="round"/>
      <!-- Forward primer -->
      <rect x="10" y="28" width="18" height="4" rx="1.5" fill="#4CAF50"/>
      <polygon points="30,30 27,27 27,33" fill="#4CAF50"/>
      <!-- Match/mismatch ticks for forward -->
      <line x1="13" y1="23" x2="13" y2="27" stroke="#4CAF50" stroke-width="1.5"/>
      <line x1="17" y1="23" x2="17" y2="27" stroke="#4CAF50" stroke-width="1.5"/>
      <line x1="21" y1="23" x2="21" y2="27" stroke="#E8451E" stroke-width="1.5"/>
      <line x1="25" y1="23" x2="25" y2="27" stroke="#4CAF50" stroke-width="1.5"/>
      <line x1="29" y1="23" x2="29" y2="27" stroke="#4CAF50" stroke-width="1.5"/>
      <!-- Reverse primer -->
      <polygon points="58,30 61,27 61,33" fill="#E8451E"/>
      <rect x="60" y="28" width="18" height="4" rx="1.5" fill="#E8451E"/>
      <!-- Match/mismatch ticks for reverse -->
      <line x1="61" y1="23" x2="61" y2="27" stroke="#4CAF50" stroke-width="1.5"/>
      <line x1="65" y1="23" x2="65" y2="27" stroke="#4CAF50" stroke-width="1.5"/>
      <line x1="69" y1="23" x2="69" y2="27" stroke="#E8451E" stroke-width="1.5"/>
      <line x1="73" y1="23" x2="73" y2="27" stroke="#4CAF50" stroke-width="1.5"/>
      <line x1="77" y1="23" x2="77" y2="27" stroke="#4CAF50" stroke-width="1.5"/>
      <!-- Magnifying glass -->
      <circle cx="65" cy="57" r="13" fill="none" stroke="#555" stroke-width="2.5"/>
      <line x1="74" y1="66" x2="83" y2="75" stroke="#555" stroke-width="3" stroke-linecap="round"/>
      <!-- Checkmark inside -->
      <polyline points="59,57 63,62 72,52" fill="none" stroke="#4CAF50" stroke-width="2.5" stroke-linecap="round" stroke-linejoin="round"/>
    </svg>
    """

    _ICON_MONITOR = """
    <svg viewBox="0 0 120 80" width="120" height="80" xmlns="http://www.w3.org/2000/svg">
      <!-- Phylogenetic tree -->
      <line x1="6" y1="40" x2="22" y2="40" stroke="#555" stroke-width="2.5"/>
      <line x1="22" y1="40" x2="22" y2="18" stroke="#555" stroke-width="2.5"/>
      <line x1="22" y1="40" x2="22" y2="62" stroke="#555" stroke-width="2.5"/>
      <!-- Upper branch -->
      <line x1="22" y1="18" x2="50" y2="18" stroke="#555" stroke-width="2.5"/>
      <line x1="50" y1="18" x2="50" y2="10" stroke="#555" stroke-width="2.5"/>
      <line x1="50" y1="18" x2="50" y2="26" stroke="#555" stroke-width="2.5"/>
      <line x1="50" y1="10" x2="78" y2="10" stroke="#555" stroke-width="2.5"/>
      <line x1="50" y1="26" x2="78" y2="26" stroke="#555" stroke-width="2.5"/>
      <!-- Lower branch -->
      <line x1="22" y1="62" x2="50" y2="62" stroke="#555" stroke-width="2.5"/>
      <line x1="50" y1="62" x2="50" y2="54" stroke="#555" stroke-width="2.5"/>
      <line x1="50" y1="62" x2="50" y2="70" stroke="#555" stroke-width="2.5"/>
      <line x1="50" y1="54" x2="78" y2="54" stroke="#555" stroke-width="2.5"/>
      <line x1="50" y1="70" x2="78" y2="70" stroke="#555" stroke-width="2.5"/>
      <!-- Leaf dots -->
      <circle cx="81" cy="10" r="3" fill="#2E86AB"/>
      <circle cx="81" cy="26" r="3" fill="#2E86AB"/>
      <circle cx="81" cy="54" r="3" fill="#2E86AB"/>
      <circle cx="81" cy="70" r="3" fill="#2E86AB"/>
      <!-- Time checkpoint lines (vertical dashed) -->
      <line x1="38" y1="4" x2="38" y2="76" stroke="#4CAF50" stroke-width="1.5" stroke-dasharray="3,3" opacity="0.6"/>
      <line x1="66" y1="4" x2="66" y2="76" stroke="#4CAF50" stroke-width="1.5" stroke-dasharray="3,3" opacity="0.6"/>
      <!-- Checkmarks in circles at checkpoint 1 -->
      <circle cx="38" cy="18" r="5" fill="white" stroke="#4CAF50" stroke-width="1.5"/>
      <polyline points="35,18 37,20 41,15" fill="none" stroke="#4CAF50" stroke-width="1.8" stroke-linecap="round" stroke-linejoin="round"/>
      <circle cx="38" cy="62" r="5" fill="white" stroke="#4CAF50" stroke-width="1.5"/>
      <polyline points="35,62 37,64 41,59" fill="none" stroke="#4CAF50" stroke-width="1.8" stroke-linecap="round" stroke-linejoin="round"/>
      <!-- Checkmarks in circles at checkpoint 2 -->
      <circle cx="66" cy="10" r="5" fill="white" stroke="#4CAF50" stroke-width="1.5"/>
      <polyline points="63,10 65,12 69,7" fill="none" stroke="#4CAF50" stroke-width="1.8" stroke-linecap="round" stroke-linejoin="round"/>
      <circle cx="66" cy="26" r="5" fill="white" stroke="#4CAF50" stroke-width="1.5"/>
      <polyline points="63,26 65,28 69,23" fill="none" stroke="#4CAF50" stroke-width="1.8" stroke-linecap="round" stroke-linejoin="round"/>
      <circle cx="66" cy="54" r="5" fill="white" stroke="#4CAF50" stroke-width="1.5"/>
      <polyline points="63,54 65,56 69,51" fill="none" stroke="#4CAF50" stroke-width="1.8" stroke-linecap="round" stroke-linejoin="round"/>
      <circle cx="66" cy="70" r="5" fill="white" stroke="#E8451E" stroke-width="1.5"/>
      <line x1="63" y1="67" x2="69" y2="73" stroke="#E8451E" stroke-width="1.8" stroke-linecap="round"/>
      <line x1="69" y1="67" x2="63" y2="73" stroke="#E8451E" stroke-width="1.8" stroke-linecap="round"/>
    </svg>
    """

    # Action buttons
    col1, col2, col3 = st.columns(3)

    with col1:
        st.markdown(
            f'<div style="text-align: center; padding: 1em 0;">{_ICON_DESIGN}</div>',
            unsafe_allow_html=True,
        )
        if st.button(
            "Design new primer sets",
            use_container_width=True,
            type="primary",
        ):
            st.session_state.workflow = "design"
            _clear_pipeline_state()
            _navigate("select_target")
            st.rerun()
        st.caption(
            "Generate optimized primer sets for your target pathogen sequences. "
            "Supports singleplex, multiplex, and probe design modes."
        )

    with col2:
        st.markdown(
            f'<div style="text-align: center; padding: 1em 0;">{_ICON_EVALUATE}</div>',
            unsafe_allow_html=True,
        )
        if st.button(
            "Evaluate existing primer sets",
            use_container_width=True,
            type="primary",
        ):
            st.session_state.workflow = "evaluate"
            _clear_pipeline_state()
            _navigate("select_target")
            st.rerun()
        st.caption(
            "Assess the performance of your existing primers against target sequences. "
            "Upload a primer set or paste individual sequences."
        )

    with col3:
        st.markdown(
            f'<div style="text-align: center; padding: 1em 0;">{_ICON_MONITOR}</div>',
            unsafe_allow_html=True,
        )
        if st.button(
            "Monitor primer performance",
            use_container_width=True,
            type="primary",
        ):
            st.session_state.workflow = "monitor"
            _clear_pipeline_state()
            _navigate("select_target")
            st.rerun()
        st.caption(
            "Track the validity of existing primer designs against newly observed "
            "pathogen sequences."
        )

    st.divider()

    # Contact
    st.markdown(
        "<div style='text-align: center; color: gray; font-size: 0.9em;'>"
        "<strong>Contact</strong><br>"
        "S. Chan Baek (<a href='mailto:baekseun@broadinstitute.org'>baekseun@broadinstitute.org</a>) · "
        "Kenneth B Hsu (<a href='mailto:khsu@broadinstitute.org'>khsu@broadinstitute.org</a>)"
        "</div>",
        unsafe_allow_html=True,
    )


# ---------------------------------------------------------------------------
# Tab 1: Files
# ---------------------------------------------------------------------------

def _tab_files():
    st.header("FASTA Files")
    st.write("Upload FASTA files to `target_seqs/original/` for use as targets, cross-reactivity, or host references.")

    # Upload
    uploaded = st.file_uploader(
        "Upload FASTA file(s)",
        type=["fa", "fasta", "fna"],
        accept_multiple_files=True,
        key="fasta_upload",
    )
    if uploaded:
        for f in uploaded:
            p = Path(f.name)
            if p.suffix.lower() not in (".fa", ".fasta", ".fna"):
                st.error(f"Unsupported extension: {p.suffix}")
                continue
            # Normalize to .fa extension
            dest = FASTA_DIR / (p.stem + ".fa")
            dest.write_bytes(f.getvalue())
            # Remove any leftover files with other FASTA extensions
            for ext in (".fasta", ".fna"):
                old = FASTA_DIR / (p.stem + ext)
                if old.exists():
                    old.unlink()
            st.success(f"Saved {dest.name}")

    st.divider()

    # File listing
    fastas = _available_fasta()
    if not fastas:
        st.info("No FASTA files found. Upload files above.")
        return

    st.subheader(f"Available files ({len(fastas)})")

    for stem in fastas:
        info = _fasta_info(stem)
        col1, col2, col3, col4 = st.columns([3, 2, 2, 1])
        col1.write(f"**{stem}**")
        col2.write(f"{info['seqs']} sequences")
        col3.write(_format_size(info["size"]))
        if col4.button("Delete", key=f"del_{stem}"):
            # Find and delete the file
            for ext in (".fa", ".fasta", ".fna"):
                p = FASTA_DIR / f"{stem}{ext}"
                if p.exists():
                    p.unlink()
                    break
            st.rerun()


# ---------------------------------------------------------------------------
# Tab 2: Configuration
# ---------------------------------------------------------------------------

def _tab_configuration():
    st.header("Pipeline Configuration")

    fastas = _available_fasta()

    mode = st.radio(
        "Design mode",
        ["Singleplex", "Multiplex", "Evaluate"],
        key="mode",
        horizontal=True,
    )

    probe = st.checkbox("Enable probe design", key="probe_enabled")

    st.divider()

    if mode == "Singleplex":
        st.subheader("Singleplex settings")
        st.multiselect(
            "TARGETS (on-target sequences to design primers for)",
            options=fastas,
            key="targets",
        )
        st.multiselect(
            "CROSS (cross-reactivity sequences to avoid)",
            options=fastas,
            key="cross",
        )
        st.multiselect(
            "HOST (host genome sequences to avoid)",
            options=fastas,
            key="host",
        )

    elif mode == "Multiplex":
        st.subheader("Multiplex settings")
        st.multiselect(
            "PANEL (targets for multiplex panel)",
            options=fastas,
            key="panel",
        )
        st.multiselect(
            "HOST (host genome sequences to avoid)",
            options=fastas,
            key="host_multiplex",
        )

    elif mode == "Evaluate":
        st.subheader("Evaluate existing primers")

        eval_method = st.radio(
            "Input method",
            ["Paste sequences", "Upload primer FASTA"],
            key="eval_method",
            horizontal=True,
        )

        if eval_method == "Paste sequences":
            st.text_input("Forward primer sequence", key="eval_for")
            st.text_input("Reverse primer sequence", key="eval_rev")
        else:
            pset_upload = st.file_uploader(
                "Upload primer set FASTA",
                type=["fa", "fasta"],
                key="eval_pset_upload",
            )
            if pset_upload:
                pset_dir = PROJECT_ROOT / "evaluate"
                pset_dir.mkdir(parents=True, exist_ok=True)
                pset_path = pset_dir / pset_upload.name
                pset_path.write_bytes(pset_upload.getvalue())
                st.session_state["eval_pset_path"] = str(pset_path)
                st.success(f"Saved {pset_upload.name}")

        st.selectbox(
            "Evaluation target (must match an uploaded FASTA)",
            options=fastas if fastas else ["(no files available)"],
            key="eval_target",
        )


# ---------------------------------------------------------------------------
# Tab 3: Parameters
# ---------------------------------------------------------------------------

def _tab_parameters():
    _init_params()
    p = st.session_state.params

    if st.button("Reset to defaults"):
        st.session_state.params = dict(DEFAULT_PARAMS)
        st.rerun()

    # --- Primer generation ---
    with st.expander("Primer generation", expanded=True):
        p["DESIGN_WINDOW"] = st.number_input(
            "Alignment window size (bp)", value=int(p["DESIGN_WINDOW"]),
            min_value=50, step=50, key="p_design_window",
            help="Window size for dividing the MSA when selecting representative sequences.",
        )
        p["MAX_PRIMER_CANDIDATES"] = st.number_input(
            "MAX_PRIMER_CANDIDATES", value=int(p["MAX_PRIMER_CANDIDATES"]),
            min_value=100, step=1000, key="p_max_cand",
        )
        c1, c2 = st.columns(2)
        p["TM_MIN"] = c1.number_input(
            "TM_MIN", value=float(p["TM_MIN"]), step=0.5, key="p_tm_min",
        )
        p["TM_MAX"] = c2.number_input(
            "TM_MAX", value=float(p["TM_MAX"]), step=0.5, key="p_tm_max",
        )
        if p["TM_MIN"] >= p["TM_MAX"]:
            st.warning("TM_MIN should be less than TM_MAX")

        p["GC_MAX"] = st.number_input(
            "GC_MAX (%)", value=float(p["GC_MAX"]), step=1.0, key="p_gc_max",
        )
        p["DG_MIN"] = st.number_input(
            "DG_MIN (kcal/mol)", value=float(p["DG_MIN"]), step=0.5, key="p_dg_min",
        )

        c3, c4 = st.columns(2)
        p["PRIMER_LEN_MIN"] = c3.number_input(
            "PRIMER_LEN_MIN", value=int(p["PRIMER_LEN_MIN"]),
            min_value=15, max_value=35, key="p_len_min",
        )
        p["PRIMER_LEN_MAX"] = c4.number_input(
            "PRIMER_LEN_MAX", value=int(p["PRIMER_LEN_MAX"]),
            min_value=15, max_value=35, key="p_len_max",
        )
        if p["PRIMER_LEN_MIN"] > p["PRIMER_LEN_MAX"]:
            st.warning("PRIMER_LEN_MIN should be <= PRIMER_LEN_MAX")

        p["TILING_STEP"] = st.number_input(
            "TILING_STEP", value=int(p["TILING_STEP"]),
            min_value=1, max_value=50, key="p_tiling",
        )

    # --- Probe generation ---
    probe_enabled = st.session_state.get("probe_enabled", False)
    if probe_enabled:
        with st.expander("Probe generation", expanded=True):
            c5, c6 = st.columns(2)
            p["PROBE_LEN_MIN"] = c5.number_input(
                "PROBE_LEN_MIN", value=int(p["PROBE_LEN_MIN"]),
                min_value=15, max_value=40, key="p_probe_len_min",
            )
            p["PROBE_LEN_MAX"] = c6.number_input(
                "PROBE_LEN_MAX", value=int(p["PROBE_LEN_MAX"]),
                min_value=15, max_value=40, key="p_probe_len_max",
            )
            if p["PROBE_LEN_MIN"] > p["PROBE_LEN_MAX"]:
                st.warning("PROBE_LEN_MIN should be <= PROBE_LEN_MAX")

            c7, c8 = st.columns(2)
            p["PROBE_TM_MIN"] = c7.number_input(
                "PROBE_TM_MIN", value=float(p["PROBE_TM_MIN"]),
                step=0.5, key="p_probe_tm_min",
            )
            p["PROBE_TM_MAX"] = c8.number_input(
                "PROBE_TM_MAX", value=float(p["PROBE_TM_MAX"]),
                step=0.5, key="p_probe_tm_max",
            )
            if p["PROBE_TM_MIN"] >= p["PROBE_TM_MAX"]:
                st.warning("PROBE_TM_MIN should be less than PROBE_TM_MAX")

            p["PROBE_HOMOPOLYMER_MAX"] = st.number_input(
                "PROBE_HOMOPOLYMER_MAX", value=int(p["PROBE_HOMOPOLYMER_MAX"]),
                min_value=1, max_value=10, key="p_probe_hp",
            )
            p["PROBE_AVOID_5PRIME_G"] = st.checkbox(
                "PROBE_AVOID_5PRIME_G", value=bool(p["PROBE_AVOID_5PRIME_G"]),
                key="p_probe_5g",
            )

        with st.expander("Probe mode configuration", expanded=False):
            p["PROBE_MAX_MISMATCHES"] = st.number_input(
                "PROBE_MAX_MISMATCHES", value=int(p["PROBE_MAX_MISMATCHES"]),
                min_value=0, max_value=10, key="p_probe_mm",
            )
            p["PROBE_AMPLICON_BUFFER"] = st.number_input(
                "PROBE_AMPLICON_BUFFER (nt)", value=int(p["PROBE_AMPLICON_BUFFER"]),
                min_value=0, max_value=100, key="p_probe_buf",
            )
            c9, c10 = st.columns(2)
            p["MIN_PROBES_PER_PAIR"] = c9.number_input(
                "MIN_PROBES_PER_PAIR", value=int(p["MIN_PROBES_PER_PAIR"]),
                min_value=0, max_value=20, key="p_probe_min_pp",
            )
            p["MAX_PROBES_PER_PAIR"] = c10.number_input(
                "MAX_PROBES_PER_PAIR", value=int(p["MAX_PROBES_PER_PAIR"]),
                min_value=1, max_value=50, key="p_probe_max_pp",
            )

    # --- Input preparation ---
    with st.expander("Input preparation (amplicon lengths)", expanded=False):
        c11, c12 = st.columns(2)
        p["AMPLEN_MIN"] = c11.number_input(
            "AMPLEN_MIN", value=int(p["AMPLEN_MIN"]),
            min_value=30, step=10, key="p_amp_min",
        )
        p["AMPLEN_MAX"] = c12.number_input(
            "AMPLEN_MAX", value=int(p["AMPLEN_MAX"]),
            min_value=30, step=10, key="p_amp_max",
        )
        if p["AMPLEN_MIN"] >= p["AMPLEN_MAX"]:
            st.warning("AMPLEN_MIN should be less than AMPLEN_MAX")

        c13, c14 = st.columns(2)
        p["OFFLEN_MIN"] = c13.number_input(
            "OFFLEN_MIN", value=int(p["OFFLEN_MIN"]),
            min_value=30, step=10, key="p_off_min",
        )
        p["OFFLEN_MAX"] = c14.number_input(
            "OFFLEN_MAX", value=int(p["OFFLEN_MAX"]),
            min_value=30, step=100, key="p_off_max",
        )

    # --- Evaluation ---
    with st.expander("Evaluation", expanded=False):
        p["NUM_TOP_SENSITIVITY"] = st.number_input(
            "NUM_TOP_SENSITIVITY", value=int(p["NUM_TOP_SENSITIVITY"]),
            min_value=10, step=10, key="p_num_top",
        )


# ---------------------------------------------------------------------------
# Tab 4: Run
# ---------------------------------------------------------------------------

def _build_command() -> list[str]:
    """Build the snakemake command list from current session state."""
    mode = st.session_state.get("mode", "Singleplex")
    probe = st.session_state.get("probe_enabled", False)
    cores = st.session_state.get("cores", 1)
    dry_run = st.session_state.get("dry_run", False)

    snakemake_bin = _find_tool("snakemake") or "snakemake"
    cmd = [snakemake_bin, "-s", "Snakefile", "--cores", str(cores)]

    config_args: list[str] = []
    if mode == "Multiplex":
        config_args.append("multiplex=1")
    if probe:
        config_args.append("probe=1")
    if mode == "Evaluate":
        config_args.append("evaluate=1")
        eval_method = st.session_state.get("eval_method", "Paste sequences")
        if eval_method == "Paste sequences":
            fwd = st.session_state.get("eval_for", "")
            rev = st.session_state.get("eval_rev", "")
            if fwd:
                config_args.append(f"for={fwd}")
            if rev:
                config_args.append(f"rev={rev}")
        else:
            pset_path = st.session_state.get("eval_pset_path", "")
            if pset_path:
                config_args.append(f"pset={pset_path}")

    if config_args:
        cmd += ["--config"] + config_args

    if dry_run:
        cmd.append("--dry-run")

    return cmd


def _preflight_checks() -> list[str]:
    """Return a list of preflight error messages. Empty = all OK."""
    errors: list[str] = []
    mode = st.session_state.get("mode", "Singleplex")

    # Check required targets
    if mode == "Singleplex":
        if not st.session_state.get("targets"):
            errors.append("No TARGETS selected in Configuration tab.")
    elif mode == "Multiplex":
        if not st.session_state.get("panel"):
            errors.append("No PANEL targets selected in Configuration tab.")
    elif mode == "Evaluate":
        eval_method = st.session_state.get("eval_method", "Paste sequences")
        if eval_method == "Paste sequences":
            if not st.session_state.get("eval_for") or not st.session_state.get("eval_rev"):
                errors.append("Forward and reverse primer sequences are required.")
        else:
            if not st.session_state.get("eval_pset_path"):
                errors.append("Upload a primer set FASTA file.")
        if not st.session_state.get("eval_target") or st.session_state.get("eval_target") == "(no files available)":
            errors.append("Select an evaluation target FASTA.")

    # Check tools
    for tool in ("snakemake", "bowtie2", "mafft", "sam2pairwise"):
        if not _check_tool(tool):
            errors.append(f"Tool '{tool}' not found in PATH.")

    return errors


def _write_pipeline_files():
    """Write Snakefile and params.txt to project root."""
    mode = st.session_state.get("mode", "Singleplex")
    probe = st.session_state.get("probe_enabled", False)

    targets = st.session_state.get("targets", [])
    cross = st.session_state.get("cross", [])
    panel = st.session_state.get("panel", [])

    if mode == "Multiplex":
        host = st.session_state.get("host_multiplex", [])
    else:
        host = st.session_state.get("host", [])

    # For evaluate mode, set TARGETS to the eval target so the template
    # doesn't fail the safety check
    if mode == "Evaluate":
        eval_target = st.session_state.get("eval_target", "")
        if eval_target and eval_target != "(no files available)":
            targets = [eval_target]

    run_id = datetime.now().strftime("%Y%m%d_%H%M%S")
    st.session_state.run_id = run_id

    snakefile_content = build_snakefile(
        targets=targets,
        cross=cross,
        host=host,
        panel=panel,
        run_id=run_id,
    )
    (PROJECT_ROOT / "Snakefile").write_text(snakefile_content)

    _init_params()
    params_content = build_params_txt(st.session_state.params)
    (PROJECT_ROOT / "params.txt").write_text(params_content)


# Pipeline rule steps grouped for progress display
_PIPELINE_RULES_DESIGN = [
    ("make_MSA", "Building multiple sequence alignment"),
    ("pick_representative_seqs", "Selecting representative sequences"),
    ("generate_primers", "Generating primer candidates"),
    ("build_index", "Building bowtie2 index"),
    ("align", "Aligning primers to targets"),
    ("parse_map", "Parsing alignment results"),
    ("process_map", "Processing alignment maps"),
    ("prepare_input", "Preparing ML model input"),
    ("evaluate", "Evaluating with ML model"),
    ("rescue_evaluate", "Rescue evaluation"),
    ("filter_primer_list", "Filtering and ranking primers"),
    ("check_coverage", "Checking primer coverage"),
]

_PIPELINE_RULES_EVALUATE = [
    ("build_index", "Building bowtie2 index"),
    ("align", "Aligning primers to targets"),
    ("parse_map", "Parsing alignment results"),
    ("process_map", "Processing alignment maps"),
    ("prepare_input", "Preparing ML model input"),
    ("evaluate_pset", "Evaluating primer set"),
]


def _get_pipeline_rules() -> list[tuple[str, str]]:
    mode = st.session_state.get("mode", "Singleplex")
    if mode == "Evaluate":
        return _PIPELINE_RULES_EVALUATE
    return _PIPELINE_RULES_DESIGN


def _render_pipeline_progress(rules: list[tuple[str, str]], completed_rules: set[str],
                               current_rule: str, detail: str = ""):
    """Render pipeline progress like the fetch UI."""
    lines = []
    current_found = False
    for rule_name, label in rules:
        if rule_name in completed_rules:
            lines.append(f"✅  {label}")
        elif rule_name == current_rule:
            lines.append(f"🔵  {label}...")
            current_found = True
        elif current_found:
            lines.append(f"⚪  {label}")
        else:
            lines.append(f"⚪  {label}")
    return "\n".join(lines)


def _save_run_config():
    """Save all run settings to a JSON file alongside the final output."""
    run_id = st.session_state.get("run_id", "")
    if not run_id:
        return

    workflow = st.session_state.get("workflow", "design")
    mode = st.session_state.get("mode", "Singleplex")

    config = {
        "run_id": run_id,
        "timestamp": datetime.now().isoformat(),
        "workflow": workflow,
        "mode": mode,
        "targets": st.session_state.get("targets", []),
        "cross_reactivity": st.session_state.get("cross", []),
        "host": st.session_state.get("host", []),
        "probe_enabled": st.session_state.get("probe_enabled", False),
        "cores": st.session_state.get("cores", 1),
    }

    # Fetch settings (target)
    fetch_prefix = "fetch"
    fetch_config = {}
    virus_sel = st.session_state.get(f"{fetch_prefix}_virus_select")
    if virus_sel:
        fetch_config["virus"] = virus_sel
    for key_suffix in ("target_name", "nuc_completeness", "segment",
                       "min_seq_length", "max_seq_length",
                       "min_release_date", "max_release_date",
                       "geo_location", "host",
                       "limit", "min_collection_date", "max_collection_date",
                       "max_ambiguous_chars", "refseq_only"):
        val = st.session_state.get(f"{fetch_prefix}_{key_suffix}")
        if val is not None and val != "" and val is not False:
            fetch_config[key_suffix] = val
    if fetch_config:
        config["fetch_settings"] = fetch_config

    # Evaluate-specific
    if mode == "Evaluate":
        config["eval_method"] = st.session_state.get("eval_method", "")
        config["eval_target"] = st.session_state.get("eval_target", "")
        if st.session_state.get("eval_method") == "Paste sequences":
            config["eval_forward"] = st.session_state.get("eval_for", "")
            config["eval_reverse"] = st.session_state.get("eval_rev", "")
        else:
            config["eval_pset_path"] = st.session_state.get("eval_pset_path", "")

    # Parameters
    config["parameters"] = dict(st.session_state.get("params", {}))

    # Determine output directory
    if mode == "Evaluate":
        out_dir = PROJECT_ROOT / "evaluate" / run_id
    else:
        out_dir = PROJECT_ROOT / "final" / run_id
    out_dir.mkdir(parents=True, exist_ok=True)

    # Custom JSON serializer for dates
    def _json_default(obj):
        if isinstance(obj, (date, datetime)):
            return obj.isoformat()
        return str(obj)

    config_path = out_dir / "run_config.json"
    config_path.write_text(json.dumps(config, indent=2, default=_json_default, ensure_ascii=False))


def _tab_run():
    st.header("Run Pipeline")

    # Controls
    c1, c2 = st.columns(2)
    max_cpu = os.cpu_count() or 1
    c1.slider("CPU cores", min_value=1, max_value=max_cpu, value=min(4, max_cpu), key="cores")
    c2.checkbox("Dry run (plan only, no execution)", key="dry_run")

    # Preflight
    errors = _preflight_checks()
    if errors:
        for e in errors:
            st.error(e)

    rules = _get_pipeline_rules()

    # Clear previous result when entering Run page for first time
    if st.session_state.get("_run_page_fresh", True):
        st.session_state.pipeline_return_code = None
        st.session_state.pipeline_log = ""
        st.session_state.pipeline_completed_rules = set()
        st.session_state._run_page_fresh = False

    # Run / Stop buttons
    col_run, col_stop = st.columns(2)
    running = st.session_state.get("pipeline_running", False)

    with col_run:
        run_clicked = st.button("Run", disabled=running or bool(errors), type="primary")
    with col_stop:
        stop_clicked = st.button("Stop", disabled=not running, type="secondary")

    # Handle Stop
    if stop_clicked:
        pid = st.session_state.get("pipeline_pid")
        if pid:
            try:
                os.killpg(os.getpgid(pid), signal.SIGTERM)
            except (ProcessLookupError, PermissionError):
                try:
                    os.kill(pid, signal.SIGTERM)
                except (ProcessLookupError, PermissionError):
                    pass
        st.session_state.pipeline_running = False
        st.session_state.pipeline_return_code = -1
        st.rerun()

    # Handle Run click — set state and rerun so buttons update first
    if run_clicked:
        st.session_state.pipeline_running = True
        st.session_state.pipeline_should_start = True
        st.session_state.pipeline_return_code = None
        st.session_state.pipeline_log = ""
        st.session_state.pipeline_completed_rules = set()
        st.rerun()

    # Execute pipeline (triggered after rerun from Run click)
    if st.session_state.get("pipeline_should_start"):
        st.session_state.pipeline_should_start = False

        _write_pipeline_files()

        # Unlock snakemake directory in case of leftover locks
        snakemake_bin = _find_tool("snakemake") or "snakemake"
        subprocess.run(
            [snakemake_bin, "-s", "Snakefile", "--unlock", "--cores", "1"],
            capture_output=True, cwd=str(PROJECT_ROOT),
        )

        env = os.environ.copy()
        env["PYTHONUNBUFFERED"] = "1"
        env_bin = str(Path(sys.executable).parent)
        if env_bin not in env.get("PATH", ""):
            env["PATH"] = env_bin + os.pathsep + env.get("PATH", "")

        cmd = _build_command()

        progress_area = st.empty()
        status_area = st.empty()

        completed_rules: set[str] = set()
        current_rule = rules[0][0] if rules else ""
        log = ""
        start_time = time.time()

        def _update_progress():
            lines = []
            for rule_name, label in rules:
                if rule_name in completed_rules:
                    lines.append(f"✅ &nbsp; {label}")
                elif rule_name == current_rule:
                    lines.append(f"⏳ &nbsp; **{label}** ...")
                else:
                    lines.append(f"⚪ &nbsp; {label}")
            elapsed = int(time.time() - start_time)
            mins, secs = divmod(elapsed, 60)
            progress_area.markdown("<br>".join(lines), unsafe_allow_html=True)
            status_area.caption(f"Elapsed: {mins}m {secs}s")

        _update_progress()

        proc = subprocess.Popen(
            cmd,
            stdout=subprocess.PIPE,
            stderr=subprocess.STDOUT,
            cwd=str(PROJECT_ROOT),
            env=env,
            text=True,
        )
        st.session_state.pipeline_pid = proc.pid

        while proc.poll() is None:
            line = proc.stdout.readline()
            if line:
                log += line
                stripped = line.strip()
                if stripped.startswith("rule ") or stripped.startswith("localrule "):
                    rule_name = stripped.split("rule ")[1].rstrip(":")
                    if current_rule and current_rule != rule_name:
                        completed_rules.add(current_rule)
                    # Mark all rules before the current one as completed (skipped by snakemake)
                    for rn, _ in rules:
                        if rn == rule_name:
                            break
                        completed_rules.add(rn)
                    current_rule = rule_name
                    _update_progress()
                elif "Finished job" in stripped:
                    _update_progress()

        # Read remaining output
        remaining = proc.stdout.read()
        if remaining:
            log += remaining

        # Mark last rule as complete
        if current_rule:
            completed_rules.add(current_rule)

        rc = proc.returncode
        st.session_state.pipeline_running = False
        st.session_state.pipeline_return_code = rc
        st.session_state.pipeline_log = log
        st.session_state.pipeline_completed_rules = completed_rules

        elapsed = int(time.time() - start_time)
        mins, secs = divmod(elapsed, 60)

        # Final progress
        lines = []
        for rule_name, label in rules:
            if rule_name in completed_rules:
                lines.append(f"✅ &nbsp; {label}")
            else:
                lines.append(f"⚪ &nbsp; {label}")
        progress_area.markdown("<br>".join(lines), unsafe_allow_html=True)
        status_area.caption(f"Completed in {mins}m {secs}s")

        if rc == 0:
            _save_run_config()
            st.success("Pipeline finished successfully!")
        else:
            st.error(f"Pipeline failed with exit code {rc}.")
            with st.expander("Show error log"):
                last_lines = "\n".join(log.splitlines()[-50:])
                st.code(last_lines, language="text")

    # Show previous run result if not currently running
    elif not running and st.session_state.get("pipeline_return_code") is not None:
        rc = st.session_state.pipeline_return_code
        completed_rules = st.session_state.get("pipeline_completed_rules", set())

        lines = []
        for rule_name, label in rules:
            if rule_name in completed_rules:
                lines.append(f"✅ &nbsp; {label}")
            else:
                lines.append(f"⚪ &nbsp; {label}")
        st.markdown("<br>".join(lines), unsafe_allow_html=True)

        if rc == 0:
            st.success("Pipeline finished successfully!")
        else:
            st.error(f"Pipeline failed with exit code {rc}.")
            log = st.session_state.get("pipeline_log", "")
            if log:
                with st.expander("Show error log"):
                    last_lines = "\n".join(log.splitlines()[-50:])
                    st.code(last_lines, language="text")


# ---------------------------------------------------------------------------
# Tab 5: Results
# ---------------------------------------------------------------------------

def _tab_results():
    st.header("Results")

    st.markdown(
        "For help interpreting output columns, see the "
        "[Output Interpretation Guide]"
        "(https://github.com/broadinstitute/qprimer_designer/blob/main/docs/output_interpretation_guide.md)."
    )

    # --- Final CSV files ---
    st.subheader("Final output files")

    if FINAL_DIR.exists():
        csvs = sorted(FINAL_DIR.glob("**/*.csv"), reverse=True)
    else:
        csvs = []

    if csvs:
        # Show relative path from FINAL_DIR for each CSV
        csv_labels = [str(c.relative_to(FINAL_DIR)) for c in csvs]

        selected_csv_label = st.selectbox(
            "Select CSV to view",
            options=csv_labels,
            key="result_csv",
        )
        if selected_csv_label:
            import pandas as pd

            csv_path = FINAL_DIR / selected_csv_label
            try:
                df = pd.read_csv(csv_path)
                st.write(f"**{len(df)} primer pairs**")

                # Summary metrics
                mc1, mc2, mc3 = st.columns(3)
                mc1.metric("Pairs", len(df))
                if "coverage" in df.columns:
                    mc2.metric("Avg coverage", f"{df['coverage'].mean():.2f}")
                if "score" in df.columns:
                    mc3.metric("Best score", f"{df['score'].max():.4f}")

                st.dataframe(df, use_container_width=True)

                # Download button
                st.download_button(
                    "Download CSV",
                    data=csv_path.read_bytes(),
                    file_name=csv_path.name,
                    mime="text/csv",
                )
            except Exception as exc:
                st.error(f"Error reading {selected_csv_label}: {exc}")
    else:
        st.info("No CSV results in `final/` yet. Run the pipeline first.")

    st.divider()

    # --- Evaluate reports ---
    st.subheader("Evaluation reports")

    if EVALUATE_DIR.exists():
        xlsx_files = sorted(EVALUATE_DIR.glob("**/*.xlsx"), reverse=True)
    else:
        xlsx_files = []

    if xlsx_files:
        import pandas as pd

        # Show run timestamp and filename
        xlsx_labels = []
        for xf in xlsx_files:
            # Path is evaluate/{run_id}/{pset_name}/{file}.xlsx
            rel = xf.relative_to(EVALUATE_DIR)
            xlsx_labels.append(str(rel))

        selected_xlsx_label = st.selectbox(
            "Select report to view",
            options=xlsx_labels,
            key="result_xlsx",
        )
        if selected_xlsx_label:
            xf = EVALUATE_DIR / selected_xlsx_label
            try:
                # Show summary sheet
                df_summary = pd.read_excel(xf, sheet_name="summary", header=None).astype(str)
                st.write("**Summary**")
                st.dataframe(df_summary, use_container_width=True, hide_index=True)

                # Show detail sheet
                df_detail = pd.read_excel(xf, sheet_name="detail")
                st.write(f"**Detail** — {len(df_detail)} target alignments")
                st.dataframe(df_detail, use_container_width=True)
            except Exception as exc:
                st.error(f"Error reading {selected_xlsx_label}: {exc}")

            st.download_button(
                f"Download {xf.name}",
                data=xf.read_bytes(),
                file_name=xf.name,
                mime="application/vnd.openxmlformats-officedocument.spreadsheetml.sheet",
                key=f"dl_{selected_xlsx_label}",
            )
    else:
        st.info("No evaluation reports found.")

    st.divider()

    # --- Cleanup ---
    st.subheader("Cleanup")
    if st.button("Delete all pipeline outputs", type="secondary"):
        result = subprocess.run(
            ["snakemake", "-s", "Snakefile", "--delete-all-output", "--cores", "1"],
            capture_output=True,
            text=True,
            cwd=str(PROJECT_ROOT),
        )
        if result.returncode == 0:
            st.success("Pipeline outputs deleted.")
        else:
            st.error(f"Cleanup failed:\n{result.stderr}")


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def _page_select_target():
    """Common page for selecting pathogen target sequences."""
    _render_workflow_progress("select_target")

    workflow = st.session_state.get("workflow", "design")
    workflow_label = {"design": "Design", "evaluate": "Evaluate", "monitor": "Monitor"}[workflow]

    st.header(f"{workflow_label} — Select Target Sequences")
    st.markdown(
        "Which pathogen do you want to detect? "
        "Select target sequences from the list below, or add new ones."
    )

    st.divider()

    # --- Select from existing FASTA files ---
    fastas = _available_fasta()

    if fastas:
        # Build display labels with sequence count and file size
        labels = []
        for stem in fastas:
            info = _fasta_info(stem)
            labels.append(f"{stem}  ({info['seqs']} seqs, {_format_size(info['size'])})")

        selected_label = st.selectbox(
            "Select a target sequence",
            options=[""] + labels,
            format_func=lambda x: "Choose a pathogen..." if x == "" else x,
            key="target_select",
        )

        if selected_label:
            # Extract stem name (everything before the first double-space)
            selected_stem = selected_label.split("  (")[0]
            st.session_state.targets = [selected_stem]
        else:
            st.session_state.targets = []
    else:
        st.info("No FASTA files found. Upload or fetch sequences below.")

    st.divider()

    # --- Add new sequences ---
    st.subheader("Don't see your pathogen?")

    add_col1, add_col2 = st.columns(2)

    with add_col1:
        st.markdown("**Upload FASTA file**")
        uploaded = st.file_uploader(
            "Upload FASTA file(s)",
            type=["fa", "fasta", "fna"],
            accept_multiple_files=True,
            key="fasta_upload",
            label_visibility="collapsed",
        )
        if uploaded:
            for f in uploaded:
                p = Path(f.name)
                if p.suffix.lower() not in (".fa", ".fasta", ".fna"):
                    st.error(f"Unsupported extension: {p.suffix}")
                    continue
                dest = FASTA_DIR / (p.stem + ".fa")
                dest.write_bytes(f.getvalue())
                for ext in (".fasta", ".fna"):
                    old = FASTA_DIR / (p.stem + ext)
                    if old.exists():
                        old.unlink()
                st.success(f"Saved {dest.name}")
            st.rerun()

    with add_col2:
        st.markdown("**Fetch from NCBI Virus (using gget)**")
        _render_fetch_ui("fetch")

    st.divider()

    # --- Continue button ---
    selected = st.session_state.get("targets", [])
    if st.button(
        "Continue →",
        disabled=not selected,
        type="primary",
        use_container_width=True,
    ):
        _navigate("select_offtarget")
        st.rerun()

    if not selected:
        st.caption("Select a target sequence to continue.")


def _render_fetch_ui(prefix: str):
    """Render fetch-from-NCBI UI. prefix is used to namespace session state keys."""
    _OTHER_OPTION = "Other (enter TaxID manually)"
    virus_options = [f"{name}  (TaxID: {tid})" for name, tid in VIRUS_PRESETS]
    virus_options.append(_OTHER_OPTION)

    selected_virus = st.selectbox(
        "Virus",
        options=virus_options,
        index=None,
        placeholder="Type a virus name to search...",
        key=f"{prefix}_virus_select",
    )

    taxid = ""
    auto_name = ""
    if selected_virus == _OTHER_OPTION:
        taxid = st.text_input(
            "Taxonomy ID",
            placeholder="e.g., 2697049",
            key=f"{prefix}_taxid_manual",
        )
    elif selected_virus:
        taxid = selected_virus.split("TaxID: ")[1].rstrip(")")
        auto_name = selected_virus.split("  (TaxID:")[0]

    auto_target = auto_name
    if not auto_name and taxid:
        for name, tid in VIRUS_PRESETS:
            if str(tid) == taxid.strip():
                auto_target = name
                break

    prev_auto = st.session_state.get(f"_prev_auto_{prefix}", "")
    if auto_target != prev_auto:
        st.session_state[f"_prev_auto_{prefix}"] = auto_target
        sanitized = re.sub(r"[^\w\-]", "_", auto_target)
        sanitized = re.sub(r"_+", "_", sanitized).strip("_")
        st.session_state[f"{prefix}_target_name"] = sanitized

    target_name = st.text_input(
        "Target name",
        placeholder="e.g., SARS-CoV-2",
        key=f"{prefix}_target_name",
        help="Name for the output FASTA file.",
    )

    nuc_completeness = st.selectbox(
        "Nucleotide completeness",
        options=["complete", "partial"],
        key=f"{prefix}_nuc_completeness",
    )

    segment = st.text_input(
        "Gene segment",
        placeholder="e.g., 1, 2, 3 (leave blank for non-segmented)",
        key=f"{prefix}_segment",
    )

    col_len1, col_len2 = st.columns(2)
    with col_len1:
        st.number_input(
            "Min sequence length (bp)", min_value=0, value=None, step=100,
            key=f"{prefix}_min_seq_length",
        )
    with col_len2:
        st.number_input(
            "Max sequence length (bp)", min_value=0, value=None, step=100,
            key=f"{prefix}_max_seq_length",
        )

    col_d1, col_d2 = st.columns(2)
    with col_d1:
        st.date_input("Min release date", value=None, key=f"{prefix}_min_release_date")
    with col_d2:
        st.date_input("Max release date", value=None, key=f"{prefix}_max_release_date")

    st.text_input("Geographic location", placeholder="e.g., USA, Africa, Nigeria",
                   key=f"{prefix}_geo_location")
    st.text_input("Host", placeholder="e.g., human, Homo sapiens", key=f"{prefix}_host")

    with st.expander("Additional parameters"):
        st.number_input("Max sequences", min_value=1, max_value=100000, value=None,
                        step=100, key=f"{prefix}_limit")
        col_c1, col_c2 = st.columns(2)
        with col_c1:
            st.date_input("Min collection date", value=None, key=f"{prefix}_min_collection_date")
        with col_c2:
            st.date_input("Max collection date", value=None, key=f"{prefix}_max_collection_date")
        st.number_input("Max ambiguous characters", min_value=0, value=10, step=1,
                        key=f"{prefix}_max_ambiguous_chars")
        st.checkbox("RefSeq only", key=f"{prefix}_refseq_only")

    if st.button("Fetch sequences", disabled=not taxid or not target_name, key=f"{prefix}_fetch_btn"):
        gget_bin = shutil.which("gget")
        if not gget_bin:
            env_bin = Path(sys.executable).parent
            candidate = env_bin / "gget"
            if candidate.exists():
                gget_bin = str(candidate)
        if not gget_bin:
            st.error("'gget' is not installed. Run: `pip install gget`")
        else:
            gget_tmp = FASTA_DIR / f".gget_tmp_{target_name}"
            gget_tmp.mkdir(parents=True, exist_ok=True)

            cmd = [gget_bin, "virus", str(taxid).strip(), "--out", str(gget_tmp)]

            nuc = st.session_state.get(f"{prefix}_nuc_completeness", "complete")
            if nuc:
                cmd.extend(["--nuc_completeness", nuc])
            seg = st.session_state.get(f"{prefix}_segment", "").strip()
            if seg:
                cmd.extend(["--segment", seg])
            min_len = st.session_state.get(f"{prefix}_min_seq_length")
            if min_len:
                cmd.extend(["--min_seq_length", str(min_len)])
            max_len = st.session_state.get(f"{prefix}_max_seq_length")
            if max_len:
                cmd.extend(["--max_seq_length", str(max_len)])
            min_rel = st.session_state.get(f"{prefix}_min_release_date")
            if min_rel:
                cmd.extend(["--min_release_date", str(min_rel)])
            max_rel = st.session_state.get(f"{prefix}_max_release_date")
            if max_rel:
                cmd.extend(["--max_release_date", str(max_rel)])
            geo = st.session_state.get(f"{prefix}_geo_location", "").strip()
            if geo:
                cmd.extend(["--geographic_location", geo])
            host_val = st.session_state.get(f"{prefix}_host", "").strip()
            if host_val:
                cmd.extend(["--host", host_val])
            if st.session_state.get(f"{prefix}_refseq_only"):
                cmd.append("--refseq_only")
            min_col = st.session_state.get(f"{prefix}_min_collection_date")
            if min_col:
                cmd.extend(["--min_collection_date", str(min_col)])
            max_col = st.session_state.get(f"{prefix}_max_collection_date")
            if max_col:
                cmd.extend(["--max_collection_date", str(max_col)])
            max_amb = st.session_state.get(f"{prefix}_max_ambiguous_chars")
            if max_amb is not None:
                cmd.extend(["--max_ambiguous_chars", str(max_amb)])

            _FETCH_STEPS = [
                ("STEP 1", "Validating input"),
                ("STEP 2", "Checking optimized pathways"),
                ("STEP 3", "Fetching metadata from NCBI"),
                ("STEP 4", "Applying metadata filters"),
                ("STEP 5", "Downloading sequences"),
                ("STEP 6", "Applying sequence filters"),
                ("STEP 7", "Saving output files"),
                ("STEP 8", "Fetching GenBank metadata"),
            ]

            progress_area = st.empty()
            detail_area = st.empty()
            status_area = st.empty()

            def _render_progress(current_step, detail=""):
                lines = []
                for i, (_, label) in enumerate(_FETCH_STEPS):
                    if i < current_step:
                        lines.append(f"✅  {label}")
                    elif i == current_step:
                        lines.append(f"🔵  {label}...")
                    else:
                        lines.append(f"⚪  {label}")
                progress_area.markdown("\n".join(lines))
                if detail:
                    detail_area.caption(detail)

            _render_progress(-1)

            try:
                env = os.environ.copy()
                env["PYTHONUNBUFFERED"] = "1"
                proc = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
                                        text=True, env=env)
                current_step = -1
                detail = ""
                start_time = time.time()
                timeout = 1800

                while proc.poll() is None:
                    if time.time() - start_time > timeout:
                        proc.kill()
                        st.error("Fetch timed out after 30 minutes.")
                        break
                    line = proc.stdout.readline()
                    if line:
                        for i, (key, _) in enumerate(_FETCH_STEPS):
                            if key in line:
                                current_step = i
                                detail = ""
                                break
                        if "retrieved" in line and "records" in line:
                            detail = line.split("INFO - ")[-1].strip() if "INFO - " in line else ""
                        elif "Downloading sequences for" in line:
                            detail = line.split("INFO - ")[-1].strip() if "INFO - " in line else ""
                        elif "Received" in line and "sequences" in line:
                            detail = line.split("INFO - ")[-1].strip() if "INFO - " in line else ""
                        elif "PROCESS COMPLETED" in line:
                            current_step = len(_FETCH_STEPS)
                        _render_progress(current_step, detail)
                    elapsed = int(time.time() - start_time)
                    mins, secs = divmod(elapsed, 60)
                    status_area.caption(f"Elapsed: {mins}m {secs}s")

                remaining = proc.stdout.read()
                if remaining and "PROCESS COMPLETED" in remaining:
                    current_step = len(_FETCH_STEPS)
                    _render_progress(current_step)

                elapsed = int(time.time() - start_time)
                mins, secs = divmod(elapsed, 60)
                status_area.caption(f"Completed in {mins}m {secs}s")

                if proc.returncode != 0:
                    _render_progress(current_step, "")
                    st.error(f"gget failed (exit code {proc.returncode})")
                else:
                    _render_progress(len(_FETCH_STEPS))
                    fasta_files = list(gget_tmp.glob("*.fa")) + list(gget_tmp.glob("*.fasta"))
                    if fasta_files:
                        dest = FASTA_DIR / f"{target_name}.fa"
                        shutil.move(str(fasta_files[0]), str(dest))
                        from qprimer_designer.adapt_cli import _deduplicate_fasta
                        n_deduped = _deduplicate_fasta(dest)
                        with open(dest) as f:
                            n_seqs = sum(1 for line in f if line.startswith(">"))
                        msg = f"Downloaded {n_seqs} unique sequences → `{dest.name}`"
                        if n_deduped > 0:
                            msg += f" ({n_deduped} duplicates removed)"
                        st.success(msg)
                        detail_area.empty()
                    else:
                        st.warning("No FASTA output found. The query may have returned no results.")
            except Exception as exc:
                st.error(f"Fetch error: {exc}")
            finally:
                shutil.rmtree(gget_tmp, ignore_errors=True)


def _page_select_offtarget():
    """Page for selecting off-target (cross-reactivity) and host sequences."""
    _render_workflow_progress("select_offtarget")

    workflow = st.session_state.get("workflow", "design")
    workflow_label = {"design": "Design", "evaluate": "Evaluate", "monitor": "Monitor"}[workflow]

    st.header(f"{workflow_label} — Off-Target & Host Sequences")
    st.markdown(
        "Optionally specify sequences to check for cross-reactivity or host background. "
        "You can skip this step if not needed."
    )

    st.divider()

    fastas = _available_fasta()
    # Exclude already-selected targets from the options
    selected_targets = st.session_state.get("targets", [])
    offtarget_options = [s for s in fastas if s not in selected_targets]

    # --- Cross-reactivity sequences ---
    st.subheader("Cross-reactivity sequences")
    st.caption("Other pathogens you want to avoid detecting (off-target amplification)")

    if offtarget_options:
        cross_labels = []
        for stem in offtarget_options:
            info = _fasta_info(stem)
            cross_labels.append(f"{stem}  ({info['seqs']} seqs, {_format_size(info['size'])})")

        selected_cross = st.multiselect(
            "Select cross-reactivity sequences",
            options=cross_labels,
            default=[],
            key="cross_select",
            placeholder="Choose pathogens to avoid...",
        )
        st.session_state.cross = [label.split("  (")[0] for label in selected_cross]
    else:
        st.info("No additional FASTA files available.")
        st.session_state.cross = []

    # Add new cross-reactivity sequences
    with st.expander("Add new cross-reactivity sequences"):
        add_cross_col1, add_cross_col2 = st.columns(2)

        with add_cross_col1:
            st.markdown("**Upload FASTA file**")
            uploaded_cross = st.file_uploader(
                "Upload FASTA",
                type=["fa", "fasta", "fna"],
                accept_multiple_files=True,
                key="cross_fasta_upload",
                label_visibility="collapsed",
            )
            if uploaded_cross:
                for f in uploaded_cross:
                    p = Path(f.name)
                    if p.suffix.lower() not in (".fa", ".fasta", ".fna"):
                        st.error(f"Unsupported extension: {p.suffix}")
                        continue
                    dest = FASTA_DIR / (p.stem + ".fa")
                    dest.write_bytes(f.getvalue())
                    for ext in (".fasta", ".fna"):
                        old = FASTA_DIR / (p.stem + ext)
                        if old.exists():
                            old.unlink()
                    st.success(f"Saved {dest.name}")
                st.rerun()

        with add_cross_col2:
            st.markdown("**Fetch from NCBI Virus (using gget)**")
            _render_fetch_ui("cross")

    st.divider()

    # --- Host sequences ---
    st.subheader("Host sequences")
    st.caption("Host genome to check for non-specific amplification (e.g., human genome)")

    if offtarget_options:
        host_labels = []
        for stem in offtarget_options:
            info = _fasta_info(stem)
            host_labels.append(f"{stem}  ({info['seqs']} seqs, {_format_size(info['size'])})")

        selected_host = st.multiselect(
            "Select host sequences",
            options=host_labels,
            default=[],
            key="host_select",
            placeholder="Choose host genome...",
        )
        st.session_state.host = [label.split("  (")[0] for label in selected_host]
    else:
        st.info("No FASTA files available.")
        st.session_state.host = []

    # Upload host sequences
    with st.expander("Upload host sequences"):
        uploaded_host = st.file_uploader(
            "Upload FASTA",
            type=["fa", "fasta", "fna"],
            accept_multiple_files=True,
            key="host_fasta_upload",
            label_visibility="collapsed",
        )
        if uploaded_host:
            for f in uploaded_host:
                p = Path(f.name)
                if p.suffix.lower() not in (".fa", ".fasta", ".fna"):
                    st.error(f"Unsupported extension: {p.suffix}")
                    continue
                dest = FASTA_DIR / (p.stem + ".fa")
                dest.write_bytes(f.getvalue())
                for ext in (".fasta", ".fna"):
                    old = FASTA_DIR / (p.stem + ext)
                    if old.exists():
                        old.unlink()
                st.success(f"Saved {dest.name}")
            st.rerun()

    st.divider()

    # --- Continue / Skip ---
    if st.button("Continue →", type="primary", use_container_width=True):
        _navigate(f"{workflow}_config")
        st.rerun()

    st.caption("You can skip this step — cross-reactivity and host checks are optional.")


def _config_design():
    """Design-specific configuration."""
    st.header("Design — Configuration")

    # Summary of selections from previous steps
    targets = st.session_state.get("targets", [])
    cross = st.session_state.get("cross", [])
    host = st.session_state.get("host", [])
    st.markdown(f"**Target:** {', '.join(targets) if targets else 'none'}")
    if cross:
        st.markdown(f"**Cross-reactivity:** {', '.join(cross)}")
    if host:
        st.markdown(f"**Host:** {', '.join(host)}")

    st.divider()

    # Design mode
    st.subheader("Design mode")
    col_single, col_multi = st.columns(2)
    with col_single:
        singleplex = st.button(
            "Singleplex",
            use_container_width=True,
            type="primary" if st.session_state.get("design_mode", "Singleplex") == "Singleplex" else "secondary",
            key="btn_singleplex",
        )
        if singleplex:
            st.session_state.design_mode = "Singleplex"
            st.rerun()
        st.caption("Design primers for a single target pathogen.")
    with col_multi:
        st.button(
            "Multiplex (coming soon)",
            use_container_width=True,
            disabled=True,
            key="btn_multiplex",
        )
        st.caption("Design a multiplexed panel for multiple targets simultaneously.")

    # Set internal mode for pipeline
    if st.session_state.get("design_mode", "Singleplex") == "Singleplex":
        st.session_state.mode = "Singleplex"

    st.divider()

    # Probe design
    st.checkbox("Enable probe design", key="probe_enabled")

    # Parameters
    _tab_parameters()

    st.divider()

    if st.button("Continue →", type="primary", use_container_width=True):
        st.session_state._run_page_fresh = True
        _navigate("run")
        st.rerun()


def _config_evaluate():
    """Evaluate-specific configuration."""
    st.header("Evaluate — Configuration")

    # Summary of selections from previous steps
    targets = st.session_state.get("targets", [])
    cross = st.session_state.get("cross", [])
    host = st.session_state.get("host", [])
    st.markdown(f"**Target:** {', '.join(targets) if targets else 'none'}")
    if cross:
        st.markdown(f"**Cross-reactivity:** {', '.join(cross)}")
    if host:
        st.markdown(f"**Host:** {', '.join(host)}")

    st.divider()

    # Set internal mode
    st.session_state.mode = "Evaluate"

    # Primer input
    st.subheader("Primer sequences to evaluate")

    eval_method = st.radio(
        "Input method",
        ["Paste sequences", "Upload primer FASTA"],
        key="eval_method",
        horizontal=True,
    )

    if eval_method == "Paste sequences":
        st.text_input("Forward primer sequence (5'→3')", key="eval_for",
                       placeholder="e.g., ATGCGATCGATCGATCG")
        st.text_input("Reverse primer sequence (5'→3')", key="eval_rev",
                       placeholder="e.g., TAGCTAGCTAGCTAGCT")
    else:
        pset_upload = st.file_uploader(
            "Upload primer set FASTA",
            type=["fa", "fasta"],
            key="eval_pset_upload",
        )
        if pset_upload:
            pset_dir = PROJECT_ROOT / "evaluate"
            pset_dir.mkdir(parents=True, exist_ok=True)
            pset_path = pset_dir / pset_upload.name
            pset_path.write_bytes(pset_upload.getvalue())
            st.session_state["eval_pset_path"] = str(pset_path)
            st.success(f"Saved {pset_upload.name}")

    st.divider()

    # Parameters
    _tab_parameters()

    st.divider()

    if st.button("Continue →", type="primary", use_container_width=True):
        st.session_state._run_page_fresh = True
        _navigate("run")
        st.rerun()


def _config_monitor():
    """Monitor-specific configuration."""
    st.header("Monitor — Configuration")

    # Summary of selections from previous steps
    targets = st.session_state.get("targets", [])
    cross = st.session_state.get("cross", [])
    host = st.session_state.get("host", [])
    st.markdown(f"**Target:** {', '.join(targets) if targets else 'none'}")
    if cross:
        st.markdown(f"**Cross-reactivity:** {', '.join(cross)}")
    if host:
        st.markdown(f"**Host:** {', '.join(host)}")

    st.divider()

    # Set internal mode
    st.session_state.mode = "Evaluate"

    # Primer input
    st.subheader("Primer sequences to monitor")

    eval_method = st.radio(
        "Input method",
        ["Paste sequences", "Upload primer FASTA"],
        key="eval_method",
        horizontal=True,
    )

    if eval_method == "Paste sequences":
        st.text_input("Forward primer sequence (5'→3')", key="eval_for",
                       placeholder="e.g., ATGCGATCGATCGATCG")
        st.text_input("Reverse primer sequence (5'→3')", key="eval_rev",
                       placeholder="e.g., TAGCTAGCTAGCTAGCT")
    else:
        pset_upload = st.file_uploader(
            "Upload primer set FASTA",
            type=["fa", "fasta"],
            key="eval_pset_upload",
        )
        if pset_upload:
            pset_dir = PROJECT_ROOT / "evaluate"
            pset_dir.mkdir(parents=True, exist_ok=True)
            pset_path = pset_dir / pset_upload.name
            pset_path.write_bytes(pset_upload.getvalue())
            st.session_state["eval_pset_path"] = str(pset_path)
            st.success(f"Saved {pset_upload.name}")

    st.divider()

    # Parameters
    _tab_parameters()

    st.divider()

    if st.button("Continue →", type="primary", use_container_width=True):
        st.session_state._run_page_fresh = True
        _navigate("run")
        st.rerun()


def _page_design():
    """Design workflow — config page."""
    _render_workflow_progress("config")
    _config_design()


def _page_evaluate():
    """Evaluate workflow — config page."""
    _render_workflow_progress("config")
    _config_evaluate()


def _page_monitor():
    """Monitor workflow — config page."""
    _render_workflow_progress("config")
    _config_monitor()


def _page_run():
    """Run pipeline page."""
    _render_workflow_progress("run")
    _tab_run()

    # Show results button after pipeline completion
    rc = st.session_state.get("pipeline_return_code")
    if rc is not None and not st.session_state.get("pipeline_running"):
        st.divider()
        if st.button("View Results →", type="primary", use_container_width=True):
            _navigate("results")
            st.rerun()


def _page_results():
    """Results page."""
    _render_workflow_progress("results")
    _tab_results()


def main():
    _render_sidebar()

    page = _current_page()

    if page == "home":
        _page_home()
    elif page == "select_target":
        _page_select_target()
    elif page == "select_offtarget":
        _page_select_offtarget()
    elif page in ("design", "design_config"):
        _page_design()
    elif page in ("evaluate", "evaluate_config"):
        _page_evaluate()
    elif page in ("monitor", "monitor_config"):
        _page_monitor()
    elif page == "run":
        _page_run()
    elif page == "results":
        _page_results()
    else:
        _page_home()


if __name__ == "__main__":
    main()
