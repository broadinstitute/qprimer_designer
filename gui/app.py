"""Streamlit GUI for the qPrimer Designer pipeline."""

import os
import select
import shutil
import signal
import subprocess
import sys
import time
from pathlib import Path

import streamlit as st

# ---------------------------------------------------------------------------
# Resolve project root (parent of gui/)
# ---------------------------------------------------------------------------
PROJECT_ROOT = Path(__file__).resolve().parent.parent
FASTA_DIR = PROJECT_ROOT / "target_seqs" / "original"
FINAL_DIR = PROJECT_ROOT / "final"
EVALUATE_DIR = PROJECT_ROOT / "evaluate"

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
# Helpers
# ---------------------------------------------------------------------------

def _available_fasta() -> list[str]:
    """Return sorted list of FASTA stem names in target_seqs/original/."""
    if not FASTA_DIR.exists():
        return []
    return sorted(
        p.stem for p in FASTA_DIR.iterdir()
        if p.suffix in (".fa", ".fasta", ".fna") and p.is_file()
    )


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


def _check_tool(name: str) -> bool:
    return shutil.which(name) is not None


# ---------------------------------------------------------------------------
# Sidebar
# ---------------------------------------------------------------------------

def _render_sidebar():
    with st.sidebar:
        st.title("qPrimer Designer")
        st.caption(f"Working directory: `{PROJECT_ROOT}`")

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

        st.divider()

        # Pipeline status indicator
        if st.session_state.get("pipeline_running"):
            st.markdown(":orange[Pipeline running...]")
        elif st.session_state.get("pipeline_return_code") is not None:
            rc = st.session_state.pipeline_return_code
            if rc == 0:
                st.markdown(":green[Pipeline finished successfully]")
            else:
                st.markdown(f":red[Pipeline failed (exit code {rc})]")


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
            # Normalize extension to .fa so the Snakefile can find it
            dest = FASTA_DIR / (Path(f.name).stem + ".fa")
            dest.write_bytes(f.getvalue())
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
    st.header("Pipeline Parameters")

    _init_params()
    p = st.session_state.params

    # Load from existing file
    col_load, col_reset = st.columns(2)
    with col_load:
        if st.button("Load from existing params.txt"):
            params_path = PROJECT_ROOT / "params.txt"
            if params_path.exists():
                loaded = parse_params(str(params_path))
                for k, v in loaded.items():
                    p[k] = v
                st.success("Loaded params.txt")
            else:
                st.warning("No params.txt found in project root")
    with col_reset:
        if st.button("Reset to defaults"):
            st.session_state.params = dict(DEFAULT_PARAMS)
            st.rerun()

    st.divider()

    # --- Representative sequence selection ---
    with st.expander("Representative sequence selection", expanded=False):
        p["DESIGN_WINDOW"] = st.number_input(
            "DESIGN_WINDOW", value=int(p["DESIGN_WINDOW"]),
            min_value=50, step=50, key="p_design_window",
        )

    # --- Primer generation ---
    with st.expander("Primer generation", expanded=True):
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

    cmd = ["snakemake", "-s", "Snakefile", "--cores", str(cores)]

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

    snakefile_content = build_snakefile(
        targets=targets,
        cross=cross,
        host=host,
        panel=panel,
    )
    (PROJECT_ROOT / "Snakefile").write_text(snakefile_content)

    _init_params()
    params_content = build_params_txt(st.session_state.params)
    (PROJECT_ROOT / "params.txt").write_text(params_content)


def _tab_run():
    st.header("Run Pipeline")

    # Controls
    c1, c2, c3 = st.columns(3)
    max_cpu = os.cpu_count() or 1
    c1.slider("CPU cores", min_value=1, max_value=max_cpu, value=min(4, max_cpu), key="cores")
    c2.checkbox("Dry run (plan only, no execution)", key="dry_run")
    c3.checkbox("Show live output", value=True, key="show_live_output")

    # Command preview
    cmd = _build_command()
    st.code(" ".join(cmd), language="bash")

    # Preflight
    errors = _preflight_checks()
    if errors:
        for e in errors:
            st.error(e)

    # Run / Stop buttons
    col_run, col_stop = st.columns(2)

    running = st.session_state.get("pipeline_running", False)

    with col_run:
        if st.button("Run", disabled=running or bool(errors), type="primary"):
            _write_pipeline_files()
            st.session_state.pipeline_log = ""
            st.session_state.pipeline_return_code = None
            st.session_state.pipeline_running = True

            env = os.environ.copy()
            env["PYTHONUNBUFFERED"] = "1"

            proc = subprocess.Popen(
                cmd,
                stdout=subprocess.PIPE,
                stderr=subprocess.STDOUT,
                cwd=str(PROJECT_ROOT),
                env=env,
                text=True,
            )
            st.session_state.pipeline_pid = proc.pid
            st.session_state._pipeline_proc = proc
            st.rerun()

    with col_stop:
        if st.button("Stop", disabled=not running):
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

    # Log output area
    st.subheader("Pipeline output")
    log_area = st.empty()
    show_live = st.session_state.get("show_live_output", True)

    if running and hasattr(st.session_state, "_pipeline_proc"):
        proc = st.session_state._pipeline_proc
        log = st.session_state.get("pipeline_log", "")

        # Non-blocking read of available output
        fd = proc.stdout.fileno()
        ready, _, _ = select.select([fd], [], [], 0.1)
        if ready:
            chunk = os.read(fd, 8192).decode("utf-8", errors="replace")
            if chunk:
                log += chunk
                st.session_state.pipeline_log = log

        # Check if process finished
        rc = proc.poll()
        if rc is not None:
            # Read any remaining output
            remaining = proc.stdout.read()
            if remaining:
                log += remaining
                st.session_state.pipeline_log = log
            st.session_state.pipeline_running = False
            st.session_state.pipeline_return_code = rc
            st.rerun()

        if show_live:
            log_area.code(log if log else "(waiting for output...)", language="text")
        else:
            log_area.info("Pipeline is running. Enable **Show live output** to see logs.")

        # Auto-refresh while running
        if st.session_state.get("pipeline_running"):
            time.sleep(1)
            st.rerun()
    else:
        log = st.session_state.get("pipeline_log", "")
        if log:
            log_area.code(log, language="text")
        else:
            log_area.info("Click Run to start the pipeline.")

    # Show inline results after successful completion
    rc = st.session_state.get("pipeline_return_code")
    if rc is not None and not st.session_state.get("pipeline_running"):
        st.divider()
        if rc == 0:
            st.success("Pipeline finished successfully!")
            mode = st.session_state.get("mode", "Singleplex")
            if mode == "Evaluate":
                xlsx_files = sorted(EVALUATE_DIR.glob("**/*.xlsx")) if EVALUATE_DIR.exists() else []
                if xlsx_files:
                    import pandas as pd

                    st.subheader("Evaluation results preview")
                    for xf in xlsx_files:
                        try:
                            summary_df = pd.read_excel(xf, sheet_name=0)
                            st.write(f"**{xf.name}**")
                            st.dataframe(summary_df, use_container_width=True)
                        except Exception as exc:
                            st.warning(f"Could not read {xf.name}: {exc}")
                    st.caption("See the **Results** tab for full details and downloads.")
            elif FINAL_DIR.exists():
                csvs = sorted(FINAL_DIR.glob("*.csv"))
                if csvs:
                    import pandas as pd

                    st.subheader("Results preview")
                    for csv_path in csvs:
                        try:
                            df = pd.read_csv(csv_path)
                            st.write(f"**{csv_path.name}** — {len(df)} primer pairs")
                            st.dataframe(df.head(10), use_container_width=True)
                        except Exception as exc:
                            st.warning(f"Could not read {csv_path.name}: {exc}")
                    st.caption("See the **Results** tab for full details and downloads.")
        else:
            st.error(f"Pipeline failed with exit code {rc}. Check the log above.")


# ---------------------------------------------------------------------------
# Tab 5: Results
# ---------------------------------------------------------------------------

def _tab_results():
    st.header("Results")

    # --- Final CSV files ---
    st.subheader("Final output files")

    if FINAL_DIR.exists():
        csvs = sorted(FINAL_DIR.glob("*.csv"))
    else:
        csvs = []

    if csvs:
        selected_csv = st.selectbox(
            "Select CSV to view",
            options=[c.name for c in csvs],
            key="result_csv",
        )
        if selected_csv:
            import pandas as pd

            csv_path = FINAL_DIR / selected_csv
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
                    file_name=selected_csv,
                    mime="text/csv",
                )
            except Exception as exc:
                st.error(f"Error reading {selected_csv}: {exc}")
    else:
        st.info("No CSV results in `final/` yet. Run the pipeline first.")

    st.divider()

    # --- Evaluate reports ---
    st.subheader("Evaluation reports")

    if EVALUATE_DIR.exists():
        xlsx_files = sorted(EVALUATE_DIR.glob("**/*.xlsx"))
    else:
        xlsx_files = []

    if xlsx_files:
        import pandas as pd

        selected_xlsx = st.selectbox(
            "Select report to view",
            options=[xf.name for xf in xlsx_files],
            key="result_xlsx",
        )
        if selected_xlsx:
            # Find full path (name may not be unique across subdirs)
            xlsx_path = next(xf for xf in xlsx_files if xf.name == selected_xlsx)
            try:
                sheet_names = pd.ExcelFile(xlsx_path).sheet_names

                for sheet in sheet_names:
                    df = pd.read_excel(xlsx_path, sheet_name=sheet)
                    st.write(f"**{sheet}** ({len(df)} rows)")
                    st.dataframe(df, use_container_width=True)

                st.download_button(
                    "Download XLSX",
                    data=xlsx_path.read_bytes(),
                    file_name=selected_xlsx,
                    mime="application/vnd.openxmlformats-officedocument.spreadsheetml.sheet",
                    key=f"dl_{selected_xlsx}",
                )
            except Exception as exc:
                st.error(f"Error reading {selected_xlsx}: {exc}")
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

def main():
    _render_sidebar()

    tab_files, tab_config, tab_params, tab_run, tab_results = st.tabs(
        ["Files", "Configuration", "Parameters", "Run", "Results"]
    )

    with tab_files:
        _tab_files()
    with tab_config:
        _tab_configuration()
    with tab_params:
        _tab_parameters()
    with tab_run:
        _tab_run()
    with tab_results:
        _tab_results()


if __name__ == "__main__":
    main()
