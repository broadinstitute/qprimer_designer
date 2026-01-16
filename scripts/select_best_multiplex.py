#!/usr/bin/env python3
"""
select_best_multiplex.py

Reads multiple per-target design CSVs, selects the top-performing rows per target,
then computes an additional score `total_deltaG` based on pairwise ViennaRNA RNAduplex
ΔG interactions among primers in the *selected multiplex pool*.

Multiplex pool definition:
  The set of selected primer pairs (rows) across all targets in PANEL.
  Each row contributes two primers: pseq_f and pseq_r.

Selection rule (per input CSV / target):
  1) Maximize sco_target
  2) Tie-break: minimize worst off-target score = max(sco_<X>) across off-target columns
  3) Tie-break: minimize sum of off-target scores = sum(sco_<X>)

`total_deltaG` computation (per selected row i):
  Let Pi = {pseq_f_i, pseq_r_i}
  Let Ppool = all primers (pseq_f and pseq_r) across all selected rows
  total_deltaG_i = sum( ΔG(p, q) for p in Pi for q in (Ppool / Pi) )

Notes:
- More negative ΔG implies stronger (worse) dimerization tendency.
- We exclude within-row interactions (i.e., pseq_f_i vs pseq_r_i) from total_deltaG
  because that is typically represented by your existing Dimer_dG column. If you want
  to include within-row too, see INCLUDE_SELF_PAIR below.

Usage:
  python select_best_multiplex.py final/A.csv final/B.csv ... --out final.csv --RNAduplex /path/to/RNAduplex
"""

from __future__ import annotations

import argparse
import os
import re
import subprocess
import sys
from typing import Dict, List, Tuple

import pandas as pd


# Include self-pairing of primer delta G in calculation
INCLUDE_SELF_PAIR = True


def get_deltaG_vienna(seq1: str, seq2: str, RNAduplex: str) -> float:
    """Compute primer–primer dimer ΔG using ViennaRNA RNAduplex (DNA mode)."""
    inp = f"{seq1}\n{seq2}\n"
    res = subprocess.run(
        [RNAduplex, "--noconv", "--paramFile=DNA"],
        input=inp,
        capture_output=True,
        text=True,
        check=True,
    )
    deltaG = re.search(r"\(\s*([-+]?\d*\.?\d+)\s*\)", res.stdout)
    if not deltaG:
        raise ValueError(f"Could not parse ΔG from RNAduplex output:\n{res.stdout}")
    return float(deltaG.group(1))


def infer_target_name(path: str) -> str:
    """Infer target name from filename like final/A.csv -> 'A'."""
    base = os.path.basename(path)
    for ext in [".csv.gz", ".tsv.gz", ".csv", ".tsv"]:
        if base.endswith(ext):
            base = base[: -len(ext)]
            break
    return base


def find_off_score_columns(df: pd.DataFrame) -> List[str]:
    """Return sco_* columns except sco_target."""
    cols = [c for c in df.columns if c.startswith("sco_")]
    return [c for c in cols if c != "sco_target"]


def select_top_rows(df: pd.DataFrame, top_n: int) -> pd.DataFrame:
    if "sco_target" not in df.columns:
        raise ValueError("Input CSV is missing required column 'sco_target'.")

    df = df.copy()
    df["sco_target"] = pd.to_numeric(df["sco_target"], errors="coerce")

    off_cols = find_off_score_columns(df)
    if off_cols:
        for c in off_cols:
            df[c] = pd.to_numeric(df[c], errors="coerce")
        df["_off_worst"] = df[off_cols].max(axis=1, skipna=True)
        df["_off_sum"] = df[off_cols].sum(axis=1, skipna=True)
    else:
        df["_off_worst"] = 0.0
        df["_off_sum"] = 0.0

    df = df.dropna(subset=["sco_target"])
    if df.empty:
        raise ValueError("No valid rows after parsing 'sco_target' as numeric.")

    df = df.sort_values(
        by=["sco_target", "_off_worst", "_off_sum"],
        ascending=[False, True, True],
        kind="mergesort",
    )

    out = df.head(top_n).drop(columns=["_off_worst", "_off_sum"])
    return out


def compute_total_deltaG_for_pool(
    combined: pd.DataFrame, RNAduplex: str, f_col: str = "pseq_f", r_col: str = "pseq_r"
) -> pd.Series:
    """
    Given a DataFrame of selected rows (the multiplex pool),
    compute total_deltaG per row as described in the module docstring.
    """
    if f_col not in combined.columns or r_col not in combined.columns:
        raise ValueError(f"Missing required primer sequence columns '{f_col}' and/or '{r_col}'.")

    # Normalize sequences (strip whitespace, uppercase)
    f_seqs = combined[f_col].astype(str).str.strip().str.upper().tolist()
    r_seqs = combined[r_col].astype(str).str.strip().str.upper().tolist()

    # Build primer list per row and a flat pool
    row_primers: List[Tuple[str, str]] = list(zip(f_seqs, r_seqs))
    pool_primers: List[str] = []
    for f, r in row_primers:
        pool_primers.extend([f, r])

    # Cache ΔG computations because many pairs repeat (and ΔG is symmetric for our purposes)
    cache: Dict[Tuple[str, str], float] = {}

    def dg(a: str, b: str) -> float:
        key = (a, b) if a <= b else (b, a)
        if key not in cache:
            cache[key] = get_deltaG_vienna(key[0], key[1], RNAduplex)
        return cache[key]

    totals: List[float] = []
    for i, (fi, ri) in enumerate(row_primers):
        # Primers belonging to this row
        my = [fi, ri]

        # Build list of "other" primers in pool
        if INCLUDE_SELF_PAIR:
            others = pool_primers[:]  # includes self; we'll handle exclusions below
        else:
            # Exclude both primers from this row (two occurrences) from the pool list
            # Do it by index: primers are [f0,r0,f1,r1,...]
            others = []
            for j, p in enumerate(pool_primers):
                if j in (2 * i, 2 * i + 1):
                    continue
                others.append(p)

        total = 0.0

        # Pairwise interactions of this pair's primers against all other primers
        for p in my:
            for q in others:
                # If INCLUDE_SELF_PAIR is True, prevent p vs itself (identical primer) from being counted
                if not INCLUDE_SELF_PAIR:
                    pass
                else:
                    if p == q:
                        continue
                total += dg(p, q)

        # Optionally include within-row forward-vs-reverse interaction
        if INCLUDE_SELF_PAIR:
            total += dg(fi, ri)

        totals.append(total)

    return pd.Series(totals, index=combined.index, name="total_deltaG")


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(
        description="Select top multiplex primer/probe sets per target, compute total_deltaG across pool, and concatenate."
    )
    p.add_argument(
        "csvs",
        nargs="+",
        help="Input per-target CSV files (e.g., final/A.csv final/B.csv ...).",
    )
    p.add_argument(
        "--out",
        required=True,
        help="Output combined CSV path (e.g., final.csv).",
    )
    p.add_argument(
        "--top-n",
        type=int,
        default=1,
        help="Number of top rows to keep per target (default: 1).",
    )
    p.add_argument(
        "--program",
        dest="program_path",
        required=True,
        help="Path to ViennaRNA RNAduplex executable (DNA mode is used via --paramFile=DNA).",
    )
    p.add_argument(
        "--target-from",
        choices=["filename", "column"],
        default="filename",
        help="How to determine the target label (default: filename).",
    )
    p.add_argument(
        "--target-column",
        default="target",
        help="If --target-from=column, read target label from this column (default: target).",
    )
    p.add_argument(
        "--f-col",
        default="pseq_f",
        help="Column name for forward primer sequence (default: pseq_f).",
    )
    p.add_argument(
        "--r-col",
        default="pseq_r",
        help="Column name for reverse primer sequence (default: pseq_r).",
    )
    return p.parse_args()


def main() -> int:
    args = parse_args()

    viennaDir = args.program_path

    if args.top_n < 1:
        raise ValueError("--top-n must be >= 1")

    selected_frames: List[pd.DataFrame] = []

    RNAduplex = f"{viennaDir}/src/bin/RNAduplex"
    if not os.path.exists(RNAduplex):
        raise FileNotFoundError(
            f"RNAduplex not found at: {RNAduplex}\n"
            "Check VIENNA_DIRECTORY in params.txt."
        )

    for path in args.csvs:
        df = pd.read_csv(path)

        if args.target_from == "column":
            if args.target_column not in df.columns:
                raise ValueError(
                    f"{path}: --target-from=column but missing column '{args.target_column}'."
                )
            for tgt, sub in df.groupby(args.target_column, dropna=False):
                sub_sel = select_top_rows(sub, args.top_n).copy()
                sub_sel["target"] = str(tgt)
                selected_frames.append(sub_sel)
        else:
            tgt = infer_target_name(path)
            sel = select_top_rows(df, args.top_n).copy()
            sel["target"] = tgt
            selected_frames.append(sel)

    combined = pd.concat(selected_frames, ignore_index=True)

    # Compute total_deltaG across the multiplex pool of selected primers
    combined["total_deltaG"] = compute_total_deltaG_for_pool(
        combined, RNAduplex=RNAduplex, f_col=args.f_col, r_col=args.r_col
    )

    # Put 'target' first for readability
    cols = combined.columns.tolist()
    if "target" in cols:
        cols = ["target"] + [c for c in cols if c != "target"]
        combined = combined[cols]

    combined.to_csv(args.out, index=False)
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except Exception as e:
        print(f"ERROR: {e}", file=sys.stderr)
        raise
