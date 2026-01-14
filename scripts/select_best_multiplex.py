#!/usr/bin/env python3
"""
select_best_multiplex.py

Reads multiple per-target design CSVs, selects the top-performing rows per target,
and writes a combined CSV with an added 'target' column.

Selection rule (per input CSV / target):
  1) Maximize sco_target
  2) Tie-break: minimize worst off-target score = max(sco_<X>) across off-target columns
  3) Tie-break: minimize sum of off-target scores = sum(sco_<X>)

By default selects 3 row per target (use --top-n to select more).
"""

from __future__ import annotations

import argparse
import os
import sys
from typing import List, Tuple

import pandas as pd


def infer_target_name(path: str) -> str:
    """Infer target name from filename like final/A.csv -> 'A'."""
    base = os.path.basename(path)
    # strip compression if any
    for ext in [".csv.gz", ".tsv.gz", ".csv", ".tsv"]:
        if base.endswith(ext):
            base = base[: -len(ext)]
            break
    return base


def find_off_score_columns(df: pd.DataFrame) -> List[str]:
    """Return sco_* columns except sco_target."""
    cols = [c for c in df.columns if c.startswith("sco_")]
    off = [c for c in cols if c != "sco_target"]
    return off


def select_top_rows(df: pd.DataFrame, top_n: int) -> pd.DataFrame:
    if "sco_target" not in df.columns:
        raise ValueError("Input CSV is missing required column 'sco_target'.")

    off_cols = find_off_score_columns(df)

    # Ensure numeric scoring columns where possible
    df = df.copy()
    df["sco_target"] = pd.to_numeric(df["sco_target"], errors="coerce")

    if off_cols:
        for c in off_cols:
            df[c] = pd.to_numeric(df[c], errors="coerce")
        df["_off_worst"] = df[off_cols].max(axis=1, skipna=True)
        df["_off_sum"] = df[off_cols].sum(axis=1, skipna=True)
    else:
        # No off-target scores provided; treat as no penalty.
        df["_off_worst"] = 0.0
        df["_off_sum"] = 0.0

    # Drop rows with missing sco_target (cannot rank them)
    df = df.dropna(subset=["sco_target"])

    if df.empty:
        raise ValueError("No valid rows after parsing 'sco_target' as numeric.")

    # Sort: sco_target desc, off_worst asc, off_sum asc
    df = df.sort_values(
        by=["sco_target", "_off_worst", "_off_sum"],
        ascending=[False, True, True],
        kind="mergesort",  # stable sort
    )

    # Select top N
    out = df.head(top_n).drop(columns=["_off_worst", "_off_sum"])
    return out


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(
        description="Select top multiplex primer/probe sets per target and concatenate."
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
        default=3,
        help="Number of top rows to keep per target (default: 3).",
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
    return p.parse_args()


def main() -> int:
    args = parse_args()

    if args.top_n < 1:
        raise ValueError("--top-n must be >= 1")

    selected_frames: List[pd.DataFrame] = []

    for path in args.csvs:
        df = pd.read_csv(path)

        if args.target_from == "column":
            if args.target_column not in df.columns:
                raise ValueError(
                    f"{path}: --target-from=column but missing column '{args.target_column}'."
                )
            # If multiple targets exist in one file, split and select per target value.
            for tgt, sub in df.groupby(args.target_column, dropna=False):
                sub_sel = select_top_rows(sub, args.top_n)
                sub_sel = sub_sel.copy()
                sub_sel["target"] = str(tgt)
                selected_frames.append(sub_sel)
        else:
            tgt = infer_target_name(path)
            sel = select_top_rows(df, args.top_n)
            sel = sel.copy()
            sel["target"] = tgt
            selected_frames.append(sel)

    combined = pd.concat(selected_frames, ignore_index=True)

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
