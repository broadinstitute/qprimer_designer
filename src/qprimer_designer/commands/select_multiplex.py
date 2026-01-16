"""Select best multiplex primer sets."""

import argparse
import os
import sys
from typing import List

import pandas as pd


def register(subparsers):
    """Register the select-multiplex subcommand."""
    parser = subparsers.add_parser(
        "select-multiplex",
        help="Select top multiplex primer sets per target",
        description="""
Read multiple per-target design CSVs, select the top-performing rows per target,
and write a combined CSV with an added 'target' column.

Selection rule (per input CSV / target):
  1) Maximize sco_target
  2) Tie-break: minimize worst off-target score = max(sco_<X>) across off-target columns
  3) Tie-break: minimize sum of off-target scores = sum(sco_<X>)
""",
    )
    parser.add_argument("csvs", nargs="+", help="Input per-target CSV files")
    parser.add_argument("--out", required=True, help="Output combined CSV path")
    parser.add_argument("--top-n", type=int, default=3, help="Number of top rows to keep per target (default: 3)")
    parser.add_argument(
        "--target-from",
        choices=["filename", "column"],
        default="filename",
        help="How to determine the target label (default: filename)",
    )
    parser.add_argument(
        "--target-column",
        default="target",
        help="If --target-from=column, read target label from this column (default: target)",
    )
    parser.set_defaults(func=run)


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
    """Select top N rows based on scoring criteria."""
    if "sco_target" not in df.columns:
        raise ValueError("Input CSV is missing required column 'sco_target'.")

    off_cols = find_off_score_columns(df)

    df = df.copy()
    df["sco_target"] = pd.to_numeric(df["sco_target"], errors="coerce")

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

    # Sort: sco_target desc, off_worst asc, off_sum asc
    df = df.sort_values(
        by=["sco_target", "_off_worst", "_off_sum"],
        ascending=[False, True, True],
        kind="mergesort",
    )

    out = df.head(top_n).drop(columns=["_off_worst", "_off_sum"])
    return out


def run(args):
    """Run the select-multiplex command."""
    if args.top_n < 1:
        raise ValueError("--top-n must be >= 1")

    selected_frames: List[pd.DataFrame] = []

    for path in args.csvs:
        df = pd.read_csv(path)

        if args.target_from == "column":
            if args.target_column not in df.columns:
                raise ValueError(f"{path}: --target-from=column but missing column '{args.target_column}'.")
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
    print(f"Wrote {len(combined)} rows to {args.out}")
