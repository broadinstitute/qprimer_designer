"""Forecast future viral variant frequencies and select representative sequences.

This command forecasts which viral lineages are likely to be most prevalent
in N days and writes a FASTA of representative sequences from those lineages.
The output FASTA can be used as the reference for primer evaluation, producing
a "future robustness" score for each primer pair.
"""

from __future__ import annotations

import random
import sys
from pathlib import Path

import pandas as pd
from Bio import SeqIO

from qprimer_designer.external.evofr import (
    _CLADE_COLUMNS,
    _ACCESSION_COLUMNS,
    _find_column,
    metadata_to_clade_counts,
    run_mlr_forecast,
)
from qprimer_designer.utils.params import parse_params


def _load_fasta_by_accession(fasta_path: Path) -> dict[str, SeqIO.SeqRecord]:
    """Load FASTA records indexed by accession (first whitespace-delimited token)."""
    records: dict[str, SeqIO.SeqRecord] = {}
    for rec in SeqIO.parse(str(fasta_path), "fasta"):
        acc = rec.id.split()[0]
        records[acc] = rec
    return records


def _sample_sequences(
    accessions: list[str],
    records: dict[str, SeqIO.SeqRecord],
    n: int,
    rng: random.Random,
) -> list[SeqIO.SeqRecord]:
    """Sample up to n records whose accessions appear in records."""
    available = [a for a in accessions if a in records]
    chosen = rng.sample(available, min(n, len(available)))
    return [records[a] for a in chosen]


def run(args) -> None:
    """Execute the forecast-variants command."""
    sequences_path = Path(args.sequences)
    metadata_path = Path(args.metadata)
    out_sequences = Path(args.out_sequences)
    out_forecast = Path(args.out_forecast)
    horizon_days = int(args.horizon)
    min_sequences = int(args.min_sequences)
    top_n = int(args.top_n)
    top_n_per_clade = int(getattr(args, "top_n_per_clade", 20))
    location = getattr(args, "location", None) or None

    out_sequences.parent.mkdir(parents=True, exist_ok=True)
    out_forecast.parent.mkdir(parents=True, exist_ok=True)

    # --- Load metadata ---
    print(f"Loading metadata from {metadata_path}...")
    try:
        sep = "\t" if str(metadata_path).endswith(".tsv") else ","
        metadata_df = pd.read_csv(metadata_path, sep=sep, low_memory=False)
    except Exception as e:
        print(f"Error reading metadata: {e}", file=sys.stderr)
        sys.exit(1)
    print(f"  {len(metadata_df)} rows in metadata")

    # --- Convert to clade counts ---
    try:
        clade_counts = metadata_to_clade_counts(metadata_df, location_filter=location)
    except ValueError as e:
        print(f"Error processing metadata: {e}", file=sys.stderr)
        sys.exit(1)

    n_clades = clade_counts["clade"].nunique()
    n_total = int(clade_counts["count"].sum())
    date_min = clade_counts["date"].min()
    date_max = clade_counts["date"].max()
    print(f"  {n_clades} clades, {n_total} sequences, {date_min} – {date_max}")

    # --- Run MLR forecast ---
    print(f"\nFitting MLR model (horizon = {horizon_days} days)...")
    try:
        forecast_df = run_mlr_forecast(
            clade_counts,
            horizon_days=horizon_days,
            min_sequences=min_sequences,
        )
    except ValueError as e:
        print(f"Error: {e}", file=sys.stderr)
        sys.exit(1)

    print(f"\nVariant frequency forecast (+{horizon_days} days):")
    print(f"  {'Clade':<30} {'Current%':>9} {'Future%':>9} {'Growth':>9}")
    print(f"  {'-'*30} {'-'*9} {'-'*9} {'-'*9}")
    for _, row in forecast_df.head(min(top_n * 2, 15)).iterrows():
        print(
            f"  {row['clade']:<30} "
            f"{row['current_freq']*100:>8.1f}% "
            f"{row['predicted_freq']*100:>8.1f}% "
            f"{row['growth_advantage']:>8.2f}x"
        )

    # Save forecast table
    forecast_df.to_csv(str(out_forecast), sep="\t", index=False)
    print(f"\nForecast table saved: {out_forecast}")

    # --- Select top N clades ---
    top_clades = forecast_df.head(top_n)["clade"].tolist()
    print(f"\nTop {top_n} clades for primer testing: {', '.join(top_clades)}")

    # --- Load sequences ---
    print(f"\nLoading sequences from {sequences_path}...")
    if not sequences_path.exists():
        print(f"Error: sequences file not found: {sequences_path}", file=sys.stderr)
        sys.exit(1)
    try:
        all_records = _load_fasta_by_accession(sequences_path)
    except Exception as e:
        print(f"Error reading sequences: {e}", file=sys.stderr)
        sys.exit(1)
    print(f"  {len(all_records)} sequences loaded")

    # --- Build accession → clade mapping ---
    acc_col = _find_column(metadata_df, _ACCESSION_COLUMNS)
    if acc_col is None:
        print(
            "Warning: no accession column found in metadata; "
            "using first column as accession.",
            file=sys.stderr,
        )
        acc_col = metadata_df.columns[0]

    clade_col = _find_column(metadata_df, _CLADE_COLUMNS)
    if clade_col is None:
        print("Error: no clade column in metadata after forecast step.", file=sys.stderr)
        sys.exit(1)

    acc_to_clade: dict[str, str] = dict(zip(
        metadata_df[acc_col].astype(str).str.strip(),
        metadata_df[clade_col].fillna("Unknown").astype(str).str.strip(),
    ))

    clade_to_accs: dict[str, list[str]] = {}
    for acc, clade in acc_to_clade.items():
        clade_to_accs.setdefault(clade, []).append(acc)

    # --- Sample sequences from predicted top clades ---
    rng = random.Random(42)
    selected: list[SeqIO.SeqRecord] = []
    for clade in top_clades:
        accs = clade_to_accs.get(clade, [])
        sampled = _sample_sequences(accs, all_records, top_n_per_clade, rng)
        selected.extend(sampled)
        print(f"  {clade:<30}: {len(sampled):>4} sequences "
              f"(of {len(accs)} in clade, {len(all_records)} total)")

    if not selected:
        print(
            "\nError: no sequences could be matched to predicted clades.",
            file=sys.stderr,
        )
        print(
            "Check that accession IDs in metadata match FASTA header IDs.\n"
            "Example FASTA header: >MW123456.1 Homo sapiens...\n"
            "Example metadata accession: MW123456.1",
            file=sys.stderr,
        )
        sys.exit(1)

    # Write output FASTA
    SeqIO.write(selected, str(out_sequences), "fasta")
    print(f"\nFuture variant FASTA: {out_sequences} ({len(selected)} sequences)")
    print(
        f"\nNext step: evaluate your primer set against {out_sequences.name} "
        "using 'adapt evaluate' or 'adapt forecast --pset <your_primers.fa>'"
    )


def register(subparsers) -> None:
    """Register the forecast-variants subcommand with qprimer."""
    p = subparsers.add_parser(
        "forecast-variants",
        help="Forecast future variant frequencies and select representative sequences",
        description=(
            "Fit a multinomial logistic regression model to sequence metadata to "
            "forecast which viral lineages will be most prevalent in the future, "
            "then select representative sequences from those lineages. The output "
            "FASTA can be used as a target for primer evaluation to measure future "
            "robustness."
        ),
        formatter_class=__import__("argparse").RawDescriptionHelpFormatter,
        epilog="""
Examples:
  qprimer forecast-variants \\
    --sequences target_seqs/original/sars_cov_2.fa \\
    --metadata target_seqs/original/sars_cov_2_metadata.csv \\
    --horizon 30 --top-n 5 \\
    --out-sequences forecast/future_variants.fa \\
    --out-forecast forecast/forecast.tsv

Metadata file format (CSV or TSV):
  Required columns (flexible naming):
    accession / sequence_id    — matches FASTA header IDs
    date / collection_date     — sequence collection date (YYYY-MM-DD)
    clade / pango_lineage      — lineage or clade label
  Optional:
    country / geographic_region — used with --location filter

The metadata file is automatically saved alongside the FASTA when using
'adapt fetch'. For manually uploaded sequences, create a TSV with the
required columns.
""",
    )
    p.add_argument(
        "--sequences",
        required=True,
        help="FASTA file of currently circulating sequences",
    )
    p.add_argument(
        "--metadata",
        required=True,
        help=(
            "Metadata CSV/TSV with accession, date, and clade columns. "
            "Generated automatically by 'adapt fetch' from gget output."
        ),
    )
    p.add_argument(
        "--horizon",
        type=int,
        default=30,
        help="Days ahead to forecast (default: 30)",
    )
    p.add_argument(
        "--min-sequences",
        type=int,
        default=10,
        dest="min_sequences",
        help="Minimum sequences per clade to include in model (default: 10)",
    )
    p.add_argument(
        "--top-n",
        type=int,
        default=5,
        dest="top_n",
        help="Number of top predicted clades to use for the output FASTA (default: 5)",
    )
    p.add_argument(
        "--top-n-per-clade",
        type=int,
        default=20,
        dest="top_n_per_clade",
        help="Maximum sequences to sample per clade (default: 20)",
    )
    p.add_argument(
        "--location",
        help="Filter to sequences from this location substring, e.g. 'USA' (default: global)",
    )
    p.add_argument(
        "--out-sequences",
        required=True,
        dest="out_sequences",
        help="Output FASTA of sequences from predicted future clades",
    )
    p.add_argument(
        "--out-forecast",
        required=True,
        dest="out_forecast",
        help="Output TSV with clade frequency forecast table",
    )
    p.add_argument(
        "--params",
        default="params.txt",
        dest="param_file",
        help="Parameters file for FORECAST_* defaults (default: params.txt)",
    )
    p.set_defaults(func=run)
