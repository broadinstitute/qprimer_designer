"""Select best multiplex primer sets."""

import argparse
import os
import sys
import tempfile
from typing import List

import pandas as pd

from qprimer_designer.utils import parse_params
from qprimer_designer.external import build_index, align_primers, compute_batch_dimer_dg


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
    parser.add_argument("--out", default=None, help="Output combined CSV path (default: auto-generated from target names)")
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
    parser.add_argument(
        "--params",
        dest="param_file",
        default=None,
        help="Parameters file (for DG_MIN threshold in cross-probe hybridization check)",
    )
    parser.set_defaults(func=run)


def infer_target_name(path: str) -> str:
    """Infer target name from filename like final/A.csv -> 'A'."""
    base = os.path.basename(path)
    # Original code: case-sensitive extension matching (e.g., ".CSV" not handled).
    # for ext in [".csv.gz", ".tsv.gz", ".csv", ".tsv"]:
    #     if base.endswith(ext):
    #         base = base[: -len(ext)]
    #         break
    #
    # Fixed: use case-insensitive extension matching.
    base_lower = base.lower()
    for ext in [".csv.gz", ".tsv.gz", ".csv", ".tsv"]:
        if base_lower.endswith(ext):
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


def _parse_sam_cross_hits(sam_path, source_target):
    """Parse SAM file and return probes that align to other targets' amplicons.

    Args:
        sam_path: Path to SAM output file
        source_target: The target that these probes belong to

    Returns:
        List of (probe_name, hit_target, hit_amplicon_name) tuples for cross-target hits
    """
    cross_hits = []
    if not os.path.exists(sam_path) or os.path.getsize(sam_path) == 0:
        return cross_hits

    with open(sam_path) as f:
        for line in f:
            if line.startswith("@"):
                continue
            fields = line.strip().split("\t")
            if len(fields) < 3:
                continue

            probe_name = fields[0]
            flag = int(fields[1])
            ref_name = fields[2]

            # Skip unmapped reads
            if flag & 4:
                continue

            # ref_name format: {target}_{pair_idx}
            # Extract target name (everything before the last underscore+digit group)
            parts = ref_name.rsplit("_", 1)
            if len(parts) < 2:
                continue
            hit_target = parts[0]

            # Only keep cross-target hits
            if hit_target != source_target:
                cross_hits.append((probe_name, hit_target, ref_name))

    return cross_hits


def check_cross_probe_hybridization(combined, dg_threshold):
    """Check for cross-amplicon probe hybridization and remove invalid probes.

    For each target's probes, check if they hybridize to amplicons from other
    targets in the multiplex panel. Uses bowtie2 for fast alignment screening
    followed by RNAduplex for thermodynamic confirmation.

    Args:
        combined: DataFrame with columns including 'target', 'valid_probes',
                  'valid_probe_sequences', and 'amplicon_seq'
        dg_threshold: ΔG threshold (kcal/mol); probes with dG < threshold are invalid

    Returns:
        Updated DataFrame with cross-hybridizing probes removed
    """
    combined = combined.copy()
    targets = combined['target'].unique().tolist()

    if len(targets) < 2:
        print("Only one target in panel, skipping cross-probe hybridization check")
        return combined

    with tempfile.TemporaryDirectory() as tmpdir:
        # 1. Write combined amplicon FASTA (all targets)
        amplicon_fasta = os.path.join(tmpdir, "amplicons.fa")
        amplicon_count = 0
        with open(amplicon_fasta, "w") as f:
            for _, row in combined.iterrows():
                target = row['target']
                amplicon = row.get('amplicon_seq', '')
                if pd.isna(amplicon) or not amplicon:
                    continue
                # Use row index as unique identifier
                seq_name = f"{target}_{amplicon_count}"
                f.write(f">{seq_name}\n{amplicon}\n")
                amplicon_count += 1

        if amplicon_count == 0:
            print("No amplicon sequences found, skipping cross-probe hybridization check")
            return combined

        # 2. Build bowtie2 index
        index_prefix = os.path.join(tmpdir, "amplicon_idx")
        print(f"Building amplicon index from {amplicon_count} sequences...")
        build_index(amplicon_fasta, index_prefix)

        # 3. Per target: write probe FASTA, align, parse cross-target hits
        all_cross_hits = []  # (probe_name, probe_seq, source_target, hit_target, hit_amplicon)

        for target in targets:
            target_rows = combined[combined['target'] == target]

            # Collect all unique probes for this target
            probe_names = []
            probe_seqs = {}
            for _, row in target_rows.iterrows():
                names = str(row.get('valid_probes', ''))
                sequences = str(row.get('valid_probe_sequences', ''))
                if not names or pd.isna(names):
                    continue
                name_list = [n.strip() for n in names.split(',') if n.strip()]
                seq_list = [s.strip() for s in sequences.split(',') if s.strip()]
                for name, seq in zip(name_list, seq_list):
                    if name not in probe_seqs:
                        probe_seqs[name] = seq
                        probe_names.append(name)

            if not probe_names:
                continue

            # Write probe FASTA
            probe_fasta = os.path.join(tmpdir, f"{target}_probes.fa")
            with open(probe_fasta, "w") as f:
                for name in probe_names:
                    f.write(f">{name}\n{probe_seqs[name]}\n")

            # Align probes against combined amplicon index
            sam_path = os.path.join(tmpdir, f"{target}_probes.sam")
            align_primers(index_prefix, probe_fasta, sam_path)

            # Parse cross-target hits
            cross_hits = _parse_sam_cross_hits(sam_path, target)
            for probe_name, hit_target, hit_amplicon in cross_hits:
                if probe_name in probe_seqs:
                    all_cross_hits.append((probe_name, probe_seqs[probe_name], target, hit_target, hit_amplicon))

        if not all_cross_hits:
            print("No cross-target probe alignments found")
            return combined

        print(f"Found {len(all_cross_hits)} cross-target probe alignment(s), checking thermodynamic stability...")

        # 4. Get amplicon sequences for RNAduplex verification
        # Build map from amplicon FASTA name -> sequence
        amplicon_map = {}
        idx = 0
        for _, row in combined.iterrows():
            target = row['target']
            amplicon = row.get('amplicon_seq', '')
            if pd.isna(amplicon) or not amplicon:
                continue
            seq_name = f"{target}_{idx}"
            amplicon_map[seq_name] = amplicon
            idx += 1

        # 5. Batch RNAduplex on (probe, cross-amplicon) pairs
        dg_pairs = []
        hit_info = []
        for probe_name, probe_seq, source_target, hit_target, hit_amplicon in all_cross_hits:
            amplicon_seq = amplicon_map.get(hit_amplicon, "")
            if not amplicon_seq:
                continue
            dg_pairs.append((probe_seq, amplicon_seq))
            hit_info.append((probe_name, source_target, hit_target))

        if not dg_pairs:
            print("No valid probe-amplicon pairs for RNAduplex check")
            return combined

        dg_values = compute_batch_dimer_dg(dg_pairs)

        # 6. Flag probes with dG < threshold
        invalid_probes = set()  # Set of (probe_name, source_target) tuples
        for (probe_name, source_target, hit_target), dg in zip(hit_info, dg_values):
            if dg < dg_threshold:
                invalid_probes.add((probe_name, source_target))
                print(f"  Removing probe {probe_name} (target {source_target}): "
                      f"cross-hybridizes to {hit_target} amplicon (dG={dg:.1f} kcal/mol)")

        if not invalid_probes:
            print("No probes failed thermodynamic cross-hybridization check")
            return combined

        print(f"Removing {len(invalid_probes)} probe(s) due to cross-amplicon hybridization")

        # 7. Update valid_probes and valid_probe_sequences columns
        for idx, row in combined.iterrows():
            target = row['target']
            names_str = str(row.get('valid_probes', ''))
            seqs_str = str(row.get('valid_probe_sequences', ''))

            if not names_str or pd.isna(names_str):
                continue

            name_list = [n.strip() for n in names_str.split(',') if n.strip()]
            seq_list = [s.strip() for s in seqs_str.split(',') if s.strip()]

            filtered_names = []
            filtered_seqs = []
            for name, seq in zip(name_list, seq_list):
                if (name, target) not in invalid_probes:
                    filtered_names.append(name)
                    filtered_seqs.append(seq)

            combined.at[idx, 'valid_probes'] = ','.join(filtered_names)
            combined.at[idx, 'valid_probe_sequences'] = ','.join(filtered_seqs)

    return combined


def run(args):
    """Run the select-multiplex command."""
    if args.top_n < 1:
        raise ValueError("--top-n must be >= 1")

    # Original code: empty csvs list raised unclear pandas error "No objects to concatenate".
    # Fixed: provide a clear error message for empty input.
    if not args.csvs:
        raise ValueError("No input CSV files provided. Please specify at least one CSV file.")

    selected_frames: List[pd.DataFrame] = []

    for path in args.csvs:
        # Original code: relied on pandas error for missing files.
        # Fixed: check file existence with a clear error message.
        if not os.path.isfile(path):
            raise FileNotFoundError(f"Input file not found: {path}")

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

    # Cross-amplicon probe hybridization check (only in probe mode)
    has_probe_cols = (
        'amplicon_seq' in combined.columns
        and 'valid_probes' in combined.columns
        and 'valid_probe_sequences' in combined.columns
    )
    if has_probe_cols:
        # Determine dG threshold from params file
        dg_threshold = -6.0  # default
        if args.param_file:
            params = parse_params(args.param_file)
            dg_threshold = float(params.get("DG_MIN", -6))
        print(f"Running cross-amplicon probe hybridization check (dG threshold: {dg_threshold})...")
        combined = check_cross_probe_hybridization(combined, dg_threshold)
    else:
        print("No probe columns found, skipping cross-amplicon hybridization check")

    out_path = args.out
    if out_path is None:
        # Auto-generate filename from target names
        targets = combined["target"].unique().tolist()
        targets.sort()
        out_path = "_".join(targets) + ".csv"

    combined.to_csv(out_path, index=False)
    print(f"Wrote {len(combined)} rows to {out_path}")
