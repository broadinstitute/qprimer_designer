"""Filter top primers from evaluation results."""

import argparse
import ast

import pandas as pd
from Bio import SeqIO

from qprimer_designer.utils import parse_params
from qprimer_designer.external import compute_dimer_dg


def load_probe_data(probe_mapping_path, probe_seqs_path):
    """Load probe mapping table and sequences."""
    mapping_df = pd.read_csv(probe_mapping_path)
    probe_seqs = {rec.id: str(rec.seq) for rec in SeqIO.parse(probe_seqs_path, "fasta")}
    return mapping_df, probe_seqs


def find_valid_probes_for_pair(
    pair_full_rows,
    probe_mapping_df,
    probe_seqs_dict,
    primer_seqs_dict,
    min_dg,
    buffer
):
    """
    Find probes compatible with a primer pair across ALL its target mappings.

    A probe is only valid if it:
    1. Maps to ALL targets the primer pair covers
    2. Falls within amplicon region for ALL those targets
    3. Has good dimer ΔG with both primers

    Args:
        pair_full_rows: List of rows from .full file for this pair (each with different target sets)
        probe_mapping_df: DataFrame with probe mappings
        probe_seqs_dict: Dict {probe_name: sequence}
        primer_seqs_dict: Dict {primer_name: sequence}
        min_dg: Minimum ΔG threshold for dimers
        buffer: Buffer distance from primer ends (nt)

    Returns:
        List of valid probe names
    """
    # Get primer sequences (same for all rows)
    first_row = pair_full_rows[0]
    pname_f = first_row['pname_f']
    pname_r = first_row['pname_r']
    fwd_seq = primer_seqs_dict[pname_f]
    rev_seq = primer_seqs_dict[pname_r]

    # Collect all target-position mappings for this primer pair
    target_positions = {}  # {target_id: [(fwd_start, prod_len), ...]}

    for pair_row in pair_full_rows:
        targets = ast.literal_eval(pair_row['targets'])
        starts = ast.literal_eval(pair_row['starts'])
        prod_len = pair_row['prod_len']

        for target_id, fwd_start in zip(targets, starts):
            if target_id not in target_positions:
                target_positions[target_id] = []
            target_positions[target_id].append((fwd_start, prod_len))

    all_targets = set(target_positions.keys())

    # Find all candidate probes that map to at least one target
    candidate_probes = set()
    for target_id in all_targets:
        probes_on_target = probe_mapping_df[probe_mapping_df['target_id'] == target_id]
        for _, probe_row in probes_on_target.iterrows():
            candidate_probes.add(probe_row['probe_name'])

    valid_probes = []

    # Check each candidate probe
    for probe_name in candidate_probes:
        if probe_name not in probe_seqs_dict:
            continue

        probe_seq = probe_seqs_dict[probe_name]

        # Check dimer formation with primers (once, doesn't depend on target)
        try:
            dg_fwd = compute_dimer_dg(probe_seq, fwd_seq)
            dg_rev = compute_dimer_dg(probe_seq, rev_seq)

            if dg_fwd <= min_dg or dg_rev <= min_dg:
                continue
        except Exception as e:
            print(f"Warning: Failed to compute dimer for {probe_name}: {e}")
            continue

        # Check if probe is valid for ALL targets
        valid_for_all_targets = True

        for target_id in all_targets:
            # Get probe mappings for this target
            probe_on_this_target = probe_mapping_df[
                (probe_mapping_df['target_id'] == target_id) &
                (probe_mapping_df['probe_name'] == probe_name)
            ]

            # Probe must map to this target
            if probe_on_this_target.empty:
                valid_for_all_targets = False
                break

            # Get probe position on this target
            probe_start = probe_on_this_target.iloc[0]['start_pos']
            probe_len = len(probe_seq)
            probe_end = probe_start + probe_len

            # Check if probe is within amplicon for at least one position on this target
            # (A target may appear multiple times with different positions)
            target_valid = False
            for fwd_start, prod_len in target_positions[target_id]:
                amplicon_start = fwd_start + buffer
                amplicon_end = fwd_start + prod_len - buffer

                if amplicon_start <= probe_start and probe_end <= amplicon_end:
                    target_valid = True
                    break

            if not target_valid:
                valid_for_all_targets = False
                break

        if valid_for_all_targets:
            valid_probes.append(probe_name)

    return valid_probes


def register(subparsers):
    """Register the filter subcommand."""
    parser = subparsers.add_parser(
        "filter",
        help="Filter top primers from evaluation results",
        description="""
Select top-performing primers based on evaluation results and write
their sequences to a FASTA file.
""",
    )
    parser.add_argument("--scores", required=True, help="Evaluation CSV (output of evaluate)")
    parser.add_argument("--init", required=True, help="Initial primer FASTA")
    parser.add_argument("--out", required=True, help="Output FASTA of passed primers")
    parser.add_argument("--params", required=True, help="Parameters file")
    parser.add_argument("--probe-mapping", help="Probe mapping CSV (for probe mode)")
    parser.add_argument("--probe-seqs", help="Probe FASTA (for probe mode)")
    parser.add_argument("--probe-out", help="Output CSV for probe assignments")
    parser.set_defaults(func=run)


def run(args):
    """Run filter command with optional probe mode."""
    params = parse_params(args.params)
    num_select = int(params.get("NUM_TOP_SENSITIVITY", 100))
    min_dg = float(params.get("DG_MIN", -6))

    # Load evaluation results
    res = pd.read_csv(args.scores)
    if res.empty:
        open(args.out, "w").close()
        csv_out = args.out.replace(".fa", "_pairs.csv")
        open(csv_out, "w").close()
        if args.probe_out:
            open(args.probe_out, "w").close()
        return

    # Load primer sequences
    primer_seqs = {rec.id: str(rec.seq) for rec in SeqIO.parse(args.init, "fasta")}

    # PROBE MODE
    if args.probe_mapping and args.probe_seqs:
        print("Probe mode enabled - filtering by probe compatibility...")

        # Load probe data
        probe_mapping_df, probe_seqs_dict = load_probe_data(
            args.probe_mapping, args.probe_seqs
        )

        # Parameters
        buffer = int(params.get("PROBE_AMPLICON_BUFFER", 20))
        min_probes = int(params.get("MIN_PROBES_PER_PAIR", 1))

        # Select top N pairs FIRST
        top_candidates = res.iloc[:num_select].copy()

        # Load .full file and process row-by-row for memory efficiency
        eval_full_path = f"{args.scores}.full"

        # Build lookup of top pairs we care about
        top_pairs_set = set(zip(top_candidates['pname_f'], top_candidates['pname_r']))

        # Read .full file in chunks and extract ALL rows for top pairs
        # Each pair may have multiple rows with different target sets
        pair_full_data = {}
        for chunk in pd.read_csv(eval_full_path, chunksize=10000):
            for _, row in chunk.iterrows():
                pair_key = (row['pname_f'], row['pname_r'])
                if pair_key in top_pairs_set:
                    if pair_key not in pair_full_data:
                        pair_full_data[pair_key] = []
                    pair_full_data[pair_key].append(row)

        print(f"Checking top {len(top_candidates)} pairs for probe compatibility...")

        # Check probe compatibility for top N pairs
        probe_assignments = []
        valid_pairs = []
        valid_probes_per_pair = {}

        for _, row in top_candidates.iterrows():
            pname_f = row['pname_f']
            pname_r = row['pname_r']
            pair_key = (pname_f, pname_r)

            # Get detailed data from .full file
            if pair_key not in pair_full_data:
                continue

            # Find valid probes for this pair
            valid_probes = find_valid_probes_for_pair(
                pair_full_data[pair_key],
                probe_mapping_df,
                probe_seqs_dict,
                primer_seqs,
                min_dg,
                buffer
            )

            # Keep pair if it has enough valid probes
            if len(valid_probes) >= min_probes:
                valid_pairs.append(pair_key)
                valid_probes_per_pair[pair_key] = valid_probes

                # Record probe assignments
                for probe in valid_probes:
                    probe_assignments.append({
                        'pname_f': pname_f,
                        'pname_r': pname_r,
                        'probe_name': probe,
                        'probe_seq': probe_seqs_dict[probe]
                    })

        print(f"Found {len(valid_pairs)} of top {len(top_candidates)} pairs with compatible probes")

        # Filter to only pairs with valid probes
        if valid_pairs:
            valid_pairs_set = set(valid_pairs)
            top = top_candidates[
                top_candidates.apply(lambda r: (r['pname_f'], r['pname_r']) in valid_pairs_set, axis=1)
            ]
        else:
            print("WARNING: No primer pairs in top N have compatible probes!")
            top = top_candidates.iloc[:0]  # Empty DataFrame

        # Add primer sequences and probe count to top DataFrame
        if not top.empty:
            top['pseq_f'] = top['pname_f'].map(primer_seqs)
            top['pseq_r'] = top['pname_r'].map(primer_seqs)
            top['num_probes'] = top.apply(
                lambda r: len(valid_probes_per_pair.get((r['pname_f'], r['pname_r']), [])),
                axis=1
            )

        # Write probe assignments CSV
        if args.probe_out:
            probe_df = pd.DataFrame(probe_assignments)
            probe_df.to_csv(args.probe_out, index=False)
            print(f"Wrote {len(probe_df)} probe assignments to {args.probe_out}")

            # Write probe FASTA
            probe_fasta_out = args.probe_out.replace(".csv", ".fa")
            unique_probes = probe_df[['probe_name', 'probe_seq']].drop_duplicates()
            with open(probe_fasta_out, "w") as out:
                for _, row in unique_probes.iterrows():
                    out.write(f">{row['probe_name']}\n{row['probe_seq']}\n")
            print(f"Wrote {len(unique_probes)} unique probes to {probe_fasta_out}")

    else:
        # STANDARD MODE (no probes)
        print("Standard mode - filtering by primer dimer...")

        # Select top N pairs FIRST
        top_candidates = res.iloc[:num_select].copy()

        # Compute primer dimer ΔG for each pair
        print(f"Computing primer dimer ΔG for top {len(top_candidates)} pairs...")
        valid_pairs = []

        for _, row in top_candidates.iterrows():
            pname_f = row['pname_f']
            pname_r = row['pname_r']

            # Get primer sequences
            fwd_seq = primer_seqs.get(pname_f)
            rev_seq = primer_seqs.get(pname_r)

            if fwd_seq is None or rev_seq is None:
                continue

            # Compute primer dimer ΔG
            try:
                dimer_dg = compute_dimer_dg(fwd_seq, rev_seq)
                if dimer_dg > min_dg:
                    valid_pairs.append((pname_f, pname_r, dimer_dg))
            except Exception as e:
                print(f"Warning: Failed to compute dimer for {pname_f}/{pname_r}: {e}")
                continue

        print(f"Found {len(valid_pairs)} of top {len(top_candidates)} pairs with primer_dimer_dg > {min_dg}")

        # Filter to only pairs with valid primer dimer
        if valid_pairs:
            valid_pairs_set = set((pf, pr) for pf, pr, _ in valid_pairs)
            top = top_candidates[
                top_candidates.apply(lambda r: (r['pname_f'], r['pname_r']) in valid_pairs_set, axis=1)
            ]

            # Add primer sequences and dimer ΔG
            dimer_dict = {(pf, pr): dg for pf, pr, dg in valid_pairs}
            top['pseq_f'] = top['pname_f'].map(primer_seqs)
            top['pseq_r'] = top['pname_r'].map(primer_seqs)
            top['primer_dimer_dg'] = top.apply(
                lambda r: dimer_dict.get((r['pname_f'], r['pname_r'])),
                axis=1
            )
        else:
            print("WARNING: No primer pairs in top N have primer_dimer_dg > DG_MIN!")
            top = top_candidates.iloc[:0]  # Empty DataFrame

    # Write primer pairs CSV (always)
    csv_out = args.out.replace(".fa", ".csv")
    if not top.empty:
        top.to_csv(csv_out, index=False)
        print(f"Wrote {len(top)} primer pairs to {csv_out}")
    else:
        open(csv_out, "w").close()
        print(f"No valid pairs - created empty {csv_out}")

    # Collect unique primer names from top pairs
    primer_names = set(top["pname_f"]) | set(top["pname_r"]) if not top.empty else set()

    # Write primer FASTA (always)
    with open(args.out, "w") as out:
        for pname in primer_names:
            if pname not in primer_seqs:
                raise KeyError(f"Primer {pname} not found in {args.init}")
            out.write(f">{pname}\n{primer_seqs[pname]}\n")

    print(f"Wrote {len(primer_names)} primers to {args.out}")
