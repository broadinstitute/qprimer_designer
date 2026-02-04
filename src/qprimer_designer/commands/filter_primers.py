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
    pair_full_data,
    probe_mapping_df,
    probe_seqs_dict,
    primer_seqs_dict,
    min_dg,
    buffer
):
    """
    Find probes compatible with a primer pair.

    Args:
        pair_full_data: Row from .full file with targets, starts, pname_f, pname_r, etc.
        probe_mapping_df: DataFrame with probe mappings
        probe_seqs_dict: Dict {probe_name: sequence}
        primer_seqs_dict: Dict {primer_name: sequence}
        min_dg: Minimum ΔG threshold for dimers
        buffer: Buffer distance from primer ends (nt)

    Returns:
        List of valid probe names
    """
    valid_probes = []

    # Parse target lists from .full file
    targets = ast.literal_eval(pair_full_data['targets'])
    starts = ast.literal_eval(pair_full_data['starts'])

    # Get primer sequences and amplicon length
    pname_f = pair_full_data['pname_f']
    pname_r = pair_full_data['pname_r']
    fwd_seq = primer_seqs_dict[pname_f]
    rev_seq = primer_seqs_dict[pname_r]
    prod_len = pair_full_data['prod_len']

    # Check each target this pair covers
    for target_id, fwd_start in zip(targets, starts):
        # Find probes mapped to this target
        probes_on_target = probe_mapping_df[probe_mapping_df['target_id'] == target_id]

        for _, probe_row in probes_on_target.iterrows():
            probe_name = probe_row['probe_name']

            # Skip if we don't have this probe sequence
            if probe_name not in probe_seqs_dict:
                continue

            probe_seq = probe_seqs_dict[probe_name]
            probe_start = probe_row['start_pos']
            probe_len = len(probe_seq)
            probe_end = probe_start + probe_len

            # Check position: must be ≥20nt from both primer ends
            amplicon_start = fwd_start + buffer
            amplicon_end = fwd_start + prod_len - buffer

            if not (amplicon_start <= probe_start and probe_end <= amplicon_end):
                continue

            # Check dimer formation with both primers
            try:
                dg_fwd = compute_dimer_dg(probe_seq, fwd_seq)
                dg_rev = compute_dimer_dg(probe_seq, rev_seq)

                # Both dimers must be above threshold
                if dg_fwd > min_dg and dg_rev > min_dg:
                    valid_probes.append(probe_name)
            except Exception as e:
                print(f"Warning: Failed to compute dimer for {probe_name}: {e}")
                continue

    return list(set(valid_probes))  # Deduplicate


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

    # Load evaluation results
    res = pd.read_csv(args.scores)
    if res.empty:
        open(args.out, "w").close()
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
        min_dg = float(params.get("DG_MIN", -6))
        buffer = int(params.get("PROBE_AMPLICON_BUFFER", 20))
        min_probes = int(params.get("MIN_PROBES_PER_PAIR", 1))

        # Select top N pairs FIRST
        top_candidates = res.iloc[:num_select]

        # Load .full file and process row-by-row for memory efficiency
        eval_full_path = f"{args.scores}.full"

        # Build lookup of top pairs we care about
        top_pairs_set = set(zip(top_candidates['pname_f'], top_candidates['pname_r']))

        # Read .full file in chunks and extract only rows for top pairs
        pair_full_data = {}
        for chunk in pd.read_csv(eval_full_path, chunksize=10000):
            for _, row in chunk.iterrows():
                pair_key = (row['pname_f'], row['pname_r'])
                if pair_key in top_pairs_set:
                    pair_full_data[pair_key] = row

        print(f"Checking top {len(top_candidates)} pairs for probe compatibility...")

        # Check probe compatibility for top N pairs
        probe_assignments = []
        valid_pairs = []

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

        # Write probe assignments (already filtered to final pairs)
        if args.probe_out:
            probe_df = pd.DataFrame(probe_assignments)
            probe_df.to_csv(args.probe_out, index=False)
            print(f"Wrote {len(probe_df)} probe assignments to {args.probe_out}")

    else:
        # STANDARD MODE (no probes)
        top = res.iloc[:num_select]

    # Collect unique primer names from top pairs
    primer_names = set(top["pname_f"]) | set(top["pname_r"])

    # Write primer FASTA
    with open(args.out, "w") as out:
        for pname in primer_names:
            if pname not in primer_seqs:
                raise KeyError(f"Primer {pname} not found in {args.init}")
            out.write(f">{pname}\n{primer_seqs[pname]}\n")

    print(f"Wrote {len(primer_names)} primers to {args.out}")
