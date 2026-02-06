"""Build final primer design output."""

import argparse
import ast
import os
import time

import pandas as pd
from Bio import SeqIO

from qprimer_designer.utils import parse_params
from qprimer_designer.external import compute_dimer_dg


def load_probe_data(probe_mapping_paths, probe_seqs_path):
    """Load probe mapping tables and sequences.

    Args:
        probe_mapping_paths: List of paths to probe mapping CSVs
        probe_seqs_path: Path to probe FASTA file

    Returns:
        Tuple of (list of mapping DataFrames, dict of probe sequences)
    """
    mapping_dfs = [pd.read_csv(path) for path in probe_mapping_paths]
    probe_seqs = {rec.id: str(rec.seq) for rec in SeqIO.parse(probe_seqs_path, "fasta")}
    return mapping_dfs, probe_seqs


def find_valid_probes_offtarget(
    pair_key,
    primer_pair_probes,
    offtarget_full_data_list,
    offtarget_probe_mappings,
    probe_seqs_dict,
    buffer
):
    """
    Find probes that do NOT fall within amplicon regions in off-targets.

    A probe is valid if it does NOT fall within amplicon regions for ANY off-target.

    Args:
        pair_key: Tuple (pname_f, pname_r)
        primer_pair_probes: List of probe names assigned to this primer pair
        offtarget_full_data_list: List of DataFrames (.full files for each off-target)
        offtarget_probe_mappings: List of DataFrames (probe mappings for each off-target)
        probe_seqs_dict: Dict {probe_name: sequence}
        buffer: Buffer distance from primer ends (nt)

    Returns:
        List of valid probe names (those NOT in off-target amplicons)
    """
    valid_probes = []

    for probe_name in primer_pair_probes:
        if probe_name not in probe_seqs_dict:
            continue

        probe_seq = probe_seqs_dict[probe_name]
        probe_len = len(probe_seq)

        # Check if this probe falls within amplicons in ANY off-target
        probe_invalid = False

        for offtarget_full, offtarget_mapping in zip(offtarget_full_data_list, offtarget_probe_mappings):
            # Get rows for this primer pair (either primer can be f or r)
            pair_rows = offtarget_full[
                (offtarget_full["pname_f"].isin(pair_key)) &
                (offtarget_full["pname_r"].isin(pair_key))
            ]

            if pair_rows.empty:
                continue

            # Collect all target-position mappings for this pair in this off-target
            target_positions = {}
            for _, row in pair_rows.iterrows():
                targets = ast.literal_eval(row['targets'])
                starts = ast.literal_eval(row['starts'])
                prod_len = row['prod_len']

                for target_id, fwd_start in zip(targets, starts):
                    if target_id not in target_positions:
                        target_positions[target_id] = []
                    target_positions[target_id].append((fwd_start, prod_len))

            # Check if probe maps to this off-target
            probe_on_offtarget = offtarget_mapping[
                offtarget_mapping['probe_name'] == probe_name
            ]

            if probe_on_offtarget.empty:
                continue

            # Check each mapping of this probe
            for _, probe_row in probe_on_offtarget.iterrows():
                target_id = probe_row['target_id']
                probe_start = probe_row['start_pos']
                probe_end = probe_start + probe_len

                # Check if this target has amplicon positions
                if target_id not in target_positions:
                    continue

                # Check if probe falls within any amplicon region
                for fwd_start, prod_len in target_positions[target_id]:
                    amplicon_start = fwd_start + buffer
                    amplicon_end = fwd_start + prod_len - buffer

                    # If probe is within amplicon, it's INVALID for off-target
                    if amplicon_start <= probe_start and probe_end <= amplicon_end:
                        probe_invalid = True
                        break

                if probe_invalid:
                    break

            if probe_invalid:
                break

        # Only include probes that did NOT fall within any off-target amplicon
        if not probe_invalid:
            valid_probes.append(probe_name)

    return valid_probes


def register(subparsers):
    """Register the build-output subcommand."""
    parser = subparsers.add_parser(
        "build-output",
        help="Build final primer designs",
        description="""
Combine ON-target and OFF-target evaluation results into a final
primer-pair ranking, applying additional filtering based on
primer-primer dimerization (ΔG).
""",
    )
    parser.add_argument("--on", dest="eval_on", required=True, help="ON-target evaluation CSV")
    parser.add_argument("--off", dest="eval_off", nargs="+", required=True, help="OFF-target evaluation CSVs")
    parser.add_argument("--fa", dest="primers", required=True, help="Primer FASTA")
    parser.add_argument("--out", dest="output", required=True, help="Output CSV")
    parser.add_argument("--name", required=True, help="Target name")
    parser.add_argument("--params", dest="param_file", required=True, help="Parameters file")
    parser.add_argument("--probe-mapping-on", help="ON-target probe mapping CSV (for probe mode)")
    parser.add_argument("--probe-mapping-off", nargs="+", help="OFF-target probe mapping CSVs (for probe mode)")
    parser.add_argument("--probe-seqs", help="Probe FASTA (for probe mode)")
    parser.set_defaults(func=run)


def run(args):
    """Run the build-output command."""
    params = parse_params(args.param_file)
    min_dg = float(params.get("DG_MIN", -6))
    num_select = int(params.get("NUM_TOP_SENSITIVITY", 100))

    print(f"Building final output for {args.name}...")
    start_time = time.time()

    # Load ON-target evaluation
    teval = pd.read_csv(args.eval_on, index_col=[0, 1]).iloc[:num_select].copy()
    teval.columns = ["cov_target", "act_target", "sco_target"]

    # Filter teval to only include pairs that exist in the filtered primer CSV
    filt_csv_path = args.primers.replace(".fa", ".csv")
    filt_pairs = pd.read_csv(filt_csv_path)
    valid_pairs = set(zip(filt_pairs['pname_f'], filt_pairs['pname_r']))
    teval = teval[teval.index.isin(valid_pairs)].copy()
    print(f"Filtered to {len(teval)} pairs from {filt_csv_path}")

    merged = teval.copy()

    # Merge OFF-target evaluations
    for off_path in args.eval_off:
        oname = os.path.basename(off_path).split(".")[1]
        oeval = pd.read_csv(off_path)

        tmp = pd.DataFrame(
            index=teval.index,
            columns=[f"cov_{oname}", f"act_{oname}", f"sco_{oname}"],
        )

        for pair in teval.index:
            sub = oeval[(oeval["pname_f"].isin(pair)) & (oeval["pname_r"].isin(pair))]
            if sub.empty:
                tmp.loc[pair] = 0
            else:
                tmp.loc[pair, f"cov_{oname}"] = sub["coverage"].sum()
                tmp.loc[pair, f"act_{oname}"] = sub["activity"].max()
                tmp.loc[pair, f"sco_{oname}"] = sub["score"].sum()

        merged = merged.join(tmp)

    primer_seqs = {s.id: str(s.seq) for s in SeqIO.parse(args.primers, "fasta")}

    for fname, rname in merged.index:
        fseq = primer_seqs[fname]
        rseq = primer_seqs[rname]
        merged.loc[(fname, rname), "pseq_f"] = fseq
        merged.loc[(fname, rname), "pseq_r"] = rseq

    final = merged.copy()

    # PROBE MODE: Filter probes by off-target amplicon overlap
    if args.probe_mapping_on and args.probe_mapping_off and args.probe_seqs:
        print("Probe mode enabled - filtering probes by off-target amplicons...")

        # Load probe data
        probe_seqs_dict = {rec.id: str(rec.seq) for rec in SeqIO.parse(args.probe_seqs, "fasta")}

        # Get probe names from filtered primer pairs CSV (top 5 probes per pair)
        pair_to_probes = {}
        for _, row in filt_pairs.iterrows():
            pair_key = (row['pname_f'], row['pname_r'])
            if 'probe_names' in row and pd.notna(row['probe_names']) and row['probe_names']:
                # Split comma-separated probe names
                probe_names = [p.strip() for p in row['probe_names'].split(',') if p.strip()]
                pair_to_probes[pair_key] = probe_names
            else:
                pair_to_probes[pair_key] = []

        # Load off-target probe mappings
        offtarget_probe_mappings = [pd.read_csv(path) for path in args.probe_mapping_off]

        # Load off-target .full files and filter by classifier >= 0.5
        buffer = 1
        offtarget_full_data_list = []
        for off_path in args.eval_off:
            full_path = f"{off_path}.full"
            if os.path.exists(full_path):
                df = pd.read_csv(full_path)
                # Filter to high-confidence predictions only
                df = df[df.get('classifier', pd.Series([1.0] * len(df))) >= 0.5]
                offtarget_full_data_list.append(df)
            else:
                print(f"Warning: .full file not found: {full_path}")
                offtarget_full_data_list.append(pd.DataFrame())

        # Find valid probes for each primer pair
        valid_probes_dict = {}
        for pair_key in final.index:
            if pair_key not in pair_to_probes:
                valid_probes_dict[pair_key] = []
                continue

            primer_pair_probes = pair_to_probes[pair_key]

            valid_probes = find_valid_probes_offtarget(
                pair_key,
                primer_pair_probes,
                offtarget_full_data_list,
                offtarget_probe_mappings,
                probe_seqs_dict,
                buffer
            )

            valid_probes_dict[pair_key] = valid_probes

        # Add probe information to final DataFrame
        final['valid_probes'] = final.index.map(lambda p: ','.join(valid_probes_dict.get(p, [])))
        final['valid_probe_sequences'] = final.index.map(
            lambda p: ','.join([probe_seqs_dict[pname] for pname in valid_probes_dict.get(p, []) if pname in probe_seqs_dict])
        )

        print(f"Probe filtering complete. Pairs with valid probes: {(final['valid_probe_sequences'].str.len() > 0).sum()}")

    # Sort by sum of off-target scores (ascending - lower is better)
    offtarget_score_cols = [col for col in final.columns if col.startswith('sco_') and col != 'sco_target']
    if offtarget_score_cols:
        final['offtarget_score_sum'] = final[offtarget_score_cols].sum(axis=1)
        final = final.sort_values('offtarget_score_sum')
        print(f"Sorted by sum of off-target scores")


    # Write output
    final.round(3).to_csv(args.output)

    runtime = time.time() - start_time
    print(f"Wrote {len(final)} primer pairs to {args.output} ({runtime:.1f} sec)")
