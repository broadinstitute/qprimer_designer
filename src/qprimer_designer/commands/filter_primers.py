"""Filter top primers from evaluation results."""

import argparse
import ast
from collections import defaultdict

import pandas as pd
from Bio import SeqIO

from qprimer_designer.utils import parse_params
from qprimer_designer.external import compute_batch_dimer_dg


def load_probe_data(probe_mapping_path, probe_seqs_path):
    """Load probe mapping table and sequences, with pre-built indexes."""
    mapping_df = pd.read_csv(probe_mapping_path)
    probe_seqs = {rec.id: str(rec.seq) for rec in SeqIO.parse(probe_seqs_path, "fasta")}

    # Build dict indexes for O(1) lookups instead of repeated DataFrame filtering
    probes_by_target = defaultdict(set)
    probe_positions = {}  # (target_id, probe_name) -> start_pos
    probe_mismatches = {}  # (target_id, probe_name) -> mismatches
    probe_indels = {}  # (target_id, probe_name) -> indels
    for row in mapping_df.itertuples(index=False):
        tid = row.target_id
        pname = row.probe_name
        probes_by_target[tid].add(pname)
        key = (tid, pname)
        if key not in probe_positions:
            probe_positions[key] = row.start_pos
            probe_mismatches[key] = row.mismatches
            probe_indels[key] = getattr(row, 'indels', 0)

    return dict(probes_by_target), probe_positions, probe_mismatches, probe_indels, probe_seqs


def _find_position_valid_probes(
    pair_full_rows,
    probes_by_target,
    probe_positions,
    probe_mismatches,
    probe_indels,
    probe_seqs_dict,
    buffer,
    max_mismatches,
):
    """
    Find probes that pass mapping, position, mismatch, and indel checks.

    Does NOT check dimer ΔG — that is handled via batching in the caller.

    Args:
        pair_full_rows: List of rows from .full file for this pair
        probes_by_target: Dict {target_id: set of probe_names}
        probe_positions: Dict {(target_id, probe_name): start_pos}
        probe_mismatches: Dict {(target_id, probe_name): mismatches}
        probe_indels: Dict {(target_id, probe_name): indels}
        probe_seqs_dict: Dict {probe_name: sequence}
        buffer: Buffer distance from primer ends (nt)
        max_mismatches: Maximum allowed mismatches for probe alignment

    Returns:
        List of probe names passing position, mismatch, and indel checks
    """
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

    # Find candidate probes that map to at least one target
    candidate_probes = set()
    for target_id in all_targets:
        candidate_probes.update(probes_by_target.get(target_id, set()))

    passing = []

    for probe_name in candidate_probes:
        if probe_name not in probe_seqs_dict:
            continue

        probe_len = len(probe_seqs_dict[probe_name])
        valid_for_all = True

        for target_id in all_targets:
            pos_key = (target_id, probe_name)
            if pos_key not in probe_positions:
                valid_for_all = False
                break

            # Check mismatch and indel counts
            mm = probe_mismatches.get(pos_key, 0)
            indel = probe_indels.get(pos_key, 0)
            if mm > max_mismatches or indel > 0:
                valid_for_all = False
                break

            probe_start = probe_positions[pos_key]
            probe_end = probe_start + probe_len

            target_valid = False
            for fwd_start, prod_len in target_positions[target_id]:
                amp_start = fwd_start + buffer
                amp_end = fwd_start + prod_len - buffer
                if amp_start <= probe_start and probe_end <= amp_end:
                    target_valid = True
                    break

            if not target_valid:
                valid_for_all = False
                break

        if valid_for_all:
            passing.append(probe_name)

    return passing


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
    try:
        res = pd.read_csv(args.scores)
    except pd.errors.EmptyDataError:
        print(f"Warning: Scores file {args.scores} is empty, writing empty output")
        open(args.out, "w").close()
        csv_out = args.out.replace(".fa", ".csv")
        open(csv_out, "w").close()
        if args.probe_out:
            open(args.probe_out, "w").close()
        return
    if res.empty:
        print(f"Warning: Scores file {args.scores} has no data rows, writing empty output")
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

        # Load probe data with pre-built indexes
        probes_by_target, probe_positions, probe_mismatches, probe_indels, probe_seqs_dict = load_probe_data(
            args.probe_mapping, args.probe_seqs
        )

        # Parameters
        buffer = int(params.get("PROBE_AMPLICON_BUFFER", 20))
        min_probes = int(params.get("MIN_PROBES_PER_PAIR", 1))
        max_probes = int(params.get("MAX_PROBES_PER_PAIR", 5))
        max_mismatches = int(params.get("PROBE_MAX_MISMATCHES", 2))

        # Select top N pairs FIRST
        top_candidates = res.iloc[:num_select].copy()

        # Load .full file and process row-by-row for memory efficiency
        eval_full_path = f"{args.scores}.full"

        # Build lookup of top pairs we care about
        top_pairs_set = set(zip(top_candidates['pname_f'], top_candidates['pname_r']))

        # Read .full file in chunks and extract ALL rows for top pairs
        # Only keep rows with classifier >= 0.5 (high confidence predictions)
        pair_full_data = {}
        for chunk in pd.read_csv(eval_full_path, chunksize=10000):
            for _, row in chunk.iterrows():
                pair_key = (row['pname_f'], row['pname_r'])
                if row.get('classifier', 1.0) < 0.5:
                    continue
                if pair_key in top_pairs_set:
                    if pair_key not in pair_full_data:
                        pair_full_data[pair_key] = []
                    pair_full_data[pair_key].append(row)

        print(f"Checking top {len(top_candidates)} pairs for probe compatibility...")

        # Phase 1: position filtering (cheap) — collect all position-passing
        # probes per pair, and accumulate dimer pairs for batch computation
        pair_candidates = []  # [(pair_key, [probe_names...])]
        all_dimer_pairs = []  # flat list of (probe_seq, primer_seq) for batch

        for _, row in top_candidates.iterrows():
            pname_f = row['pname_f']
            pname_r = row['pname_r']
            pair_key = (pname_f, pname_r)

            if pair_key not in pair_full_data:
                continue

            candidates = _find_position_valid_probes(
                pair_full_data[pair_key],
                probes_by_target,
                probe_positions,
                probe_mismatches,
                probe_indels,
                probe_seqs_dict,
                buffer,
                max_mismatches,
            )

            if not candidates:
                continue

            pair_candidates.append((pair_key, candidates))
            fwd_seq = primer_seqs[pname_f]
            rev_seq = primer_seqs[pname_r]
            for probe_name in candidates:
                probe_seq = probe_seqs_dict[probe_name]
                all_dimer_pairs.append((probe_seq, fwd_seq))
                all_dimer_pairs.append((probe_seq, rev_seq))

        # Phase 2: batch dimer computation (single RNAduplex process)
        all_dg = []
        if all_dimer_pairs:
            print(f"Computing {len(all_dimer_pairs)} dimer ΔG values in batch...")
            try:
                all_dg = compute_batch_dimer_dg(all_dimer_pairs)
            except Exception as e:
                print(f"Warning: Batch dimer computation failed: {e}")
                all_dg = [None] * len(all_dimer_pairs)

        # Phase 3: distribute results and filter by ΔG threshold
        valid_pairs = []
        valid_probes_per_pair = {}
        dg_idx = 0

        for pair_key, candidates in pair_candidates:
            pair_valid_probes = []
            for probe_name in candidates:
                dg_fwd = all_dg[dg_idx]
                dg_rev = all_dg[dg_idx + 1]
                dg_idx += 2
                if dg_fwd is not None and dg_rev is not None:
                    if dg_fwd > min_dg and dg_rev > min_dg:
                        pair_valid_probes.append(probe_name)

            if len(pair_valid_probes) >= min_probes:
                valid_pairs.append(pair_key)
                valid_probes_per_pair[pair_key] = pair_valid_probes

        print(f"Found {len(valid_pairs)} of top {len(top_candidates)} pairs with compatible probes")

        # Filter to only pairs with valid probes
        if valid_pairs:
            valid_pairs_set = set(valid_pairs)
            top = top_candidates[
                top_candidates.apply(lambda r: (r['pname_f'], r['pname_r']) in valid_pairs_set, axis=1)
            ].copy()
        else:
            print("WARNING: No primer pairs in top N have compatible probes!")
            top = top_candidates.iloc[:0].copy()  # Empty DataFrame

        # Add primer sequences, probe count, and probe names to top DataFrame
        if not top.empty:
            top['pseq_f'] = top['pname_f'].map(primer_seqs)
            top['pseq_r'] = top['pname_r'].map(primer_seqs)
            top['num_probes'] = top.apply(
                lambda r: len(valid_probes_per_pair.get((r['pname_f'], r['pname_r']), [])),
                axis=1
            )
            top['probe_names'] = top.apply(
                lambda r: ','.join(valid_probes_per_pair.get((r['pname_f'], r['pname_r']), [])[:max_probes]),
                axis=1
            )

        # Write probe assignments (limit to first N probes per pair)
        if args.probe_out:
            probe_assignments = []
            for pair_key in valid_pairs:
                pname_f, pname_r = pair_key
                pair_probes = valid_probes_per_pair[pair_key]
                for probe_name in pair_probes[:max_probes]:
                    probe_assignments.append({
                        'pname_f': pname_f,
                        'pname_r': pname_r,
                        'probe_name': probe_name,
                        'probe_seq': probe_seqs_dict[probe_name]
                    })

            probe_df = pd.DataFrame(probe_assignments)
            probe_df.to_csv(args.probe_out, index=False)
            print(f"Wrote {len(probe_df)} probe assignments to {args.probe_out} (max {max_probes} per pair)")

            # Write probe FASTA with unique probes from limited set
            probe_fasta_out = args.probe_out.replace(".csv", ".fa")
            unique_probes = probe_df.drop_duplicates(subset=['probe_name'])

            with open(probe_fasta_out, "w") as out:
                for _, row in unique_probes.iterrows():
                    out.write(f">{row['probe_name']}\n{row['probe_seq']}\n")
            print(f"Wrote {len(unique_probes)} unique probes to {probe_fasta_out}")

    else:
        # STANDARD MODE (no probes)
        print("Standard mode - filtering by primer dimer...")

        # Select top N pairs FIRST
        top_candidates = res.iloc[:num_select].copy()

        # Collect all sequence pairs for batch dimer computation
        pairs_to_check = []
        pair_info = []
        for _, row in top_candidates.iterrows():
            pname_f = row['pname_f']
            pname_r = row['pname_r']
            fwd_seq = primer_seqs.get(pname_f)
            rev_seq = primer_seqs.get(pname_r)
            if fwd_seq is None or rev_seq is None:
                continue
            pairs_to_check.append((fwd_seq, rev_seq))
            pair_info.append((pname_f, pname_r))

        # Single batch RNAduplex call for all pairs
        print(f"Computing primer dimer ΔG for {len(pairs_to_check)} pairs in batch...")
        try:
            dg_values = compute_batch_dimer_dg(pairs_to_check)
        except Exception as e:
            print(f"Warning: Batch dimer computation failed: {e}")
            dg_values = [None] * len(pairs_to_check)

        valid_pairs = []
        for (pname_f, pname_r), dg in zip(pair_info, dg_values):
            if dg is not None and dg > min_dg:
                valid_pairs.append((pname_f, pname_r, dg))

        print(f"Found {len(valid_pairs)} of {len(pairs_to_check)} pairs with primer_dimer_dg > {min_dg}")

        # Filter to only pairs with valid primer dimer
        if valid_pairs:
            valid_pairs_set = set((pf, pr) for pf, pr, _ in valid_pairs)
            top = top_candidates[
                top_candidates.apply(lambda r: (r['pname_f'], r['pname_r']) in valid_pairs_set, axis=1)
            ].copy()

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
            top = top_candidates.iloc[:0].copy()  # Empty DataFrame

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
