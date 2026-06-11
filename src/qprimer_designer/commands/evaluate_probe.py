"""Evaluate probes by direct wobble-aware matching against reference sequences.

Replaces bowtie2-based probe alignment for evaluate mode. For each primer pair,
extracts the amplicon region from each target sequence and slides the probe
across it using wobble-weighted mismatch counting.
"""

import ast
from pathlib import Path

import pandas as pd
from Bio import SeqIO

from qprimer_designer.utils.probe import build_match_string, slide_probe_match, wobble_mismatch_count
from qprimer_designer.utils.sequences import complement_dna, reverse_complement_dna


def register(subparsers):
    """Register the evaluate-probe subcommand."""
    parser = subparsers.add_parser(
        "evaluate-probe",
        help="Evaluate probes by direct wobble-aware matching",
        description=(
            "Match probes against amplicon regions extracted from reference "
            "sequences using wobble-aware mismatch counting (G-T/A-G = 0.20 "
            "penalty). Replaces bowtie2-based probe alignment."
        ),
    )
    parser.add_argument("--probe-fa", dest="probe_fa", required=True,
                        help="Probe FASTA file")
    parser.add_argument("--eval-full", dest="eval_full", required=True,
                        help="Path to .eval.full file (has amplicon coordinates)")
    parser.add_argument("--ref", required=True,
                        help="Reference FASTA file")
    parser.add_argument("--out", required=True,
                        help="Output probe mapping CSV")
    parser.add_argument("--max-mismatches", dest="max_mismatches", type=float,
                        default=3.0,
                        help="Max wobble-adjusted mismatches (default: 3.0)")
    parser.add_argument("--max-indels", dest="max_indels", type=int, default=0,
                        help="Max indels allowed (default: 0)")
    parser.set_defaults(func=run)


def run(args):
    """Run evaluate-probe: match probes to amplicon regions."""
    probe_fa = Path(args.probe_fa)
    eval_full_path = Path(args.eval_full)
    ref_path = Path(args.ref)
    out_path = Path(args.out)

    if not probe_fa.exists() or probe_fa.stat().st_size == 0:
        print("No probe FASTA provided or empty. Writing empty output.")
        out_path.parent.mkdir(parents=True, exist_ok=True)
        pd.DataFrame(columns=[
            'probe_name', 'probe_seq', 'target_id', 'start_pos',
            'orientation', 'mismatches', 'indels', 'tseq', 'match'
        ]).to_csv(out_path, index=False)
        return

    if not eval_full_path.exists() or eval_full_path.stat().st_size == 0:
        print("No .eval.full file found. Writing empty output.")
        out_path.parent.mkdir(parents=True, exist_ok=True)
        pd.DataFrame(columns=[
            'probe_name', 'probe_seq', 'target_id', 'start_pos',
            'orientation', 'mismatches', 'indels', 'tseq', 'match'
        ]).to_csv(out_path, index=False)
        return

    # Load probes with strand info from FASTA description
    probes = {}
    probe_strands = {}
    for rec in SeqIO.parse(str(probe_fa), "fasta"):
        probes[rec.id] = str(rec.seq).upper()
        # Parse strand from description (e.g. "probe_r0_1 strand=fwd")
        strand = 'fwd'  # default
        desc = rec.description
        if 'strand=' in desc:
            strand = desc.split('strand=')[1].split()[0]
        probe_strands[rec.id] = strand
    print(f"  Loaded {len(probes)} probe(s)")

    # Load reference sequences
    ref_seqs = {}
    for rec in SeqIO.parse(str(ref_path), "fasta"):
        ref_seqs[rec.id] = str(rec.seq).upper()
    print(f"  Loaded {len(ref_seqs)} reference sequence(s)")

    # Parse eval.full to get amplicon regions per target
    eval_df = pd.read_csv(eval_full_path)
    print(f"  Loaded {len(eval_df)} primer pair(s) from eval.full")

    max_mm = args.max_mismatches
    max_indels = args.max_indels
    mappings = []

    for _, row in eval_df.iterrows():
        # Parse targets and starts lists
        try:
            targets = ast.literal_eval(str(row['targets']))
            starts = ast.literal_eval(str(row['starts']))
        except (ValueError, SyntaxError):
            continue

        prod_len = int(row['prod_len'])

        for target_id, amp_start in zip(targets, starts):
            if target_id not in ref_seqs:
                continue

            ref_seq = ref_seqs[target_id]
            amp_start = int(amp_start)
            amp_end = min(amp_start + prod_len, len(ref_seq))
            amplicon = ref_seq[amp_start:amp_end]

            if len(amplicon) < 10:
                continue

            for probe_name, probe_seq in probes.items():
                # Test both orientations: TaqMan probes bind single-
                # stranded DNA after denaturation, so either strand works.
                # MSA strand may differ from genome strand at the actual
                # amplicon position.
                probe_rc = reverse_complement_dna(probe_seq)
                for match_seq, orientation in [
                    (probe_seq, '+'),
                    (probe_rc, '-'),
                ]:
                    probe_len = len(match_seq)
                    for i in range(len(amplicon) - probe_len + 1):
                        window = amplicon[i:i + probe_len]
                        mm, indels = wobble_mismatch_count(
                            match_seq.upper(), window.upper())
                        if mm <= max_mm and indels <= max_indels:
                            abs_pos = amp_start + i
                            window_comp = complement_dna(window.upper())
                            match_str = build_match_string(
                                match_seq.upper(), window_comp)
                            mappings.append({
                                'probe_name': probe_name,
                                'probe_seq': probe_seq,
                                'target_id': target_id,
                                'start_pos': abs_pos,
                                'orientation': orientation,
                                'mismatches': round(mm, 2),
                                'indels': indels,
                                'tseq': window_comp,
                                'match': match_str,
                            })

    # Deduplicate: same probe + target + position + orientation
    df = pd.DataFrame(mappings)
    if not df.empty:
        df = df.drop_duplicates(
            subset=['probe_name', 'target_id', 'start_pos', 'orientation'])

    out_path.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(out_path, index=False)

    n_targets_hit = df['target_id'].nunique() if not df.empty else 0
    print(f"  Probe mapping: {len(df)} hits across {n_targets_hit} target(s)")
    print(f"  Saved: {out_path}")
