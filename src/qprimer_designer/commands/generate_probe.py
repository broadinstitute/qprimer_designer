"""Generate candidate probe sequences from target sequences."""

import argparse
import random
import time
from collections import defaultdict

import pandas as pd
from Bio import SeqIO
from Bio.SeqUtils import gc_fraction

from qprimer_designer.utils import get_tm, parse_params, get_probe_params, reverse_complement_dna, has_homopolymer, sanitize_iupac
from qprimer_designer.external import compute_self_dimer_dg


def register(subparsers):
    """Register the generate-probe subcommand."""
    parser = subparsers.add_parser(
        "generate-probe",
        help="Generate candidate probe sequences from target sequences",
        description="""
Generate candidate probe sequences from target sequences by tiling across
all possible lengths, then filtering by melting temperature (Tm), GC content,
self-dimer free energy (ΔG), homopolymer runs, and 5' nucleotide constraints.
""",
    )
    parser.add_argument("--in", dest="target_seqs", required=True,
                       help="Input FASTA of target sequences")
    parser.add_argument("--out", dest="probe_seqs", required=True,
                       help="Output FASTA for probes")
    parser.add_argument("--params", dest="param_file", required=True,
                       help="Parameters file (params.txt)")
    parser.add_argument("--name", required=True,
                       help="Name prefix for probe IDs")
    parser.set_defaults(func=run)


def generate_probes(
    target_seqs,
    len_min: int,
    len_max: int,
    max_tm: float,
    min_tm: float,
    max_gc: float,
    min_dg: float,
    homopolymer_max: int,
    avoid_5prime_G: bool,
    max_num: int,
):
    """Generate and filter probe candidates across multiple target sequences."""
    # Generate all candidate probes with step=1 from both strands
    probes = {}
    for target_seq in target_seqs:
        seq_len = len(target_seq)

        # Generate from forward strand
        for probe_len in range(len_min, len_max + 1):
            for i in range(0, seq_len - probe_len + 1):
                seq = target_seq[i : i + probe_len]
                if "N" not in seq:
                    probes[seq] = i

        # Generate from reverse complement strand
        target_seq_rc = reverse_complement_dna(target_seq)
        for probe_len in range(len_min, len_max + 1):
            for i in range(0, seq_len - probe_len + 1):
                seq = target_seq_rc[i : i + probe_len]
                if "N" not in seq:
                    # Convert RC position to original sequence coordinates
                    # Position in original = len - position_in_rc - probe_len
                    original_pos = seq_len - i - probe_len
                    probes[seq] = original_pos

    print(f">> Probes with unique sequence: {len(probes)}")

    # Filter probes
    filtered = {}
    features = defaultdict(dict)

    for pseq in probes:
        # Check 5' G constraint
        if avoid_5prime_G and pseq[0] == 'G':
            continue

        # Check homopolymer
        if has_homopolymer(pseq, homopolymer_max):
            continue

        # Calculate features
        tm = get_tm(pseq)
        gc = gc_fraction(pseq)

        # Filter by Tm and GC
        if min_tm <= tm <= max_tm and gc <= max_gc:
            dG = compute_self_dimer_dg(pseq)

            # Only add features if all filters pass
            if dG >= min_dg:
                features[pseq]["Tm"] = round(tm, 1)
                features[pseq]["GC"] = round(gc, 2)
                features[pseq]["len"] = len(pseq)
                features[pseq]["position"] = probes[pseq]
                features[pseq]["dG"] = round(dG, 1)
                filtered[pseq] = probes[pseq]

    print(f">> Probes filtered: {len(filtered)}")

    # Subsample if needed
    if len(filtered) > max_num:
        selected_seqs = random.sample(list(filtered.keys()), max_num)
        filtered = {seq: filtered[seq] for seq in selected_seqs}
        print(f">> Randomly subsampled to: {len(filtered)}")

    return filtered, features


def run(args):
    """Run the generate-probe command."""
    params = parse_params(args.param_file)
    probe_params = get_probe_params(params)

    len_min = probe_params["len_min"]
    len_max = probe_params["len_max"]
    min_tm = probe_params["min_tm"]
    max_tm = probe_params["max_tm"]
    max_gc = probe_params["max_gc"]
    min_dg = probe_params["min_dg"]
    homopolymer_max = probe_params["homopolymer_max"]
    avoid_5prime_G = probe_params["avoid_5prime_g"]
    max_num = probe_params["max_num"]

    target_seqs = [sanitize_iupac(str(s.seq)) for s in SeqIO.parse(args.target_seqs, "fasta")]

    print(f"Generating probes from {args.target_seqs}...")
    print(f"Probe lengths: {len_min}-{len_max}")
    print(f"Tm range: {min_tm}-{max_tm}")
    start_time = time.time()

    filtered, features = generate_probes(
        target_seqs, len_min, len_max,
        max_tm, min_tm, max_gc, min_dg,
        homopolymer_max, avoid_5prime_G, max_num,
    )

    # Write output FASTA
    probe_list = list(filtered.keys())
    with open(args.probe_seqs, "w") as fout:
        for i, seq in enumerate(probe_list):
            pname = f"{args.name}_{i+1}_probe"
            features[seq]["pname"] = pname
            fout.write(f">{pname}\n{seq}\n")

    # Write features CSV
    if features:
        features_df = pd.DataFrame(features).T
        features_df = features_df.reset_index(names="pseq")[["pname", "pseq", "len", "position", "Tm", "GC", "dG"]]
    else:
        features_df = pd.DataFrame(columns=["pname", "pseq", "len", "position", "Tm", "GC", "dG"])

    fname = args.probe_seqs.replace(".fa", ".feat")
    features_df.to_csv(fname, index=False)

    runtime = time.time() - start_time
    print(f"Wrote {len(probe_list)} probes to {args.probe_seqs} ({runtime:.1f} sec)")
