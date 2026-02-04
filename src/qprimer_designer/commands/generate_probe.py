"""Generate candidate probe sequences from target sequences."""

import argparse
import random
import time
from collections import defaultdict

import pandas as pd
from Bio import SeqIO
from Bio.SeqUtils import gc_fraction

from qprimer_designer.utils import get_tm, parse_params
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


def has_homopolymer(seq: str, max_len: int) -> bool:
    """Check if sequence has homopolymer run > max_len."""
    count = 1
    for i in range(1, len(seq)):
        if seq[i] == seq[i-1]:
            count += 1
            if count > max_len:
                return True
        else:
            count = 1
    return False


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
    # Generate all candidate probes with step=1
    probes = {}
    for target_seq in target_seqs:
        for probe_len in range(len_min, len_max + 1):
            for i in range(0, len(target_seq) - probe_len + 1):
                seq = target_seq[i : i + probe_len]
                if "N" not in seq:
                    probes[seq] = i

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

        features[pseq]["Tm"] = round(tm, 1)
        features[pseq]["GC"] = gc
        features[pseq]["len"] = len(pseq)

        # Filter by Tm and GC
        if min_tm <= tm <= max_tm and gc <= max_gc:
            dG = compute_self_dimer_dg(pseq)
            features[pseq]["dG"] = round(dG, 1)

            if dG >= min_dg:
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

    # Probe-specific parameters
    len_min = int(params.get("PROBE_LEN_MIN", 24))
    len_max = int(params.get("PROBE_LEN_MAX", 28))
    min_tm = float(params.get("PROBE_TM_MIN", 65))
    max_tm = float(params.get("PROBE_TM_MAX", 70))
    homopolymer_max = int(params.get("PROBE_HOMOPOLYMER_MAX", 3))
    avoid_5prime_G = params.get("PROBE_AVOID_5PRIME_G", "True").lower() in ("true", "1", "yes")

    # Shared parameters
    max_gc = float(params.get("GC_MAX", 60))
    min_dg = float(params.get("DG_MIN", -6))
    max_num = int(params.get("MAX_PRIMER_CANDIDATES", 10000))

    target_seqs = [str(s.seq) for s in SeqIO.parse(args.target_seqs, "fasta")]

    print(f"Generating probes from {args.target_seqs}...")
    print(f"Probe length range: {len_min}-{len_max}")
    print(f"Tm range: {min_tm}-{max_tm}")
    print(f"Using step=1 for comprehensive tiling")
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
    features_df = pd.DataFrame(features).T
    if features_df.empty:
        features_df = pd.DataFrame(columns=["pname", "pseq", "len", "Tm", "GC", "dG"])
    else:
        features_df = features_df.reset_index(names="pseq")[["pname", "pseq", "len", "Tm", "GC", "dG"]]

    fname = args.probe_seqs.replace(".fa", ".feat")
    features_df.to_csv(fname, index=False)

    runtime = time.time() - start_time
    print(f"Wrote {len(probe_list)} probes to {args.probe_seqs} ({runtime:.1f} sec)")
