"""Generate candidate PCR primers from target sequences."""

import argparse
import bisect
import random
import time
from collections import defaultdict

import pandas as pd
from Bio import SeqIO
from Bio.SeqUtils import gc_fraction

from qprimer_designer.utils import reverse_complement_dna, get_tm, parse_params, get_primer_params
from qprimer_designer.external import compute_self_dimer_dg


def register(subparsers):
    """Register the generate subcommand."""
    parser = subparsers.add_parser(
        "generate",
        help="Generate candidate primers from target sequences",
        description="""
Generate candidate forward and reverse PCR primers from target sequences
by tiling, then filtering by melting temperature (Tm), GC content, and
self-dimer free energy (ΔG).
""",
    )
    parser.add_argument("--in", dest="target_seqs", required=True, 
                        help="Input FASTA of target sequences")
    parser.add_argument("--out", dest="primer_seqs", required=True, 
                        help="Output FASTA for primers")
    parser.add_argument("--params", dest="param_file", required=True, 
                        help="Parameters file (params.txt)")
    parser.add_argument("--name", required=True, 
                        help="Name prefix for primer IDs")
    parser.set_defaults(func=run)


def generate_primers_single(target_seq: str, step: int, primer_lens: list[int], min_amp_len: int):
    """Generate forward and reverse primer candidates from a single target."""
    target_seq_rc = reverse_complement_dna(target_seq)

    forps, revps = {}, {}
    max_len = max(primer_lens)

    for i in range(0, len(target_seq) - max_len - min_amp_len, step):
        for primer_len in primer_lens:
            # Check if this length fits at this position
            if i + primer_len <= len(target_seq):
                fseq = target_seq[i : i + primer_len]
                if "N" not in fseq:
                    forps[fseq] = i

                rseq = target_seq_rc[i : i + primer_len]
                if "N" not in rseq:
                    revps[rseq] = len(target_seq) - i

    return forps, revps


def count_primer_pairs(fors: dict, revs: dict, min_amp_len: int, max_amp_len: int) -> int:
    """Count possible amplicons based on forward start and reverse end positions."""
    sts = sorted(fors.values())
    ens = sorted(revs.values())

    count = 0
    for st in sts:
        left = bisect.bisect_left(ens, st + min_amp_len)
        right = bisect.bisect_right(ens, st + max_amp_len)
        count += right - left

    return count


def generate_primers_multi(
    target_seqs,
    step: int,
    primer_lens: list[int],
    min_amp_len: int,
    max_amp_len: int,
    max_tm: float,
    min_tm: float,
    max_gc: float,
    min_dg: float,
):
    """Generate and filter primer candidates across multiple target sequences."""
    forps, revps = {}, {}
    for tseq in target_seqs:
        flist, rlist = generate_primers_single(tseq, step, primer_lens, min_amp_len)
        forps.update(flist)
        revps.update(rlist)

    uniq_pairs = count_primer_pairs(forps, revps, min_amp_len, max_amp_len)
    print(f">> Primers with a unique sequence: {len(forps)} forwards, {len(revps)} reverses, {uniq_pairs} pairs")

    for_filt, rev_filt = {}, {}
    features = defaultdict(dict)

    for unfilt, filt in zip([forps, revps], [for_filt, rev_filt]):
        for pseq in unfilt:
            # Calculate features
            tm = get_tm(pseq)
            gc = gc_fraction(pseq)

            # Filter by Tm and GC
            if min_tm <= tm <= max_tm and gc <= max_gc / 100.0:
                dG = compute_self_dimer_dg(pseq)

                # Only add features if all filters pass
                if dG >= min_dg:
                    features[pseq]["Tm"] = round(tm, 1)
                    features[pseq]["GC"] = round(gc, 2)
                    features[pseq]["len"] = len(pseq)
                    features[pseq]["dG"] = round(dG, 1)
                    filt[pseq] = unfilt[pseq]

    valid_pairs = count_primer_pairs(for_filt, rev_filt, min_amp_len, max_amp_len)
    print(f">> Primers filtered: {len(for_filt)} forwards, {len(rev_filt)} reverses, {valid_pairs} pairs")

    return for_filt, rev_filt, features


def run(args):
    """Run the generate command."""
    params = parse_params(args.param_file)
    primer_params = get_primer_params(params)

    max_num = primer_params["max_num"]
    step = primer_params["step"]
    primer_lens = primer_params["primer_lens"]
    min_amp_len = primer_params["min_amp_len"]
    max_amp_len = primer_params["max_amp_len"]
    max_tm = primer_params["max_tm"]
    min_tm = primer_params["min_tm"]
    max_gc = primer_params["max_gc"]
    min_dg = primer_params["min_dg"]

    target_seqs = [str(s.seq) for s in SeqIO.parse(args.target_seqs, "fasta")]

    print(f"Generating primers from {args.target_seqs}...")
    print(f"Primer lengths: {primer_lens}")
    print(f"Tm range: {min_tm}-{max_tm}")
    start_time = time.time()

    for_filt, rev_filt, features = generate_primers_multi(
        target_seqs, step, primer_lens, min_amp_len, max_amp_len,
        max_tm, min_tm, max_gc, min_dg,
    )

    forwards = list(for_filt.keys())
    reverses = list(rev_filt.keys())

    # Randomly subsample to keep output bounded
    if len(forwards) > max_num // 2:
        forwards = random.sample(forwards, max_num // 2)
    if len(reverses) > max_num // 2:
        reverses = random.sample(reverses, max_num // 2)

    # Write output FASTA
    with open(args.primer_seqs, "w") as fout:
        for i, seq in enumerate(forwards):
            pname = f"{args.name}_{i+1}_f"
            features[seq]["pname"] = pname
            features[seq]["forrev"] = "f"
            fout.write(f">{pname}\n{seq}\n")

        for i, seq in enumerate(reverses):
            pname = f"{args.name}_{i+1}_r"
            features[seq]["pname"] = pname
            features[seq]["forrev"] = "r"
            fout.write(f">{pname}\n{seq}\n")

    # Write features CSV
    if features:
        features_df = pd.DataFrame(features).T
        features_df = features_df.reset_index(names="pseq")[["pname", "pseq", "forrev", "len", "Tm", "GC", "dG"]]
    else:
        features_df = pd.DataFrame(columns=["pname", "pseq", "forrev", "len", "Tm", "GC", "dG"])

    fname = args.primer_seqs.replace(".fa", ".feat")
    features_df.to_csv(fname, index=False)

    runtime = time.time() - start_time
    print(f"Wrote {len(forwards)+len(reverses)} primers to {args.primer_seqs} ({runtime:.1f} sec)")
