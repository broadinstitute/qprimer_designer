"""Generate candidate PCR primers from target sequences."""

import argparse
import bisect
import random
import time
from collections import defaultdict

import pandas as pd
from Bio import SeqIO
from Bio.SeqUtils import gc_fraction

from qprimer_designer.utils import reverse_complement_dna, get_tm, parse_params
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
    parser.add_argument("--in", dest="target_seqs", required=True, help="Input FASTA of target sequences")
    parser.add_argument("--out", dest="primer_seqs", required=True, help="Output FASTA for primers")
    parser.add_argument("--params", dest="param_file", required=True, help="Parameters file (params.txt)")
    parser.add_argument("--name", required=True, help="Name prefix for primer IDs")
    parser.set_defaults(func=run)


def generate_primers_single(target_seq: str, step: int, primer_len: int, min_amp_len: int):
    """Generate forward and reverse primer candidates from a single target."""
    target_seq_rc = reverse_complement_dna(target_seq)

    forps, revps = {}, {}
    for i in range(0, len(target_seq) - primer_len - min_amp_len, step):
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
    primer_len: int,
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
        flist, rlist = generate_primers_single(tseq, step, primer_len, min_amp_len)
        forps.update(flist)
        revps.update(rlist)

    uniq_pairs = count_primer_pairs(forps, revps, min_amp_len, max_amp_len)
    print(f">> Primers with a unique sequence: {len(forps)} forwards, {len(revps)} reverses, {uniq_pairs} pairs")

    for_filt, rev_filt = {}, {}
    features = defaultdict(dict)

    for unfilt, filt in zip([forps, revps], [for_filt, rev_filt]):
        for pseq in unfilt:
            tm = get_tm(pseq)
            gc = gc_fraction(pseq)

            features[pseq]["Tm"] = round(tm, 1)
            features[pseq]["GC"] = gc
            features[pseq]["len"] = len(pseq)

            if min_tm <= tm <= max_tm and gc <= max_gc / 100.0:
                dG = compute_self_dimer_dg(pseq)
                features[pseq]["dG"] = round(dG, 1)

                if dG >= min_dg:
                    filt[pseq] = unfilt[pseq]

    valid_pairs = count_primer_pairs(for_filt, rev_filt, min_amp_len, max_amp_len)
    print(f">> Primers filtered: {len(for_filt)} forwards, {len(rev_filt)} reverses, {valid_pairs} pairs")

    return for_filt, rev_filt, features


def run(args):
    """Run the generate command."""
    params = parse_params(args.param_file)

    max_num = int(params.get("MAX_PRIMER_CANDIDATES", 10000))
    step = int(params.get("TILING_STEP", 1))
    primer_len = int(params.get("PRIMER_LEN", 20))
    min_amp_len = int(params.get("AMPLEN_MIN", 60))
    max_amp_len = int(params.get("AMPLEN_MAX", 200))
    max_tm = float(params.get("TM_MAX", 60))
    min_tm = float(params.get("TM_MIN", 55))
    max_gc = float(params.get("GC_MAX", 60))
    min_dg = float(params.get("DG_MIN", -8))

    target_seqs = [str(s.seq) for s in SeqIO.parse(args.target_seqs, "fasta")]

    print(f"Generating primers from {args.target_seqs}...")
    start_time = time.time()

    for_filt, rev_filt, features = generate_primers_multi(
        target_seqs, step, primer_len, min_amp_len, max_amp_len,
        max_tm, min_tm, max_gc, min_dg,
    )

    forwards = list(for_filt.keys())
    reverses = list(rev_filt.keys())

    # Randomly subsample to keep output bounded
    if len(forwards) > max_num // 2:
        forwards = random.sample(forwards, max_num // 2)
    if len(reverses) > max_num // 2:
        reverses = random.sample(reverses, max_num // 2)

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

    features_df = pd.DataFrame(features).T
    if features_df.empty:
        features_df = pd.DataFrame(columns=["pname", "pseq", "forrev", "len", "Tm", "GC", "dG"])
    else:
        features_df = features_df.reset_index(names="pseq")[["pname", "pseq", "forrev", "len", "Tm", "GC", "dG"]]

    fname = args.primer_seqs.replace(".fa", ".feat")
    features_df.to_csv(fname, index=False)

    runtime = time.time() - start_time
    print(f"Wrote {len(forwards)+len(reverses)} primers to {args.primer_seqs} ({runtime:.1f} sec)")
