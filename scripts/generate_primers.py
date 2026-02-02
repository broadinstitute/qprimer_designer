#!/usr/bin/env python
# coding: utf-8

"""
generate_primers.py

Purpose
-------
Generate candidate forward and reverse PCR primers from one or more
target sequences by tiling, then filtering by:

- Primer length
- Melting temperature (Tm; Primer3-like conditions)
- GC upper bound
- Self-dimer ΔG (ViennaRNA RNAduplex, DNA parameters)

Outputs
-------
1) FASTA file of primer sequences (forward + reverse)
2) Feature table (.feat) with primer metadata:
   pname, pseq, forrev, len, Tm, GC, dG

Key assumptions
---------------
- Primer sequences containing 'N' are discarded.
- Primer sequences are treated as UNIQUE by sequence string:
  if the same sequence occurs multiple times in the target,
  only one instance is retained.
- GC content is stored as a FRACTION (0–1), not percent.
"""

import argparse
import bisect
import time
import subprocess
import re
import random
import os

import pandas as pd
from Bio.SeqUtils import MeltingTemp, gc_fraction
from Bio import Seq, SeqIO
from collections import defaultdict


# ------------------------------------------------------------
# Helper functions
# ------------------------------------------------------------

def reverse_complement_dna(seq: str) -> str:
    """Return reverse-complement of a DNA sequence."""
    return str(Seq.Seq(seq).reverse_complement())


def get_Tm(seq: str) -> float:
    """Primer melting temperature (Primer3-like ionic conditions)."""
    return MeltingTemp.Tm_NN(seq, Na=50, Mg=1.5, dNTPs=0.6)


def get_deltaG_vienna(seq1: str, seq2: str, RNAduplex: str) -> float:
    """
    Compute duplex free energy (ΔG) using ViennaRNA RNAduplex
    in DNA parameter mode.
    """
    inp = f"{seq1}\n{seq2}\n"
    res = subprocess.run(
        [RNAduplex, "--noconv", "--paramFile=DNA"],
        input=inp,
        capture_output=True,
        text=True,
        check=True,
    )
    deltaG = re.search(r"\(\s*([-+]?\d*\.?\d+)\s*\)", res.stdout)
    if not deltaG:
        raise ValueError(
            f"Could not parse ΔG from RNAduplex output:\n{res.stdout}\n{seq1}"
        )
    return float(deltaG.group(1))


# ------------------------------------------------------------
# Primer generation logic
# ------------------------------------------------------------

def generate_primers_single(
    targetSeq: str,
    step: int,
    primerLen: int,
    minAmpLen: int,
):
    """
    Generate forward and reverse primer candidates from a single target.

    Notes
    -----
    - Primers are keyed by sequence string, enforcing uniqueness.
    - Forward primers store start positions.
    - Reverse primers store approximate end positions on the + strand.
    """
    targetSeq_rc = reverse_complement_dna(targetSeq)

    forps, revps = {}, {}
    for i in range(0, len(targetSeq) - primerLen - minAmpLen, step):
        fseq = targetSeq[i : i + primerLen]
        if "N" not in fseq:
            forps[fseq] = i

        rseq = targetSeq_rc[i : i + primerLen]
        if "N" not in rseq:
            revps[rseq] = len(targetSeq) - i

    return forps, revps


def count_primer_pairs(
    fors: dict,
    revs: dict,
    minAmpLen: int,
    maxAmpLen: int,
) -> int:
    """
    Count possible amplicons based on forward start and reverse end positions.
    """
    sts = sorted(fors.values())
    ens = sorted(revs.values())

    count = 0
    for st in sts:
        left = bisect.bisect_left(ens, st + minAmpLen)
        right = bisect.bisect_right(ens, st + maxAmpLen)
        count += (right - left)

    return count


def generate_primers_multi(
    targetSeqs,
    step: int,
    primerLen: int,
    minAmpLen: int,
    maxAmpLen: int,
    maxTm: float,
    minTm: float,
    maxGc: float,
    minDg: float,
    viennaDir: str,
):
    """
    Generate and filter primer candidates across multiple target sequences.
    """
    RNAduplex = f"{viennaDir}/src/bin/RNAduplex"
    if not os.path.exists(RNAduplex):
        raise FileNotFoundError(
            f"RNAduplex not found at: {RNAduplex}\n"
            "Check VIENNA_DIRECTORY in params.txt."
        )

    forps, revps = {}, {}
    for tseq in targetSeqs:
        flist, rlist = generate_primers_single(
            tseq, step, primerLen, minAmpLen
        )
        forps.update(flist)
        revps.update(rlist)

    uniqPairs = count_primer_pairs(
        forps, revps, minAmpLen, maxAmpLen
    )
    print(
        f">> Primers with a unique sequence: "
        f"{len(forps)} forwards, {len(revps)} reverses, {uniqPairs} pairs"
    )

    forFilt, revFilt = {}, {}
    features = defaultdict(dict)

    for unfilt, filt in zip([forps, revps], [forFilt, revFilt]):
        for pseq in unfilt:
            tm = get_Tm(pseq)
            gc = gc_fraction(pseq)  # fraction (0–1)

            features[pseq]["Tm"] = round(tm, 1)
            features[pseq]["GC"] = gc
            features[pseq]["len"] = len(pseq)

            if minTm <= tm <= maxTm and gc <= maxGc / 100.0:
                dG = get_deltaG_vienna(pseq, pseq, RNAduplex)
                features[pseq]["dG"] = round(dG, 1)

                if dG >= minDg:
                    filt[pseq] = unfilt[pseq]

    validPairs = count_primer_pairs(
        forFilt, revFilt, minAmpLen, maxAmpLen
    )
    print(
        f">> Primers filtered: "
        f"{len(forFilt)} forwards, {len(revFilt)} reverses, {validPairs} pairs"
    )

    return forFilt, revFilt, features


# ------------------------------------------------------------
# Parameter parsing
# ------------------------------------------------------------

def parse_params(paramFile: str):
    """
    Parse params.txt.

    Numeric values are parsed as float first, then cast where needed.
    """
    params = {}
    for l in open(paramFile):
        if "=" in l:
            name, value = l.split("=", 1)
            try:
                params[name.strip()] = float(value.strip())
            except ValueError:
                params[name.strip()] = value.strip()

    #viennaDir = params["VIENNA_DIRECTORY"]
    maxNum = int(params["MAX_PRIMER_CANDIDATES"])
    step = int(params["TILING_STEP"])
    primerLen = int(params["PRIMER_LEN"])
    minAmpLen = int(params["AMPLEN_MIN"])
    maxAmpLen = int(params["AMPLEN_MAX"])
    maxTm = float(params["TM_MAX"])
    minTm = float(params["TM_MIN"])
    maxGc = float(params["GC_MAX"])
    minDg = float(params["DG_MIN"])

    return (
        #viennaDir,
        maxNum,
        step,
        primerLen,
        minAmpLen,
        maxAmpLen,
        maxTm,
        minTm,
        maxGc,
        minDg,
    )


# ------------------------------------------------------------
# Main
# ------------------------------------------------------------

def main():
    parser = argparse.ArgumentParser(
        prog="python -u generate_primers.py",
        description="Generate primer candidates",
    )
    parser.add_argument("--in", dest="target_seqs", required=True)
    parser.add_argument("--out", dest="primer_seqs", required=True)
    parser.add_argument("--params", dest="param_file", required=True)
    parser.add_argument("--program", dest="program_path", required=True)
    parser.add_argument("--name", dest="name", required=True)
    args = parser.parse_args()

    viennaDir = args.program_path

    targetSeqs = [
        str(s.seq) for s in SeqIO.parse(args.target_seqs, "fasta")
    ]

    (
        #viennaDir,
        maxNum,
        step,
        primerLen,
        minAmpLen,
        maxAmpLen,
        maxTm,
        minTm,
        maxGc,
        minDg,
    ) = parse_params(args.param_file)

    print(f"Generating primers from {args.target_seqs}...")
    startTime = time.time()

    forFilt, revFilt, features = generate_primers_multi(
        targetSeqs,
        step,
        primerLen,
        minAmpLen,
        maxAmpLen,
        maxTm,
        minTm,
        maxGc,
        minDg,
        viennaDir,
    )

    forwards = list(forFilt.keys())
    reverses = list(revFilt.keys())

    # Randomly subsample to keep output bounded
    if len(forwards) > maxNum // 2:
        forwards = random.sample(forwards, maxNum // 2)
    if len(reverses) > maxNum // 2:
        reverses = random.sample(reverses, maxNum // 2)

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
        features_df = pd.DataFrame(
            columns=["pname","pseq","forrev","len","Tm","GC","dG"]
        )
    else:
        features_df = (
            features_df
            .reset_index(names="pseq")
            [["pname","pseq","forrev","len","Tm","GC","dG"]]
        )

    fname = args.primer_seqs.replace(".fa", ".feat")
    features_df.to_csv(fname, index=False)

    runtime = time.time() - startTime
    print(
        f"Wrote {len(forwards)+len(reverses)} primers to "
        f"{args.primer_seqs} (%.1f sec)" % runtime
    )


if __name__ == "__main__":
    main()

