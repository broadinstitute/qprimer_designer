#!/usr/bin/env python3
# coding: utf-8

"""
prepare_features.py

Purpose
-------
Compute primer features (.feat) for a given primer FASTA,
without generating new primers.

Input
-----
- FASTA file containing primers (*_for, *_rev)

Output
------
- Feature table (.feat):
  pname, pseq, forrev, len, Tm, GC, dG
"""

import argparse
import os
import re
import subprocess
import pandas as pd
from Bio import SeqIO
from Bio.SeqUtils import MeltingTemp, gc_fraction
from Bio.Seq import Seq

def get_Tm(seq: str) -> float:
    return MeltingTemp.Tm_NN(seq, Na=50, Mg=1.5, dNTPs=0.6)


def get_deltaG_vienna(seq1: str, seq2: str, RNAduplex: str) -> float:
    inp = f"{seq1}\n{seq2}\n"
    res = subprocess.run(
        [RNAduplex, "--noconv", "--paramFile=DNA"],
        input=inp,
        capture_output=True,
        text=True,
        check=True,
    )
    m = re.search(r"\(\s*([-+]?\d*\.?\d+)\s*\)", res.stdout)
    if not m:
        raise ValueError(f"Could not parse ΔG from RNAduplex output:\n{res.stdout}")
    return float(m.group(1))

def main():
    parser = argparse.ArgumentParser(
        description="Compute primer features for evaluation primers"
    )
    parser.add_argument("--fa", required=True, help="Primer FASTA")
    parser.add_argument("--out", required=True, help="Output .feat file")
    parser.add_argument("--params", required=True)
    parser.add_argument("--program", required=True, help="ViennaRNA directory")
    args = parser.parse_args()

    RNAduplex = f"{args.program}/src/bin/RNAduplex"
    if not os.path.exists(RNAduplex):
        raise FileNotFoundError(f"RNAduplex not found: {RNAduplex}")

    records = list(SeqIO.parse(args.fa, "fasta"))
    if not records:
        raise ValueError(f"No primers found in {args.fa}")

    rows = []

    for rec in records:
        pname = rec.id
        pseq = str(rec.seq)

        if pname.endswith("_for"):
            forrev = "f"
        elif pname.endswith("_rev"):
            forrev = "r"
        else:
            raise ValueError(f"Invalid primer name: {pname}")

        tm = get_Tm(pseq)
        gc = gc_fraction(pseq)
        dg = get_deltaG_vienna(pseq, pseq, RNAduplex)

        rows.append({
            "pname": pname,
            "pseq": pseq,
            "forrev": forrev,
            "len": len(pseq),
            "Tm": round(tm, 1),
            "GC": gc,
            "dG": round(dg, 1),
        })

    df = pd.DataFrame(rows)
    df.to_csv(args.out, index=False)

    print(f"Wrote {len(df)} primer features to {args.out}")

if __name__ == "__main__":
    main()

