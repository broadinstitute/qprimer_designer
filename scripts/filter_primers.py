#!/usr/bin/env python
# coding: utf-8

"""
filter_primer.py

Purpose
-------
Select top-performing primers based on evaluation results and write
their sequences to a FASTA file.

This script:
- Reads ranked primer-pair scores (CSV)
- Selects the top N primer pairs
- Extracts unique forward + reverse primer IDs
- Writes corresponding sequences from the initial primer FASTA

Key assumptions
---------------
- Evaluation results are sorted by descending score
- Primer IDs in the score table exist in the initial FASTA
- NUM_TOP_SENSITIVITY is defined in params.txt
"""

# ------------------------------------------------------------
# Imports
# ------------------------------------------------------------

import argparse
import pandas as pd
from Bio import SeqIO


# ------------------------------------------------------------
# Parameter parsing
# ------------------------------------------------------------

def parse_params(paramFile: str) -> int:
    """
    Extract NUM_TOP_SENSITIVITY from params.txt.
    """
    params = {}
    for l in open(paramFile):
        if "=" in l:
            name, value = l.split("=", 1)
            try:
                params[name.strip()] = float(value.strip())
            except ValueError:
                params[name.strip()] = value.strip()
    return int(params["NUM_TOP_SENSITIVITY"])


# ------------------------------------------------------------
# Main
# ------------------------------------------------------------

def main():
    parser = argparse.ArgumentParser(
        description="Filter top primers from evaluation results"
    )
    parser.add_argument(
        "--scores",
        required=True,
        help="Evaluation CSV (output of evaluate_primers.py)",
    )
    parser.add_argument(
        "--init",
        required=True,
        help="Initial primer FASTA",
    )
    parser.add_argument(
        "--out",
        required=True,
        help="Output FASTA of passed primers",
    )
    parser.add_argument(
        "--params",
        required=True,
        help="Parameters file",
    )
    args = parser.parse_args()

    numSelect = parse_params(args.params)

    # Load evaluation results
    res = pd.read_csv(args.scores)
    if res.empty:
        open(args.out, "w").close()
        return

    # Load initial primer sequences
    primer_seqs = {
        rec.id: str(rec.seq)
        for rec in SeqIO.parse(args.init, "fasta")
    }

    # Select top N primer pairs
    top = res.iloc[:numSelect]

    # Collect unique primer names (forward + reverse)
    primer_names = set(top["pname_f"]) | set(top["pname_r"])

    # Write FASTA
    with open(args.out, "w") as out:
        for pname in primer_names:
            if pname not in primer_seqs:
                raise KeyError(f"Primer {pname} not found in {args.init}")
            out.write(f">{pname}\n{primer_seqs[pname]}\n")


if __name__ == "__main__":
    main()

