#!/usr/bin/env python
# coding: utf-8

"""
build_final_output.py

Purpose
-------
Combine ON-target and OFF-target evaluation results into a final
primer-pair ranking, applying additional filtering based on
primer–primer dimerization (ΔG).

Outputs
-------
Final CSV containing:
- ON-target coverage / activity / score
- OFF-target coverage / activity / score (per off-target)
- Primer sequences
- Primer–primer dimer ΔG

Key assumptions
---------------
- Target evaluation file is sorted by descending score
- OFF-target evaluation files follow the same schema
- Primer FASTA contains all primer IDs referenced in evaluations
"""


# ------------------------------------------------------------
# Imports
# ------------------------------------------------------------

import argparse
import time
import subprocess
import re
import os

import pandas as pd
from Bio import SeqIO


# ------------------------------------------------------------
# Helper functions
# ------------------------------------------------------------

def get_deltaG_vienna(seq1: str, seq2: str, RNAduplex: str) -> float:
    """Compute primer–primer dimer ΔG using ViennaRNA RNAduplex (DNA mode)."""
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
            f"Could not parse ΔG from RNAduplex output:\n{res.stdout}"
        )
    return float(deltaG.group(1))


# ------------------------------------------------------------
# Parameter parsing
# ------------------------------------------------------------

def parse_params(paramFile: str):
    """Extract required parameters from params.txt."""
    params = {}
    for l in open(paramFile):
        if "=" in l:
            name, value = l.split("=", 1)
            try:
                params[name.strip()] = float(value.strip())
            except ValueError:
                params[name.strip()] = value.strip()

    viennaDir = params["VIENNA_DIRECTORY"]
    minDg = float(params["DG_MIN"])
    numSelect = int(params["NUM_TOP_SENSITIVITY"])

    return viennaDir, minDg, numSelect


# ------------------------------------------------------------
# Main
# ------------------------------------------------------------

def main():
    parser = argparse.ArgumentParser(
        prog="python -u build_final_output.py",
        description="Build final primer designs",
    )

    # === REQUIRED BY SNAKEFILE ===
    parser.add_argument("--on", dest="eval_on", required=True)
    parser.add_argument("--off", dest="eval_off", nargs="+", required=True)
    parser.add_argument("--fa", dest="primers", required=True)
    parser.add_argument("--out", dest="output", required=True)
    parser.add_argument("--name", required=True)
    parser.add_argument("--params", dest="param_file", required=True)

    args = parser.parse_args()

    viennaDir, minDg, numSelect = parse_params(args.param_file)

    print(f"Building final output for {args.name}...")
    startTime = time.time()

    # --------------------------------------------------------
    # Load ON-target evaluation
    # --------------------------------------------------------

    teval = (
        pd.read_csv(args.eval_on, index_col=[0, 1])
        .iloc[:numSelect]
        .copy()
    )
    teval.columns = ["cov_target", "act_target", "sco_target"]
    merged = teval.copy()

    # --------------------------------------------------------
    # Merge OFF-target evaluations
    # --------------------------------------------------------

    for off_path in args.eval_off:
        oname = os.path.basename(off_path).split(".")[1]
        oeval = pd.read_csv(off_path)

        tmp = pd.DataFrame(
            index=teval.index,
            columns=[
                f"cov_{oname}",
                f"act_{oname}",
                f"sco_{oname}",
            ],
        )

        for pair in teval.index:
            sub = oeval[
                (oeval["pname_f"].isin(pair))
                & (oeval["pname_r"].isin(pair))
            ]
            if sub.empty:
                tmp.loc[pair] = 0
            else:
                tmp.loc[pair, f"cov_{oname}"] = sub["coverage"].sum()
                tmp.loc[pair, f"act_{oname}"] = sub["activity"].max()
                tmp.loc[pair, f"sco_{oname}"] = sub["score"].sum()

        merged = merged.join(tmp)

    # --------------------------------------------------------
    # Primer–primer dimer filtering
    # --------------------------------------------------------

    RNAduplex = f"{viennaDir}/src/bin/RNAduplex"
    if not os.path.exists(RNAduplex):
        raise FileNotFoundError(
            f"RNAduplex not found at: {RNAduplex}\n"
            "Check VIENNA_DIRECTORY in params.txt."
        )

    primer_seqs = {
        s.id: str(s.seq)
        for s in SeqIO.parse(args.primers, "fasta")
    }

    for (fname, rname) in merged.index:
        fseq = primer_seqs[fname]
        rseq = primer_seqs[rname]
        merged.loc[(fname, rname), "Dimer_dG"] = get_deltaG_vienna(
            fseq, rseq, RNAduplex
        )
        merged.loc[(fname, rname), "pseq_f"] = fseq
        merged.loc[(fname, rname), "pseq_r"] = rseq

    final = merged[merged["Dimer_dG"] > minDg]

    # --------------------------------------------------------
    # Write output
    # --------------------------------------------------------

    final.round(3).to_csv(args.output)

    runtime = time.time() - startTime
    print(
        f"Wrote {len(final)} primer pairs to {args.output} "
        f"(%.1f sec)." % runtime
    )


if __name__ == "__main__":
    main()

