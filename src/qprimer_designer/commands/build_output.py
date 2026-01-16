"""Build final primer design output."""

import argparse
import os
import time

import pandas as pd
from Bio import SeqIO

from qprimer_designer.utils import parse_params
from qprimer_designer.external import compute_dimer_dg


def register(subparsers):
    """Register the build-output subcommand."""
    parser = subparsers.add_parser(
        "build-output",
        help="Build final primer designs",
        description="""
Combine ON-target and OFF-target evaluation results into a final
primer-pair ranking, applying additional filtering based on
primer-primer dimerization (ΔG).
""",
    )
    parser.add_argument("--on", dest="eval_on", required=True, help="ON-target evaluation CSV")
    parser.add_argument("--off", dest="eval_off", nargs="+", required=True, help="OFF-target evaluation CSVs")
    parser.add_argument("--fa", dest="primers", required=True, help="Primer FASTA")
    parser.add_argument("--out", dest="output", required=True, help="Output CSV")
    parser.add_argument("--name", required=True, help="Target name")
    parser.add_argument("--params", dest="param_file", required=True, help="Parameters file")
    parser.set_defaults(func=run)


def run(args):
    """Run the build-output command."""
    params = parse_params(args.param_file)
    min_dg = float(params.get("DG_MIN", -8))
    num_select = int(params.get("NUM_TOP_SENSITIVITY", 100))

    print(f"Building final output for {args.name}...")
    start_time = time.time()

    # Load ON-target evaluation
    teval = pd.read_csv(args.eval_on, index_col=[0, 1]).iloc[:num_select].copy()
    teval.columns = ["cov_target", "act_target", "sco_target"]
    merged = teval.copy()

    # Merge OFF-target evaluations
    for off_path in args.eval_off:
        oname = os.path.basename(off_path).split(".")[1]
        oeval = pd.read_csv(off_path)

        tmp = pd.DataFrame(
            index=teval.index,
            columns=[f"cov_{oname}", f"act_{oname}", f"sco_{oname}"],
        )

        for pair in teval.index:
            sub = oeval[(oeval["pname_f"].isin(pair)) & (oeval["pname_r"].isin(pair))]
            if sub.empty:
                tmp.loc[pair] = 0
            else:
                tmp.loc[pair, f"cov_{oname}"] = sub["coverage"].sum()
                tmp.loc[pair, f"act_{oname}"] = sub["activity"].max()
                tmp.loc[pair, f"sco_{oname}"] = sub["score"].sum()

        merged = merged.join(tmp)

    # Primer-primer dimer filtering
    primer_seqs = {s.id: str(s.seq) for s in SeqIO.parse(args.primers, "fasta")}

    for fname, rname in merged.index:
        fseq = primer_seqs[fname]
        rseq = primer_seqs[rname]
        merged.loc[(fname, rname), "Dimer_dG"] = compute_dimer_dg(fseq, rseq)
        merged.loc[(fname, rname), "pseq_f"] = fseq
        merged.loc[(fname, rname), "pseq_r"] = rseq

    final = merged[merged["Dimer_dG"] > min_dg]

    # Write output
    final.round(3).to_csv(args.output)

    runtime = time.time() - start_time
    print(f"Wrote {len(final)} primer pairs to {args.output} ({runtime:.1f} sec)")
