"""Filter top primers from evaluation results."""

import argparse

import pandas as pd
from Bio import SeqIO

from qprimer_designer.utils import parse_params


def register(subparsers):
    """Register the filter subcommand."""
    parser = subparsers.add_parser(
        "filter",
        help="Filter top primers from evaluation results",
        description="""
Select top-performing primers based on evaluation results and write
their sequences to a FASTA file.
""",
    )
    parser.add_argument("--scores", required=True, help="Evaluation CSV (output of evaluate)")
    parser.add_argument("--init", required=True, help="Initial primer FASTA")
    parser.add_argument("--out", required=True, help="Output FASTA of passed primers")
    parser.add_argument("--params", required=True, help="Parameters file")
    parser.set_defaults(func=run)


def run(args):
    """Run the filter command."""
    params = parse_params(args.params)
    num_select = int(params.get("NUM_TOP_SENSITIVITY", 100))

    # Load evaluation results
    res = pd.read_csv(args.scores)
    if res.empty:
        open(args.out, "w").close()
        return

    # Load initial primer sequences
    primer_seqs = {rec.id: str(rec.seq) for rec in SeqIO.parse(args.init, "fasta")}

    # Select top N primer pairs
    top = res.iloc[:num_select]

    # Collect unique primer names (forward + reverse)
    primer_names = set(top["pname_f"]) | set(top["pname_r"])

    # Write FASTA
    with open(args.out, "w") as out:
        for pname in primer_names:
            if pname not in primer_seqs:
                raise KeyError(f"Primer {pname} not found in {args.init}")
            out.write(f">{pname}\n{primer_seqs[pname]}\n")

    print(f"Wrote {len(primer_names)} primers to {args.out}")
