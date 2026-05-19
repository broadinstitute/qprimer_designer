"""Compute primer features from existing primer sequences."""

import argparse

import pandas as pd
from Bio import SeqIO
from Bio.SeqUtils import gc_fraction

from qprimer_designer.utils import get_tm
from qprimer_designer.external import compute_batch_dimer_dg


def register(subparsers):
    """Register the prepare-features subcommand."""
    parser = subparsers.add_parser(
        "prepare-features",
        help="Compute features for existing primers",
        description="""
Compute primer features (Tm, GC%, dG) for a given primer FASTA file.
This is useful when you want to use your own primers instead of
generating them with the 'generate' command.
""",
    )
    parser.add_argument("--fa", required=True, help="Input primer FASTA file")
    parser.add_argument("--out", required=True, help="Output feature table (.feat)")
    parser.set_defaults(func=run)


def infer_primer_orientation(pname: str) -> str:
    """
    Infer primer orientation from name.

    Supports:
    - *_f or *_for → forward
    - *_r or *_rev → reverse
    """
    pname_lower = pname.lower()
    if pname_lower.endswith("_f") or pname_lower.endswith("_for"):
        return "f"
    elif pname_lower.endswith("_r") or pname_lower.endswith("_rev"):
        return "r"
    elif pname_lower.endswith("_p") or pname_lower.endswith("_pro"):
        return "p"
    else:
        raise ValueError(
            f"Cannot infer primer orientation from name: {pname}\n"
            "Expected suffix: _f, _for, _r, _rev, _p, or _pro"
        )


def run(args):
    """Run the prepare-features command."""
    records = list(SeqIO.parse(args.fa, "fasta"))
    if not records:
        raise ValueError(f"No primers found in {args.fa}")

    print(f"Computing features for {len(records)} primers from {args.fa}...")

    # Compute all features, batch dG in a single subprocess
    primer_data = []
    for rec in records:
        pname = rec.id
        pseq = str(rec.seq)
        forrev = infer_primer_orientation(pname)
        tm = get_tm(pseq)
        gc = gc_fraction(pseq)
        primer_data.append((pname, pseq, forrev, tm, gc))

    pairs = [(pseq, pseq) for _, pseq, _, _, _ in primer_data]
    dg_values = compute_batch_dimer_dg(pairs)

    rows = []
    for (pname, pseq, forrev, tm, gc), dg in zip(primer_data, dg_values):
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
