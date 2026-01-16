"""Prepare input for ML evaluation from mapped alignments."""

import argparse
import sys
import time
from collections import defaultdict
from pathlib import Path

import numpy as np
import pandas as pd
from Bio import SeqIO
from pandas.errors import EmptyDataError

from qprimer_designer.utils import parse_params, get_tm, fast_reverse_complement


def register(subparsers):
    """Register the prepare-input subcommand."""
    parser = subparsers.add_parser(
        "prepare-input",
        help="Prepare input for ML evaluation",
        description="""
Convert mapped primer-target alignments into a feature table suitable
for ML evaluation. Pairs forward and reverse primers into candidate
amplicons and computes amplicon-level features.
""",
    )
    parser.add_argument("--in", dest="mapped", required=True, help="Mapped alignments file")
    parser.add_argument("--out", dest="ml_input", required=True, help="Output CSV for ML input")
    parser.add_argument("--ref", dest="reference", required=True, help="Reference FASTA")
    parser.add_argument("--params", dest="param_file", required=True, help="Parameters file")
    parser.add_argument("--reftype", dest="reftype", required=True, choices=["on", "off"], help="on-target or off-target")
    parser.add_argument("--features", dest="pri_features", required=True, help="Primer features CSV")
    parser.add_argument("--prev", default="", help="Previous evaluation file (for off-target restriction)")
    parser.set_defaults(func=run)


def _ensure_list(x):
    """Ensure aggregated values are always lists."""
    if isinstance(x, list):
        return x
    if pd.isna(x):
        return []
    return [x]


def run(args):
    """Run the prepare-input command."""
    params = parse_params(args.param_file)
    primer_len = int(params.get("PRIMER_LEN", 20))
    min_amp_len = int(params.get("AMPLEN_MIN", 60))
    max_amp_len = int(params.get("AMPLEN_MAX", 200))
    min_off_len = int(params.get("OFFLEN_MIN", 60))
    max_off_len = int(params.get("OFFLEN_MAX", 2000))
    num_select = int(params.get("NUM_TOP_SENSITIVITY", 100))

    print(f"Preparing ML input from {args.mapped}...")
    start_time = time.time()
    nlines = 0

    # Load reference sequences
    tarseqs = {s.id: str(s.seq) for s in SeqIO.parse(args.reference, "fasta")}

    # Primer feature table
    feats = pd.read_csv(args.pri_features, index_col=0)

    # Previous evaluation (OFF-target restriction)
    valid_pairs = None
    if args.prev:
        try:
            teval = pd.read_csv(args.prev).iloc[:num_select]
            arr = teval[['pname_f', 'pname_r']].to_numpy()
            comb = np.vstack([
                np.column_stack([arr[:, 0], arr[:, 0]]),
                np.column_stack([arr[:, 0], arr[:, 1]]),
                np.column_stack([arr[:, 1], arr[:, 0]]),
                np.column_stack([arr[:, 1], arr[:, 1]]),
            ])
            valid_pairs = pd.DataFrame(comb, columns=['pname_f', 'pname_r']).drop_duplicates(ignore_index=True)
        except (FileNotFoundError, EmptyDataError):
            valid_pairs = None

    # Load mapped alignments
    cols = ['pname', 'orientation', 'tname', 'start', 'pseq', 'tseq', 'match']
    try:
        raw = pd.read_table(args.mapped, sep='\t', names=cols)
        if raw.empty:
            Path(args.ml_input).touch()
            sys.exit()
    except EmptyDataError:
        Path(args.ml_input).touch()
        sys.exit()

    raw['orientation'] = raw['orientation'] % 256

    maptbl = (
        raw.groupby(['pname', 'pseq', 'tseq', 'orientation'], sort=False)
        .agg(
            tnames=('tname', list),
            starts=('start', list),
            match=('match', 'first'),
        )
        .reset_index()
    )

    maptbl['indel'] = maptbl['pseq'].str.count('-') + maptbl['tseq'].str.count('-')
    maptbl['mm'] = maptbl['pseq'].str.len() - maptbl['match'].str.count(r'\|') - maptbl['indel']
    maptbl = maptbl.join(feats[['forrev', 'len', 'Tm', 'GC']], on='pname')

    drop_cols = ['orientation', 'forrev', 'match']

    if args.reftype == 'on':
        fors = maptbl[(maptbl['orientation'] == 0) & (maptbl['forrev'] == 'f')].drop(columns=drop_cols)
        revs = maptbl[(maptbl['orientation'] == 16) & (maptbl['forrev'] == 'r')].drop(columns=drop_cols)
        minl, maxl = min_amp_len, max_amp_len
        lfunc = max
    else:
        fors = maptbl[maptbl['orientation'] == 0].drop(columns=drop_cols)
        revs = maptbl[maptbl['orientation'] == 16].drop(columns=drop_cols)
        minl, maxl = min_off_len, max_off_len
        lfunc = min

    revs['pseq'] = revs['pseq'].apply(fast_reverse_complement)
    revs['tseq'] = revs['tseq'].apply(fast_reverse_complement)

    revs_meta = revs.reset_index().rename(columns={'index': 'r_id'})
    rev_index = defaultdict(list)

    for r in revs_meta.itertuples(index=False):
        for t, st in zip(_ensure_list(r.tnames), _ensure_list(r.starts)):
            rev_index[t].append((r.r_id, int(st)))

    if args.reftype == 'off' and valid_pairs is not None:
        allowed_r_by_f = valid_pairs.groupby('pname_f')['pname_r'].agg(set).to_dict()
        rname_by_id = revs_meta.set_index('r_id')['pname']

    allpairs = []

    for f in fors.itertuples():
        f_id = f.Index
        fname = f.pname

        tnamesf = _ensure_list(f.tnames)
        startsf = _ensure_list(f.starts)
        if not tnamesf:
            continue

        targets_by_r = defaultdict(set)
        amplens_by_r = defaultdict(set)
        tms_by_r = defaultdict(list)

        for t_f, st_f in zip(tnamesf, startsf):
            cand = rev_index.get(t_f)
            if not cand:
                continue

            st_f = int(st_f)
            for r_id, st_r in cand:
                if args.reftype == 'off' and valid_pairs is not None:
                    if rname_by_id[r_id] not in allowed_r_by_f.get(fname, set()):
                        continue

                ampseq = tarseqs[t_f][st_f - 1: st_r - 1 + primer_len]
                if minl <= len(ampseq) <= maxl:
                    targets_by_r[r_id].add(t_f)
                    amplens_by_r[r_id].add(len(ampseq))
                    tms_by_r[r_id].append(get_tm(ampseq))

        if not targets_by_r:
            continue

        rev_ids = list(targets_by_r.keys())
        revsub = revs.loc[rev_ids].copy()
        revsub['targets'] = [sorted(targets_by_r[r]) for r in revsub.index]
        revsub['prod_len'] = [lfunc(amplens_by_r[r]) for r in revsub.index]
        revsub['prod_Tm'] = [round(np.mean(tms_by_r[r]), 1) for r in revsub.index]

        forsub = fors.loc[[f_id]].copy().drop(['tnames', 'starts'], axis=1)
        pairs = forsub.merge(
            revsub.drop(['tnames', 'starts'], axis=1),
            how='cross',
            suffixes=('_f', '_r'),
        )

        nlines += len(pairs)
        allpairs.append(pairs)

    if not allpairs:
        Path(args.ml_input).touch()
    else:
        pd.concat(allpairs, ignore_index=True).to_csv(args.ml_input, index=False)

    runtime = time.time() - start_time
    print(f"Wrote {nlines} lines to {args.ml_input} ({runtime:.1f} sec)")
