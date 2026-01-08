#!/usr/bin/env python
# coding: utf-8

"""
prepare_input.py

Purpose
-------
Convert mapped primer–target alignments into a feature table suitable
for ML evaluation. This script pairs forward and reverse primers into
candidate amplicons and computes amplicon-level features.

Key assumptions
---------------
- Alignment coordinates (`start`) are 1-based.
- Orientation 0  = forward strand
- Orientation 16 = reverse strand
- Primer feature table is indexed by primer name (pname).
- For ON-target runs, only (F,R) primer orientations are allowed.
- For OFF-target runs, allowed primer pairs are restricted by the prev input.
"""

import argparse
import time
import pandas as pd
import numpy as np
import os
import sys
from Bio.SeqUtils import MeltingTemp
from Bio import SeqIO
from collections import defaultdict
from pathlib import Path
from pandas.errors import EmptyDataError


def get_Tm(seq):
    """Primer3-like melting temperature."""
    return MeltingTemp.Tm_NN(seq, Na=50, Mg=1.5, dNTPs=0.6)


def _ensure_list(x):
    """Ensure aggregated values are always lists."""
    if isinstance(x, list):
        return x
    if pd.isna(x):
        return []
    return [x]


def parse_params(paramFile):
    """Parse params.txt with minimal assumptions."""
    params = {}
    for l in open(paramFile):
        if '=' in l:
            name, value = l.split('=', 1)
            try:
                params[name.strip()] = float(value.strip())
            except ValueError:
                params[name.strip()] = value.strip()

    primerLen = int(params['PRIMER_LEN'])
    minAmpLen = int(params['AMPLEN_MIN'])
    maxAmpLen = int(params['AMPLEN_MAX'])
    minOffLen = int(params['OFFLEN_MIN'])
    maxOffLen = int(params['OFFLEN_MAX'])
    numSelect = int(params['NUM_TOP_SENSITIVITY'])

    return primerLen, minAmpLen, maxAmpLen, minOffLen, maxOffLen, numSelect


# Fast reverse-complement translation
rev_table = str.maketrans({
    'A':'T','T':'A','C':'G','G':'C',
    'a':'t','t':'a','c':'g','g':'c',
    '-':'-','N':'N','n':'n'
})


def main():
    parser = argparse.ArgumentParser(
        prog='python -u prepare_input.py',
        description='Prepare input for ML'
    )

    # === Snakefile-compatible arguments ===
    parser.add_argument('--in', dest='mapped', required=True)
    parser.add_argument('--out', dest='ml_input', required=True)
    parser.add_argument('--ref', dest='reference', required=True)
    parser.add_argument('--params', dest='param_file', required=True)
    parser.add_argument('--reftype', dest='reftype', required=True, choices=['on','off'])
    parser.add_argument('--features', dest='pri_features', required=True)
    parser.add_argument('--prev', dest='prev', default='')

    args = parser.parse_args()

    primerLen, minAmpLen, maxAmpLen, minOffLen, maxOffLen, numSelect = parse_params(args.param_file)

    print(f'Preparing ML input from {args.mapped}...')
    startTime = time.time()
    nlines = 0

    # Load reference sequences
    tarseqs = {s.id: str(s.seq) for s in SeqIO.parse(args.reference, 'fasta')}

    # Primer feature table
    feats = pd.read_csv(args.pri_features, index_col=0)

    # === Previous evaluation (OFF-target restriction) ===
    validPairs = None
    if args.prev:
        try:
            teval = pd.read_csv(args.prev).iloc[:numSelect]
            arr = teval[['pname_f', 'pname_r']].to_numpy()
            comb = np.vstack([
                np.column_stack([arr[:,0], arr[:,0]]),
                np.column_stack([arr[:,0], arr[:,1]]),
                np.column_stack([arr[:,1], arr[:,0]]),
                np.column_stack([arr[:,1], arr[:,1]])
            ])
            validPairs = (
                pd.DataFrame(comb, columns=['pname_f','pname_r'])
                .drop_duplicates(ignore_index=True)
            )
        except (FileNotFoundError, EmptyDataError):
            validPairs = None

    # === Load mapped alignments ===
    cols = ['pname','orientation','tname','start','pseq','tseq','match']
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
        raw.groupby(['pname','pseq','tseq','orientation'], sort=False)
           .agg(
               tnames=('tname', list),
               starts=('start', list),
               match=('match','first')
           )
           .reset_index()
    )

    maptbl['indel'] = maptbl['pseq'].str.count('-') + maptbl['tseq'].str.count('-')
    maptbl['mm'] = (
        maptbl['pseq'].str.len()
        - maptbl['match'].str.count(r'\|')
        - maptbl['indel']
    )

    maptbl = maptbl.join(feats[['forrev','len','Tm','GC']], on='pname')

    dropCols = ['orientation','forrev','match']

    if args.reftype == 'on':
        fors = maptbl[(maptbl['orientation']==0) & (maptbl['forrev']=='f')].drop(columns=dropCols)
        revs = maptbl[(maptbl['orientation']==16) & (maptbl['forrev']=='r')].drop(columns=dropCols)
        minl, maxl = minAmpLen, maxAmpLen
        lfunc = max
    else:
        fors = maptbl[maptbl['orientation']==0].drop(columns=dropCols)
        revs = maptbl[maptbl['orientation']==16].drop(columns=dropCols)
        minl, maxl = minOffLen, maxOffLen
        lfunc = min

    revs['pseq'] = revs['pseq'].str.translate(rev_table).str[::-1]
    revs['tseq'] = revs['tseq'].str.translate(rev_table).str[::-1]

    revs_meta = revs.reset_index().rename(columns={'index':'r_id'})
    rev_index = defaultdict(list)

    for r in revs_meta.itertuples(index=False):
        for t, st in zip(_ensure_list(r.tnames), _ensure_list(r.starts)):
            rev_index[t].append((r.r_id, int(st)))

    if args.reftype == 'off' and validPairs is not None:
        allowed_r_by_f = validPairs.groupby('pname_f')['pname_r'].agg(set).to_dict()
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
                if args.reftype == 'off' and validPairs is not None:
                    if rname_by_id[r_id] not in allowed_r_by_f.get(fname, set()):
                        continue

                # Coordinates assumed 1-based
                ampseq = tarseqs[t_f][st_f-1 : st_r-1 + primerLen]
                if minl <= len(ampseq) <= maxl:
                    targets_by_r[r_id].add(t_f)
                    amplens_by_r[r_id].add(len(ampseq))
                    tms_by_r[r_id].append(get_Tm(ampseq))

        if not targets_by_r:
            continue

        rev_ids = list(targets_by_r.keys())
        revsub = revs.loc[rev_ids].copy()
        revsub['targets'] = [sorted(targets_by_r[r]) for r in revsub.index]
        revsub['prod_len'] = [lfunc(amplens_by_r[r]) for r in revsub.index]
        revsub['prod_Tm'] = [round(np.mean(tms_by_r[r]),1) for r in revsub.index]

        forsub = fors.loc[[f_id]].copy().drop(['tnames','starts'], axis=1)
        pairs = forsub.merge(
            revsub.drop(['tnames','starts'], axis=1),
            how='cross',
            suffixes=('_f','_r')
        )

        nlines += len(pairs)
        allpairs.append(pairs)

    if not allpairs:
        Path(args.ml_input).touch()
    else:
        pd.concat(allpairs, ignore_index=True).to_csv(args.ml_input, index=False)

    runtime = time.time() - startTime
    print(f'Wrote {nlines} lines to {args.ml_input} (%.1f sec)' % runtime)


if __name__ == '__main__':
    main()

