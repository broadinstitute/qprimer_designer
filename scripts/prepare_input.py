#!/usr/bin/env python
# coding: utf-8


import argparse
import time
import pandas as pd
import numpy as np
import os
import sys
from Bio.SeqUtils import MeltingTemp, gc_fraction
from Bio import Seq, SeqIO
from collections import defaultdict
from pathlib import Path
from pandas.errors import EmptyDataError

def get_Tm(seq):
    return MeltingTemp.Tm_NN(seq, Na=50, Mg=1.5, dNTPs=.6) # Primer3 condition

def _ensure_list(x):
    if isinstance(x, list):
        return x
    if pd.isna(x):
        return []
    return [x]


def parse_params(paramFile):
    params = {}
    for l in open(paramFile):
        if '=' in l:
            name, value = l.split('=')
            try:
                params[name.strip()] = int(value.strip())
            except ValueError:
                params[name.strip()] = value.strip()
    primerLen = params['PRIMER_LEN']
    minAmpLen = params['AMPLEN_MIN']
    maxAmpLen = params['AMPLEN_MAX']
    minOffLen = params['OFFLEN_MIN']
    maxOffLen = params['OFFLEN_MAX']
    numSelect = params['NUM_TOP_SENSITIVITY']
    return primerLen, minAmpLen, maxAmpLen, minOffLen, maxOffLen, numSelect


rev_table = str.maketrans({
    'A':'T','T':'A','C':'G','G':'C',
    'a':'t','t':'a','c':'g','g':'c',
    '-':'-','N':'N','n':'n' })

def main():
    parser = argparse.ArgumentParser(prog='python -u prepare_input.py', description='Prepare input for ML')
    parser.add_argument('--in', dest='mapped', required=True, help='Parsed alignments')
    parser.add_argument('--out', dest='ml_input', required=True, help='Input for ML')
    parser.add_argument('--ref', dest='reference', required=True, help='Target FASTA file')
    parser.add_argument('--params', dest='param_file', required=True, help='Parameters file')
    parser.add_argument('--target', dest='target', required=True, help='on or off')
    parser.add_argument('--features', dest='pri_features', required=True, help='Table of primer features')
    parser.add_argument('--teval', dest='teval', default='-', help='.result file from evaluation of sensitivity')
    args = parser.parse_args()
    primerLen, minAmpLen, maxAmpLen, minOffLen, maxOffLen, numSelect = parse_params(args.param_file)
    
    print(f'Preparing ML input from {args.mapped}...')
    startTime = time.time()
    nlines = 0
    tarseqs = { s.id:str(s.seq) for s in SeqIO.parse(args.reference, 'fasta') }
    feats = pd.read_csv(args.pri_features, index_col=0) # index is 'pname'
    try:
        teval = pd.read_csv(args.teval).iloc[:numSelect]
        arr = teval[['pname_f', 'pname_r']].to_numpy()
        comb = np.vstack([
            np.column_stack([arr[:, 0], arr[:, 0]]),  # (F, F)
            np.column_stack([arr[:, 0], arr[:, 1]]),  # (F, R)
            np.column_stack([arr[:, 1], arr[:, 0]]),  # (R, F)
            np.column_stack([arr[:, 1], arr[:, 1]]),  # (R, R)
        ])
        validPairs = (pd.DataFrame(comb, columns=['pname_f', 'pname_r']).drop_duplicates(ignore_index=True))
    except FileNotFoundError:
        pass

    cols = ['pname','orientation','tname','start','pseq','tseq','match']
    try:    
        raw = pd.read_table(args.mapped, sep='\t', names=cols)
        if raw.empty:
            Path(args.ml_input).touch()
            print(f'No alignments in {args.mapped}.')
            sys.exit()
    except EmptyDataError:
        Path(args.ml_input).touch()
        print(f'No alignments in {args.mapped}.')
        sys.exit()
    
    raw['orientation'] = raw['orientation'] % 256
    maptbl = (raw.groupby(['pname','pseq','tseq','orientation'], sort=False)
                .agg(
                    tnames=('tname', list),
                    starts=('start', list),
                    match=('match', 'first')
                ).reset_index())
    
    maptbl['indel'] = maptbl['pseq'].str.count('-') + maptbl['tseq'].str.count('-')
    maptbl['mm'] = maptbl['pseq'].str.len() - maptbl['match'].str.count(r'\|') - maptbl['indel']
    maptbl = maptbl.join(feats[['forrev','len','Tm','GC']], on='pname')
    dropCols = ['orientation', 'forrev', 'match'] 
    if args.target in ['on','On','ON']:
        fors = maptbl.loc[(maptbl['orientation']==0) & (maptbl['forrev']=='f')].drop(columns=dropCols).copy()
        revs = maptbl.loc[(maptbl['orientation']==16) & (maptbl['forrev']=='r')].drop(columns=dropCols).copy()
        minl = minAmpLen
        maxl = maxAmpLen
        lfunc = max
    else:
        fors = maptbl.loc[(maptbl['orientation']==0)].drop(columns=dropCols).copy()
        revs = maptbl.loc[(maptbl['orientation']==16)].drop(columns=dropCols).copy()
        minl = minOffLen
        maxl = maxOffLen
        lfunc = min
    
    revs['pseq'] = revs['pseq'].str.translate(rev_table).str[::-1] #  .apply(rev_com_enc)
    revs['tseq'] = revs['tseq'].str.translate(rev_table).str[::-1] #  .apply(rev_com_enc)
    
    revs_meta = revs.reset_index().rename(columns={"index": "r_id"})
    rev_index = defaultdict(list)
    for r in revs_meta.itertuples(index=False):
        tnames = _ensure_list(getattr(r, "tnames"))
        starts = _ensure_list(getattr(r, "starts"))
        r_id   = getattr(r, "r_id")
        for t, st in zip(tnames, starts):
            rev_index[t].append((r_id, int(st)))

    if args.target in ['off','Off','OFF']:
        allowed_r_by_f = validPairs.groupby('pname_f')['pname_r'].agg(set).to_dict()
        rname_by_id = revs_meta.set_index('r_id')['pname']

    allpairs = []
    for f in fors.itertuples():
        f_id    = f.Index
        fname   = getattr(f, 'pname')
        tnamesf = _ensure_list(getattr(f, "tnames"))
        startsf = _ensure_list(getattr(f, "starts"))
        if not tnamesf or not startsf:
            continue

        targets_by_r = defaultdict(set)
        amplens_by_r = defaultdict(set)
        tms_by_r = defaultdict(list)
        for t_f, st_f in zip(tnamesf, startsf):
            cand = rev_index.get(t_f)
            st_f = int(st_f)
            if not cand:
                continue
            
            for r_id, st_r in cand:
                if args.target in ['off','Off','OFF']:
                    if rname_by_id[r_id] not in allowed_r_by_f.get(fname):
                        continue

                ampseq = tarseqs[t_f][st_f-1:st_r-1+primerLen]
                if minl <= len(ampseq) <= maxl:
                    targets_by_r[r_id].add(t_f)
                    amplens_by_r[r_id].add(len(ampseq))
                    tms_by_r[r_id].append(get_Tm(ampseq))

        if not targets_by_r:
            continue

        rev_ids = list(targets_by_r.keys())
        revsub = revs.loc[rev_ids].copy()
        revsub['targets'] = [sorted(list(targets_by_r[r_id])) for r_id in revsub.index]
        revsub['prod_len'] = [lfunc(list(amplens_by_r[r_id])) for r_id in revsub.index]
        revsub['prod_Tm'] = [round(np.average(tms_by_r[r_id]),1) for r_id in revsub.index]
        forsub = fors.loc[[f_id]].copy().drop(['tnames','starts'], axis=1)
        pairs = forsub.merge(revsub.drop(['tnames','starts'], axis=1), how="cross", suffixes=("_f", "_r"))
        nlines += len(pairs)
        allpairs.append(pairs)
        
    allpairs = pd.concat(allpairs, ignore_index=True)
    if allpairs.empty:
        Path(args.ml_input).touch()
    else:
        allpairs.to_csv(args.ml_input, index=False)

    runtime = (time.time() - startTime)
    print(f'Wrote {nlines} lines to {args.ml_input} (%.1f sec)' % runtime)

if __name__ == '__main__':
    main()

