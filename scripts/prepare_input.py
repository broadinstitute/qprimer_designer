#!/usr/bin/env python
# coding: utf-8

# In[ ]:


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

# In[ ]:


def reverse_complement_dna(seq):
    return str(Seq.Seq(seq).reverse_complement())


def get_Tm(seq):
    return MeltingTemp.Tm_NN(seq, Na=50, Mg=1.5, dNTPs=.6) # Primer3 condition


def _ensure_list(x):
    if isinstance(x, list):
        return x
    if pd.isna(x):
        return []
    return [x]


def rev_com_enc(enc):
    revcoms = {'A':'T',
               'T':'A',
               'C':'G',
               'G':'C',
               'a':'t',
               't':'a',
               'c':'g',
               'g':'c',
               '-':'-',
               'N':'N',
               'n':'n'
              }
    return ''.join([revcoms[char] for char in enc[::-1]])


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
    return primerLen, minAmpLen, maxAmpLen, minOffLen, maxOffLen


# In[ ]:


def main():
    parser = argparse.ArgumentParser(prog='python -u prepare_input.py', description='Prepare input for ML')
    parser.add_argument('--in', dest='mapped', required=True, help='Parsed alignments')
    parser.add_argument('--out', dest='ml_input', required=True, help='Input for ML')
    parser.add_argument('--ref', dest='reference', required=True, help='Target FASTA file')
    parser.add_argument('--params', dest='param_file', required=True, help='Parameters file')
    parser.add_argument('--target', dest='target', required=True, help='on or off')
    args = parser.parse_args()
    primerLen, minAmpLen, maxAmpLen, minOffLen, maxOffLen = parse_params(args.param_file)
    
    print(f'Preparing ML input from {args.mapped}...')
    startTime = time.time()
    nlines = 0
    tarseqs = { s.id:str(s.seq) for s in SeqIO.parse(args.reference, 'fasta') }

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

    grp = raw.groupby(['pname','pseq','tseq'])
    maptbl = grp.first()
    maptbl['tnames'] = grp['tname'].apply(list)
    maptbl['starts'] = grp['start'].apply(list)
    maptbl['orientation'] = maptbl['orientation'].apply(lambda x: x % 256)
    maptbl.reset_index(inplace=True)
    maptbl['forrev'] = maptbl['pname'].apply(lambda x: x.split('_')[-1])
    maptbl['pgap'] = maptbl['pseq'].apply(lambda x: x.count('-'))
    maptbl['tgap'] = maptbl['tseq'].apply(lambda x: x.count('-'))
    maptbl['indel'] = maptbl['pgap'] + maptbl['tgap']
    maptbl['mm'] = maptbl['pseq'].apply(len) - maptbl['match'].apply(lambda x: x.count('|')) - maptbl['indel']
    maptbl['pseq_raw'] = maptbl['pseq'].apply(lambda x: x.replace('-',''))
    maptbl['len'] = maptbl['pseq_raw'].apply(len)
    maptbl['Tm'] = maptbl['pseq_raw'].apply(get_Tm)
    maptbl['GC'] = maptbl['pseq_raw'].apply(gc_fraction)

    subcols = ['pname','pseq','pseq_raw','tseq','orientation','mm','indel','len','Tm','GC','tnames','starts']
    if args.target in ['on','On','ON']:
        fors = maptbl.loc[(maptbl['orientation']==0)&(maptbl['forrev']=='f'), subcols]
        revs = maptbl.loc[(maptbl['orientation']==16)&(maptbl['forrev']=='r'), subcols]
        minl = minAmpLen
        maxl = maxAmpLen
        lfunc = max
    else:
        fors = maptbl.loc[(maptbl['orientation']==0), subcols]
        revs = maptbl.loc[(maptbl['orientation']==16), subcols]
        minl = minOffLen
        maxl = maxOffLen
        lfunc = min
        
    revs['pseq'] = revs['pseq'].apply(rev_com_enc)
    revs['tseq'] = revs['tseq'].apply(rev_com_enc)

    revs_meta = revs.reset_index().rename(columns={"index": "r_id"})
    rev_index = defaultdict(list)
    for r in revs_meta.itertuples(index=False):
        tnames = _ensure_list(getattr(r, "tnames"))
        starts = _ensure_list(getattr(r, "starts"))
        r_id   = getattr(r, "r_id")
        for t, st in zip(tnames, starts):
            rev_index[t].append((r_id, int(st)))

    header_flag = True
    for f in fors.itertuples():
        f_id    = f.Index
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
        revsub['prod_Tm'] = [np.average(tms_by_r[r_id]) for r_id in revsub.index]
        forsub = fors.loc[[f_id]].copy().drop(['tnames','starts'],axis=1)
        pairs = forsub.merge(revsub.drop(['tnames','starts'],axis=1), how="cross", suffixes=("_f", "_r"))
        nlines += len(pairs)
        if os.path.exists(args.ml_input) and header_flag:
            mode = 'w'
        else:
            mode = 'a'
        pairs.to_csv(args.ml_input, mode=mode, index=False, header=header_flag)
        header_flag = False
    
    if not os.path.exists(args.ml_input):
        Path(args.ml_input).touch()

    runtime = (time.time() - startTime)
    print(f'Wrote {nlines} lines to {args.ml_input} (%.1f sec)' % runtime)

if __name__ == '__main__':
    main()

