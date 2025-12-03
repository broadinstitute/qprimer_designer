#!/usr/bin/env python
# coding: utf-8

import argparse
import time
import subprocess
import re
import os
import pandas as pd
from Bio import SeqIO
from functools import reduce


def get_deltaG_vienna(seq1, seq2, RNAduplex):
    inp = f'{seq1}\n{seq2}\n'
    res = subprocess.run([RNAduplex, '--noconv', '--paramFile=DNA'], input=inp,
                         capture_output=True, text=True, check=True)
    deltaG = re.search(r'\(([-+]?\d*\.?\d+)\)', res.stdout)
    if not deltaG:
        raise ValueError(f'Could not parse ΔG from RNAduplex output:\n{res.stdout}')
    return float(deltaG.group(1))  


def parse_params(paramFile):
    params = {}
    for l in open(paramFile):
        if '=' in l:
            name, value = l.split('=')
            try:
                params[name.strip()] = int(value.strip())
            except ValueError:
                params[name.strip()] = value.strip()
    viennaDir = params['VIENNA_DIRECTORY']
    minDg = params['DG_MIN']
    numSelect = params['NUM_TOP_SENSITIVITY']
    return viennaDir, minDg, numSelect


def main():
    parser = argparse.ArgumentParser(prog='python -u build_final_output.py', 
                                     description='Build final primer designs')
    parser.add_argument('--name', dest='name', required=True, help='Pathogen name')
    parser.add_argument('--teval', dest='tar_eval', required=True, help='.result file from sensitivity evaluation')
    parser.add_argument('--oeval', dest='off_eval', nargs='+', required=True, help='List of specificity evaluation')
    parser.add_argument('--primers', dest='primers', required=True, help='Primer FASTA file')
    parser.add_argument('--params', dest='param_file', required=True, help='Parameter file')
    parser.add_argument('--out', dest='output', required=True, help='Final output')
    args = parser.parse_args()
    viennaDir, minDg, numSelect = parse_params(args.param_file)
    
    print(f'Building final output for {args.name}...')
    startTime = time.time()
    
    teval = pd.read_csv(args.tar_eval, index_col=[0,1]).iloc[:numSelect]
    teval.columns = ['cov_target', 'act_target', 'sco_target']
    merged = teval.copy()

    for oeval in args.off_eval:
        oname = oeval.split('/')[1].split('.')[1]
        oeval = pd.read_csv(oeval)
        tmp = pd.DataFrame(index=teval.index, columns=[f'cov_{oname}', f'act_{oname}', f'sco_{oname}'])
        for pair in teval.index:
            sub = oeval[(oeval['pname_f'].isin(set(pair)))&(oeval['pname_r'].isin(set(pair)))]
            if sub.empty:
                tmp.loc[pair] = 0
            else:
                tmp.loc[pair, f'cov_{oname}'] = sub['coverage'].sum()
                tmp.loc[pair, f'act_{oname}'] = sub['activity'].max()
                tmp.loc[pair, f'sco_{oname}'] = sub['score'].sum()
        merged = merged.join(tmp)


    RNAduplex = f'{viennaDir}/src/bin/RNAduplex'
    priseqs = { s.id:str(s.seq) for s in SeqIO.parse(args.primers, 'fasta') }
    for pair in merged.index:
        fname, rname = pair
        fseq = priseqs[fname]
        rseq = priseqs[rname]
        deltaG = get_deltaG_vienna(fseq, rseq, RNAduplex)
        merged.loc[pair, 'Dimer_dG'] = deltaG
        merged.loc[pair, 'pseq_f'] = fseq
        merged.loc[pair, 'pseq_r'] = rseq
    final = merged[merged['Dimer_dG'] > minDg]
    final.round(3).to_csv(args.output)
    
    runtime = (time.time() - startTime)
    print(f'Wrote {len(final)} primer pairs to {args.output} (%.1f sec).' % runtime)
    
if __name__ == '__main__':
    main()

