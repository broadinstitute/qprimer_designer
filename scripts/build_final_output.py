#!/usr/bin/env python
# coding: utf-8

# In[ ]:


import argparse
import time
import subprocess
import re
import os
import pandas as pd
from Bio import SeqIO
from functools import reduce


# In[ ]:


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


# In[ ]:


def main():
    parser = argparse.ArgumentParser(prog='python -u build_final_output.py', 
                                     description='Build final primer designs')
    parser.add_argument('--name', dest='name', required=True, help='Pathogen name')
    parser.add_argument('--resdir', dest='resdir', required=True, help='Directory ML results were saved')
    parser.add_argument('--refs', dest='refs', nargs='+', required=True, help='List of evalauted references')
    parser.add_argument('--primers', dest='primers', required=True, help='Primer FASTA file')
    parser.add_argument('--params', dest='param_file', required=True, help='Parameter file')
    parser.add_argument('--out', dest='output', required=True, help='Final output')
    args = parser.parse_args()
    viennaDir, minDg, numSelect = parse_params(args.param_file)
    
    refname = ', '.join(args.refs)
    print(f'Building final output for {args.name} evaluated across {refname}...')
    startTime = time.time()
    
    fname = f'{args.resdir}/{args.name}.{args.name}.result'
    merged  = pd.read_csv(fname).iloc[:numSelect]
    merged = merged[['pname_f', 'pname_r', 'coverage', 'activity', 'score']]
    merged.columns = ['pname_f', 'pname_r', f'cov_{args.name}', f'act_{args.name}', f'sco_{args.name}']
    for ref in args.refs:
        filename = f'{args.resdir}/{args.name}.{ref}.result'
        if ref == args.name or os.path.getsize(filename) == 0:
            continue
        res = pd.read_csv(filename)
        res.columns = ['pname_f', 'pname_r', f'cov_{ref}', f'act_{ref}', f'sco_{ref}']
        merged = pd.merge(merged, res,  on=['pname_f','pname_r'], how='left').fillna(0)
    merged = merged.set_index(['pname_f','pname_r'])
    
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
    final.to_csv(args.output)
    
    runtime = (time.time() - startTime)
    print(f'Wrote {len(final)} primer pairs to {args.output} (%.1f sec).' % runtime)
    
if __name__ == '__main__':
    main()


# In[ ]:




