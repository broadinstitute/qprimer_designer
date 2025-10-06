#!/usr/bin/env python
# coding: utf-8

# In[18]:


import argparse
import bisect
import time
import subprocess
import re
import random
from Bio.SeqUtils import MeltingTemp, gc_fraction
from Bio import Seq, SeqIO


# In[19]:


def reverse_complement_dna(seq):
    return str(Seq.Seq(seq).reverse_complement())


def get_Tm(seq):
    return MeltingTemp.Tm_NN(seq, Na=50, Mg=1.5, dNTPs=.6) # Primer3 condition


def get_GC_percentage(seq):
    return gc_fraction(seq)*100


def get_deltaG_vienna(seq1, seq2, RNAduplex):
    inp = f'{seq1}\n{seq2}\n'
    res = subprocess.run([RNAduplex, '--noconv', '--paramFile=DNA'], input=inp,
                         capture_output=True, text=True, check=True)
    deltaG = re.search(r'\(\s*([-+]?\d*\.?\d+)\s*\)', res.stdout)
    if not deltaG:
        raise ValueError(f'Could not parse ΔG from RNAduplex output:\n{res.stdout}\n{seq1}')
    return float(deltaG.group(1))  


def generate_primers_single(targetSeq, step, primerLen, minAmpLen):   
    targetSeq_rc = reverse_complement_dna(targetSeq)
    
    forps, revps = {}, {}
    for i in range(0, len(targetSeq)-primerLen-minAmpLen, step):
        fseq = targetSeq[i:i+primerLen]
        if 'N' not in fseq:
            forStart = i
            forps[fseq] = forStart
        rseq = targetSeq_rc[i:i+primerLen]
        if 'N' not in rseq:
            revEnd = len(targetSeq)-i
            revps[rseq] = revEnd
    return forps, revps


def generate_primers_multi(targetSeqs, step, primerLen, minAmpLen, maxAmpLen, 
                           maxTm, minTm, maxGc, minDg, viennaDir):
    RNAduplex = f'{viennaDir}/src/bin/RNAduplex'
    
    forps, revps = {}, {}
    for tseq in targetSeqs:
        flist, rlist = generate_primers_single(tseq, step, primerLen, minAmpLen)
        forps.update(flist)
        revps.update(rlist)
    uniqPairs = count_primer_pairs(forps, revps, minAmpLen, maxAmpLen)
    print(f'>> Primers with a unique sequence: {len(forps)} forwards, {len(revps)} reverses, {uniqPairs} pairs')
    
    forFilt, revFilt = {}, {}
    for unfilt, filt in zip([forps,revps], [forFilt,revFilt]):
        for pseq in unfilt:
            if minTm <= get_Tm(pseq) <= maxTm and get_GC_percentage(pseq) <= maxGc:
                if get_deltaG_vienna(pseq, pseq, RNAduplex) >= minDg:
                    filt[pseq] = unfilt[pseq]           
    validPairs = count_primer_pairs(forFilt, revFilt, minAmpLen, maxAmpLen)
    print(f'>> Primers filtered: {len(forFilt)} forwards, {len(revFilt)} reverses, {validPairs} pairs')
    return forFilt, revFilt


def count_primer_pairs(fors, revs, minAmpLen, maxAmpLen):
    sts = sorted(fors.values())
    ens_sorted = sorted(revs.values())
    count = 0
    for st in sts:
        left = bisect.bisect_left(ens_sorted, st+minAmpLen)
        right = bisect.bisect_right(ens_sorted, st+maxAmpLen)
        count += (right-left)
    return count


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
    maxNum = params['MAX_PRIMER_CANDIDATES']
    step = params['TILING_STEP']
    primerLen = params['PRIMER_LEN']
    minAmpLen = params['AMPLEN_MIN']
    maxAmpLen = params['AMPLEN_MAX']
    maxTm = params['TM_MAX']
    minTm = params['TM_MIN']
    maxGc = params['GC_MAX']
    minDg = params['DG_MIN']
    return viennaDir, maxNum, step, primerLen, minAmpLen, maxAmpLen, maxTm, minTm, maxGc, minDg


# In[23]:


def main():
    parser = argparse.ArgumentParser(prog='python -u generate_primers.py', 
                                     description='Generate primer candidates')
    parser.add_argument('--in', dest='target_seqs', required=True, help='FASTA file of target sequences')
    parser.add_argument('--out', dest='primer_seqs', required=True, help='FASTA file of primers')
    parser.add_argument('--params', dest='param_file', required=True, help='Parameters file')
    parser.add_argument('--name', dest='name', required=True, help='Pathogen name')
    args = parser.parse_args()

    targetSeqs = [ str(s.seq) for s in SeqIO.parse(args.target_seqs, 'fasta') ]
    viennaDir, maxNum, step, primerLen, minAmpLen, maxAmpLen, \
    maxTm, minTm, maxGc, minDg = parse_params(args.param_file)
    
    print(f'Generating primers from {args.target_seqs}...')
    startTime = time.time()
    
    forFilt, revFilt = generate_primers_multi(targetSeqs, step, primerLen, minAmpLen, maxAmpLen, 
                                              maxTm, minTm, maxGc, minDg, viennaDir)
    forwards = list(forFilt.keys())
    reverses = list(revFilt.keys())
    if len(forwards) > maxNum//2:
        forwards = random.sample(forwards, maxNum//2)
    if len(reverses) > maxNum//2:
        reverses = random.sample(reverses, maxNum//2)
    with open(args.primer_seqs, 'w') as fout:
        for i, forseq in enumerate(forwards):
            fout.write(f'>{args.name}_{i+1}_f\n{forseq}\n')
        for i, revseq in enumerate(reverses):
            fout.write(f'>{args.name}_{i+1}_r\n{revseq}\n')
    
    runtime = (time.time() - startTime)
    print(f'Wrote {len(forwards)+len(reverses)} primers to {args.primer_seqs} (%.1f sec)' % runtime)

if __name__ == '__main__':
    main()


# In[ ]:




