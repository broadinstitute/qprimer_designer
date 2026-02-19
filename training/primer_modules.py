#!/usr/bin/env python
# coding: utf-8

# **Packages**
# - numpy 1.23.5
# - pandas 2.2.3
# - biopython 1.83
# - levenshtein 0.27.1

# In[2]:


import numpy as np
import pandas as pd
import os
from difflib import SequenceMatcher
from collections import defaultdict
from Bio.SeqUtils import MeltingTemp, gc_fraction
from Bio import Seq, SeqIO
from Levenshtein import distance


# In[3]:


from adapt_pcr_external_tools import *
get_ipython().run_line_magic('env', 'DATAPATH $DATAPATH')


# ---

# **Class Primer**  

# In[4]:


def reverse_complement_dna(seq):
    return str(Seq.Seq(seq).reverse_complement())

def delete(files):
    for f in files:
        _ = get_ipython().getoutput('rm -f $f')

def get_dG(ct, maxdG=1):
    key = 'ENERGY'
    res = get_ipython().getoutput('grep $key $ct # "  45  ENERGY = -2.1  For1"')
    if res: 
        return float(res[0].strip().split()[3])
    return maxdG

def eval_self_hairpin(fa, fapath, ta, numst=1):
    ta = ta + 273
    ct = '%s/hairpin.ct' % fapath
    _ = get_ipython().getoutput('$Fold $fa $ct -t $ta -m $numst -d')
    dG = get_dG(ct)
    delete([ct])
    return dG

def eval_dimer(fa1, fa2, fapath, ta, numst=1):
    ta = ta + 273
    ct = '%s/dimer.ct' % fapath
    _ = get_ipython().getoutput('$bifold $fa1 $fa2 $ct -t $ta -m $numst -d -i # dna parameter, intermolecular inx only')
    dG = get_dG(ct)
    delete([ct])
    return dG
    
def eval_binding_energy(seq1, seq2, fapath, ta, show_efn=False, numst=1):
    seq1, seq2 = sorted([seq1,seq2],key=lambda x:len(x))
    startPos, dist = min([(i,distance(seq1,seq2[i:i+len(seq1)])) for i in range(len(seq2)-len(seq1))],
                         key=lambda x:x[1])
    seq2 = reverse_complement_dna(seq2[startPos:startPos+len(seq1)])
    
    fa1 = '%s/1.fa'%fapath
    with open(fa1,'wt') as outf:
        outf.write('>1\n%s\n'%seq1)
    fa2 = '%s/2.fa'%fapath
    with open(fa2,'wt') as outf:
        outf.write('>2\n%s\n'%seq2)
    
    ta = ta + 273
    ct = '%s/dimer.ct' % fapath
    efn = '%s/efn.txt' % fapath
    _ = get_ipython().getoutput('$bifold $fa1 $fa2 $ct -t $ta -m $numst -d -i # dna parameter, intermolecular inx only')
    _ = get_ipython().getoutput('$efn2 $ct $efn -t $ta -d -p -w')
    enInfo = parse_efn2(efn, len(seq1))
    
    index = 0
    energies = {}
    for pos,en in enInfo.items():
        energies[index] = {'energy':en, 'position':len(seq1)-pos+1}
        index += 1
        
    if show_efn:
        print(open(efn,'rt').read().strip())
    delete([fa1,fa2,ct,efn])
    return pd.DataFrame(energies).T

def parse_efn2(efnFile, length):
    poses, ens = defaultdict(list), {}
    for l in open(efnFile,'rt'):
        if 'for stack of' in l or 'for closure of' in l:
            items = l.strip().split()
            etype = items[4]
            energy = float(items[2])
            pos = int(items[6].split('-')[0])
            poses[etype].append(pos)
            ens[pos] = energy
        elif 'Exterior loop' in l:
            ens[0] = float(l.strip().split()[3])
    d = {}
    i = 1
    while i<=length:
        if i in poses['stack']:
            d[i] = ens[i]; i+=1
        elif i in poses['closure']:
            if i==length:
                d[i] = ens[i]; i+=1
                continue
            nexti = min([p for p in poses['stack'] if p>i]+[i+1])
            for j in range(i,nexti):
                d[j] = ens[i]/(nexti-i)
            i = nexti
        else:
            d[i] = 0; i+=1
    d[length] += ens[0] # dG for exterior loop (end of match) is added to dg for the first bp
    return d

def eval_mismatches(seq1, seq2, show_match=False, window=1):
    mmtypes={'s':'sub', 'd':'del', 'i':'ins'}
    short, long = sorted([seq1,seq2],key=lambda x:len(x))
    startPos, dist = min([(i,distance(short,long[i:i+len(short)])) for i in range(len(long)-len(short))],
                         key=lambda x:x[1])
    long = long[startPos:startPos+len(short)]
    match1, match2, combined = parse_opcodes(SequenceMatcher(None, short, long).get_opcodes(), short, long)
    if show_match:
        print('seq1\t%s\nseq2\t%s\ncomb\t%s'%(match1,match2,combined))
    
    index = 0
    mismatches = defaultdict(list)
    for i,mm in enumerate(combined[::-1]):
        if mm in 'sdi':
            mismatches[index] = {'type':mmtypes[mm], 'position':i+1}
            index += 1
    return pd.DataFrame(mismatches).T

def parse_opcodes(opcodes, lseq, rseq):
    left, right, combined = '', '', ''
    for code, l1, l2, r1, r2 in opcodes:
        if code=='equal':
            combined += '-'*(l2-l1)
            left += lseq[l1:l2]
            right += rseq[r1:r2]
        elif code=='replace':
            combined += 's'*(l2-l1)
            left += lseq[l1:l2]
            right += rseq[r1:r2]
        elif code=='delete':
            combined += 'd'*(l2-l1)
            left += lseq[l1:l2]
            right += '-'*(l2-l1)
        else:
            combined += 'i'*(r2-r1)
            left += '-'*(r2-r1)
            right += rseq[r1:r2]
    return left, right, combined


# In[5]:


class Primer:
    '''This is an object for a single primer. 
        > __init__: pid, seq, type, fapath
        > Instance variables:
            - pid: Primer ID
            - seq: Primer sequence
            - type: 'for' or 'rev'
            - fapath: Path to save files
            - features: length, Tm, GC ratio, free energies for forming self-hairpin or dimier
            - fa: fasta file of the sequence
        > Instance functions:
            - set_features: set features
            - eval_dna_binding_energy: ...
            - eval_dna_mismatches: ...
    '''

    def __init__(self, **kwargs):
        self.pid = kwargs.get('pid',None)
        self.seq = kwargs.get('seq',None)
        self.type = kwargs.get('type',None)
        self.fapath = kwargs.get('fapath',None)
        self.features = pd.DataFrame()
        self.fa = '%s/%s.fa'%(self.fapath, self.pid)
        if not os.path.exists(self.fa):
            with open(self.fa,'wt') as outf:
                outf.write('>%s\n%s\n'%(self.pid, self.seq))
    
    def set_features(self, ta):
        length = len(self.seq)
        tm = MeltingTemp.Tm_NN(self.seq, Na=50, Mg=1.5, dNTPs=.6) # Primer3 parameters
        gcRatio = gc_fraction(self.seq) 
        dgHpn = eval_self_hairpin(self.fa, self.fapath, ta)
        dgDmr = eval_dimer(self.fa, self.fa, self.fapath, ta)
        features = { 'length':length, 'Tm':tm, 'GC ratio':gcRatio, 'dG hairpin':dgHpn, 'dG self-dimer':dgDmr }
        self.features = pd.DataFrame({self.pid:features}) # pid on columns
    
    def eval_dna_binding_energy(self, dnaid, dnaseq, ta, show_efn=False): # dnaseq should include For primer seq
        if self.type in ['rev','reverse']:
            dnaseq = reverse_complement_dna(dnaseq)
        entbl = eval_binding_energy(self.seq, dnaseq, self.fapath, ta, show_efn) 
        if not entbl.empty:
            entbl.loc[:,'primer'] = self.pid
            entbl.loc[:,'target'] = dnaid # columns: position, energy, primer, target
        return entbl.sort_values('position') # position counts from the 5' end of the primer binding site
    
    def eval_dna_mismatches(self, dnaid, dnaseq, show_match=False): # dnaseq should include For primer seq
        if self.type in ['rev','reverse']:
            dnaseq = reverse_complement_dna(dnaseq)
        mmtbl = eval_mismatches(self.seq, dnaseq, show_match)
        if not mmtbl.empty:
            mmtbl.loc[:,'primer'] = self.pid
            mmtbl.loc[:,'target'] = dnaid # columns: position, type, primer, target
        return mmtbl # position counts from the 5' end of the primer binding site


# In[6]:


# ## test
# INPUT_FILES = { 
#     'target_seq':'072024_target_design/1101_targets_edit.fa',
#     'target_mut_info':'072024_target_design/1101_mutagenesis_infos_edit.csv',
#     'primer_pair_info':'072024_target_design/1101_primer_pairs.csv',
#     'target_map':'012025_biomark/data/0130_map_targets.xlsx', 
#     'primer_map':'012025_biomark/data/0130_map_primers.xlsx',
#     'IFC_map':'012025_biomark/misc/IFC_mapping_9696.csv'}

# TARGET_SEQS = {s.id:str(s.seq) for s in SeqIO.parse(INPUT_FILES['target_seq'],'fasta')}
# PRIMER_PAIR_INFO = pd.read_csv(INPUT_FILES['primer_pair_info'],index_col=0)
# PRIMER_PAIR_INFO.index = ['P%i'%i for i in range(1,361)]
# PRIMER_PAIR_INFO.head(2)


# In[7]:


# fapath = '012025_biomark/misc'
# pairid = 'P1'
# row = PRIMER_PAIR_INFO.loc[pairid]

# forid = 'For%i'%row['Forward_ID']
# forseq = row['Forward_Seq']

# revid = 'Rev%i'%row['Reverse_ID']
# revseq = row['Reverse_Seq']

# forp = Primer(pid=forid, seq=forseq, type='for', fapath=fapath)
# revp = Primer(pid=revid, seq=revseq, type='rev', fapath=fapath)


# In[8]:


# ta = 55
# revp.set_features(ta)
# revp.features


# In[11]:


# m1 = TARGET_SEQS['M1']
# df = forp.eval_dna_mismatches('M1',m1,show_match=True)
# df


# In[10]:


# df = forp.eval_dna_binding_energy('M2',m2,ta,show_efn=True)
# df


# ---
# **Class Primer Pair**

# In[10]:


def get_product(targetSeq, forseq, revseq): # target should include the forward primer seq
    targetRc = reverse_complement_dna(targetSeq)
    tlen = len(targetSeq)
    dists = []
    for i in range(tlen-len(forseq)):
        ford = distance(forseq, targetSeq[i:i+len(forseq)])
        seqRc = reverse_complement_dna(targetSeq[i+len(forseq):])
        for j in range(len(seqRc)-len(revseq)):
            revd = distance(revseq, seqRc[j:j+len(revseq)])
            dists.append((ford+revd,i,tlen-j))
    dmin, fst, ren = min(dists,key=lambda x:x[0])
    return targetSeq[fst:ren], fst, ren

def calc_ta(forseq, revseq, prodseq):
    fortm = MeltingTemp.Tm_NN(forseq, Na=50, Mg=1.5, dNTPs=.6)
    revtm = MeltingTemp.Tm_NN(revseq, Na=50, Mg=1.5, dNTPs=.6)
    prodtm = MeltingTemp.Tm_NN(prodseq, Na=50, Mg=1.5, dNTPs=.6)
    return .3*min(fortm, revtm) + .7*prodtm - 14.9


# In[1]:


class PrimerPair:
    '''This is an object for a primer pair. 
        > __init__: pairid, forward, reverse
        > Instance variables:
            - pairid: Primer pair ID
            - forward: Forward primer object
            - reverse: Reverse primer object
            - fapath: Path to save files
            - features_each: Features of each pri
        > Instance functions:
            - set_features: set features
            - eval_target_binding: ...
    '''
    def __init__(self, **kwargs):
        self.pairid = kwargs.get('pairid',None)
        self.forward = kwargs.get('forward',None)
        self.reverse = kwargs.get('reverse',None)
        self.fapath = self.forward.fapath
        self.features = pd.DataFrame()
        
    def set_features(self, ta):
        self.forward.set_features(ta)
        self.reverse.set_features(ta)
        dgDmr = eval_dimer(self.forward.fa, self.reverse.fa, self.fapath, ta)
        combined = self.forward.features.join(self.reverse.features)
        combined.loc[:,self.pairid] = combined.mean(axis=1)
        combined.loc['dG dimer',self.pairid] = dgDmr
        self.features = combined      
        
    def eval_target_binding(self, targetId, targetSeq, ta, show_match=False, off_target=False):
        prodSeq, st, en = get_product(targetSeq, self.forward.seq, self.reverse.seq)
        
        tm = MeltingTemp.Tm_NN(prodSeq, Na=50, Mg=1.5, dNTPs=.6)
        taRyc = calc_ta(self.forward.seq, self.reverse.seq, prodSeq)
        features = { 'start':st, 'end':en, 'product':prodSeq, 'length':len(prodSeq), 
                     'Tm':tm, 'GC ratio':gc_fraction(prodSeq) }
        features = pd.DataFrame({self.pairid:features}) # pairid on columns
        
        ensFor = self.forward.eval_dna_binding_energy(targetId, prodSeq, ta)
        ensRev = self.reverse.eval_dna_binding_energy(targetId, prodSeq, ta)
        ensAll = pd.concat([ensFor,ensRev],axis=0,ignore_index=True)
        if not ensAll.empty:
            ensAll.loc[:,'start'] = st
            ensAll.loc[:,'end'] = en
            ensAll.loc[:,'pair'] = self.pairid
        
        mmsFor = self.forward.eval_dna_mismatches(targetId, prodSeq, show_match)
        mmsRev = self.reverse.eval_dna_mismatches(targetId, prodSeq, show_match)
        mmsAll = pd.concat([mmsFor,mmsRev],axis=0,ignore_index=True)
        if not mmsAll.empty:
            mmsAll.loc[:,'start'] = st
            mmsAll.loc[:,'end'] = en
            mmsAll.loc[:,'pair'] = self.pairid

        ensOffAll = pd.DataFrame()
        if off_target:
            nblock = 'N'*20
            offId = '%s_off'%targetId
            offSeq = targetSeq[:st]+nblock+targetSeq[st+20:en-20]+nblock+targetSeq[en:]
            
            ensOffFor = self.forward.eval_dna_binding_energy(offId, offSeq, ta)
            ensOffRev = self.reverse.eval_dna_binding_energy(offId, offSeq, ta)
            ensOffAll = pd.concat([ensOffFor,ensOffRev],axis=0,ignore_index=True)
    
            print('Off target sites')
            mmsFor = self.forward.eval_dna_mismatches(offId, offSeq, show_match=True)
            mmsRev = self.reverse.eval_dna_mismatches(offId, offSeq, show_match=True)

        return features, ensAll, mmsAll, ensOffAll


# In[2]:


# fapath = '012025_biomark/misc'
# pairid = 'P1'
# row = PRIMER_PAIR_INFO.loc[pairid]

# forid = 'For%i'%row['Forward_ID']
# forseq = row['Forward_Seq']

# revid = 'Rev%i'%row['Reverse_ID']
# revseq = row['Reverse_Seq']

# forp = Primer(pid=forid, seq=forseq, type='for', fapath=fapath)
# revp = Primer(pid=revid, seq=revseq, type='rev', fapath=fapath)
# pair = PrimerPair(pairid=pairid, forward=forp, reverse=revp)


# In[3]:


# ta = 55
# pair.set_features(ta)
# pair.features


# In[14]:


# targetId = 'M15'
# targetSeq = TARGET_SEQS[targetId]
# features, ens, mms, ensoff = pair.eval_target_binding(targetId,targetSeq,ta,show_match=True)


# In[15]:


# features


# In[17]:


# mms


# In[16]:


# ens


# In[ ]:




