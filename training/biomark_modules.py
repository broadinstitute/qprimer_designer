#!/usr/bin/env python
# coding: utf-8

# **Packages**
# - numpy 1.23.5
# - pandas 2.2.3

# In[3]:


import numpy as np
import pandas as pd
import os
import warnings
from collections import defaultdict
from io import StringIO


# In[8]:


THRES = .1
CT_MAX = 40.0
CT_CONTAM = 36
CT_MAXDIFF = 4
print('Parameters:\n - Ct threshold: %.1f\n - Ct max: %.1f\n - Ct contamination: %.1f\n\
 - Ct max difference: %.1f'%(THRES, CT_MAX, CT_CONTAM, CT_MAXDIFF))


# In[12]:


def parse_reaction_maps(tmap, pmap):
    index = 0
    reactions = {}
    trep = defaultdict(int)
    rrep = defaultdict(int)
    for trow in tmap.index:
        for tcol in tmap.columns:
            twell = 'A%s'%str(len(tmap.columns)*trow+tcol+1).zfill(2)
            target = tmap.loc[trow,tcol]
            trep[target] += 1
            prep = defaultdict(int)
            for prow in pmap.index:
                for pcol in pmap.columns:
                    pwell = 'S%s'%str(len(pmap.columns)*prow+pcol+1).zfill(2)
                    primer = pmap.loc[prow,pcol]
                    prep[primer] += 1
                    rrep[(target,primer)] += 1
                    chamber = '%s-%s'%(pwell,twell)
                    reactions[index] = { 'chamber':chamber, 'assay':twell, 'sample':pwell, 
                                         'target':target, 'primer':primer, 'reaction':'%s-%s'%(target,primer),
                                         'trep':trep[target], 'prep':prep[primer], 'rep':rrep[(target,primer)] }
                    index += 1
    chip = '%s.%s'%(pwell[1:],twell[1:])
    reactions = pd.DataFrame(reactions).T
    targets = set(reactions['target'].tolist())
    primers = set(reactions['primer'].tolist())
    return chip, targets, primers, reactions

def parse_biomark_rawdata(resultf):
    rawdata = {}
    chunks = open(resultf,'rt').read().split('\n\n')
    for chunk in chunks:
        label = chunk.split('\n')[0].strip()
        if label.startswith('Bkgd') or label.startswith('Raw'):
            label = '_'.join(np.array(label.split())[[-1,0]]) # ex. "Raw Data for ROX" -> "rox_raw"
            datatbl = pd.read_csv(StringIO('\n'.join(chunk.split('\n')[1:])),index_col=0)
            datatbl = datatbl.loc[:,~datatbl.columns.str.contains('^Unnamed')]
            datatbl.columns = [int(i) for i in datatbl.columns]
            rawdata[label] = datatbl
    print('Raw data table list: %s\n'%', '.join(sorted(rawdata.keys())))
    return rawdata

def remove_noise(row):
    vs = row.values[::-1]
    newvs = [vs[0]]
    lmin = vs[0]
    for i in range(1,len(vs)):
        if vs[i]<lmin:
            lmin = vs[i]
        newvs.append(lmin)
    return pd.Series({i+1:v for i,v in enumerate(newvs[::-1])})

def min_max_scale(row):
    row = remove_noise(row)
    minv = row[[1,2]].mean()
    maxv = row[[39,40]].mean()
    if maxv-minv<=.1:
        return (row-minv).clip(lower=0)
    return (row-minv).clip(lower=0) * maxv / (maxv-minv)

def qc_outliers(sigType, option, sigs, rntbl, nsamples, nassays, stdcnt=3, initCycle=1):
    print('\tQC %s %s signals...'%(option,sigType))
    initSigs = sigs[initCycle]
    if option in ['high','High']:
        ub = initSigs.mean() + stdcnt*initSigs.std()
        chs = list(sigs[initSigs>ub].index)
    elif option in ['low','Low']:
        lb = initSigs.mean() - stdcnt*initSigs.std()
        chs = list(sigs[initSigs<lb].index)
    else:
        print('option parameter is not valid'); return [],[]

    rns = rntbl[rntbl['chamber'].isin(chs)]    
    assaysToExclude, samplesToExclude = [], []
    assaysToPrint, samplesToPrint = [], []
    for (awell,target),grp in rns.groupby(['assay','target']):
        if len(grp)>nsamples/3:
            assaysToPrint.append('%s %s (%i)'%(target,awell,len(grp)))
            assaysToExclude.extend(grp.loc[grp['assay']==awell,'chamber'].tolist())
    for (swell,primer),grp in rns.groupby(['sample','primer']):
        if len(grp)>nassays/3:
            samplesToPrint.append('%s %s (%i)'%(primer,swell,len(grp)))
            samplesToExclude.extend(grp.loc[grp['sample']==swell,'chamber'].tolist())
    
    if assaysToPrint or samplesToPrint:
        print('\t %i chambers from %i assays and %i samples show too %s %s signals'\
          %(len(rns),len(set(rns['assay'])),len(set(rns['sample'])), option, sigType))
        if assaysToPrint:
            print('\t Targets:', ', '.join(assaysToPrint))
        if samplesToPrint:
            print('\t Primers:', ', '.join(samplesToPrint))
    return assaysToExclude, samplesToExclude

def qc_controls(sigtbl, rntbl, controls): 
    targetPcs, targetNcs, primerPcs, primerNcs = controls
    rns = rntbl.set_index('chamber').copy()
    rns = rns.join(pd.DataFrame({'ct':sigtbl.apply(get_ct,axis=1)}), how='inner')
    
    noPrimer, noTarget, contamPrimer, contamTarget = [], [], [], []
    for (primer,swell),grp in rns[rns['target'].isin(targetPcs)].groupby(['primer','sample']):
        if primer not in primerNcs and (grp['ct']>=CT_MAX).all():
            noPrimer.append((primer,swell))   
    for (target,awell),grp in rns[rns['primer'].isin(primerPcs)].groupby(['target','assay']):
        if target not in targetNcs and (grp['ct']>=CT_MAX).all():
            noTarget.append((target,awell))  
    for (primer,swell),grp in rns.groupby(['primer','sample']):
        ncs = grp[grp['target'].isin(targetNcs)]
        els = grp[~grp['target'].isin(targetNcs)]
        if (ncs['ct']<CT_CONTAM).all() and ncs['ct'].mean()-els['ct'].mean()<CT_MAXDIFF:
            contamPrimer.append((primer,swell))
    for (target,awell),grp in rns.groupby(['target','assay']):
        ncs = grp[grp['primer'].isin(primerNcs)]
        els = grp[~grp['primer'].isin(primerNcs)]
        if (ncs['ct']<CT_CONTAM).all() and ncs['ct'].mean()-els['ct'].mean()<CT_MAXDIFF:
            contamTarget.append((target,awell))
    
    if noPrimer:
        print('\tPrimer presence check: %i wells show no amplification with target PCs'%len(noPrimer))
        print('\t',' '.join(['%s (%s)'%(p,w) for p,w in noPrimer]))
    if contamPrimer:
        print('\tPrimer contamination check: %i wells show amplification with target NCs'%len(contamPrimer)) 
        print('\t',' '.join(['%s (%s)'%(p,w) for p,w in contamPrimer]))
    if noTarget:
        print('\tTarget presence check: %i wells show no amplification with primer PCs'%len(noTarget)) 
        print('\t',' '.join(['%s (%s)'%(t,w) for t,w in noTarget]))
    if contamTarget:
        print('\tTarget contamination check: %i targets show amplification with primer NCs'%len(contamTarget))
        print('\t',' '.join(['%s (%s)'%(t,w) for t,w in contamTarget]))
    return noPrimer, contamPrimer, noTarget, contamTarget

def get_ct(row):
    lowCys = [cyc for cyc,sig in row.items() if sig<THRES]
    if not lowCys:
        print('Error: all signals are stronger than threshold.')
        return -1
    return max(lowCys)

def get_ct_cont(row):
    for cyc,sig in row.items():
        if cyc==CT_MAX and sig<THRES:
            return CT_MAX
        if sig<THRES and row[cyc+1]>=THRES:
            return cyc + (THRES-sig)/(row[cyc+1]-sig)
    print('Error: all signals are stronger than threshold.')
    return -1


# In[18]:


class Biomark:
    def __init__(self, **kwargs):
        self.inputs = kwargs.get('inputs',None)
        self.params = kwargs.get('params',None)
        tmap = pd.ExcelFile(self.inputs['tmapf']).parse(self.params['tset'],header=None) 
        pmap = pd.ExcelFile(self.inputs['pmapf']).parse(self.params['pset'],header=None) 
        self.chip, self.targets, self.primers, self.reactions = parse_reaction_maps(tmap, pmap)
        self.rawdata = parse_biomark_rawdata(self.inputs['resultf'])
        
        self.controls = []
        self.normsigs = pd.DataFrame()
        self.adjsigs = pd.DataFrame()
        self.exclude = {}
        self.output = pd.DataFrame()
        
    def qc_raw_sigs(self, cond, option):
        dye = self.params[cond]
        sigs = self.rawdata['%s_Raw'%dye] / self.rawdata['%s_Bkgd'%dye]
        nsamples, nassays = [int(n) for n in self.chip.split('.')]
        return qc_outliers(dye, option, sigs, self.reactions, nsamples, nassays)
    
    def set_norm_sigs(self, initCycles=[1,2,3], exclude_high=False, exclude_low=True):
        print('Normalizing signals...')
        print('discarding high or low signals: %s, %s'%(exclude_high, exclude_low))
        exclude = defaultdict(list)
        for cond in ['dyeRef','dyeAmp']:
            for option in ['high','low']:
                for label,excl in zip(['assay','sample'],self.qc_raw_sigs(cond, option)):
                    if excl:
                        qc = '%s %s %s'%(label,option,self.params[cond])
                        self.reactions.loc[self.reactions['chamber'].isin(excl),'QC_signals'] = qc
                        exclude[option].extend(excl)
        ampsigs = self.rawdata['%s_Raw'%self.params['dyeAmp']] / self.rawdata['%s_Bkgd'%self.params['dyeAmp']]
        refsigs = self.rawdata['%s_Raw'%self.params['dyeRef']] / self.rawdata['%s_Bkgd'%self.params['dyeRef']]
        normsigs = (ampsigs / refsigs)
        maxv = np.percentile(normsigs[40], 95)
        normsigs = (normsigs / maxv).clip(upper=1)
        self.normsigs = normsigs.apply(min_max_scale,axis=1)
        if exclude_high and exclude['high']:
            self.normsigs.drop(set(exclude['high']), inplace=True)
            print('%i reactions with high signals were excluded.'%len(set(exclude['high'])))
        if exclude_low and exclude['low']:
            self.normsigs.drop(set(exclude['low']), inplace=True)
            print('%i reactions with low signals were excluded.'%len(set(exclude['low'])))
        print('Normalizing signals completed.\n')
    
    def qc_norm_sigs(self):
        texcel = pd.ExcelFile(self.inputs['tmapf'])
        pexcel = pd.ExcelFile(self.inputs['pmapf'])
        if not ('controls' in texcel.sheet_names and 'controls' in pexcel.sheet_names):
            print('No information about controls.')
            return [], [], [], []
        targetPcs, targetNcs = [ t.split(';') for t in texcel.parse('controls')['targets'] ]
        primerPcs, primerNcs = [ p.split(';') for p in pexcel.parse('controls')['primers'] ]
        self.controls = [targetPcs, targetNcs, primerPcs, primerNcs]
        #print('> Target controls: %s (P), %s (N)'%(', '.join(targetPcs),', '.join(targetNcs)))
        #print('> Primer controls: %s (P); %s (N)'%(', '.join(primerPcs),', '.join(primerNcs)))
        return qc_controls(self.normsigs, self.reactions, self.controls)
           
    def adjust_basal_sig(self, exclude_lowq=False): # adjust basal signal
        print('Checking target/primer absence or contamination...')
        print('discarding low-quality reactions: %s'%(exclude_lowq))
        noPrimer, contamPrimer, noTarget, contamTarget = self.qc_norm_sigs()
        targetPcs, targetNcs, primerPcs, primerNcs = self.controls
        
        exclude = []
        for label,ls in zip(['no','contaminated'],[noPrimer,contamPrimer]):
            if ls:
                swells = list(zip(*ls))[1]
                self.reactions.loc[self.reactions['sample'].isin(swells),'QC_controls'] = '%s primer'%label
                exclude.extend(self.reactions.loc[self.reactions['sample'].isin(swells),'chamber'].tolist())
        for label,ls in zip(['no','contaminated'],[noTarget,contamTarget]):
            if ls:
                awells = list(zip(*ls))[1]
                self.reactions.loc[self.reactions['assay'].isin(awells),'QC_controls'] = '%s target'%label
                exclude.extend(self.reactions.loc[self.reactions['assay'].isin(awells),'chamber'].tolist())
        
        if exclude_lowq:
            self.adjsigs = self.normsigs.copy().drop(exclude)
            print('\t%i reactions were excluded.'%len(exclude))
        else:
            self.adjsigs = self.normsigs.copy()
        print('Filtering low-quality reactions completed.\n')
        
    def set_cts(self):        
        print('Calculating readouts...')
        targetPcs, targetNcs, primerPcs, primerNcs = self.controls
        rns = self.reactions.set_index('chamber').reindex(self.adjsigs.index)
        rns['ct'] = self.adjsigs.apply(get_ct_cont,axis=1)
        ncCts = {swell:grp['ct'].mean() for swell,grp in rns[rns['target'].isin(targetNcs)].groupby('sample')}
        pcCts = {awell:grp['ct'].mean() for awell,grp in rns[rns['primer'].isin(primerPcs)].groupby('assay')}
        rns['nc_ct'] = rns['sample'].apply(lambda x:ncCts[x])
        rns['pc_ct'] = rns['assay'].apply(lambda x:pcCts[x])
        rns['score'] = ((rns['nc_ct']-rns['ct'])/(rns['nc_ct']-rns['pc_ct'])).clip(lower=0)
        
        cts = rns.pivot_table(values='ct', index='reaction', columns='rep')
        scores = rns.pivot_table(values='score', index='reaction', columns='rep')
        self.output = cts.join(scores,lsuffix='_ct',rsuffix='_score',how='inner')
        self.output['mean'] = scores.mean(axis=1)
        self.output['std'] = scores.std(axis=1)
        print('Calculating Ct values and PCR scores completed.\n')


# In[21]:


# DATAPATH = '/home/jupyter/ADAPT_PCR_share/safe/biomark'
# INPUTS = {
#     'tmapf':'%s/0130_map_targets.xlsx'%DATAPATH,
#     'pmapf':'%s/0130_map_primers.xlsx'%DATAPATH }
# PARAMS = {'dyeRef':'ROX', 'dyeAmp':'EvaGreen'}


# In[22]:


# output = '2025-01-16_target-set1_primer-set11.csv'
# _, ts, ps = output.replace('.csv','').split('_')
# tset = ts.split('-')[1]
# pset = ps.split('-')[1]
# inputs = INPUTS.copy()
# inputs.update({'resultf':'%s/data/%s'%(DATAPATH,output)})
# params = PARAMS.copy()
# params.update({'tset':tset, 'pset':pset})
# biomark = Biomark(inputs=inputs, params=params)
# biomark.set_norm_sigs(exclude_high=False, exclude_low=True)
# biomark.adjust_basal_sig(exclude_lowq=True)
# biomark.set_cts()


# In[ ]:





# ---

# In[8]:


# class Biomark:
#     def __init__(self, **kwargs):
#         self.inputs = kwargs.get('inputs',None)
#         self.params = kwargs.get('params',None)
#         tmap = pd.ExcelFile(self.inputs['tmapf']).parse(self.params['tset'],header=None) 
#         pmap = pd.ExcelFile(self.inputs['pmapf']).parse(self.params['pset'],header=None) 
#         self.chip, self.targets, self.primers, self.reactions = parse_reaction_maps(tmap, pmap)
#         self.rawdata = parse_biomark_rawdata(self.inputs['resultf'])

#         self.normsigs = pd.DataFrame()
#         self.adjsigs = pd.DataFrame()
#         self.exclude = {}
#         self.cts = pd.DataFrame()
        
#     def qc_raw_sigs(self, cond, option):
#         dye = self.params[cond]
#         sigs = self.rawdata['%s_Raw'%dye] / self.rawdata['%s_Bkgd'%dye]
#         nsamples, nassays = [int(n) for n in self.chip.split('.')]
#         return qc_outliers(dye, option, sigs, self.reactions, nsamples, nassays)
    
#     def set_norm_sigs(self, initCycles=[1,2,3], exclude_high=False, exclude_low=True):
#         print('Setting normalized signals and discarding abnormal signals...')
#         print('excluding high signals: %s, excluding low signals: %s'%(exclude_high, exclude_low))
#         exclude = defaultdict(list)
#         for cond in ['dyeRef','dyeAmp']:
#             for option in ['high','low']:
#                 for label,excl in zip(['assay','sample'],self.qc_raw_sigs(cond, option)):
#                     if excl:
#                         qc = '%s %s %s'%(label,option,self.params[cond])
#                         self.reactions.loc[self.reactions['chamber'].isin(excl),'QC_signals'] = qc
#                         exclude[option].extend(excl)
#         ampsigs = self.rawdata['%s_Raw'%self.params['dyeAmp']] / self.rawdata['%s_Bkgd'%self.params['dyeAmp']]
#         refsigs = self.rawdata['%s_Raw'%self.params['dyeRef']] / self.rawdata['%s_Bkgd'%self.params['dyeRef']]
#         normsigs = (ampsigs / refsigs)
#         maxv = np.percentile(normsigs[40], 95)
#         normsigs = (normsigs / maxv).clip(upper=1)
#         self.normsigs = normsigs.apply(min_max_scale,axis=1)
#         if exclude_high and exclude['high']:
#             self.normsigs.drop(set(exclude['high']), inplace=True)
#             print('%i reactions with high signals were excluded.'%len(set(exclude['high'])))
#         if exclude_low and exclude['low']:
#             self.normsigs.drop(set(exclude['low']), inplace=True)
#             print('%i reactions with low signals were excluded.'%len(set(exclude['low'])))
#         print('Normalizing signals completed.\n')
    
#     def qc_norm_sigs(self):
#         texcel = pd.ExcelFile(self.inputs['tmapf'])
#         pexcel = pd.ExcelFile(self.inputs['pmapf'])
#         if not ('controls' in texcel.sheet_names and 'controls' in pexcel.sheet_names):
#             print('No information about controls.')
#             return [], [], [], []
#         targetPcs, targetNcs = [ t.split(';') for t in texcel.parse('controls')['targets'] ]
#         primerPcs, primerNcs = [ p.split(';') for p in pexcel.parse('controls')['primers'] ]
#         controls = [targetPcs, targetNcs, primerPcs, primerNcs]
#         #print('> Target controls: %s (P), %s (N)'%(', '.join(targetPcs),', '.join(targetNcs)))
#         #print('> Primer controls: %s (P); %s (N)'%(', '.join(primerPcs),', '.join(primerNcs)))
#         return qc_controls(self.normsigs, self.reactions, controls)
           
#     def adjust_basal_sig(self, exclude_lowq=False): # adjust basal signal
#         print('Checking quality of the normalized signals...')
#         print('excluding no/contaminated targets/primers: %s'%exclude_lowq)
#         noPrimer, contamPrimer, noTarget, contamTarget, controls = self.qc_norm_sigs()
#         targetPcs, targetNcs, primerPcs, primerNcs = controls
        
#         exclude = []
#         for label,ls in zip(['no','contaminated'],[noPrimer,contamPrimer]):
#             if ls:
#                 swells = list(zip(*ls))[1]
#                 self.reactions.loc[self.reactions['sample'].isin(swells),'QC_controls'] = '%s primer'%label
#                 exclude.extend(self.reactions.loc[self.reactions['sample'].isin(swells),'chamber'].tolist())
#         for label,ls in zip(['no','contaminated'],[noTarget,contamTarget]):
#             if ls:
#                 awells = list(zip(*ls))[1]
#                 self.reactions.loc[self.reactions['assay'].isin(awells),'QC_controls'] = '%s target'%label
#                 exclude.extend(self.reactions.loc[self.reactions['assay'].isin(awells),'chamber'].tolist())
        
#         if exclude_lowq:
#             adjsigs = self.normsigs.copy().drop(exclude)
#             print('\t%i reactions were excluded.'%len(exclude))
#         else:
#             adjsigs = self.normsigs.copy()
#         rns = self.reactions.set_index('chamber').reindex(adjsigs.index)
#         rns['ct'] = adjsigs.apply(get_ct,axis=1)
#         basalCts = {swell:grp['ct'].mean() for swell,grp in rns[rns['target'].isin(targetNcs)].groupby('sample')}
#         rns['basal_ct'] = rns['sample'].apply(lambda x:basalCts[x])
# #         rawcts = pd.DataFrame({'ct':self.normsigs.apply(get_ct,axis=1)})
# #         rns = self.reactions.set_index('chamber').join(rawcts,how='inner')
# #         ctcts = rns.pivot_table(index='sample',columns='target',values='ct',aggfunc='mean')[targetNcs].min(axis=1)
# #         rns['conct'] = rns['primer'].apply(lambda x:ctcts[x])
#         deltaCt = (rns['ct']-rns['basal_ct']).clip(upper=0)
#         self.adjsigs = (adjsigs.T * (1-deltaCt.apply(lambda x:2**x))).T
#         print('Adjusting signals completed.\n')
        
#     def set_cts(self):
#         print('Extracting ct values...')
#         cts = pd.DataFrame({'ct_raw':self.normsigs.apply(get_ct,axis=1),
#                             'ct_adj':self.adjsigs.apply(get_ct,axis=1)}).dropna()

#         joined = self.reactions.set_index('chamber').join(cts,how='inner')
#         print('\tcts from %i reactions are getting calculated.'%len(joined))
#         pivot = joined.pivot_table(values='ct_raw', index='reaction', columns='rep')
#         self.cts_raw = pivot.copy()
#         self.cts_raw['mean'] = pivot.mean(axis=1)
#         self.cts_raw['std'] = pivot.std(axis=1)
        
#         pivot = joined.pivot_table(values='ct_adj', index='reaction', columns='rep')
#         self.cts_adj = pivot.copy()
#         self.cts_adj['mean'] = pivot.mean(axis=1)
#         self.cts_adj['std'] = pivot.std(axis=1)
#         print('Calculating ct values completed.\n')


# In[5]:


# tmapf = '/home/jupyter/ADAPT_PCR_share/safe/biomark/Target_Map.xlsx'
# pmapf = '/home/jupyter/ADAPT_PCR_share/safe/biomark/Primer_Map.xlsx'
# resultf = '/home/jupyter/ADAPT_PCR_share/safe/biomark/2025-01-08_Serial_dilution_biomark.csv'

# tset = 'rep1'
# pset = 'rep1'
# dyeRef = 'ROX'
# dyeAmp = 'EvaGreen'

# inputs = {'tmapf':tmapf, 'pmapf':pmapf, 'resultf':resultf}
# params = {'tset':tset, 'pset':pset, 'dyeRef':dyeRef, 'dyeAmp':dyeAmp}

# pilot = Biomark(inputs=inputs, params=params)
# pilot.set_norm_sigs(exclude_high=False, exclude_low=True)
# pilot.adjust_basal_sig(exclude_lowq=True)
# pilot.set_cts()


# In[40]:


# tmapf = '012025_biomark/data/0130_map_targets.xlsx'
# pmapf = '012025_biomark/data/0130_map_primers.xlsx'
# resultf = '012025_biomark/data/2025-01-16_target-set1_primer-set11.csv'

# tset = 'set1'
# pset = 'set11'
# dyeRef = 'ROX'
# dyeAmp = 'EvaGreen'

# inputs = {'tmapf':tmapf, 'pmapf':pmapf, 'resultf':resultf}
# params = {'tset':tset, 'pset':pset, 'dyeRef':dyeRef, 'dyeAmp':dyeAmp}

# pilot = Biomark(inputs=inputs, params=params)
# pilot.set_norm_sigs(exclude_high=False, exclude_low=True)
# pilot.adjust_basal_sig(exclude_lowq=True)
# pilot.set_cts()


# In[ ]:





# In[ ]:





# In[ ]:





# In[ ]:





# In[ ]:





# In[ ]:


'''This is an object for a biomark X run. 
    > __init__: inputs, params
    > Instance variables:
        - inputs: a dictionary of files below
            + tmapf: an excel file with target configuration on IFC
            + pmapf: an excel file with primer configuration on IFC
            + resultf: the raw csv file from Biomark X run
        - params: a dictionary of parameters below
            + tset: target set name (should be one of sheet names of tmapf)
            + pset: primer set name (should be one of sheet names of pmapf)
            + dyeRef: reference dye (ex. ROX)
            + dyeAmp: amplification dye (ex. EvaGreen)
        - rawdata: tables of raw signals 
        - reactions: reaction conditions for each chamber
        - results:
            + rawsigs:
            + normsigs:
    > Instance functions:
        - calculate_signal: get (evagreen_raw - evagreen_bkgd) / (rox_raw - rox_bkgd)
        - normalize_min: normalize min value of signals
        - get_norm_sig: get normalized signals
        - qc_no_rox: 
'''

