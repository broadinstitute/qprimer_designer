#!/usr/bin/env python
# coding: utf-8

# **Packages**
# - RNAstructure: version 6.5 (released on Jun 14, 2024)
# - primer3: libprimer3 release 2.6.1 

# In[5]:


TOOLPATH = '/home/jupyter/ADAPT_PCR_share/safe/tools' 
DATAPATH = '%s/RNAstructure/data_tables' % TOOLPATH
Fold = '%s/RNAstructure/exe/Fold' % TOOLPATH
bifold = '%s/RNAstructure/exe/bifold' % TOOLPATH
efn2 = '%s/RNAstructure/exe/efn2' % TOOLPATH
oligotm = '%s/primer3/src/oligotm' % TOOLPATH
primer3 = '%s/primer3/src/primer3_core' % TOOLPATH
RNAduplex = '%s/ViennaRNA-2.7.0/src/bin/RNAduplex' % TOOLPATH
bowtie2 = '%s/bowtie2-2.5.3-linux-x86_64/bowtie2' % TOOLPATH
sam2pairwise = '%s/sam2pairwise/src/sam2pairwise' % TOOLPATH


# In[6]:


def show_available_tools():
    print('Available tools: Fold, bifold, efn2, oligotm, primer3, RNAduplex, bowtie2, sam2pairwise')
show_available_tools()


# In[ ]:




