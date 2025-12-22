#!/usr/bin/env python
# coding: utf-8

# In[ ]:


import matplotlib as mpl
import matplotlib.font_manager as fm
import matplotlib.pyplot as plt

font_files = fm.findSystemFonts(fontpaths='/home/jupyter/ADAPT_PCR_share/safe/resources/Helvetica')
for font_file in font_files:
    fm.fontManager.addfont(font_file)
mpl.rcParams['font.sans-serif'] = 'Helvetica'
mpl.rcParams['axes.spines.right'] = 'off'
mpl.rcParams['axes.spines.top'] = 'off'
mpl.rcParams['figure.figsize'] = (2.5,2.5)
mpl.rcParams['axes.labelsize']: '10'
mpl.rcParams['xtick.labelsize']: '10'
mpl.rcParams['ytick.labelsize']: '10'

