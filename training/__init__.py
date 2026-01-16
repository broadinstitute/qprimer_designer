#!/usr/bin/env python
# coding: utf-8
"""Training module matplotlib configuration."""

import os
import matplotlib as mpl
import matplotlib.font_manager as fm
import matplotlib.pyplot as plt

# Load custom fonts if QPRIMER_FONT_PATH is set, otherwise use system sans-serif
font_path = os.environ.get('QPRIMER_FONT_PATH')
if font_path and os.path.isdir(font_path):
    font_files = fm.findSystemFonts(fontpaths=font_path)
    for font_file in font_files:
        fm.fontManager.addfont(font_file)
    mpl.rcParams['font.sans-serif'] = 'Helvetica'
else:
    # Use system default sans-serif font
    mpl.rcParams['font.sans-serif'] = ['DejaVu Sans', 'Helvetica', 'Arial', 'sans-serif']

mpl.rcParams['axes.spines.right'] = False
mpl.rcParams['axes.spines.top'] = False
mpl.rcParams['figure.figsize'] = (2.5, 2.5)
mpl.rcParams['axes.labelsize'] = 10
mpl.rcParams['xtick.labelsize'] = 10
mpl.rcParams['ytick.labelsize'] = 10

