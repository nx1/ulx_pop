# -*- coding: utf-8 -*-



import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import sys
import os

sys.path.append("..") # Adds higher directory to python modules path.
os.chdir('..')
from auxil import load_df_a
import matplotlib.colors as mcolors


df_a = load_df_a(transient_only=True)


import matplotlib
plt.rcParams.update({'font.size': 8})
matplotlib.rcParams['mathtext.fontset'] = 'stix'
matplotlib.rcParams['font.family'] = 'STIXGeneral'


hist_color = '#440256'

plt.figure(figsize=(5.2,3.9))
plt.subplot(2,2,3)
plt.ylim(0, 1)
plt.xlim(0, 45)
plt.xticks(np.arange(0, 50, 5))
plt.yticks(np.arange(0, 1.2, 0.2))
plt.hist2d(df_a['dincl'], df_a['ratio'], bins=[45,45], cmap='viridis')
plt.xlabel('Precessional angle $\Delta i$')
plt.ylabel('Ratio')


plt.subplot(2,8,13)
plt.xticks([])
plt.yticks([])
plt.ylim(0, 1)
plt.hist(df_a['ratio'],bins=45, orientation="horizontal", color=hist_color)

plt.subplot(8,2,7)
plt.xlim(0, 45)
plt.xticks([])
plt.yticks([])
plt.hist(df_a['dincl'],bins=45, color=hist_color)
hspace = 0.018
wspace = 0.018
plt.subplots_adjust(hspace=hspace, wspace=wspace)

plt.savefig('../reports/figures/2dhist_dincl_ratio.png', format='png', dpi=1000)
plt.savefig('../reports/figures/2dhist_dincl_ratio.eps', format='eps')
plt.savefig('../reports/figures/2dhist_dincl_ratio.pdf', format='pdf')
