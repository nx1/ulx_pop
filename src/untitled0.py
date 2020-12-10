# -*- coding: utf-8 -*-
"""
Created on Thu Nov 19 11:38:55 2020

@author: norma
"""

"""
Panel 1 : Unbeamed emission split NS/BH
Panel 2 : Beamed Emission split NS/BH
Panel 3 : Beamed + duty_cycle (observed) split NS/BH
Panel 4 : Beamed + duty_cycle + precession (observed) split NS/BH
Panel 5 : Beamed + duty_cycle + prcession (observed) combined NS/BH different %_BH
"""
import numpy as np
import matplotlib.pyplot as plt
import sqlite3
import subprocess
from uuid import uuid4

from tqdm import tqdm

import populations
from constants import R_sol


# XLF settings
N = 500

bin_min = 38
bin_max = 44
bin_width = 0.25
bins = np.arange(bin_min, bin_max, bin_width)
nbins = len(bins)
bin_centers = 0.5 * (bins[:-1] + bins[1:])



print('loading csv')
df = populations.startrack_v2_mt_1_all()
pop = populations.Population(df)

cols = ['idum_run', 'iidd_old', 't', 'dt', 'K_a', 'mttype', 'log_Lx_iso', 'log_Lx1', 'theta_half_deg', 'lmxrb']

df2 = pop.df[cols]


df2


gb = pop.gb_sys(df2)
gb

from tqdm import tqdm
tqdm.pandas()


def sample(df):
    samp = np.random.choice(df['log_Lx_iso'], size=10000)
    return samp


potato = gb.progress_apply(sample)
val = np.stack(potato.values, axis=0)
val = val.T


hists = [np.histogram(val[i], bins=pop.bins)[0] for i in range(len(val))]

mean = np.mean(hists, axis=0)
std = np.std(hists, axis=0)

plt.errorbar(bin_centers, mean, yerr=std, )


