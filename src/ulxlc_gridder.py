# -*- coding: utf-8 -*-
"""
Created on Wed Dec  9 15:47:42 2020

@author: norma
"""
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import ctypes

import populations

from ulxlc import ULXLC


save_every = 1000
N_mc = 10000    # Number of mc iterations
bh_ratio = 0.5  # Black hole ratio
N = 500         # Sample Size
Z = 0.02        # Metallicity to use

columns_to_use = ['Lx1', 'theta_half_deg']


# Load population
df = populations.startrack_v2_mt_1_all(nrows=10000)
pop = populations.Population(df)        # 
pop.filter_df_ulx_by_Z(Z)               # Filter Z
idx = pop.sample_ulxs(bh_ratio, size=N) # Sampled ids
df_sampled = pop.df.loc[idx]            # Sampled Rows


Lx = (ctypes.c_double * N)(*df_sampled.Lx1.values)
thetas = (ctypes.c_double * N)(*df_sampled.theta_half_deg.values)

ulxlc = ULXLC(100, 1.0)
ulxlc.grid_incl_dincl(Lx, thetas)


ulxlc.
res = pd.DataFrame()
res['N_mc'] = N_mc
res['N'] = N
res['bh_ratio'] = bh_ratio
res['Z'] = Z

res['Lx'] = Lx
res['theta'] = theta

res['dincl'] = dincl
res['incl'] = incl
res['lc_zero_incl_max_flux'] = 
res['lc_max_flux'] = 
res['lc_min_flux'] =
res['lc_ulx_lim'] = 
res['lc_flux_scaling_constant'] = 
res['classification'] = classification


# Initialize ulxlc with 100 arr size 1.0 stepsize

