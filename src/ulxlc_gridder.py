# -*- coding: utf-8 -*-
"""
Created on Wed Dec  9 15:47:42 2020

@author: norma
"""
import numpy as np
import matplotlib.pyplot as plt

import populations

from ulxlc import ULXLC


N_mc = 10000    # Number of mc iterations
bh_ratio = 0.5  # Black hole ratio
N = 500         # Sample Size
save_every = 1000
Z = 0.002       # Metallicity to use

# Load population
df = populations.startrack_v2_mt_1_all()
pop = populations.Population(df)
pop.filter_df_ulx_by_Z(Z)
idx = pop.sample_ulxs(bh_ratio, size=N) # Sampled ids
df_sampled = pop.df.loc[idx]            # Sampled Rows


ulxlc = ULXLC(100, 1.0) #Initialize ulxlc with 100 arr size 1.0 stepsize





        
