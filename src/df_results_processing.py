#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Sep  4 13:19:18 2019

@author: nk7g14

This file is used to run through the simulation outputs obtained from
df_a_analysis.py which are in the form of several csv files

each one containing for a given ratio of Black holes to neutron stars
the number of alive, dead and transient ULXs
"""

import pandas as pd
import glob
import matplotlib.pyplot as plt
import numpy as np

csv_files = glob.glob('../data/interim/sims/*.csv')

results = pd.DataFrame()

for file in csv_files:
    print(file)
    res = pd.read_csv(file)
    results = pd.concat([results, res])
results = results.drop('Unnamed: 0', axis=1)


states = ['Alive', 'Dead', 'Trans']

res=[]
for bh in results['BH_ratio'].unique():
    for state in states:
        cut = results[results['BH_ratio']==bh]
        
        percentiles = np.percentile(cut[state], [15.9, 50, 84.1])
        #1 sig: [15.9, 50, 84.1]
        #2 sig: [2.3, 50, 95.45]
        mean = percentiles[1]
        lower = mean - percentiles[0]
        upper = percentiles[2] - mean
        res.append([bh, state, mean/500*100, lower/500*100, upper/500*100])
        
        
        # corner.corner(cut[state], quantiles=[0.159, 0.5, 0.841], bins=30, show_titles=True)
        # plt.title(round(bh,2))



plt.figure(figsize=(4,3))
plt.rcParams.update({'font.size': 8})
df_res = pd.DataFrame(res)
df_res.columns=['BH_ratio', 'state', 'mean', 'lower', 'upper']
df_res_alive = df_res[df_res['state'] == 'Alive']
df_res_dead = df_res[df_res['state'] == 'Dead']
df_res_trans = df_res[df_res['state'] == 'Trans']

plt.xlabel('$P_{BH}/P_{NS}$')
plt.ylabel('% of systems')

plt.plot(df_res_alive['BH_ratio'], df_res_alive['mean'],
              label = 'Alive', c='#9bc53d')

plt.plot(df_res_dead['BH_ratio'], df_res_dead['mean'],
              label = 'Dead', c='#e55934')

plt.plot(df_res_trans['BH_ratio'], df_res_trans['mean'],
              label = 'Trans', c='#3454d1')

plt.fill_between(df_res_alive['BH_ratio'], 
                 df_res_alive['mean']-df_res_alive['lower'], 
                 df_res_alive['mean']+df_res_alive['upper'],
                 alpha=0.8, color='#cccccc', label='Alive', hatch='/////')
                 

plt.fill_between(df_res_dead['BH_ratio'], 
                 df_res_dead['mean']-df_res_dead['lower'], 
                 df_res_dead['mean']+df_res_dead['upper'],
                 alpha=0.8, color='#000000', label='Dead', hatch='+++++')
                 

plt.fill_between(df_res_trans['BH_ratio'], 
                 df_res_trans['mean']-df_res_trans['lower'], 
                 df_res_trans['mean']+df_res_trans['upper'],
                 alpha=0.8, color='#666666', label='Transient', hatch='-----')
   
plt.legend(loc='right')
plt.tight_layout()

plt.savefig('../reports/figures/ADT_BHNS.png', format='png', dpi=1000)
plt.savefig('../reports/figures/ADT_BHNS.eps', format='eps')
