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

def read_csvs(path):
    results = pd.DataFrame()
    files = glob.glob(path+'*.csv')
    for file in files:
        res = pd.read_csv(file)
        results = pd.concat([results, res])
    results = results.drop('Unnamed: 0', axis=1)
    results.columns = ['Alive', 'BH_ratio', 'Dead', 'Trans']
    return results

def process_results(results):
    res=[]
    states = ['Alive', 'Dead', 'Trans']
    for bh in results['BH_ratio'].unique():
        for state in states:
            cut = results[results['BH_ratio']==bh]
            # percentiles = np.percentile(cut[state], [15.9, 50, 84.1])
            std_dev = np.std(cut[state])
            #1 sig: [15.9, 50, 84.1]
            #2 sig: [2.3, 50, 95.45]
            mean = np.mean(cut[state])
            lower = mean - std_dev
            upper = mean + std_dev
            print(bh, state)
            print('mean:', mean)
            print('std:', std_dev)
            
            res.append([bh, state, mean/500*100, lower/500*100, upper/500*100])
            
            
            # corner.corner(cut[state], quantiles=[0.159, 0.5, 0.841], bins=30, show_titles=True)
            # plt.title(round(bh,2))
    df_res = pd.DataFrame(res)
    df_res.columns=['BH_ratio', 'state', 'mean', 'lower', 'upper']
    return df_res

def plot_results(df_res, Z):
    plt.figure(figsize=(4,3))
    plt.rcParams.update({'font.size': 8})
    
    
    df_res_alive = df_res[df_res['state'] == 'Alive']
    df_res_dead = df_res[df_res['state'] == 'Dead']
    df_res_trans = df_res[df_res['state'] == 'Trans']
    
    # plt.title('Metallicity Z = {}'.format(Z))
    plt.xlabel('$P_{BH}/P_{NS}$')
    plt.ylabel('% of systems')
    plt.ylim(0, 100)

    plt.plot(df_res_alive['BH_ratio'], df_res_alive['mean'],
                  label = 'Alive', c='#9bc53d')

    plt.plot(df_res_dead['BH_ratio'], df_res_dead['mean'],
                  label = 'Dead', c='#e55934')

    plt.plot(df_res_trans['BH_ratio'], df_res_trans['mean'],
                  label = 'Trans', c='#3454d1')

    plt.fill_between(df_res_alive['BH_ratio'],df_res_alive['lower'],df_res_alive['upper'],
                     alpha=0.8, color='#cccccc', label='Alive', hatch='/////')


    plt.fill_between(df_res_dead['BH_ratio'],df_res_dead['lower'],df_res_dead['upper'],
                     alpha=0.8, color='#000000', label='Dead', hatch='+++++')


    plt.fill_between(df_res_trans['BH_ratio'],df_res_trans['lower'],df_res_trans['upper'],
                     alpha=0.8, color='#666666', label='Transient', hatch='-----')

    plt.legend(loc='right')
    plt.tight_layout()
    plt.savefig('../reports/figures/ADT_BHNS_Z_{}.png'.format(Z), format='png', dpi=1000)
    plt.savefig('../reports/figures/ADT_BHNS_Z_{}.eps'.format(Z), format='eps')

results_Z_all = read_csvs('../data/interim/sims/')
results_Z_002 = read_csvs('../data/interim/sims_with_metallicity/0.02/')
results_Z_0002 = read_csvs('../data/interim/sims_with_metallicity/0.002/')
results_Z_00002 = read_csvs('../data/interim/sims_with_metallicity/0.0002/')

df_res_Z_all = process_results(results_Z_all)
df_res_Z_002 = process_results(results_Z_002)
df_res_Z_0002 = process_results(results_Z_0002)
df_res_Z_00002 = process_results(results_Z_00002)

plot_results(df_res_Z_all, 'all')
plot_results(df_res_Z_002, '0.02')
plot_results(df_res_Z_0002, '0.002')
plot_results(df_res_Z_00002, '0.0002')

print('Number of runs for Z=All : {}'.format(len(results_Z_all)/len(results_Z_all['BH_ratio'].unique())))
print('Number of runs for Z=All : {}'.format(len(results_Z_002)/len(results_Z_002['BH_ratio'].unique())))
print('Number of runs for Z=All : {}'.format(len(results_Z_0002)/len(results_Z_0002['BH_ratio'].unique())))
print('Number of runs for Z=All : {}'.format(len(results_Z_00002)/len(results_Z_00002['BH_ratio'].unique())))



