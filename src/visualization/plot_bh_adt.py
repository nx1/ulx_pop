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
from scipy.optimize import curve_fit


def read_csvs(path):
    results = pd.DataFrame()
    files = glob.glob(path+'*.csv')
    for file in files:
        res = pd.read_csv(file)
        results = pd.concat([results, res])
    results = results.drop('Unnamed: 0', axis=1)
    results.columns = ['Alive', 'BH_ratio', 'Dead', 'Trans']
    return results

def process_ADT_results(results):
    """Process the alive, dead and transient results as a function of the total
    population.
    """
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
#            print(bh, state)
#            print('mean:', mean)
#            print('std:', std_dev)
            
            res.append([bh, state, mean/500*100, std_dev/500*100])
            
            
            # corner.corner(cut[state], quantiles=[0.159, 0.5, 0.841], bins=30, show_titles=True)
            # plt.title(round(bh,2))
    df_res = pd.DataFrame(res)
    df_res.columns=['BH_ratio', 'state', 'mean', 'error']
    return df_res


def process_AT_results(results):
    """Process only Alive & Transient results
    """
    res=[]
    states = ['Alive', 'Trans']
    for bh in results['BH_ratio'].unique():
        for state in states:
            cut = results[results['BH_ratio'] == bh]
            state_avg_number_of_systems = np.mean(cut[state])
            alive_and_dead_avg_number_of_systems = np.mean(cut['Alive'] + cut['Trans'])
            
            state_std_dev = np.std(cut[state])

            res.append([bh,
                        state,
                        state_avg_number_of_systems/alive_and_dead_avg_number_of_systems*100,
                        state_std_dev/alive_and_dead_avg_number_of_systems*100])
    df_res = pd.DataFrame(res)
    df_res.columns=['BH_ratio', 'state', 'mean', 'error']
    return df_res


def plot_ADT(df_res, Z, axarr, arr_row, arr_col):
#    axarr[arr].figure(figsize=(4,3))
#    axarr[arr].rcParams.update({'font.size': 8})
#    axarr[arr].title(str(Z))
    alive_color = '#3E9651'
    transient_color = '#DA7C30'
    dead_color = '#396AB1'
    
    alive_fmt = '-'
    dead_fmt =  '--'
    transient_fmt = ':'

    
    df_res_alive = df_res[df_res['state'] == 'Alive']
    df_res_dead = df_res[df_res['state'] == 'Dead']
    df_res_trans = df_res[df_res['state'] == 'Trans']
    
    axarr[arr_row, arr_col].text(x=0.1, y=5, s=f'Z = {Z}')
    axarr[arr_row, arr_col].set_xlabel('$\%_{BH}$')
    axarr[arr_row, arr_col].set_ylabel('% of systems')
    axarr[arr_row, arr_col].set_ylim(0, 100)
    
    axarr[arr_row, arr_col].plot(df_res_alive['BH_ratio'], df_res_alive['mean'],
         label='Alive', color='black',
         linestyle=alive_fmt, linewidth=0.8)
    
    axarr[arr_row, arr_col].plot(df_res_dead['BH_ratio'], df_res_dead['mean'],
         label='Dead', color='black',
         linestyle=dead_fmt, linewidth=0.8)
    
    axarr[arr_row, arr_col].plot(df_res_trans['BH_ratio'], df_res_trans['mean'],
         label='Transient', color='black',
         linestyle=transient_fmt, linewidth=0.8)
    
    
#    axarr[arr_row, arr_col].errorbar(df_res_alive['BH_ratio'], df_res_alive['mean'],
#         yerr=df_res_alive['error'], capsize=1.0, color='black',
#         fmt='none', elinewidth=0.25)
#    
#    axarr[arr_row, arr_col].errorbar(df_res_dead['BH_ratio'], df_res_dead['mean'],
#         yerr=df_res_dead['error'], capsize=1.0, color='black',
#         fmt='none', elinewidth=0.25)
#        
#    axarr[arr_row, arr_col].errorbar(df_res_trans['BH_ratio'], df_res_trans['mean'],
#         yerr=df_res_trans['error'], capsize=1.0, color='black',
#         fmt='none', elinewidth=0.25)
    
    axarr[arr_row, arr_col].legend(loc='right')
#    plt.savefig('../reports/figures/ADT_BHNS_Z_{}.png'.format(Z), format='png', dpi=1000)
#    plt.savefig('../reports/figures/ADT_BHNS_Z_{}.eps'.format(Z), format='eps')
#    plt.savefig('../reports/figures/ADT_BHNS_Z_{}.pdf'.format(Z), format='pdf')
#    


def plot_AT(df_res, Z, axarr, arr_row, arr_col):
    alive_color = '#3E9651'
    transient_color = '#DA7C30'
    
    #Code for these figures are in /reports/investigations.ipynb
    PERCENT_ALIVE_EARNSHAW = 0.8148148148148148 * 100
    PERCENT_ALIVE_EARNSHAW_ERROR = 0.12162617354124897 * 100
    PERCENT_ALIVE_EARNSHAW_UPPER = PERCENT_ALIVE_EARNSHAW + PERCENT_ALIVE_EARNSHAW_ERROR
    PERCENT_ALIVE_EARNSHAW_LOWER = PERCENT_ALIVE_EARNSHAW - PERCENT_ALIVE_EARNSHAW_ERROR
    
    PERCENT_TRANS_EARNSHAW = 0.18518518518518517 * 100
    PERCENT_TRANS_EARNSHAW_ERROR = 0.03840635908519479 * 100
    PERCENT_TRANS_EARNSHAW_UPPER = PERCENT_TRANS_EARNSHAW + PERCENT_TRANS_EARNSHAW_ERROR
    PERCENT_TRANS_EARNSHAW_LOWER = PERCENT_TRANS_EARNSHAW - PERCENT_TRANS_EARNSHAW_ERROR
    
    #Filtering dfs by state
    df_res_alive = df_res[df_res['state'] == 'Alive']
    df_res_trans = df_res[df_res['state'] == 'Trans']
    
    #Curve Fitting
    def func(x, m, c):
        return m*x + c
    
    def func_inv(y, m, c):
        return (y - c) / m

    res_alive = curve_fit(func, df_res_alive['BH_ratio'], df_res_alive['mean'], sigma=df_res_alive['error'])
    res_trans = curve_fit(func, df_res_trans['BH_ratio'], df_res_trans['mean'], sigma=df_res_trans['error'])
    
    m_alive = res_alive[0][0]
    c_alive = res_alive[0][1]
    
    m_trans = res_trans[0][0]
    c_trans = res_trans[0][1]
    
    x = np.arange(0, 1, 0.001)
    y_alive = [func(i, m_alive, c_alive) for i in x]
    y_trans = [func(i, m_trans, c_trans) for i in x]
    
    # Plotting
    
    # Plot errorbars
    axarr[arr_row, arr_col].errorbar(df_res_alive['BH_ratio'], df_res_alive['mean'],
         yerr=df_res_alive['error'], label='Alive', capsize=1.0, color=alive_color,
         fmt='none', elinewidth=0.8)
    
    axarr[arr_row, arr_col].errorbar(df_res_trans['BH_ratio'], df_res_trans['mean'],
         yerr=df_res_trans['error'], label='Transient', capsize=1.0, color=transient_color,
         fmt='none', elinewidth=0.8)
    
    # Plot Line of best fits
    axarr[arr_row, arr_col].plot(x, y_alive, color='black', linewidth=0.8)
    axarr[arr_row, arr_col].plot(x, y_trans, color='grey', linewidth=0.8)
    
    # Plot Labels and annotations
    axarr[arr_row, arr_col].set_ylim(0, 100)
    
    axarr[arr_row, arr_col].set_xlabel('$\%_{BH}$')
    axarr[arr_row, arr_col].set_ylabel('% of systems')

    axarr[arr_row, arr_col].text(x=0, y=5, s=f'Z = {Z}', fontsize=6)
    axarr[arr_row, arr_col].text(x=0, y=10, s=f'alive bf: y = {m_alive:.2f}x + {c_alive:.2f}', fontsize=6)
    axarr[arr_row, arr_col].text(x=0, y=15, s=f'trans bf: y = {m_trans:.2f}x + {c_trans:.2f}', fontsize=6)
    
#    axarr[arr_row, arr_col].text(x=0, y=PERCENT_ALIVE_EARNSHAW+2,
#         s=f'Earnshaw % alive: {PERCENT_ALIVE_EARNSHAW:.2f}+-{PERCENT_ALIVE_EARNSHAW_ERROR:.2f}',
#         fontsize=6)
#    axarr[arr_row, arr_col].text(x=0, y=PERCENT_TRANS_EARNSHAW+2,
#         s=f'Earnshaw % trans: {PERCENT_TRANS_EARNSHAW:.2f}+-{PERCENT_TRANS_EARNSHAW_ERROR:.2f}',
#         fontsize=6)
    
    
    # X values interp
    alive_interp_x = func_inv(PERCENT_ALIVE_EARNSHAW, m_alive, c_alive)
    alive_upper_interp_x = func_inv(PERCENT_ALIVE_EARNSHAW_UPPER, m_alive, c_alive)
    alive_lower_interp_x = func_inv(PERCENT_ALIVE_EARNSHAW_LOWER, m_alive, c_alive)
    
    trans_interp_x = func_inv(PERCENT_TRANS_EARNSHAW, m_trans, c_trans)
    trans_upper_interp_x = func_inv(PERCENT_TRANS_EARNSHAW_UPPER, m_trans, c_trans)
    trans_lower_interp_x = func_inv(PERCENT_TRANS_EARNSHAW_LOWER, m_trans, c_trans)
    
    

    
#    axarr[arr_row, arr_col].hlines(PERCENT_ALIVE_EARNSHAW, xmin=0, xmax=alive_interp_x,
#         color='r', linewidth=0.8, linestyle='--')
    axarr[arr_row, arr_col].hlines(PERCENT_ALIVE_EARNSHAW_UPPER, xmin=0, xmax=alive_upper_interp_x,
         color='b', linewidth=0.8, linestyle='--')
    axarr[arr_row, arr_col].hlines(PERCENT_ALIVE_EARNSHAW_LOWER, xmin=0, xmax=alive_lower_interp_x,
         color='b', linewidth=0.8, linestyle='--')
    
#    axarr[arr_row, arr_col].hlines(PERCENT_TRANS_EARNSHAW, xmin=0, xmax=trans_interp_x,
#         color='r', linewidth=0.8, linestyle='--')
#    axarr[arr_row, arr_col].hlines(PERCENT_TRANS_EARNSHAW_UPPER, xmin=0, xmax=trans_upper_interp_x,
#         color='b', linewidth=0.8, linestyle='--')
#    axarr[arr_row, arr_col].hlines(PERCENT_TRANS_EARNSHAW_LOWER, xmin=0, xmax=trans_lower_interp_x,
#         color='b', linewidth=0.8, linestyle='--')
    
    
    axarr[arr_row, arr_col].vlines(alive_upper_interp_x, ymin=0, ymax=PERCENT_ALIVE_EARNSHAW_UPPER,
         color='b', linewidth=0.8, linestyle='--')
    axarr[arr_row, arr_col].vlines(alive_lower_interp_x, ymin=0, ymax=PERCENT_ALIVE_EARNSHAW_LOWER,
         color='b', linewidth=0.8, linestyle='--')
    
    axarr[arr_row, arr_col].text(x=alive_upper_interp_x, y=2,
         s=f'{alive_upper_interp_x:.2f}', c='red')
    axarr[arr_row, arr_col].text(x=alive_lower_interp_x, y=2,
         s=f'{alive_lower_interp_x:.2f}', c='red')
    
    
    axarr[arr_row, arr_col].legend(loc='right')
    
    
    # Old plotting method
    """
    axarr[arr_row, arr_col].plot(df_res_alive['BH_ratio'], df_res_alive['mean'],
         c=alive_color)

    axarr[arr_row, arr_col].plot(df_res_trans['BH_ratio'], df_res_trans['mean'],
         c=transient_color)
                  
    axarr[arr_row, arr_col].fill_between(df_res_alive['BH_ratio'],df_res_alive['lower'],df_res_alive['upper'],
                     alpha=0.8, color=alive_color, label='Alive')


    axarr[arr_row, arr_col].fill_between(df_res_trans['BH_ratio'],df_res_trans['lower'],df_res_trans['upper'],
                     alpha=0.8, color=transient_color, label='Transient')
    """

def plot_all_ADT():
    import matplotlib
    fontsize = 14
    matplotlib.rcParams['mathtext.fontset'] = 'stix'
    matplotlib.rcParams['font.family'] = 'STIXGeneral'
    fig, axarr = plt.subplots(2,2, figsize=(6,6))
    plt.gcf().subplots_adjust(bottom=0.15)
    plot_ADT(df_adt_Z_all, 'all', axarr, 0, 0)
    plot_ADT(df_adt_Z_002, '0.02', axarr, 0, 1)
    plot_ADT(df_adt_Z_0002, '0.002', axarr, 1, 0)
    plot_ADT(df_adt_Z_00002, '0.0002', axarr, 1, 1)
    plt.tight_layout()
    plt.savefig('../reports/figures/ADT_BHNS_array.png', format='png', dpi=1000)
    plt.savefig('../reports/figures/ADT_BHNS_array.eps', format='eps')
    plt.savefig('../reports/figures/ADT_BHNS_array.pdf', format='pdf')
    
def plot_all_AT():
    import matplotlib
    fontsize = 14
    matplotlib.rcParams['mathtext.fontset'] = 'stix'
    matplotlib.rcParams['font.family'] = 'STIXGeneral'
    fig, axarr = plt.subplots(2,2, figsize=(6,6))
    plt.subplots_adjust(wspace=0.295)
    plot_AT(df_at_Z_all, 'all', axarr, 0, 0)
    plot_AT(df_at_Z_002, '0.02', axarr, 0, 1)
    plot_AT(df_at_Z_0002, '0.002', axarr, 1, 0)
    plot_AT(df_at_Z_00002, '0.0002', axarr, 1, 1)
    plt.savefig('../reports/figures/AT_BHNS_array.png', format='png', dpi=1000)
    plt.savefig('../reports/figures/AT_BHNS_array.eps', format='eps')
    plt.savefig('../reports/figures/AT_BHNS_array.pdf', format='pdf')
    

results_Z_all = read_csvs('../data/interim/sims/')
results_Z_002 = read_csvs('../data/interim/sims_with_metallicity/0.02/')
results_Z_0002 = read_csvs('../data/interim/sims_with_metallicity/0.002/')
results_Z_00002 = read_csvs('../data/interim/sims_with_metallicity/0.0002/')

df_adt_Z_all = process_ADT_results(results_Z_all)
df_adt_Z_002 = process_ADT_results(results_Z_002)
df_adt_Z_0002 = process_ADT_results(results_Z_0002)
df_adt_Z_00002 = process_ADT_results(results_Z_00002)

df_at_Z_all = process_AT_results(results_Z_all)
df_at_Z_002 = process_AT_results(results_Z_002)
df_at_Z_0002 = process_AT_results(results_Z_0002)
df_at_Z_00002 = process_AT_results(results_Z_00002)

#plot_all_ADT()
plot_all_AT()





"""
plot_ADT_stackplot(df_adt_Z_all, 'all')
plot_ADT_stackplot(df_adt_Z_002, '0.02')
plot_ADT_stackplot(df_adt_Z_0002, '0.002')
plot_ADT_stackplot(df_adt_Z_00002, '0.0002')

plot_AT_stackplot(df_at_Z_all, 'all')
plot_AT_stackplot(df_at_Z_002, '0.02')
plot_AT_stackplot(df_at_Z_0002, '0.002')
plot_AT_stackplot(df_at_Z_00002, '0.0002')

"""
print('Number of runs for Z=All : {}'.format(len(results_Z_all)/len(results_Z_all['BH_ratio'].unique())))
print('Number of runs for Z=All : {}'.format(len(results_Z_002)/len(results_Z_002['BH_ratio'].unique())))
print('Number of runs for Z=All : {}'.format(len(results_Z_0002)/len(results_Z_0002['BH_ratio'].unique())))
print('Number of runs for Z=All : {}'.format(len(results_Z_00002)/len(results_Z_00002['BH_ratio'].unique())))






df_res_alive = df_at_Z_all[df_at_Z_all['state'] == 'Alive']





plt.plot()
