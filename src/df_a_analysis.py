#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Aug 14 16:00:10 2019
@author: nk7g14

This script is used to simulate sampling from a population of ULXs of a given
black hole to neutron star ratio and evaluating the number of ULXs that are
either always visible, sometimes visible or never visible.

It does this by querying the simulation results obtained and saved in the large
dataframe df_a_full.csv (432mb) which contains the simulation outputs from
ulxlc.
"""
import pandas as pd
import numpy as np
import time
import uuid

from auxil import load_systems_dataframe

def ChooseBHNS(BH_RATIO, df_systems):
    """
    Selects a bh or ns system based on a given ratio. 0 = no black holes
    """
    r = np.random.random()
    bh = df_systems[df_systems['is_bh']==1].index.values
    ns = df_systems[df_systems['is_bh']==0].index.values

    if r > BH_RATIO:
        choice = np.random.choice(ns)
        kind = 0
    if r < BH_RATIO:
        choice = np.random.choice(bh)
        kind = 1 #kind = 1 is bh
    return [choice, kind]

def evaluate(ratio, results):
    if ratio==1:
        results['alive']+=1
    if ratio == 0:
        results['dead']+=1
    else:
        results['transient']+=1

#Import csv files
df_a = pd.read_csv('../data/processed/df_a_full.csv')
df_systems = load_systems_dataframe(ulx_only=True, beamed=False, half_opening_l_45=False)


#specify testing parameters
metallicities = [0.02, 0.002, 0.0002]*10
black_hole_ratios = list(np.arange(0, 1.05, 0.05))*100




for Z in metallicities:
    all_results = []
    df_systems_subset = df_systems[df_systems['Z'] == Z]
    
    #converting to numpy arrays for faster calculations
    black_holes = df_systems_subset[df_systems_subset['is_bh']==1].index
    neutron_stars = df_systems_subset[df_systems_subset['is_bh']==0].index
    
    df_a_bh = df_a.loc[df_a['system_num'].isin(black_holes)]
    df_a_ns = df_a.loc[df_a['system_num'].isin(neutron_stars)]
    
    bh_system_numbers = np.array(df_a_bh['system_num'])
    ns_system_numbers = np.array(df_a_ns['system_num'])
    
    ratios = np.array(df_a['ratio'])
    
    
    
    for bh_ratio in black_hole_ratios:
        print('black hole ratio:', bh_ratio, 'Z:', Z)
        results = {'bh_ratio':bh_ratio,
                   'alive':0,
                   'dead':0,
                   'transient':0}
        
        N = 500
        chosen_systems = [ChooseBHNS(bh_ratio, df_systems_subset) for i in range(N)]
        
        t0 = time.time()
        for system, kind in chosen_systems:
            # print('system_number:', system, 'type:', kind)
            if kind == 0:
                try:
                    index = np.random.choice(np.where(ns_system_numbers==system)[0])
                    ratio = ratios[index]
                    evaluate(ratio, results)
                except ValueError:
                    #system not simulated, opening angle > 45
                    results['alive'] += 1
            if kind == 1:
                try:
                    index = np.random.choice(np.where(bh_system_numbers==system)[0])
                    ratio = ratios[index]
                    evaluate(ratio, results)
                except ValueError:
                    #print('system_not_found')
                    results['alive'] += 1
                    
        print('time_taken:', time.time()- t0)
        print(results)
        all_results.append(results)
        print('=============')
        
    filename = str(uuid.uuid4())+'.csv'
    results_df = pd.DataFrame(all_results)
    results_df.to_csv('../data/interim/sims_with_metallicity/{}/{}'.format(Z,filename))