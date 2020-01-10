#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jan 10 15:17:45 2020

@author: nk7g14
"""

#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Aug 14 16:00:10 2019

@author: nk7g14

This file is used to calculate the number of alive/dead transient ULXs for a
given ratio of black holes to neutron stars.

It does this by querying the simulation results obtained and saved in the large
dataframe df_a_full.csv (432mb) which contains the simulation outputs from
ulxlc.
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import time
from multiprocessing import Pool


def LoadSystemsDataframe():
    """
    Load all startrack simulation outputs.
    """
    df_master = pd.read_csv('../data/processed/all_systems_df.csv') #36420 SYSTEMS
    df = df_master[df_master['Lx'] > 1E39]    #Only ULX -     992 ULXs
    # df = df[df['b'] < 1]                    #Only Beamed -  227 Beamed ULXs
    # df = df[df['theta_half_deg'] < 45]      #thetha < 45 -  151 Beamed ULXs with half opening angles < 45
    df = df.reset_index()
    df = df.drop(columns=['index', 'Unnamed: 0'])
    return df


def FilterNSBH(df):
    """
    Filter a given dataframe into two dataframes containing BH and NS
    """
    df_bh = df[df['is_bh'] == 1]
    df_ns = df[df['is_bh'] == 0]
    return df_bh, df_ns


def ChooseBHNS(BH_RATIO):
    """
    Selects a bh or ns system based on a given ratio. 0 = no black holes
    """
    r = np.random.random()
    ns = df_ns.index.values
    bh = df_bh.index.values

    if r > BH_RATIO:
        choice = np.random.choice(ns)
    if r < BH_RATIO:
        choice = np.random.choice(bh)
    return choice


    
def Alive_Dead_Transient(BH_RATIO):
    """
    Simulate choosing N number of ULXs and calculating alive/dead analysis
    returns the number of alive, transient and dead systems from the simulation.
    """
    N_alive = 0
    N_dead = 0
    N_transient = 0
    
    N = 500 #Number of ULX Draws
    
    print('Simulating {} ULXs | BH_RATIO = {}'.format(N, BH_RATIO))
    t0 = time.time()
    for system in chosen_systems:
        if system.isBeamed():
            result = system.GetSingleRow()
            if result['ratio'].values == 0:
                N_dead += 1
            elif result['ratio'].values == 1.0:
                N_alive += 1
            else:
                N_transient += 1
        else:
            N_alive +=1

    time_taken = round((time.time() - t0), 2)

    results = {'BH_ratio':BH_RATIO,
               'Alive':N_alive,
               'Dead':N_dead,
               'Trans':N_transient}

    print('Done! Time Taken: {}'.format(time_taken))
    return results




# GLOBAL VARIABLES
df_a = pd.read_csv('../data/processed/df_a_full.csv')
df_systems = LoadSystemsDataframe()
METALLICITY = 0.02
df_systems = df_systems[df_systems['Z'] == METALLICITY]

df_bh, df_ns = FilterNSBH(df_systems)

BH_RATIO = list(np.arange(0, 1.05, 0.05))*100


for bh in BH_RATIO:
    alive = 0
    dead = 0 
    trans = 0
    print('black hole ratio:', bh)
    NUMBER_OF_SYSTEMS = 500
    chosen_systems = [ChooseBHNS(bh) for i in range(NUMBER_OF_SYSTEMS)]
    t0 = time.time()
    for sys_num in chosen_systems:
        subset = df_a.loc[df_a['system_num'] == sys_num]
        if len(subset) == 0:
            alive+=1
        else:
            result = subset.sample(1)
            if result['ratio'].values == 0:
                dead += 1
            elif result['ratio'].values == 1.0:
                alive += 1
            else:
                trans += 1
        #print(row)
    time_taken = round((time.time() - t0), 2)
    print(alive, dead, trans)
    print('Done! Time Taken: {}'.format(time_taken))
    
    

if __name__ == '__main__':
    
    # p = Pool(2)
    # results = p.map(Alive_Dead_Transient, BH_RATIO)
    results = [Alive_Dead_Transient(bh) for bh in BH_RATIO]
    
    results_df = pd.DataFrame.from_dict(results)
    random_name = 'sim_'+str(np.round(np.random.random(),4))+'.csv'
    results_df.to_csv('../data/interim/sims_with_metallicity/{}/{}'.format(METALLICITY,random_name))


