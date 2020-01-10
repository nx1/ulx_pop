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


def CreateBHNSdict(df):
    """
    create dictionary with index keys corresponding to ns or bh
    """
    bh_ns_dict = {}
    df_bh, df_ns = FilterNSBH(df)
    for i in df_bh.index:
        bh_ns_dict[i] = 1

    for i in df_ns.index:
        bh_ns_dict[i] = 0
    return bh_ns_dict


def AppendIsBHtoDataFrame(df_a, bh_ns_dict):
    """
    Append 'is black hole' column to alive dataframe.
    """
    ns_bh_list = []
    for i in df_a['system_num']:
        ns_bh_list.append(bh_ns_dict[i])
    df_a['is_bh'] = ns_bh_list
    return df_a


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


def FilterNonzero(df):
    """
    Remove all simulations with always on or always off alive/dead times
    """
    df_nonzero = df[df['ratio'] != 0]
    df_nonzero = df_nonzero[df_nonzero['ratio'] != 1]
    return df_nonzero
   
    

class System:
    def __init__(self, num):
        self.num = num
        self.params = df_systems.loc[num]
        
    def isBeamed(self):
        """
        Checks if system displays beaming or not
        """
        if self.params['theta_half_deg'] >= 45:
            return False
        return True
        
    def GetSingleSimulation(self):
        """
        Fetches the output for an abritrary simulation corresponding to a
        specific inclination and range of precessional angles from 1 to 45.
        """
        self.df_a = df_a[df_a['system_num']==self.num]
        r_mcmc = np.random.choice(self.df_a['MCMC_iter'].unique())
        self.df_a = self.df_a[self.df_a['MCMC_iter'] == r_mcmc]
        
        if len(self.df_a) != 45:
            r_sim_num = np.random.choice(self.df_a['sim_num'].unique())
            self.df_a = self.df_a[self.df_a['sim_num'] == r_sim_num]
        
        if len(self.df_a) != 45:
            r_inclination = np.random.choice(self.df_a['inclination'].unique())
            self.df_a = self.df_a[self.df_a['inclination'] == r_inclination]
        assert len(self.df_a) == 45
        return self.df_a
    
    def GetSingleRow(self):
        self.df_a = df_a[df_a['system_num']==self.num]
        row = self.df_a.sample(n=1)
        return row

def Alive_Dead_Transient(BH_RATIO):
    """
    Simulate choosing N number of ULXs and calculating alive/dead analysis
    returns the number of alive, transient and dead systems from the simulation.
    """
    N_alive = 0
    N_dead = 0
    N_transient = 0
    
    N = 500 #Number of ULX Draws
    
    chosen_systems = [System(ChooseBHNS(BH_RATIO)) for i in range(N)]
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

if __name__ == '__main__':
    BH_RATIO = list(np.arange(0, 1.05, 0.05))*100
    # p = Pool(2)
    # results = p.map(Alive_Dead_Transient, BH_RATIO)
    results = [Alive_Dead_Transient(bh) for bh in BH_RATIO]
    
    results_df = pd.DataFrame.from_dict(results)
    random_name = 'sim_'+str(np.round(np.random.random(),4))+'.csv'
    results_df.to_csv('../data/interim/sims_with_metallicity/{}/{}'.format(METALLICITY,random_name))


