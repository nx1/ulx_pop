#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Aug 14 16:00:10 2019

@author: nk7g14
"""

import pandas as pd
import numpy as np
import corner
import matplotlib.pyplot as plt
import time
from tqdm import tqdm

def LoadSystemsDataframe():
    """
    Load all startrack simulation outputs.
    """
    df_master = pd.read_csv('dataframe.csv') #36420 SYSTEMS
    df = df_master[df_master['Lx'] > 1E39]  #Only ULX -     992 ULXs
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


def CreateBHNSdict(df_bh, df_ns):
    """
    create dictionary with index keys corresponding to ns or bh
    """
    bh_ns_dict = {}
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
        self.params = df_systems.iloc[num]
        
    def isBeamed(self):
        """
        Checks if system displays beaming or not
        """
        if self.params['theta_half_deg'] >= 45:
            return False
        else:
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


class Simulate:
    def __init__(self):
        self.N_alive = 0
        self.N_dead = 0
        self.N_transient = 0
        self.df_a = pd.read_csv('df_a_full.csv')
        
    def Alive_Dead_Transient(self, BH_RATIO):
        """
        Simulate choosing N number of ULXs and calculating alive/dead analysis
        returns the number of alive, transient and dead systems from the simulation.
        """
        N = 500 #Number of ULX Draws
        chosen_systems = [System(ChooseBHNS(BH_RATIO)) for i in range(N)]
        for system in chosen_systems:
            if system.isBeamed():
                result = system.GetSingleSimulation()
                if np.sum(result['ratio']) == 0:
                    self.N_dead += 1
                elif np.sum(result['ratio']) == 45:
                    self.N_alive += 1
                else:
                    self.N_transient += 1
            else:
                self.N_alive +=1



def Simulate(BH_RATIO):
    """
    Simulate choosing N number of ULXs and calculating alive/dead analysis
    returns the number of alive, transient and dead systems from the simulation.
    """
    NUMBER_OF_ULX_DRAWS = 500
    chosen_systems = [System(ChooseBHNS(BH_RATIO)) for i in range(NUMBER_OF_ULX_DRAWS)]
    
    N_alive = 0
    N_dead = 0
    N_transient = 0
    transient_results = pd.DataFrame()
    
    for system in chosen_systems:
        if system.isBeamed():
            result = system.GetSingleSimulation()
            # result = result.sort_values(by=['dincl'])
            # transient_results = pd.concat([transient_results, result])
            if np.sum(result['ratio']) == 0:
                N_dead += 1
            if np.sum(result['ratio']) == 45:
                N_alive += 1
        else:
            N_alive +=1
    N_transient = NUMBER_OF_ULX_DRAWS - N_alive - N_dead
    return [BH_RATIO, N_alive, N_dead, N_transient]

          
def MCMC_BHNS_Number():
    results = []
    
    for i in tqdm(range(100)):
        for BH_RATIO in np.arange(0,1.1,0.1):
            result = Simulate(BH_RATIO)
            results.append(result)

    results_df = pd.DataFrame.from_records(results)
    results_df.columns = ['BH_NS', 'N_alive', 'N_dead', 'N_trans']
    return results_df


def ProcessAndPlotResults(results):
    n_alives = [] 
    n_deads = []
    n_trans = []
    n_alives_std = [] 
    n_dead_std = []
    n_trans_std = []
    
    for BH_RATIO in np.arange(0,1.1,0.1):
        r_df = results_df[results_df['BH_NS']==BH_RATIO]
        
        mean = r_df.mean()
        std = r_df.std()
        
        n_alives.append(mean['N_alive'])
        n_alives_std.append(std['N_alive'])
        n_deads.append(mean['N_dead'])
        n_dead_std.append(std['N_dead'])
        n_trans.append(mean['N_trans'])
        n_trans_std.append(std['N_trans'])
    
    import matplotlib
    matplotlib.rcParams['mathtext.fontset'] = 'stix'
    matplotlib.rcParams['font.family'] = 'STIXGeneral'    
    plt.figure(figsize=(5,4)) 
    plt.errorbar(np.arange(0,1.1,0.1), n_alives, yerr=n_alives_std, label='Alive')
    plt.errorbar(np.arange(0,1.1,0.1), n_deads, yerr=n_dead_std, label='Dead')
    plt.errorbar(np.arange(0,1.1,0.1), n_trans, yerr=n_trans_std, label='Transient')
    plt.xlabel('BH/NS ratio')
    plt.ylabel('Number')
    plt.legend()
    plt.savefig('bhns_number.eps', format='eps', bbox_inches = "tight")


def plotdfHist(df, name):
    import matplotlib
    matplotlib.rcParams['mathtext.fontset'] = 'stix'
    matplotlib.rcParams['font.family'] = 'STIXGeneral'
    plt.figure(figsize=(5.5,5.5))
    plt.title(name)
    plt.xlabel('Precession angle $\Delta i$')
    plt.ylabel('Alive / Dead Ratio')
    x, y = df['dincl'].values, df['ratio'].values
    color = df['is_bh'].values
    # corner.hist2d(x, y, bins=20)
    # plt.scatter(x,y,s=0.1,c=color)
    plt.hist2d(x,y, bins=45, cmap='viridis')
    plt.colorbar()
    # plt.savefig('dincl_vs_ratio.eps', format='eps', bbox_inches = "tight")


#GLOBAL VARIABLES
df_a = pd.read_csv('df_a_full.csv')
df_systems = LoadSystemsDataframe()


df_bh, df_ns = FilterNSBH(df_systems)

bh_ns_dict = CreateBHNSdict(df_bh, df_ns)

df_a = AppendIsBHtoDataFrame(df_a, bh_ns_dict)

df_a_ns = df_a[df_a['is_bh'] == 0]
df_a_bh = df_a[df_a['is_bh'] == 1]

df_a_ns_nonzero = FilterNonzero(df_a_ns)
df_a_bh_nonzero = FilterNonzero(df_a_bh)
df_a_nonzero = FilterNonzero(df_a)


results = []
for BH_NS in np.arange(0.1, 1.0, 0.05):
    results.append(Main(BH_NS))

ProcessAndPlotResults(results)

plotdfHist(df_results, 'results')






plotdfHist(df_a_ns, 'ns')
plotdfHist(df_a_bh, 'bh')

plotdfHist(df_a_ns_nonzero, 'ns nonzero')
plotdfHist(df_a_bh_nonzero, 'bh nonzero')
