# -*- coding: utf-8 -*-
"""
Created on Thu Jan 28 13:18:40 2021

@author: norma

ULX Population Monte-Carlo Simulation.

This file contains the code to excecute the Monte-Carlo simulation for

"""

import os
import ctypes
from itertools import product
import numpy as np
import pandas as pd


import populations
from ulxlc import MC_input, MC_output, ULXLC
from tqdm import tqdm

def filename_from_mc_param(mc_param):
    fn = str(param_mc)
    fn = fn.replace(' ', '')
    fn = fn.replace('{', '')
    fn = fn.replace('}', '')
    
    fn = fn.replace('\'', '')
    fn = fn.replace(':', '-')
    #fn = fn.replace(',', '_')
    return fn
    

ulxlc = ULXLC()
N_sys = ulxlc.get_N_sys()

grid_mc = {'N_mc'      : [10000],
           'N_sys'     : [N_sys],
           'Z'         : ['all', 0.02, 0.002, 0.0002],
           'bh_ratio'  : [0.0, 0.25, 0.5, 0.75, 1.00],
           'dincl_max' : [46, 21],
           'period'    : ['P_wind_days', 'P_sup_days'],
           'duty_cycle': [1.0, 0.2]}

# In Python versions earlier than 3.6 the
# order of the grid_params may not be consistent here.
grid_params   = list(grid_mc.keys())                 # list of 'str' containing param names
grid_N_params = len(grid_mc)                         # Number of parameters
grid_iterations = sum([len(i) for i in grid_mc.values()]) # number of grid iterations
grid_product  = product(*grid_mc.values()) # Gridsearch Iterator

iters = 0 # iteration counter
Z_last = None # Last metallicity used

# Single MC iteration
for v in grid_product:
    param_mc = {} # Array for storing current iteration sim params
    
    # Populate param_dict
    for i in range(grid_N_params):
        param_name = grid_params[i]
        param_val  = v[i]
        param_mc[param_name] = param_val
        
    # Filenames for storing output
    fn = filename_from_mc_param(param_mc)
    savepath_counts = f'../data/MC/{fn},counts.csv'
    savepath_erass = f'../data/MC/{fn},erass.csv'
    print(fn)
    
    # Check if already simulated
    if os.path.isfile(savepath_counts) and  os.path.isfile(savepath_erass):
        print(f'{savepath_counts} exists')
        print(f'{savepath_erass} exists')
        iters+=1
        continue

    if Z_last != param_mc['Z']:
        # Load population
        print('Loading population')
        df = populations.startrack_v2_mt_1_all() 
        pop = populations.Population(df)
    
        # Filter population
        pop.filter_non_bh_ns()
        pop.filter_non_thermal_or_nuclear_mt()
        pop.filter_df_ulx_by_Z(param_mc['Z'])
        Z_last = param_mc['Z']
    
    # Create arrays for storing outputs
    counts = [0] * param_mc['N_mc']
    erass  = [0] * param_mc['N_mc']
    
    # Main Loop
    for i in tqdm(range(param_mc['N_mc'])):
        # Sample population
        df_sampled = pop.sample_systems(param_mc['bh_ratio'], size=param_mc['N_sys'], subset='ulx', return_df=True)
        
        # Factor in duty cycle
        if param_mc['duty_cycle'] != 1.0:
            df_sampled['rand'] = np.random.random(size=param_mc['N_sys'])
            df_sampled['Lx1'] = np.where((df_sampled['lmxrb']==1)
                                       & (df_sampled['rand']>param_mc['duty_cycle']), 0, df_sampled['Lx1'])
        
        # Create input and output arrays
        inp = MC_input.from_df(df_sampled, param_mc['dincl_max'], param_mc['period'])
        out = MC_output.initialize()
        
        # Call to C
        ulxlc.libc.sim(ctypes.byref(inp), ctypes.byref(out))
        
        # Collect output
        out.collect()
        
        # Save outputs to arrays
        counts[i] = out.res_counts
        erass[i]  = out.res_erass
        
    
    # Create DataFrames from output arrays
    df_counts = pd.DataFrame(counts)
    df_erass = pd.DataFrame(erass)
    
    df_counts.to_csv(savepath_counts, index=False)
    df_erass.to_csv(savepath_erass, index=False)
    
    print(f'iters={iters} \t {param_mc} \t fn={fn}')
    iters+=1
