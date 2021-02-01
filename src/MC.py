# -*- coding: utf-8 -*-
"""
Created on Thu Jan 28 13:18:40 2021

@author: norma

pop
Z
N_sys
N_mc
p_BH
duty_cycle
dincl_max
bh_ratio
period
"""
import ctypes
from itertools import product
import numpy as np
import pandas as pd


import populations
from ulxlc import MC_input, MC_output, ULXLC
from tqdm import tqdm



class MonteCarlo:
    def __init__(self):
        pass
    
    def set_grid(self, grid_mc):
        """
        Set parameter grid

        Parameters
        ----------
        grid_mc : TYPE
            Dictionary with parameters names (`str`) as keys and lists of
        """
        self.grid_mc = grid_mc
        
    def sim_mc(self, par):
    """
    Single MC iteration of simulation.

    Parameters
    ----------
    par : dict
        Dictionary of parameters

    Returns
    -------
    None.

    """
    pass

# mc = MonteCarlo()
# mc.set_grid(grid_mc)
# mc.set_population(pop)
# mc.run()
# mc.results()




grid_mc = {'bh_ratio'  : [0.0, 0.25, 0.5, 0.75, 1.00],
           'dincl_max' : [46, 21],
           'Z'         : ['all', 0.02, 0.002, 0.0002],
           'period'    : ['P_wind_days', 'P_sup_days']}

# In Python versions earlier than 3.6 the
# order of the grid_params may not be consistent here.
grid_params   = list(grid_mc.keys()) 
grid_N_params = len(grid_mc)
grid_product  = product(*grid_mc.values())





    

df_info = pd.DataFrame()
for v in grid_product:
    param_mc = {}
    print(v)
    
    # Populate param_dict
    for i in range(grid_N_params):
        param_name = grid_params[i]
        param_val  = v[i]
        param_mc[param_name] = param_val    
    
    print(param_mc)
    
    sim_mc(param_mc)




N_mc = 10000
N_sys = 500
bh_ratio = 0.5
dincl_max = 46
Z = 'all'
period = 'P_wind_days'


par = {'N_sys' : N_sys,
       'bh_ratio' : bh_ratio,
       'dincl_max' : dincl_max,
       'Z' : Z,
       'period' : period}


# Load population
df = populations.startrack_v2_mt_1_all(nrows=50000) 
pop = populations.Population(df)

# Filter population
pop.filter_non_bh_ns()
pop.filter_non_thermal_or_nuclear_mt()
pop.filter_df_ulx_by_Z(Z)

ulxlc = ULXLC()

# Main Loop
for i in tqdm(range(N_mc)):
    # Sample population
    df_sampled = pop.sample_systems(par['bh_ratio'], size=par['N_sys'], subset='ulx', return_df=True)
    
    # Create input and output arrays
    inp = MC_input.from_df(df_sampled, par['dincl_max'], par['period'])
    out = MC_output.initialize()
    
    # Call to C
    ulxlc.libc.sim(ctypes.byref(inp), ctypes.byref(out))
    
    # Collect output
    out.collect()
    
    
    
    