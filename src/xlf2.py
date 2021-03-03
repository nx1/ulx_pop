import os
import ctypes
from itertools import product

import numpy as np
import matplotlib.pyplot as plt
from tqdm import tqdm

import populations
from ulxlc import ULXLC
from MC import filename_from_mc_param

def filename_from_mc_param(mc_param):
    fn = str(param_mc)
    fn = fn.replace(' ', '')
    fn = fn.replace('{', '')
    fn = fn.replace('}', '')
    
    fn = fn.replace('\'', '')
    fn = fn.replace(':', '-')
    #fn = fn.replace(',', '_')
    return fn

def calc_cum_hist(data, bins):
    """Calculate cumulative histogram.
    Sometimes known as the survival function."""
    hist = np.histogram(data, bins=bins)[0]
    hist_cum = np.cumsum(hist)
    hist_surv = sum(hist) - hist_cum
    # assert(hist.max() == N_sys) # Force full sample in hist
    return hist_surv

# Initialize ULXLC
ulxlc = ULXLC()
N_sys = ulxlc.get_N_sys()

# XLF Settings
f_hist  = [np.histogram, calc_cum_hist]
bin_min = 38
bin_max = 43
bin_width = 0.25
bins = np.arange(bin_min, bin_max, bin_width)
nbins = len(bins)
bin_centers = 0.5 * (bins[:-1] + bins[1:])


# Grid to iterate over
grid_xlf = {'N_mc'      : [1000],
            'N_sys'     : [N_sys],
            'Z'         : ['all'],
            'bh_ratio'  : [0.0, 0.25, 0.5, 0.75, 1.00],
            'dincl_max' : [46, 21],
            'duty_cycle': [1.0, 0.2],
            'L': ['L_iso', 'Lx1', 'L_b', 'L_d', 'L_b_d', 'L_prec'],
            'subset': ['all', 'BH', 'NS', 'alive', 'transient']}

# In Python versions earlier than 3.6 the
# order of the grid_params may not be consistent here.
grid_params   = list(grid_xlf.keys())                 # list of 'str' containing param names
grid_N_params = len(grid_xlf)                         # Number of parameters
grid_iterations = sum([len(i) for i in grid_xlf.values()]) # number of grid iterations
grid_product  = product(*grid_xlf.values()) # Gridsearch Iterator

iters = 0 # iteration counter
Z_last = None # Last metallicity used

# Single MC iteration
for v in grid_product:
    param_mc = {} 
    # Populate param_dict
    for i in range(grid_N_params):
        param_name = grid_params[i]
        param_val  = v[i]
        param_mc[param_name] = param_val
        
    # Filenames for storing output
    fn = filename_from_mc_param(param_mc)
    savepath_hist = f'../data/XLF/{fn}.npy'

    # Check if already simulated
    if os.path.isfile(savepath_hist): 
        print(f'{savepath_hist} exists')
        iters+=1
        continue

    # Check metallicity
    if Z_last != param_mc['Z']:
        # Load population
        print('Loading population')
        df = populations.startrack_v2_mt_1_all(nrows=1000) 
        pop = populations.Population(df)
    
        # Filter population
        pop.filter_non_bh_ns()
        pop.filter_non_thermal_or_nuclear_mt()
        pop.filter_df_ulx_by_Z(param_mc['Z'])
        Z_last = param_mc['Z']

    N_mc = 100
    
    hists = np.ndarray((N_mc,nbins-1), dtype=np.int32)
    # Get L
    # L = get_L()
    for i in tqdm(range(N_mc)):
        df_sampled = pop.sample_systems(param_mc['bh_ratio'],
                                        size=param_mc['N_sys'],
                                        subset='all',
                                        return_df=True)
        L = np.random.normal(loc=1e35, scale=5e37, size=N_sys)
        L = np.log10(L)
        hists[i] = calc_cum_hist(L, bins)
    
    print(hists)
    
    np.save(savepath_hist, hists)
    hist_l = np.load(savepath_hist)
    print(hist_l)
