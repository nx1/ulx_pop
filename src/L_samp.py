"""
L_samp.py

Samples and saves luminosities for later use.

"""


import os
import ctypes
from itertools import product

import numpy as np
import matplotlib.pyplot as plt
from tqdm import tqdm

import populations
from ulxlc import ULXLC

def create_param_dict(v, grid):
    """Create param dictionary from tuple of values.
       tuple: (10000, 500, ...)
       grid:  {'N_mc':[], 'N_sys':[], ...}
        
       returns:
       dict = {'N_mc':10000, 'N_sys':500  }
    """
    param_dict = {} 
    
    # Populate param_dict
    for i in range(len(grid)):
        param_name = grid_params[i]
        param_dict[param_name] = v[i]
    return param_dict
 
    
def filename_from_param_dict(param_dict):
    fn = str(param_dict)
    fn = fn.replace(' ', '')
    fn = fn.replace('{', '')
    fn = fn.replace('}', '')
    fn = fn.replace('\'', '')
    fn = fn.replace(':', '-')
    #fn = fn.replace(',', '_') #Not good if you've got _ in param_name
    return fn

def filename(param_dict, save_subset, L):
    fn = filename_from_param_dict(param_dict)
    fn = f'{fn},save-{save_subset},L-{L}'
    return fn

def calc_Lx_prec(df_sampled, ulxlc, N):                
    """Calculate Lx prec and obtain classification."""
    c_Lx_prec = (ctypes.c_double * N)()    
    c_lc_classification = (ctypes.c_double * N)()    

    c_Lx      = (ctypes.c_double * N)(*df_sampled['Lx1'].values)
    c_thetas  = (ctypes.c_double * N)(*df_sampled['theta_half_deg'].values)
    c_incls   = (ctypes.c_double * N)(*df_sampled['inclination'].values)
    c_dincls  = (ctypes.c_double * N)(*df_sampled['dincl'].values)

    ulxlc.xlf_calc_L_prec(c_Lx_prec, c_lc_classification, c_Lx, c_thetas, c_incls, c_dincls, N)
    Lx_prec = np.ctypeslib.as_array(c_Lx_prec)
    lc_classification = np.ctypeslib.as_array(c_lc_classification)
    return Lx_prec, lc_classification


def calc_cum_hist(data, bins):
    """Calculate cumulative histogram.
    Sometimes known as the survival function."""
    hist = np.histogram(data, bins=bins)[0]
    hist_cum = np.cumsum(hist)
    hist_surv = sum(hist) - hist_cum
    # assert(hist.max() == N_sys) # Force full sample in hist
    return hist_surv

def load_population(d):
    # Load population
    print('Loading population')
    df = populations.startrack_v2_mt_1_all() #nrows=10000
    pop = populations.Population(df)

    # Filter population
    pop.filter_non_bh_ns()
    pop.filter_non_thermal_or_nuclear_mt()
    pop.filter_df_ulx_by_Z(d['Z'])
    return pop


if __name__=='__main__':
    # Initialize ULXLC
    ulxlc = ULXLC()
    N_sys = ulxlc.get_N_sys()
 
    # Grid to iterate over
    grid_xlf = {'N_mc'      : [100],
                'N_sys'     : [N_sys],
                'Z'         : ['all', 0.02, 0.002, 0.0002],
                'bh_ratio'  : [0.0, 0.25, 0.5, 0.75, 1.00],
                'dincl_max' : [46,21],
                'duty_cycle': [0.2],
                'pop_subset': ['all']} 

    # subset of sampled population to save L for
    save_subsets = ['all', 'ns', 'bh', 'lmxrb', 'alive', 'trans']

    # Luminosities to create histograms for
    grid_L = ['Lx_iso', 'Lx1', 'Lx1_vis_b', 'Lx1_vis_b_d', 'Lx1_prec', 'Lx1_prec_vis']


    # In Python versions earlier than 3.6 the
    # order of the grid_params may not be consistent here.
    grid_params   = list(grid_xlf.keys())                      # list of 'str' containing param names
    grid_N_params = len(grid_xlf)                              # Number of parameters
    grid_iterations = sum([len(i) for i in grid_xlf.values()]) # number of grid iterations

    iters = 0     # iteration counter
    Z_last = None # Last metallicity used

    # Single MC iteration
    for v in product(*grid_xlf.values()):
        # v is tuple of param values (10000,500,'all',46,21,1.0,0.2,'L_iso','all')
        # d is a dictionary of current params = {'N_sys':10000, etc..., })
        d = create_param_dict(v, grid_xlf)
        print(f'v={v}')
        print(f'd={d}')
        print(f'iters={iters}')
        print(f'Z_last={Z_last}')
        # Filenames for storing output
        fn = filename(d, save_subsets[0], grid_L[0])
        savepath = f'../data/L_samp/{fn}.npy'

        print(f'savepath{savepath}')
        
        # Check if already simulated
        if os.path.isfile(savepath): 
            print(f'{savepath} exists')
            iters+=1
            continue

        # Load population while checking if already loaded 
        if Z_last != d['Z']:
            pop = load_population(d)
            pop.df = pop.df[['Unnamed: 0', 'original_row', 't', 'dt', 'K_a', 'K_b', 'L_a', 'L_b', 'mttype', 'Lxmt', 'Lx', 'idum_run', 'iidd_old', 'Z', 'Lx_iso', 'b', 'Lx1','theta_half_deg', 'lmxrb']]
            Z_last = d['Z']
      
        # Create dictionary for storing N_mc x N_sys luminosities
        Ls = {}
        
        # populate dictionary with arrays
        for k in save_subsets:
            for L in grid_L:
                key = f'{L}-{k}'
                Ls[key] = np.ndarray((d['N_mc'], d['N_sys']), dtype=np.float)

        # Sample N_mc times
        for i in tqdm(range(d['N_mc'])):
            # sample systems
            df_sampled = pop.sample_systems(d['bh_ratio'],
                                            size=d['N_sys'],
                                            subset=d['pop_subset'],
                                            return_df=True)

            df_sampled['d']           = np.where(df_sampled['lmxrb']==1, d['duty_cycle'], 1.0)       # Duty cycle
            df_sampled['p_obs']       = df_sampled['b'] * df_sampled['d']                            # Observation probability (b*d)
            df_sampled['r']           = np.random.random(size=d['N_sys'])                            # Random number for visibility
            df_sampled['vis_d']       = df_sampled['r'] < df_sampled['d']                            # duty cycle onlu visibility roll
            df_sampled['vis_b']       = df_sampled['r'] < df_sampled['b']                            # Beaming only visbility roll
            df_sampled['vis_b_d']     = df_sampled['r'] < df_sampled['p_obs']                        # duty cycle + beaming visbility roll
            df_sampled['inclination'] = np.random.randint(0, 91, size=d['N_sys'])                    # Random Inclination
            df_sampled['dincl']       = np.random.randint(0, d['dincl_max'], size=d['N_sys'])        # Random precessional angle
            df_sampled['Lx1_vis_b']   = np.where(df_sampled['vis_b_d']==True, df_sampled['Lx1'], 0)  # Visible luminosity (b)
            df_sampled['Lx1_vis_b_d'] = np.where(df_sampled['vis_b_d']==True, df_sampled['Lx1'], 0)  # Visible luminosity (b*d) 
            
            # Calculate new luminosity from precession and get classification
            Lx_prec, lc_classification = calc_Lx_prec(df_sampled, ulxlc, d['N_sys'])
            
            # Add calculated values to df
            df_sampled['lc_classification'] = lc_classification
            df_sampled['Lx1_prec']          = Lx_prec                               # Luminosity /w precession
            df_sampled['Lx1_prec_vis']      = np.where(df_sampled['vis_d']==True,   # Luminosity /w precession + duty cycle
                                                       df_sampled['Lx1_prec'],
                                                       0)
            # iterate over all luminosities
            
            for L in grid_L:
                for k in save_subsets:
                    df = df_sampled.copy()
                    if k=='all':
                        pass
                    elif k=='ns':
                        df[L] = np.where(df['K_a']==13, df[L], 0)
                    elif k=='bh':
                        df[L] = np.where(df_sampled['K_a']==14, df_sampled[L], 0)
                    elif k=='lmxrb':
                        df[L] = np.where(df_sampled['lmxrb']==1, df_sampled[L], 0)
                    elif k=='alive':
                        df[L] = np.where(((df_sampled['lc_classification'] == 2) | (df_sampled['lc_classification']==3)), df_sampled[L], 0)
                    elif k=='trans':
                        df[L] = np.where(df_sampled['lc_classification']==1, df_sampled[L], 0)
                    else:
                        raise KeyError(f'{k} is not a valid save subset')
                    
                    L_val = df[L].values
                    key = f'{L}-{k}'
                    Ls[key][i] = L_val


        for L in grid_L:
            for k in save_subsets:
                fn = filename(d, k, L)
                savepath = f'../data/L_samp/{fn}.npy'
                key = f'{L}-{k}'
                print(d)
                print(f'saved to: {savepath}')
                np.save(savepath, Ls[key])
        iters+=1
