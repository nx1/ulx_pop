"""
xlf2.py

setup ulxlc
setup xlf settings
setup grid to run

while in the grid:
    - Load population
    - filter population
    - Sample population 10000x
    - For each sample:

	       in simulation: 'Lxmt', 'Lx', 'L_a', 'L_b','Lx_iso', 'Lx1',
		   Additional:    '
        calculate the luminosities that aren't calculated.

			- Beamed (b*d)		 : log_Lx1_vis     : Lx1_b_d         calc_Lx1_vis()
			- Precession (ulxlc) : log_Lx1_prec    : Lx1_prec        calc_Lx1_prec()
			- Precession (d)	 : log_Lx1_prec_vis
 
    # Luminosities to create histograms forr
    #         Isotropic      Beamed    Beamed (b*d)    precession      precession (d)
    grid_L = ['log_Lx_iso', 'log_Lx1', 'log_Lx1_vis', 'log_Lx1_prec', 'log_Lx1_prec_vis']

           

        for each L: #Lx, Lx1, Lx_iso, L_
            
       calc_(L)     


       - calculate d_vis
       - calculate p_obs
       - calcualte 1 - p_obs
       - calculate r
       -        
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
    df = populations.startrack_v2_mt_1_all(nrows=10000) #nrows=10000
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

    # XLF Settings
    bin_min = 38
    bin_max = 43
    bin_width = 0.25
    bins = np.arange(bin_min, bin_max, bin_width)
    nbins = len(bins)
    bin_centers = 0.5 * (bins[:-1] + bins[1:])


    # Grid to iterate over
    grid_xlf = {'N_mc'      : [11],
                'N_sys'     : [N_sys],
                'Z'         : ['all', 0.02, 0.002, 0.0002],
                'bh_ratio'  : [0.0, 0.25, 0.5, 0.75, 1.00],
                'dincl_max' : [46, 21],
                'duty_cycle': [0.2],
                'pop_subset': ['all']} 

    # subset of sampled population to create xlfs for
    xlf_subsets = ['all', 'ns', 'bh', 'lmxrb', 'alive', 'trans']

    # Luminosities to create histograms for
    #         Isotropic      Beamed    Beamed (b*d)    precession      precession (d)
    grid_L = ['log_Lx_iso', 'log_Lx1', 'log_Lx1_vis', 'log_Lx1_prec', 'log_Lx1_prec_vis']

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
        print(f'v={v}\td={d}\titers={iters}\tZ_last={Z_last}')

        # Filenames for storing output
        fn = filename_from_param_dict(d)
        
        print('SORT THIS:')
        savepath_hist = f'../data/XLF/{fn},xlf_subset-lmxrb,L-log_Lx1.npy'
        print(savepath_hist)
        
        # Check if already simulated
        if os.path.isfile(savepath_hist): 
            print(f'{savepath_hist} exists')
            iters+=1
            continue

        # Load population while checking for already loaded 
        if Z_last != d['Z']:
            pop = load_population(d)
            Z_last = d['Z']
      
        # Create array for storing N_mc x N_sys luminosities
        luminosities = {}

        # Create array for storing 10000 histograms of size nbins-1
        hists = {}
       
        for sub in xlf_subsets:
            for L in grid_L:
                key = f'xlf_subset-{sub},L-{L}'
                hists[key] =  np.ndarray((d['N_mc'],nbins-1), dtype=np.int32)
                luminosities[key] = np.ndarray((d['N_mc'], d['N_sys']), dtype=np.float)

        #hists = np.ndarray((d['N_mc'],nbins-1), dtype=np.int32)
        
        for i in tqdm(range(d['N_mc'])):
            # sample systems
            df_sampled = pop.sample_systems(d['bh_ratio'],
                                            size=d['N_sys'],
                                            subset=d['pop_subset'],
                                            return_df=True)

            df_sampled['d']           = np.where(df_sampled['lmxrb']==1, d['duty_cycle'], 1.0)      # Duty cycle
            df_sampled['p_obs']       = df_sampled['b'] * df_sampled['d']                           # Observation probability (b*d)
            df_sampled['1-p_obs']     = 1 - df_sampled['p_obs']                                     # Non-Observation prob
            df_sampled['r']           = np.random.random(size=d['N_sys'])                           # Random number for visibility
            df_sampled['vis_d']       = df_sampled['r'] < df_sampled['d']                           # duty cycle visibility roll
            df_sampled['vis_b_d']     = df_sampled['r'] < df_sampled['p_obs']                       # duty cycle + beaming visbility roll
            df_sampled['inclination'] = np.random.randint(0, 91, size=d['N_sys'])                   # Random Inclination
            df_sampled['dincl']       = np.random.randint(0, d['dincl_max'], size=d['N_sys'])       # Random precessional angle
            df_sampled['Lx1_vis']     = np.where(df_sampled['vis_b_d']==True, df_sampled['Lx1'], 0) # Visible luminosity
            df_sampled['log_Lx1_vis'] = np.log10(df_sampled['Lx1_vis']) 
            
            # Calculate new luminosity from precession and get classification
            Lx_prec, lc_classification = calc_Lx_prec(df_sampled, ulxlc, d['N_sys'])
            
            # Add calculated values to df
            df_sampled['lc_classification'] = lc_classification
            df_sampled['log_Lx1_prec']     = np.log10(Lx_prec)
            df_sampled['log_Lx1_prec_vis'] = np.where(df_sampled['vis_d']==True,
                                                      df_sampled['log_Lx1_prec'],
                                                      0)

            print(df_sampled['lc_classification'].value_counts())

            df_ns    = df_sampled[df_sampled['K_a'] == 13]
            df_bh    = df_sampled[df_sampled['K_a'] == 14]
            df_lmxrb = df_sampled[df_sampled['lmxrb']==1]
            df_alive = df_sampled[(df_sampled['lc_classification'] == 2) | (df_sampled['lc_classification']==3)]
            df_trans = df_sampled[df_sampled['lc_classification'] == 1]

            subsets = {'all'  :df_sampled,
                       'ns'   :df_ns,
                       'bh'   :df_bh,
                       'lmxrb':df_lmxrb,
                       'alive':df_alive,
                       'trans':df_trans}

            # iterate over all subsets
            for k, df in subsets.items():
                # Create histograms for each L
                for L in grid_L:
                    print(df)
                    luminosities[L][i] = df[L].values
                    print(f'L={L}') 
                    print(f'{df[L].values}')
                    print(luminosities)
                    sub_L = f'xlf_subset-{k},L-{L}' # key to store eg 'all-log_Lx1'
                    print(f'sub_L = {sub_L}')
                    print(df[['K_a', 'vis_d', 'vis_b_d', 'lc_classification', 'log_Lx1_vis', 'log_Lx1_prec', 'log_Lx1_prec_vis']])
                    hists[sub_L][i] = calc_cum_hist(df[L], bins)
                    print(hists[sub_L][i])
        
        # Save histograms for each L
        for k, hist in hists.items():
            savepath_hist = f'../data/XLF/{fn},{k}.npy'
            print('SAVEPATH_HIST=', savepath_hist)
            print(savepath_hist)
            np.save(savepath_hist, hists[k])
            iters+=1
        #hist_L = np.load(savepath_hist)
        #print(hist_L)

