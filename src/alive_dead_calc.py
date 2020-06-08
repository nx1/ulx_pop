#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Feb 14 15:38:42 2019

@author: nk7g14

This script runs through the output curves generated from the xspec model
ulxlc and performs alive/dead time analysis for each of them.

Lightcurves are stored in .txt files

The filename format goes as:
    simulation_number - dincl - inclination.txt

where dincl is the precession angle

1) Load in all lightcurves into df_dict

For each simulation_number dincl combo there is a 0 inclination system and
a random inclination system.

3) Find the 0 inclination system and set the scaling factor, c, and the limiting
value, N_lim.

c = Lx / max(curves[key]['Flux']) #Scaling factor
N_lim = 1E39 / c                  #Limiting value
"""

import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from multiprocessing import Pool
from tqdm import tqdm
import os
from pathlib import Path
import pickle

from auxil import load_systems_dataframe
from curvefunc import calc_alive_time

# =============================================================================
# Functions
# =============================================================================


def GetLx(systems_df, system_id):
    return systems_df.loc[system_id]['Lx']


def PlotCurve(filename):
    #Plots a specific curve based on filename
    system_id, dincl, inclination = split_curve_filename(filename)
    curves = find_curve_by_id_and_dincl(df_dict, system_id, dincl)
    N_lim = Normalise(curves)
    
    fig, ax = plt.subplots()

    for i in curves.keys():
        ax.plot(curves[i]['Time'], curves[i]['Flux'], label=i)

    ax.axhline(y=N_lim, color='r', linestyle='-', label='limit = '+str(N_lim))
    ax.legend()
    return plt.show()


def Normalise(curves, Lx):
    '''
    inputs :
        curves = dictionary of two curves
    Takes two curves and Normalises them based on inclination
    
    '''
    N_lim = None
    for key in curves:
        system_id, dincl, inclination = split_curve_filename(key)
        
        # print('-----------------')
        # print('Curve:', key)
        # print('Lx:', Lx)
        if inclination == '0.0':
            max_flux = max(curves[key]['Flux']) #Maximum Flux
            c = Lx / max_flux                   #Scaling factor
            N_lim = 1E39 / c                    #Limiting value
            # print('Curve Maximum:', max_flux)
            # print('Normalisation factor c:', c)
            # print('Resulting N_lim:', N_lim)
    return N_lim


def ResultsDictToPandas(r_dict):
    incls = []
    dincls_list = []
    sim_num_list = []
    df_a = pd.DataFrame.from_dict(results_dict, orient='index')

    for key, row in tqdm(df_a.iterrows()):
        
        split_key = key.split('-')
        sim_num = split_key[0]
        inclination = split_key[-1]
        dincl = split_key[1]
        
        incls.append(float(inclination))
        dincls_list.append(float(dincl))
        sim_num_list.append(int(sim_num))
    df_a['system_num'] = sim_num_list
    df_a['inclination'] = incls
    df_a['dincl'] = dincls_list

    
    df_a.columns = ['alive', 'dead','MCMC_iter', 'BH_NS', 'sim_num', 'system_num', 'inclination', 'dincl']
    df_a['ratio'] = df_a['alive']/(df_a['dead']+df_a['alive'])
    return df_a


def CalculateNormalizationLimit(zero_inclination_curve, Lx):
    '''
    Calculates the value a lightcurve for a paticular system and precession angle
    must exceed in order to be classed as a ULX.
    
    This is done by taking the maximum flux from the lightcurve at 0 inclination and
    normalizing it to it's Lx value and 1e39 erg/s
    '''
    max_flux = max(zero_inclination_curve['Flux'])
    c = Lx / max_flux
    N_lim = 1E39 / c
    return N_lim

def calc_alive_dead_curve(key):
    system_id, dincl, inclination = split_curve_filename(key)
    
    N_lim = norm_lookup[(system_id, dincl)]
    if N_lim == None:
        Alive, Dead = 1, 0
    else:
        Alive, Dead = calc_alive_time(df_dict[key], N_lim)
    return Alive, Dead, MCMC_iteration, BH_NS, simulation_number





def main():
    # =============================================================================
    # Loading curves & Systems df
    # =============================================================================
    results_dict = {}
    systems_df = load_systems_dataframe(True, False, False)
    Lx_arr = systems_df['Lx']
    
    #norm_lookup contains a lookup table containing N_lim where we have simulated
    #for each of the 991 ULXs 0 inclination curves over all precessional angles
    norm_lookup_path = Path('../data/interim/norm_lookup.pickle')
    with open(norm_lookup_path, 'rb') as handle:
        norm_lookup = pickle.load(handle)
    
    
    ulxlc_folder = 'ulxlc'
    
    pbar = tqdm(range(3, 100))
    
    """
    MCMC_ITERATION
        |____BH_NS_RATIO
                |_________SIMULATION_NUMBER
                            |________________CURVE_FILE.TXT
    """
    
    
    #Loop over all MCMC iterations
    for MCMC_iteration in pbar:
        #Set working directory to the MCMC Folder
        working_dir = '{}/curves/{}'.format(ulxlc_folder, MCMC_iteration)
        
        #Find all the subfolders within the MCMC folder corresponding to BHNS ratios
        BH_NS_folders = os.listdir(working_dir)
        #Loop over all BHNS folders
        for BH_NS in BH_NS_folders:
            #Loop over all 500 folders in each BH_NS folder
            for simulation_number in range(500):
                #for each folder we need to load all the curves and perform
                # alive/dead time analysis on them all
                
                curve_folder = '{}/{}/{}'.format(working_dir, BH_NS, simulation_number)
                #Load up all the curves in the folder
                df_dict = load_all_curves_from_path(curve_folder)
    
                pbar.set_description('%s %s %s' % (MCMC_iteration, BH_NS, simulation_number))
    
                p = Pool()
                results = p.map(calc_alive_dead_curve, df_dict.keys())
                p.close()
                
                results_dict.update(dict(zip(df_dict.keys(), results)))
    
    df_a = ResultsDictToPandas(results_dict)
    df_a.to_csv('df_a_full.csv')
    
    #Load all curves
    #Find corresponding pairs
    #Normalise
    #Calculate alive/dead time
    
    #An idea may be to simulate all the systems for 0 inclination and determine
    #their maximum values, Lx, and N_lim in advance as it would also save having to simulate
    #the 0 inclination systems more than once.
    
    #This would require simulating every systems for all dincls used 0 inclination.
    #Would take approximately 992 * 45 = 44640 simulations and they would not have to be repeated.
    
    
    def PlotHistogramResults():
        plt.figure()
        # plt.hist2d(df_a_nonzero['dincl'].values, df_a_nonzero['ratio'].values,bins=80)
        corner.hist2d(df_a_nonzero['dincl'].values, df_a_nonzero['ratio'].values,bins=20)
        plt.xlabel('Precession angle')
        plt.ylabel('alive/dead ratio')
        plt.title('151 Beamed BH ULXs, 1000 iterations per ulx in the range 0 - 45 dincl')
        plt.colorbar()
    
    # =============================================================================
    # Calculate explicit alive/dead for each system
    # =============================================================================
    results_dict = {}
    for sim_num in range(len(systems_df)):
        Lx = GetLx(sim_num) #Get Lx for system
        for dincl in dincl_list:    #Run through all precession angles
            curves = find_curve_by_id_and_dincl(df_dict, sim_num, dincl) #Two curves
            N_lim = Normalise(curves) #Find normalization limit
            for key in curves:
                Alive, Dead = calc_alive_time(df_dict[key], N_lim)
                results_dict[key] = Alive, Dead
                
        print(sim_num, '/', len(systems_df))
    
    
    PlotCurve('0-10.0-0')
    
    '''
    Things you can plot:
        FOR DIFFERENT Z AND TAGE:
            Alive Time vs dincl
            Dead Time vs dincl
            Ratio vs dincl
            beaming vs alive time
    '''
    
    def plot_curve():
        import matplotlib
        fontsize = 10
        matplotlib.rcParams['mathtext.fontset'] = 'stix'
        matplotlib.rcParams['font.family'] = 'STIXGeneral'
        plt.figure(figsize=(5.5,1.7))
        plt.xlabel('Time', fontsize=fontsize)
        plt.ylabel('Flux', fontsize=fontsize)
        plt.plot(curve['Time'], curve['Flux'], c='black', linewidth=1.0)
        plt.savefig('lightcurve.eps', format='eps', bbox_inches = "tight")
        
        



if __name__ == "__main__":
    main()