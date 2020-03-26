#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Feb 17 13:25:55 2020

@author: nk7g14

This file is used to simulate sampling ulx lightcurves obtained via ULXLC
as would be done by eRASS in intervals of six months
"""
import numpy as np
import pandas as pd
from pathlib import Path

import os
from tqdm import tqdm

from uuid import uuid4
import auxil
import curvefunc
import batchfunc

import matplotlib.pyplot as plt

def get_parameters(row):
    """Create parameters dictionary from df_a row"""
    parameters = {'period': ULXLC_PERIOD,
                    'phase': row['phase'],
                    'theta': row['theta'],
                    'inclination': row['inclination'],
                    'dincl': row['dincl'],
                    'beta': row['beta'],
                    'dopulse': 1,
                    'norm': row['norm']}
    return parameters


def get_lightcurves(row):
    """Obtain simulation lightcurve and 0 inclination normalisation curve.
    for a given row of simulation details obtained from df_a"""
    xcm_filename = 'temp.xcm'
    lc_filename = 'lc.txt'
    zero_incl_lc_filename = '0_lc.txt'
    
    parameters = get_parameters(row)
    
    batchfunc.run_ulxlc(xcm_filename, parameters, row['system_id'], lc_filename, append_to_file=False)
    parameters['inclination'] = 0
    batchfunc.run_ulxlc(xcm_filename, parameters, row['system_id'], zero_incl_lc_filename, append_to_file=False)
    
    curve = curvefunc.load_curve_file(lc_filename)
    curve_0 = curvefunc.load_curve_file(zero_incl_lc_filename)
    
    os.remove(lc_filename)
    os.remove(zero_incl_lc_filename)
    return curve, curve_0


def sample_curve(curve, curve_period, N_lim, sampling_interval, number_of_repeats):
    """Samples a lightcurve of a given period as would be done by eRASS
    (8 cycles of 6 months)
    parameters:
    -----------
    curve: lightcurve dataframe
    curve_period: period of the lightcurve
    sampling interval: how often to sample the lightcurve after selecting
    a random position somwhere in the first period.
    number_of_repeats: number of MC iterations
    """
    time_arr = np.array(curve['Time'])
    flux_arr = np.array(curve['Flux'])
    
    truth = []
    for n in tqdm(range(number_of_repeats)):
        start_time = np.random.uniform(0, 50)
        sample_times = np.array([(start_time + i*30*6)%curve_period for i in range(0,9)])
        sample_indexes = [np.argmin(np.abs(time_arr - t)) for t in sample_times]
        
        fluxes = flux_arr[sample_indexes]
        
        truth.append(list(np.where(fluxes > N_lim, True, False)))
        
    return truth
    
def truth_table_processor(truth):
    """Runs through the 8*N is ulx truth table
    to create a new dataframe with each column corresponding to an eRASS
    cycle and weather or not the source was observed as a transient ulx or not.
    thanks to based James for the help
    """
    potato = []
    for row in truth:
        transient=8*[True]
        i=0
        while row[i]==row[0] and i!=8:
            transient[i]=False
            i+=1
        potato.append(transient)
    
    df = pd.DataFrame(potato)
    return df.mean()
    
                
def plot_lightcurve(curve, N_lim, P_wind, sample_times, fluxes):
    plt.figure()
    plt.plot(curve['Time'], curve['Flux'], label='inclination:'+str(row['inclination']))
    #plt.plot(curve_0['Time'], curve_0['Flux'], label='inclination:0')
    plt.axhline(N_lim, c='r', label='ULX limit')
    plt.axvline(P_wind, c='blue', label=r'$P_{wind}$', linestyle='--', linewidth=0.8)
    plt.xlabel('Time (days)')
    plt.ylabel('Flux (normalised units)')
    plt.scatter(sample_times, fluxes, marker='x', c='r', s=12, label='sampling points')
    plt.legend()
    plt.title(r'$P_{wind}:$'+str(P_wind))
    plt.savefig(f'{SAVE_FOLDER}/lc_plots/{index}.eps')
    plt.savefig(f'{SAVE_FOLDER}/lc_plots/{index}.png', dpi=200)     
            

def make_simulation_settings_file():
    with  open(f'{SAVE_FOLDER}/settings.info', "a") as f:
        f.write(f'SIMULATION_UUID: {SIMULATION_UUID} \n')
        f.write(f'MAXIMUM_P_WIND: {MAXIMUM_P_WIND} \n')
        f.write(f'ULXLC_PERIOD: {ULXLC_PERIOD} \n')
        f.write(f'NUMBER_OF_SAMPLED_SYSTEMS: {NUMBER_OF_SAMPLED_SYSTEMS} \n')
        f.write(f'NUMBER_OF_SAMPLING_ITERATIONS: {NUMBER_OF_SAMPLING_ITERATIONS} \n')
        f.write(f'SAMPLING_INTERVAL: {SAMPLING_INTERVAL} \n')
        
        f.write(f'Z: {Z}')

if __name__ == '__main__':
    SIMULATION_UUID = str(uuid4())
    MAXIMUM_P_WIND = 365 * 4
    ULXLC_PERIOD = 50 #Period to output lightcurves from ulxlc
    NUMBER_OF_SAMPLED_SYSTEMS = 500 #How many samples we wish to create
    NUMBER_OF_SAMPLING_ITERATIONS = 100000 #How many MC sampling iterations
    SAMPLING_INTERVAL = 30*6 #days
    Z = 0.0002 #Metallicity to test
    
    print(f'RUNNING SIMULATION: {SIMULATION_UUID}')
    
    SAVE_FOLDER = Path(f'../data/interim/eROSITA_sampling_results/{SIMULATION_UUID}')


    df_systems = auxil.load_systems_dataframe(ulx_only=True,
                                          beamed=True, 
                                          half_opening_l_45=True)
    df_systems = df_systems.drop(['index'], axis=1)
    df_systems = df_systems[df_systems['P_wind_days'] < MAXIMUM_P_WIND] 
    df_systems = df_systems[df_systems['Z'] == Z]
    
    #Only select rows with system ids correspondng to our subset
    df_a = auxil.load_df_a(transient_only=False)
    df_a = df_a[df_a['system_id'].isin(df_systems.index)]
    
    # Create dataframe of systems to test
    sampled_systems = df_a.sample(n=NUMBER_OF_SAMPLED_SYSTEMS, replace=True)
    
    # Saving important information at start.
    os.makedirs(f'{SAVE_FOLDER}/lc_plots', exist_ok=True)
    df_systems.to_csv(f'{SAVE_FOLDER}/df_systems.csv')
    sampled_systems.to_csv(f'{SAVE_FOLDER}/sampled_systems.csv')
    make_simulation_settings_file()
    
    df_results = pd.DataFrame()
    
    for index, row in sampled_systems.iterrows():
        system_id = row['system_id']
        systems_row = df_systems.loc[system_id]
        P_wind = systems_row['P_wind_days']
        Lx = systems_row['Lx']
        
        curve, curve_0 = get_lightcurves(row)
        
        curve = curvefunc.scale_light_curve_period(curve, ULXLC_PERIOD, P_wind)  
        curve_0 = curvefunc.scale_light_curve_period(curve_0, ULXLC_PERIOD, P_wind)
        
        # curve = curvefunc.multiply_curve(curve, P_wind, 6*365)
        
        N_lim = curvefunc.get_lightcurve_ULX_flux_limit(curve_0, Lx)
        
        truth = sample_curve(curve,curve_period=P_wind,
                             N_lim=N_lim,
                             sampling_interval=SAMPLING_INTERVAL,
                             number_of_repeats=NUMBER_OF_SAMPLING_ITERATIONS)  
        
        results = truth_table_processor(truth)
        df_results[index] = results
        
#        plot_lightcurve(curve, N_lim, P_wind, sample_times, fluxes)
        df_results.to_csv(f'{SAVE_FOLDER}/results.csv')
