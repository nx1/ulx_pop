#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jan 27 11:24:02 2020

@author: nk7g14
"""

import glob
import curvefunc
import pandas as pd
from io import StringIO
import numpy as np
from auxil import load_systems_dataframe, load_0_inclination_N_lim_dict
import os
from uuid import uuid4
from tqdm import tqdm


def Normalise(curve, Lx):
    '''
    inputs :
        curves = dictionary of two curves
    Takes two curves and Normalises them based on inclination
    '''
    N_lim = None
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


def calc_alive_time(curve, limit):
    '''
    Calculates for a given curve the amount of time above and below a given
    limit
    '''
    if max(curve['Flux']) < limit:
        Alive = 0
        Dead = curve['Time'].iloc[-1]
        return Alive, Dead
    if min(curve['Flux']) > limit:
        Alive = curve['Time'].iloc[-1]
        Dead = 0
        return Alive, Dead

    
    curve['Above'] = curve['Flux'].map(lambda x: 1 if x >= limit else -1)

    df2 = curve[curve['Above'].diff() != 0]
    df2['time_diff'] = df2['Time'].diff()

    if len(df2) == 1:   #If the Curve never goes above the limit
        Alive = 0.0
        Dead = curve['Time'].iloc[-1]
    elif df2['Above'][0] == -1: #If the curve starts below the limit
        df_above = df2[df2['Above'] == 1]
        df_below = df2[df2['Above'] == -1]
        Dead = np.sum(df_above['time_diff']) + df_above['Time'].iloc[0]
        Alive = np.sum(df_below['time_diff'])
    elif df2['Above'][0] == 1: #If the curve starts above the limit
        df_above = df2[df2['Above'] == 1]
        df_below = df2[df2['Above'] == -1]
        Dead = np.sum(df_above['time_diff'])
        Alive = np.sum(df_below['time_diff']) + df_below['Time'].iloc[0]
    return Alive, Dead

def create_simulation_info_dict(simulation_info_list):
    info = dict([tuple(simulation_info_list[i].split(':')) for i in range(len(simulation_info_list))])
    return info

if __name__ == '__main__':
    systems_df = load_systems_dataframe(True, True, True)
    Lx_series = systems_df['Lx']
    
    limit_dict = load_0_inclination_N_lim_dict()
    curve_files = glob.glob('./new_curves/*.txt')
    done_files = list(pd.read_csv('done_files.csv')['file'])
    
    save_every = 50000
    
    results_list = []
    
    for run_number, curve_file in tqdm(enumerate(files_to_do)):
        
        with open(curve_file, 'r') as file:
            data = file.read()
            simulation_info_list = data.splitlines()[-9:]
            simulation_info = create_simulation_info_dict(simulation_info_list)
            curve = curvefunc.load_curve_file_skip_footer(StringIO(data))
        
        Lx = Lx_series[882]
        key = simulation_info['system_id'] + '-' + simulation_info['dincl'].strip()
        limit = limit_dict[key]
        alive, dead = calc_alive_time(curve, limit)
        simulation_info['alive'] = alive
        simulation_info['dead'] = dead
        simulation_info['ratio'] = np.divide(alive,dead)
        results_list.append(simulation_info)
        if run_number%save_every == 0:
            print('run_number:', run_number, 'saving...')
            df_results = pd.DataFrame(results_list)
            df_results.to_csv('new_curve_results/'+str(uuid4())+'.csv')
            results_list = []
        
    
    
        
        
