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

systems_df = load_systems_dataframe(True, True, True)
Lx_series = systems_df['Lx']

limit_dict = load_0_inclination_N_lim_dict()

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



#curve_files = glob.glob('./new_curves/*.txt')
#for curve_file in curve_files:
    

with open('./new_curves/a7232574-74cb-4240-a5b1-2e0fbf5262ca.txt', 'r') as file:
    data = file.read()
    simulation_info = data.splitlines()[-9:]
    system_id = int(simulation_info[0].split(':')[1])
    
    curve = curvefunc.load_curve_file_skip_footer(StringIO(data))
    Lx = Lx_series[882]
    key = str(system_id) + '-' + str(dincl)
    limit_dict[key]
    alive, dead = calc_alive_time(curve, limit)
    