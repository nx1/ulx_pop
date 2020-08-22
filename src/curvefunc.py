#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jan 22 17:23:33 2020

@author: nk7g14

This script contains functions related to the curves obtained from ulxlc

Lightcurves are stored in .txt files

The filename format goes as:
    simulation_number - dincl - inclination.txt
"""
import pandas as pd
import numpy as np
from pathlib import Path

import auxil


def load_curve_file(path):
    '''
    path: pathlib Path
    '''
    curve = pd.read_csv(path, delimiter=' ', header=None,
                        names=['Time', 'Time_Err', 'Flux'], skiprows=3)
    return curve

def load_curve_file_skip_footer(path):
    '''
    path: pathlib Path
    '''
    curve = pd.read_csv(path, delimiter=' ', header=None,
                        names=['Time', 'Time_Err', 'Flux'], skiprows=3,
                        skipfooter=9)
    return curve



def calc_alive_time(curve, limit):
    '''
    Calculates for a given curve the amount of time above and below a given
    limit
    '''
    if max(curve['Flux']) < limit:
        alive = 0
        dead = curve['Time'].iloc[-1]
        return alive, dead
    if min(curve['Flux']) > limit:
        alive = curve['Time'].iloc[-1]
        dead = 0
        return alive, dead

    
    curve['Above'] = curve['Flux'].map(lambda x: 1 if x >= limit else -1)

    df2 = curve[curve['Above'].diff() != 0]
    df2['time_diff'] = df2['Time'].diff()

    if len(df2) == 1:   #If the Curve never goes above the limit
        alive = 0.0
        dead = curve['Time'].iloc[-1]
    elif df2['Above'][0] == -1: #If the curve starts below the limit
        df_above = df2[df2['Above'] == 1]
        df_below = df2[df2['Above'] == -1]
        dead = np.sum(df_above['time_diff']) + df_above['Time'].iloc[0]
        alive = np.sum(df_below['time_diff'])
    elif df2['Above'][0] == 1: #If the curve starts above the limit
        df_above = df2[df2['Above'] == 1]
        df_below = df2[df2['Above'] == -1]
        dead = np.sum(df_above['time_diff'])
        alive = np.sum(df_below['time_diff']) + df_below['Time'].iloc[0]
    return alive, dead


def multiply_curve(curve, period, time_range):
    """
    Takes the first period of the curve and multiplies
    it n times to create a new curve
    period: period of the curve
    time_range: how long we want our time series to be
    """
    period_end_time = auxil.find_nearest(curve['Time'], period)
    period_end_index = curve.loc[curve['Time'] == period_end_time].index[0]
    
    new_curve = curve.copy()
    result = curve.copy()
    
    number_of_cycles = int(time_range / period) + 1
    
    for i in range(1, number_of_cycles):
        new_curve['Time'] = curve['Time'] + (i * period_end_time)
        result = pd.concat([result, new_curve])    
    return result
    
    
def scale_light_curve_period(curve, original_period, new_period):
    """Scale a lightcurve to a new period."""
    working = curve.copy()
    working['Time'] = (working['Time'] / original_period) * new_period
    return working


def get_lightcurve_ULX_flux_limit(curve_0, Lx):
    '''
    obtain ULX normalisation limit for two curves.
    curve_0: curve at 0 inclination
    Lx: maximum value of 0 inclination curve.
    '''
    zero_incl_max_flux = max(curve_0['Flux'])
    c = Lx / zero_incl_max_flux     #Curve Normalisation constant
    N_lim = 1E39 / c
    return N_lim
