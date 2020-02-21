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


def load_all_curves_from_path(path, skip_footer=False):
    '''
    path: pathlib Path
    '''
    df_dict = {}
    for p in path.glob('*.txt'):
        filename = p.stem
        if skip_footer:
            df_dict[filename] = load_curve_file_skip_footer(p)
        else:
            df_dict[filename] = load_curve_file(p)
    return df_dict


def split_curve_filename(filename):
    """splits a string of form system_id-dincl-inclination"""
    system_id = filename.split('-')[0]
    dincl = filename.split('-')[1]
    inclination = filename.split('-')[2]
    return system_id, dincl, inclination


def find_curve_by_id_and_dincl(df_dict, system_id, dincl):
    '''
    Returns the lightcurves corresponding to a
    given simulation number and dincl
    '''
    CurveDict = {}
    for i in df_dict.keys():
        curve_system_id, curve_dincl, inclination = split_curve_filename(i)
        if str(curve_system_id) == str(system_id) and str(curve_dincl) == str(dincl):
            CurveDict[i] = df_dict[i]
    return CurveDict


def get_simulation_info(systems_df, filename):
    '''
    Returns the row of simulation info for a given simulation key
    '''
    system_id, dincl, inclination = split_curve_filename(filename)
    row = systems_df.loc[system_id]
    row['dincl'] = dincl
    row['inclination'] = inclination
    return row


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
    curve['Time'] = (curve['Time'] / original_period) * new_period
    return curve



def get_lightcurve_ULX_flux_limit(curve_0, Lx):
    '''
    obtain ULX normalisation limit for two curves.
    curve_0: curve at 0 inclination
    Lx: maximum value of 0 inclination curve.
    '''
    N_lim = None
    zero_incl_max_flux = max(curve_0['Flux'])
    c = Lx / zero_incl_max_flux     #Curve Normalisation constant
    N_lim = 1E39 / c
    return N_lim

# =============================================================================
# Datasets
# =============================================================================

def load_151_systems_zero_inclination_curves():
    curves_dir = Path("../data/interim/curves/151_systems_0_inclination_curves")
    curves_path = Path(curves_dir)
    df_dict = load_all_curves_from_path(curves_path, skip_footer=True)
    return df_dict