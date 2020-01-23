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

from pathlib import Path


def load_curve_file(path):
    '''
    path: pathlib Path
    '''
    curve = pd.read_csv(path, delimiter=' ', header=None,
                        names=['Time', 'Time_Err', 'Flux'], skiprows=3)
    return curve


def load_all_curves_from_path(path):
    '''
    path: pathlib Path
    '''
    df_dict = {}
    for p in path.glob('*.txt'):
        filename = p.stem
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

# =============================================================================
# Datasets
# =============================================================================

def load_zero_inclination_curves():
    curves_dir = "../data/interim/curves/227_systems_0_inclination_curves"
    curves_path = Path(curves_dir)
    df_dict = load_all_curves_from_path(curves_path)
    return df_dict