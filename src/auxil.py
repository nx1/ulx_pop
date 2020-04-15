#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jan 15 10:57:04 2020

@author: nk7g14
"""
import pandas as pd
import numpy as np
import math
from pathlib import Path
import ast

def load_systems_dataframe(ulx_only=False, beamed=False, half_opening_l_45=False):
    systems_df_path = Path('../data/processed/all_systems_df.csv')
    df = pd.read_csv(systems_df_path)
    if ulx_only:
        df = df[df['Lx'] > 1E39]
    if beamed:
        df = df[df['b'] < 1]
    if half_opening_l_45:
        df = df[df['theta_half_deg'] < 45]
        
    df = df.drop(['Unnamed: 0'], axis=1)
    return df

def find_nearest(array,value):
    idx = np.searchsorted(array, value, side="left")
    if idx > 0 and (idx == len(array) or math.fabs(value - array[idx-1]) < math.fabs(value - array[idx])):
        return array[idx-1]
    else:
        return array[idx]

def load_0_inclination_N_lim_dict():
    with open("../data/interim/0_inclination_N_lim.txt","r") as f:
        N_lim_dict = ast.literal_eval(f.read())
    return N_lim_dict


def load_df_a(transient_only=False):
    df_a_path = Path('../data/processed/df_a.csv')
    df = pd.read_csv(df_a_path)
    if transient_only:
        df = df[(df['ratio'] < 1) & (df['ratio'] != 0)]
    df = df.drop(['Unnamed: 0'], axis=1)
    return df