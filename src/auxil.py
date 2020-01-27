#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jan 15 10:57:04 2020

@author: nk7g14
"""
import pandas as pd
from pathlib import Path
import ast

def load_systems_dataframe(ulx_only=False, beamed=False, half_opening_l_45=False):
    systems_df_path = Path('../data/processed/all_systems_df.csv')
    df = pd.read_csv(systems_df_path)
    if ulx_only:
        df = df[df['Lx'] > 1E39]
        df = df.reset_index()
    if beamed:
        df = df[df['b'] < 1]
    if half_opening_l_45:
        df = df[df['theta_half_deg'] < 45]
        
    df = df.drop(['index', 'Unnamed: 0'], axis=1)
    return df

def load_0_inclination_N_lim_dict():
    with open("../data/interim/0_inclination_N_lim.txt","r") as f:
        N_lim_dict = ast.literal_eval(f.read())
    return N_lim_dict