#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jan 15 10:57:04 2020

@author: nk7g14
"""
import pandas as pd

def load_systems_dataframe(ulx_only=False, beamed=False, half_opening_l_45=False):
    df = pd.read_csv('../data/processed/all_systems_df.csv')
    if ulx_only:
        df = df[df['Lx'] > 1E39]
    if beamed:
        df = df[df['b'] < 1]
    if half_opening_l_45:
        df = df[df['theta_half_deg'] < 45]
    df = df.reset_index()
    df = df.drop(columns=['index', 'Unnamed: 0'])
    return df