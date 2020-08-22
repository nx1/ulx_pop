#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jan 15 10:57:04 2020

@author: nk7g14
"""
import os
from pathlib import Path
from urllib.request import urlopen
from shutil import copyfileobj

import pandas as pd
import numpy as np
import math
import ast

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

def load_curve_classifications():
    curve_classifications = pd.read_csv('../data/processed/curve_classifications.csv')
    return curve_classifications

def sample_by_bh_ratio(systems_df, bh_ratio, n):
    ns_ratio = 1 - bh_ratio
    
    ns_systems = systems_df[systems_df['is_bh']==0].index
    bh_systems = systems_df[systems_df['is_bh']==1].index
    
    bh_weights = [bh_ratio/len(bh_systems)]*len(bh_systems)
    ns_weights = [ns_ratio/len(ns_systems)]*len(ns_systems)
    selected_systems = np.random.choice([*bh_systems, *ns_systems], size=n, p=[*bh_weights, *ns_weights])
    return selected_systems


def download_file(url, path):
    """
    Download file from a url

    Parameters
    ----------
    url : string
        url of item to download
    path : TYPE
        path to download file to

    """
    if os.path.isfile(path):
        print(f'{path} already exists, not downloading.')
    else:
        print(f'downloading {path}')
        with urlopen(url) as in_stream, open(path, 'wb') as out_file:
            copyfileobj(in_stream, out_file)