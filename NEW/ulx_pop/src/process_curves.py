#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Aug 12 13:57:56 2019

@author: nk7g14
"""

import pandas as pd

def LoadSystems():
    '''
    Loads all systems from STARTRACK from csv files
    '''
    df_master = pd.read_csv('dataframe.csv') #36420 SYSTEMS
    df = df_master[df_master['Lx'] > 1E39]  #Only ULX -     992 ULXs
    # df = df[df['b'] < 1]                    #Only Beamed -  227 Beamed ULXs
    # df = df[df['theta_half_deg'] < 45]      #thetha < 45 -  151 Beamed ULXs with half opening angles < 45
    df = df.reset_index()
    # df = df.drop(columns=['index', 'Unnamed: 0'])
    return df

def LoadCurve(filename):
    curve = pd.read_csv(filename, delimiter=' ', header=None,
                names=['Time', 'Time_Err', 'Flux'], skiprows=3)
    return curve

def LoadCurves(folder):
    '''
    Loads all simulation curves outputted from xspec model ulxlc
    '''
    #Importing Files
    curve_files = glob.glob('./{}/*.txt'.format(folder))
    filename = [i.split('/')[-1][:-4] for i in curve_files]
    df_dict = {}      #Dictionary for storing dataframes for lightcurves
    pbar = tqdm(curve_files)
    for filename in pbar:
        pbar.set_description('Loading curve: %s' % filename)
        filename_split = filename.split('/')[-1][:-4]
        df_dict[filename_split] = LoadCurve(filename)
    return df_dict


results_dict = {}
df = LoadSystems()
Lx_arr = df['Lx']