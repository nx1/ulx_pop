#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jan  9 11:21:17 2020

@author: nk7g14

This plot is used to illustrate the normalisation process

the two curves were created with parameters:

Curve 1:    
1    1   ulxlc      period     days     10.0000      +/-  0.0          
2    1   ulxlc      phase               0.0          +/-  0.0          
3    1   ulxlc      theta      deg      25.0000      +/-  0.0          
4    1   ulxlc      incl       deg      0.0          +/-  0.0          
5    1   ulxlc      dincl      deg      15.0000      +/-  0.0          
6    1   ulxlc      beta       c        0.200000     +/-  0.0          
7    1   ulxlc      dopulse             0            frozen
8    1   ulxlc      norm                1.00000      +/-  0.0           

Curve 2:
    
"""
import glob
import pandas as pd
import matplotlib.pyplot as plt

def LoadSystems():
    df_master = pd.read_csv('dataframe.csv')
    df = df_master[df_master['Lx'] < 1E39]
    df = df_master[df_master['b'] < 1]
    df = df.reset_index()
    df.columns
    df = df.drop(columns=['index', 'Unnamed: 0'])
    return df


def LoadCurves():
    #Importing Files
    curve_files = glob.glob('*.txt')
    print('Loading curve files...')
    df_dict={}      #Dictionary for storing dataframes for lightcurves
    for filename in curve_files:
        df_dict[filename[:-4]] = pd.read_csv(filename, delimiter=' ',
               header=None, names=['Time', 'Time_Err', 'Flux'], skiprows=3)
    print('Loading curve files...DONE!')
    return df_dict


curves = LoadCurves()

curve1 = curves['0_incl']
curve2 = curves['20_incl']

plt.plot(curve1['Time'], curve1['Flux'])
plt.plot(curve2['Time'], curve2['Flux'])