#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jul 17 14:20:34 2019

@author: nk7g14
"""
import glob
import pandas as pd
from tqdm import tqdm
import matplotlib.pyplot as plt


def LoadCurves(folder):
    '''
    Loads all simulation curves outputted from xspec model ulxlc
    '''
    #Importing Files
    curve_files = glob.glob('./{}/*.txt'.format(folder))
    df_dict = {}      #Dictionary for storing dataframes for lightcurves
    pbar = tqdm(curve_files)
    for filename in pbar:
        pbar.set_description('Loading %s' % filename)
        split1 = filename.split('/')
        split2 = split1[-1][:-4]
        df_dict[split2] = pd.read_csv(filename, delimiter=' ',
               header=None, names=['Time', 'Time_Err', 'Flux'], skiprows=3)
        
    return df_dict


def LoadCurves2(folder):
    '''
    Loads all simulation curves outputted from xspec model ulxlc
    '''
    #Importing Files
    curve_files = glob.glob('./{}/*.txt'.format(folder))
    df_dict = {}      #Dictionary for storing dataframes for lightcurves
    
    pbar = tqdm(curve_files)
    for filename in pbar:
        pbar.set_description('Loading %s' % filename)
        
        params = filename.split('/')[-1][:-4] #only return the numbers from filename
        params_split = params.split('-')
        
        sim_num = params_split[0]
        opening_angle = params_split[1]
        inclination = params_split[2]
        lcurve = pd.read_csv(filename, delimiter=' ',
               header=None, names=['Time', 'Time_Err', 'Flux'], skiprows=3)
        
        
        row = [sim_num, opening_angle, inclination, lcurve]
        df_dict[params] = row
        
    return df_dict

if __name__ == '__main__':
    # df_all_systems = pd.read_csv('dataframe.csv')
    # df_all_systems = df_all_systems.drop('Unnamed: 0', axis=1)
    
    # df_ulx = df_all_systems[df_all_systems['Lx'] > 1E39]
    # df_beamed_ulx = df_ulx[df_ulx['b'] < 1]
    
    # # curves = LoadCurves('curves0')
    # curves2 = LoadCurves2('curves0')
    
    
    
    
    

