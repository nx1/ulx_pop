#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Feb 17 10:49:27 2020

@author: nk7g14
"""
import sys
import os

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

sys.path.append("..") # Adds higher directory to python modules path.

from auxil import load_systems_dataframe


# =============================================================================
# Plotting Functions
# =============================================================================
def plot_mass_Lx(df):
    """
    Plots Mass vs Luminosity for all systems in a specified dataframe
    splits into subplots according to Z and tage and colors them based on if
    they are black holes or neutron stars.

    Parameters
    ----------
    df : Pandas Dataframe
        Pandas dataframe containing the simulation output
    """
    fig, axarr = plt.subplots(3, 10)
    for i, tage in enumerate(df['tage'].unique()):
        for j, Z in enumerate(df['Z'].unique()):
            # print(i, j, Z, tage)
            mask = (df['Z'] == Z) & (df['tage'] == tage)
            
            m = np.log10(df[mask]['m'])
            Lx = np.log10(df[mask]['Lx'])
            is_bh = df[mask]['is_bh']
            axarr[j, i].scatter(m, Lx , s=1.0, c=is_bh, cmap='coolwarm')
            axarr[0, i].set_title('t = %s' % tage)
            axarr[j, 0].set_ylabel('Z = %s' % Z)

def plot_number_evolution(df_master):
    f, axarr = plt.subplots(3,1)

    pivot = pd.pivot_table(df_master, index = ['Z','tage'], aggfunc='count', columns='is_bh')
    for i, z in enumerate(df_master['Z'].unique()):
        N_BH_arr = []
        N_NS_arr = []
        for tage in df_master['tage'].unique():
            sub = df_master[(df_master['Z'] == z) & (df_master['tage'] == tage)]
            try:
                N_BH = sub.is_bh.value_counts()[1]
            except KeyError:
                N_BH = 0
            try:
                N_NS = sub.is_bh.value_counts()[0]
            except KeyError:
                N_NS = 0
            N_BH_arr.append(N_BH)
            N_NS_arr.append(N_NS)
                    
        axarr[i].plot(df_master['tage'].unique(), N_BH_arr, label='BH', c='black')
        axarr[i].text(x=100, y=200 , s='Z = ' + str(z))
        axarr[i].plot(df_master['tage'].unique(), N_NS_arr, label='NS', c='blue')
        axarr[i].legend()
        axarr[i].set_xlabel('tage')
        axarr[i].set_ylabel('Number of systems')
            
    

def CountBHNS(df):
    N = len(df)
    N_bh = len(df[df['is_bh']==1])
    N_ns = len(df[df['is_bh']==0])
    return [N, N_bh, N_ns]


if __name__ == "__main__":
    os.chdir('..')
    df_master = load_systems_dataframe(ulx_only=False, beamed=False, half_opening_l_45=False)
    plot_mass_Lx(df_master)
    plot_number_evolution(df_master)
