# -*- coding: utf-8 -*-
"""
Created on Tue May 19 12:16:45 2020

@author: norma

populations.py

This script contains the main immutable populations used throughout this work.
"""
from pathlib import Path

import pandas as pd


def startrack():
    systems_df_path = Path('../data/processed/all_systems_df.csv')
    df = pd.read_csv(systems_df_path)
    df = df.drop(['Unnamed: 0'], axis=1)
    return df

def ulx():
    df = startrack()
    df = df[df['Lx'] > 1E39]
    return df

def ulx_beamed():
    df = ulx()
    df = df[df['b'] < 1]
    return df

def ulx_beamed_l_45():
    df = ulx_beamed()
    df = df[df['theta_half_deg'] < 45]
    return df

def ulx_beamed_l_45_P_wind_l_4_years():
    df = ulx_beamed_l_45()
    df = df[df['P_wind_days'] < 4*365]
    return df

def all():
    populations = {'startrack' :                        startrack(),
                   'ulx' :                              ulx(),
                   'ulx_beamed_l_45' :                  ulx_beamed_l_45(),
                   'ulx_beamed_l_45_P_wind_l_4_years' : ulx_beamed_l_45_P_wind_l_4_years()}
    return populations

if __name__ == "__main__":
    df_full = startrack()
    df_ulx = ulx()
    df_ulx_beamed = ulx_beamed()
    df_beamed_l_45 = ulx_beamed_l_45()
    df_beamed_l_45_P_wind_l_4_years = ulx_beamed_l_45_P_wind_l_4_years()
    dict_all = all()
