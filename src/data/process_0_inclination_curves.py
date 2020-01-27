#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jan 27 15:34:45 2020

@author: nk7g14
"""

import sys
sys.path.append("..")

from curvefunc import load_151_systems_zero_inclination_curves
from auxil import load_systems_dataframe

import os
os.chdir("..")
zero_inclination_curves = load_151_systems_zero_inclination_curves()

systems_df = load_systems_dataframe(True, True, True)
Lx_arr = systems_df['Lx']


def calc_norm_limit(curve, Lx):
    """Calculate the normalisation limit for a paticular curve
    """
    max_flux = max(curve['Flux']) #Maximum Flux
    c = Lx / max_flux                   #Scaling factor
    N_lim = 1E39 / c                    #Limiting value
    return N_lim


normalisation_limits = {}
for key, curve in zero_inclination_curves.items():
    # print(key)
    system_id = int(key.split('-')[0])
    curve = zero_inclination_curves[key]
    Lx = Lx_arr[system_id]
    N_lim = calc_norm_limit(curve, Lx)
    normalisation_limits[key] = N_lim
    
f = open("../data/interim/0_inclination_N_lim.txt","w")
f.write(str(normalisation_limits))
f.close()