#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jan 24 16:57:46 2020

@author: nk7g14

ulxlc_runner.py
"""

import numpy as np
from pathlib import Path
import itertools
from multiprocessing import Pool
import os

import batchfunc
from auxil import load_systems_dataframe


systems_df = load_systems_dataframe(True, True, True)
save_directory = Path('../data/interim/curves/MC_curves_eta_0.08_ns/gridsearch')
os.chdir(save_directory)

ulxlc_parameters = {'period': 10.0,
                'phase': 0.0,
                'theta': 10,
                'inclination': 0,
                'dincl': 12,
                'beta': 0.2,
                'dopulse': 0,
                'norm': 1.0}


indexs = systems_df.index
inclinations = np.arange(1,91)
dincls = np.arange(0,46)

def wrapper(id_dincl_tuple):
    system_id, dincl, i = id_dincl_tuple
    ulxlc_parameters['dincl'] = dincl
    filename =  f'{system_id}-{dincl}-{i}'
    xcm_n = f'{filename}.xcm'
    lc_n = f'{filename}.txt'
    batchfunc.run_ulxlc(xcm_n, ulxlc_parameters, lc_n)

with Pool(5) as p:
        p.map(wrapper, itertools.product(indexs, dincls, inclinations))
