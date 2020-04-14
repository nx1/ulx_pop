#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jan 24 16:57:46 2020

@author: nk7g14
"""

import numpy as np
from uuid import uuid4

import batchfunc
from auxil import load_systems_dataframe

systems_df = load_systems_dataframe(True, True, True)

systems_df['theta_half_deg']

inclinations = np.arange(0,91)
dincls = np.arange(0,46)

ulxlc_parameters = {'period': 10.0,
                'phase': 0.0,
                'theta': 10,
                'inclination': 0,
                'dincl': 12,
                'beta': 0.2,
                'dopulse': 0,
                'norm': 1.0}




while True: 
    sample = systems_df['theta_half_deg'].sample()
    selected_theta = sample.values[0]
    selected_id = sample.index[0]
    ulxlc_parameters['theta'] = selected_theta
    for i in inclinations:
        ulxlc_parameters['inclination'] = i
        for d in dincls:
            ulxlc_parameters['dincl'] = d

            filename =  f'{selected_id}-{d}-{i}'
            xcm_n = f'{filename}.xcm'
            lc_n = f'{filename}.txt'
            batchfunc.run_ulxlc(xcm_n, ulxlc_parameters, selected_id, lc_n)

