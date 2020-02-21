#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jan 27 13:34:01 2020

@author: nk7g14
This file is used to run ulxlc for all 151 beamed and <45 opening angle
ulxlc at 0 inclination for normalisation purposes.
"""


import numpy as np
from uuid import uuid4

import batchfunc
from auxil import load_systems_dataframe

systems_df = load_systems_dataframe(True, True, True)

systems_df['theta_half_deg']

dincls = np.arange(0,46)

ulxlc_parameters = {'period': 10.0,
                'phase': 0.0,
                'theta': 10,
                'inclination': 0,
                'dincl': 12,
                'beta': 0.2,
                'dopulse': 0,
                'norm': 1.0}

for system_id in systems_df.index:
    selected_theta = systems_df['theta_half_deg'][system_id]
    ulxlc_parameters['theta'] = selected_theta   
    for d in dincls:
        ulxlc_parameters['dincl'] = d
        xcm_n = str(uuid4())+'.xcm'
        lc_filename = str(system_id) + '-' + str(d)
        lc_n = '151_systems_0_inclination_curves/'+lc_filename+'.txt'
        batchfunc.run_ulxlc(xcm_n, ulxlc_parameters, system_id, lc_n)
        
        

            

            
            
            