# -*- coding: utf-8 -*-
"""
Created on Wed Jul  8 11:50:03 2020

@author: norma
"""
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

import curvefunc
import ulxlc


data_file = '../data/external/variable_ulxlc/tlag_f0p2_m10_mdot20_test.dat'

df = pd.read_csv(data_file, sep=' ', header=None, names=['time', 'c2', 'c3', 'c4', 'period'])

times = df['time']
periods = df['period']

t0 = times[0]
p0 = periods[0]


    
ulxlc_period = 50.0

ulxlc_parameters = {'period': ulxlc_period,
                'phase': 0.0,
                'theta': 10.0,
                'inclination': 15.0,
                'dincl': 20.0,
                'beta': 0.2,
                'dopulse': 1,
                'norm': 1.0}

lc_filename = 'lc.txt'

ulxlc.run_ulxlc('test.xcm', ulxlc_parameters, lc_filename)

curve = curvefunc.load_curve_file(lc_filename)


p_last = p0
t_start = t0

for i in range(len(periods)):
    p = periods[i]
    t = times[i]
    
    if p != p_last:
        print(f'index: {i}')
        print(f'p: {p} p_last: {p_last}')
        print(f't:{t} t_start: {t_start}')
        
        
        time_range = abs(t-t_start)
        #import pdb; pdb.set_trace()
        
        scaled_curve = curvefunc.scale_light_curve_period(curve, ulxlc_period, p)
        plt.axvline(t, c='r', linestyle='--', linewidth=0.5)
        scaled_curve = curvefunc.multiply_curve(scaled_curve, p, time_range )
        scaled_curve['Time'] = scaled_curve['Time'] + t_start
        plt.plot(scaled_curve['Time'], scaled_curve['Flux'])
        t_start = t
    p_last = p