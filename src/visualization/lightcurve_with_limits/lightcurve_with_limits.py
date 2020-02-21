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
   3    1   ulxlc      theta      deg      5.60000      +/-  0.0          
   4    1   ulxlc      incl       deg      36.0000      +/-  0.0          
   5    1   ulxlc      dincl      deg      42.0000      +/-  0.0          
   6    1   ulxlc      beta       c        0.200000     +/-  0.0          
   7    1   ulxlc      dopulse             0            frozen
   8    1   ulxlc      norm                1.00000      +/-  0.0         

Curve 2:
    
   1    1   ulxlc      period     days     10.0000      +/-  0.0          
   2    1   ulxlc      phase               0.0          +/-  0.0          
   3    1   ulxlc      theta      deg      5.60000      +/-  0.0          
   4    1   ulxlc      incl       deg      0.0      +/-  0.0          
   5    1   ulxlc      dincl      deg      42.0000      +/-  0.0          
   6    1   ulxlc      beta       c        0.200000     +/-  0.0          
   7    1   ulxlc      dopulse             0            frozen
   8    1   ulxlc      norm                1.00000      +/-  0.0       
    
"""
import glob
import pandas as pd
import matplotlib.pyplot as plt


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

curve1 = curves['ff']
curve2 = curves['ff_incl']
curve3 = curves['test']
#
curve2['Flux'] = curve2['Flux']/max(curve1['Flux']) * 1.900110e+41
curve1['Flux'] = curve1['Flux']/max(curve1['Flux']) * 1.900110e+41

import matplotlib
fontsize = 10
matplotlib.rcParams['mathtext.fontset'] = 'stix'
matplotlib.rcParams['font.family'] = 'STIXGeneral'


plt.figure(figsize=(8,3))
plt.tight_layout()
plt.gcf().subplots_adjust(bottom=0.15)

plt.plot(curve1['Time'], curve1['Flux'], label='$i=0^{\circ}$', c='grey', linestyle='--', linewidth=0.8)
plt.plot(curve2['Time'], curve2['Flux'], label='$i=36^{\circ}$', c='black', linestyle='-', linewidth=0.8)
plt.yscale('log')
plt.xlabel('Time')
plt.ylabel(r'Flux ($\mathrm{erg} \ \mathrm{s}^{-1}$)')

plt.axhline(y=max(curve1['Flux']), linestyle='--', linewidth=0.5, c='b')
plt.axhline(y=1e39, linestyle='--', linewidth=0.5, c='red')
plt.xlim(0, 60)

plt.text(x=48.3, y=0.45e39, s=r'$L_{ULX} = 1 \times 10^{39} \ \mathrm{erg} \ \mathrm{s}^{-1}$', c='r', fontsize=9)
plt.text(x=49, y=1e41, s=r'$L_x = 1.9 \times 10^{41} \ \mathrm{erg} \ \mathrm{s}^{-1}$', c='b', fontsize=9)
plt.legend()
plt.savefig('../../../reports/figures/lightcurve_w_limits.eps')
plt.savefig('../../../reports/figures/lightcurve_w_limits.png', dpi=200)

