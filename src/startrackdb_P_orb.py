# -*- coding: utf-8 -*-
"""
Created on Fri Jun  5 12:12:18 2020

@author: norma

Kepler's third law states that the orbital period T of two point masses
orbiting each other in a circular or elliptical orbit is given by:
    
    T = 2 \pi \sqrt{\frac{a^3}{GM} }
    
    Where M is the mass of the more massive body
"""
import numpy as np
import pandas as pd

import matplotlib.pyplot as plt

csv_files = ['../data/processed/startrackdb/bhdb180327_ZZ_0.02_Sal_-2.3_SS_1_BHSPIN_1_kick_6_data1_mt_1.csv',
              '../data/processed/startrackdb/bhdb180327_ZZ_0.002_Sal_-2.3_SS_1_BHSPIN_1_kick_6_data1_mt_1.csv',
              '../data/processed/startrackdb/bhdb180327_ZZ_0.0002_Sal_-2.3_SS_1_BHSPIN_1_kick_6_data1_mt_1.csv']


df = pd.read_csv(csv_files[0])

df['idum'].unique()
df['iidd'].value_counts()
df['Mb'].value_counts() #Mb is the mass of the CO
df['e']
df['Mb'].hist(bins=100)



R_sol = 6.957e8
M_sol = 2e30
G = 6.67E-11

# Kepler's third law
df['P'] = 2 * np.pi * np.sqrt((R_sol*df['a'])**3 / (G*M_sol*(df['Mb'] + df['Ma'])) )
df['P_Ma'] = 2 * np.pi * np.sqrt((R_sol*df['a'])**3 / (G*M_sol*df['Ma']) )
df['P_Mb'] = 2 * np.pi * np.sqrt((R_sol*df['a'])**3 / (G*M_sol*df['Mb']) )

df['P_days'] = df['P'] / (60*60*24)
df['P_Ma_days'] = df['P_Ma'] / (60*60*24)
df['P_Mb_days'] = df['P_Mb'] / (60*60*24)

#Townsend Charles
df['P_sup'] = 22.9 * df['P']
df['P_sup_days'] = df['P_sup'] / (60*60*24)



def func(x, a, b, c):
    return a * np.exp(-b * x) + c




df

df = df[df['mt'] == 1]
df = df[df['mttype'] == 1]

plt.figure()
plt.xlabel('Ma')
plt.ylabel('a')
plt.scatter(df['Ma'], df['a'], s=0.5)

plt.figure()
plt.xlabel('Ma')
plt.ylabel('a')
plt.scatter(df['Ma'], df['a'], s=0.5)

plt.figure()
plt.xlabel('Lx')
plt.ylabel('a')
plt.scatter(df['Lxmt'], df['a'], s=0.5)

plt.figure()
plt.xlabel('log Ma')
plt.ylabel('log a')
plt.scatter(np.log(df['Ma']), np.log(df['a']), s=0.5)


plt.figure()
plt.xlabel('a')
plt.ylabel('period')
plt.scatter(df['a'], df['P'], s=0.5)

plt.scatter(df['a']**3, df['P']**2, s=0.5)
plt.xlabel('a^3')
plt.ylabel('period^2')
