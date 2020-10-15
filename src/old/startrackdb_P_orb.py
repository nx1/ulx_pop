# -*- coding: utf-8 -*-
"""
Created on Fri Jun  5 12:12:18 2020

@author: norma

mt = 1 files have been extracted from https://universeathome.pl/universe/blackholes.php 
and saved in ../data/processed/startrackdb/*.csv

Kepler's third law states that the orbital period T of two point masses
orbiting each other in a circular or elliptical orbit is given by:
    
    T = 2 \pi \sqrt{\frac{a^3}{GM} }
    
    Where M is the mass of the more massive body
    
The columns originally provided are listed below:
t dt Ma Mb Ka Kb a e Ra Rb La Lb aspina aspinb mt mttype Lxmt Lx Vsmx Vsmy Vsmz Mzamsa Mzamsb a0 e0 idum iidd evroute

    t is the age of the system (i.e. time since ZAMS) [Myr]
    dt is the length of the actual time-step [Myr]
    a is the semi-major axis of the system [Rsun]
    e is the eccentricity of the system
    Rx is the radius of the star x [Rsun]
    Lx is the luminosity of the star x [Lsun]
    aspinx is a spin of the star x (non-zero only if x is a black hole)
    mt is 1/0 if the mass transfer is on/off
    mttype is the timescale of mass transfer (1 - nuclear, 4 - thermal, 5 - mass transfer from white dwarf)
    Lxmt is a luminosity of an accretion disk during RLOF [erg/s] (non-zero only if mass transfer onto a compact accretor is present)
    Lx is a luminosity coming from wind accretion [erg/s] (non-zero only if a donor has a significant wind and accretor is a compact object)
    Vsmx, Vsmy, Vsmz is a 3D centre of mass velocity [Rsun/day]
    Mzamsa/Mzamsb is the mass on ZAMS for star A/B [Msun]
    a0 is the separation on ZAMS [Rsun]
    e0 is the eccentricity on ZAMS
    idum, iidd - identifiers of a system (i.e. all lines with identical idum and iidd refer to the same binary)
    evroute is a symbolical description of binary evolution (see below)

Additional columns are given by:
    Z         : Metallicity
    P_orb     : Orbital period from Kepler's Third Law (using both masses)
    P_sup     : Superorbital period from Towndend & Charles (2020) P_sup = 22.1 +- 0.1 * P_orb
    P_sup_err : 1 sigma Error on P_sup
    eta       : Accretion efficiency
    mdot      : Mass accretion rate to the disk
    
    
Notes:
    - I am calculating mdot via Lxmt / L_edd
    - When selecting mt=1 there only seems to be BH systems
        email greg
    
"""
from pathlib import Path
import numpy as np
import pandas as pd

import matplotlib.pyplot as plt

import populations


def combine_csvs(csv_files):
    df = pd.DataFrame()
    for csv_file in csv_files:
        path = Path(csv_file)
        stem = path.stem
        Z = float(stem.split('_')[2])
        df_csv = pd.read_csv(csv_file)
        df_csv['Z'] = Z

    df = df.append(df_csv)
    return df



csv_files = ['../data/processed/startrackdb/bhdb180327_ZZ_0.02_Sal_-2.3_SS_1_BHSPIN_1_kick_6_data1_mt_1.csv',
             '../data/processed/startrackdb/bhdb180327_ZZ_0.002_Sal_-2.3_SS_1_BHSPIN_1_kick_6_data1_mt_1.csv',
             '../data/processed/startrackdb/bhdb180327_ZZ_0.0002_Sal_-2.3_SS_1_BHSPIN_1_kick_6_data1_mt_1.csv']


df = combine_csvs(csv_files)
#df = pd.read_csv(csv_files[0])
csv_rows = len(df)




# =============================================================================
# Additional Quantities
# =============================================================================
R_sol = 6.957e8
M_sol = 2e30
G = 6.67E-11
c = 3e8

df['P_orb'] = 2 * np.pi * np.sqrt((R_sol*df['a'])**3 / (G*M_sol*(df['Mb'] + df['Ma'])) )
df['P_sup'] = 22.1 * df['P_orb']
df['P_sup_err'] = df['P_sup'] * np.sqrt((0.1/22.1)**2)
df['eta'] = 0.1
df['mdot'] = df['Lxmt'] / (1.2e38 * df['Ma'])
# df['mdot'] = (df['Lxmt'] / (df['eta'] * c**2)) / (1.2e38 * df['Ma'])

# In days
df['P_orb_days'] = df['P_orb'] / (60*60*24)
df['P_sup_days'] = df['P_sup'] / (60*60*24)
df['P_sup_err_days'] = df['P_sup_err'] / (60*60*24)


# =============================================================================
# Exploratory
# =============================================================================
# Dataframe properties
df_unique_systems_count = df.groupby(['idum', 'iidd']).size().reset_index().rename(columns={0:'count'})
N_unique_systems = len(df_unique_systems_count['count'])
N_unique_luminosities = len(df['Lxmt'].unique())
Ka_unique = df['Ka'].value_counts()
Kb_unique = df['Kb'].value_counts()


# Mass distributions
df['Ma'].hist(bins=100)
df['Mb'].hist(bins=100)

# Spin distributions
plt.figure()
df['aspina'].hist(bins=50)

# Mass Transfer Luminosity distribution
df['Lxmt'].hist(bins=100, log=True)

# Mass vs semi-major axis
plt.figure()
plt.xlabel('Ma')
plt.ylabel('a')
plt.scatter(df['Ma'], df['a'], s=0.5)

#Log Mass vs semi-major axis
plt.figure()
plt.xlabel('log Ma')
plt.ylabel('log a')
plt.scatter(np.log(df['Ma']), np.log(df['a']), s=0.5)

# Mass Transfer Luminosity vs semi-major axis
plt.figure()
plt.xlabel('Lx')
plt.ylabel('a')
plt.scatter(df['Lxmt'], df['a'], s=0.5)

# Orbital period vs semi-major axis
plt.figure()
plt.xlabel('a')
plt.ylabel('oribtal period')
plt.scatter(df['a'], df['P_orb'], s=0.5)

# semi-major axis^3 vs oribtal period^2 (Kepler's third law)
plt.scatter(df['a']**3, df['P_orb']**2, s=0.5)
plt.xlabel('semi-major axis^3')
plt.ylabel('period^2')

# =============================================================================
# Further filtering
# =============================================================================
df = df[df['mttype'] == 1]  # 



# =============================================================================
# P_orb prediction
# =============================================================================
"""We will use a joint probability distribution on the variables
X = Ma
Y = mdot
P_orb

""" 


df_ulx = populations.ulx()

df = df[df['mdot'] < 1000]

x = np.array(df['Ma'])
y = np.array(df['mdot'])
z = np.array(df['a'])

# Mass
x_min = 0
x_max = 50
N_x_bins = 100
xbins = np.linspace(x_min, x_max, num=N_x_bins)

# mdot
y_min = 0
y_max = 1000
N_y_bins = 100
ybins = np.linspace(y_min, y_max, num=N_y_bins)

# Return the bin index in which each value of x
x_groups = np.digitize(x, xbins)
y_groups = np.digitize(y, ybins)


# groups now holds the index of the bin into which radius[i] falls
# loop through all bin indexes and select the corresponding masses
# perform your aggregation on the selected masses

m_a_agg = {}
for i in range(len(xbins)+1):
    selected_as = z[x_groups==i] #Selected semi-major axes
    m_a_agg[i] = selected_as


mdot_a_agg = {}
for i in range(len(ybins)+1):
    selected_as = z[y_groups==i]  #Selected semi-major axes
    mdot_a_agg[i] = selected_as


m = 25.9
mdot = 3.90

m_index = np.digitize(m, xbins)
mdot_index = np.digitize(mdot, xbins)

m_as = mdot_a_agg[m_index]
mdot_as = mdot_a_agg[mdot_index]

plt.hist(mdot_as)





plt.figure()
plt.xlabel('mass')
x_hist = plt.hist(x, bins=100, density=False)
plt.figure()
plt.xlabel('mdot')
y_hist = plt.hist(y, bins=100, density=False)


Ngrid = 100
grid = np.linspace(0, 2, Ngrid + 1)


H, xbins, ybins = np.histogram2d(x, z, grid)
plt.imshow(H)
#H /= np.sum(H)




df['P_orb']
df['P_orb'].describe()

df



