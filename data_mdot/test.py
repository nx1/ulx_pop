#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Feb  4 09:41:07 2019

@author: nk7g14
This program takes the raw data from the STARTRACK population synthesis code
and creates a Pandas dataframe containing many additional columns for later
use.
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import glob


df_dict = {}  #Dictionary for storing dataframes for m_dot


#CGS UNITS
eta = 0.1   #Efficiency factor for Eddington Mass
c = 3E10     #Speed of Light in cm/s
Msol = 1.989E33 #Mass of sun in g
Myr = 31557600 * 1E6 #1 Myr in Seconds
yr = 31557600 #1 yr in Seconds


a_bh = 0.998 #Black hole max spin
a_ns = 0.001 #NS average spin

G_SI = 6.67408E-11 #N.(m^2)/(kg)^2 
G = 6.674E-8 #(cm)^3 g^-1 s^-2
pi = np.pi


'''
FILE NAMES:
    eg: 'Z_0.002_tage_20.dat'
    'Z_0.002' = Metallicity = 0.002
    Metallicity range: 0.0002 - 0.02 
    'tage_20' = Age of system (time since ZAMS in Myr).
    Age range: 10 - 10000 


COLUMN NAMES:

    mdot:
        mass transfer rate (in solar masses per Myr; to the disk, not to the accretor, 
        which is lower due to the mass loss in disk wind),
    
    m:
        mass of a compact object (in solar masses; we assumed that these with a mass 
        lower than 2.5 Msun are NSs, whereas the heavier ones are BHs).
    
    idum, iidd:
        identificators of the system in the simulations. They have no physical meaning,
        but will allow to find more info about some potentially odd, or interesting system.
    
    Z:
        Metallicity | range: 0.0002 - 0.02 
    
    tage:
        Age of system (time since ZAMS in Myr). | range: 10 - 10000   
    
    Mass Accretion rates:
        mdot_gs:
            Mass Accretion rate in grams per second
        mdot_ratio:
            Eddington Ratio
        MEdd:
            Mass Accretion rate for eddington.
            defined by the equation LEdd / n c^2 where n = 1/12 for BH and 0.2 for NS
            (See Grzegorz Wiktorowicz 2017 Appendix)      
    
    Luminosities:
        LEdd:
            Eddington Luminosity in erg/s from (1.2E38 * m)
        XLsph:
            A.R King method with log contribution 
        XLsph2:
            A.R King Method without log contribution
        LXtot:
            Grzegorz Wiktorowicz Method
        Lx:
            Beamed Luminosity from Lxtot/b
            In the 2018 paper by Greg they propose the following:
                Lx = Ledd * (1 + ln(mdot))  for mdot > 1
                Lx = Ledd * mdot            for mdot < 1
        ratio:
            Ratio of LXtot/Ledd
        ratio_beamed:
            Ratio of Lx/Ledd
            
    Beaming:
    b:
        Beaming parameter set by the following:
            b = 1 for m_dot < 8.5
            b = 73/mdot**2 >= 8.5
            b = 3.2E-3 for >150
        
        Examples of b to theta (Full opening angle):
            b = 0.0032  | theta = 9.2 deg
            b = 0.01    | theta = 16.8 deg
            b = 0.022   | theta = 24 deg
            b = 0.077   | theta = 45 deg
            b = 0.1     | theta = 53 deg
            b = 0.2     | theta = 73.9 deg
            b = 0.3     | theta = 91 deg
            b = 1.0     | theta = 180 deg (No beaming)
            
    theta:
        cone opening angle in rad
        Calculated from b = 1 - cos(theta/2)
    theta_deg:
        cone opening angle in rad
    theta_half_deg:
        cone half opening angle in rad


ULXLC:
    norm - Normalisation
    period - Period of the light curve
    phase - Phase since the start of the lightcurve
    theta - Half opening angle
    incl - Inclination of system (90 is edge on)
    dincl - Presccesion angle
    beta - velocity of the outflow; given in units of the speed of light
    dopulse - 0 or 1 for full lightcurve or 1 pulse


PROCESS:
    1) Create light curve for given system for variety of dincl
        norm = 1 
        period = 10
        phase = 0
        theta = from df
        incl = random(0,90)
        dincl = [0, 5 , 10, 30, 45]
        beta = 0.5
        dopulse = 0 
       
    2) Save light curve
    3) Filename  = mdot,m,Z,tage,dincl
            
    how to export the model from XSPEC:
        after plotting the model
        >ipl
        >wdata [filename]
        
        
    1) use python to run xspec script
    2) once you get to PLT> bit use python to run the command wdata
    3) use python to repeat this for different parameters
    4) Once you have all files perform the analysis in python


'''

def ImportFiles():
    files = glob.glob('*.dat')

    for filename in files:
        # print(filename)
        df_dict[filename] = pd.read_csv(filename, delim_whitespace=True,
               header=None, names=['mdot', 'm', 'idum', 'iidd'])
    
        cut_string = filename.split(sep='_', maxsplit=4)
        df_dict[filename]['Z'] = float(cut_string[1])
        df_dict[filename]['tage'] = float(cut_string[3][:-4])
        
    
    # Concat all the datafames into one.
    df_master = pd.concat(df_dict, ignore_index=True)
    df_master = df_master.drop(['idum'], axis=1) #Drop useless columns
    df_master = df_master.drop(['iidd'], axis=1) #Drop useless columns
    df_master['is_bh'] = np.where(df_master['m'] < 2.5, 0 , 1) #Add type column
    return df_master



df_master = ImportFiles()

m = df_master['m']
mdot = df_master['mdot']


# Mass accretion rate in g/s
df_master['mdot_gs'] = mdot * (Msol/Myr)   
mdot_gs = df_master['mdot_gs']


# Eddington Luminosity in ERG/S
df_master['LEdd'] = 1.2E38 * m 
LEdd = df_master['LEdd']

# Eddington Mass, Obtained from the following prescription:
# Ledd = eta * M_dot_Eddington * c^2
# where eta = 0.2 for NS and 1/12 for bh
df_master['MEdd'] = np.where(m <= 2.5, LEdd / (0.2 * c**2), LEdd / (1/12 * c**2))
MEdd = df_master['MEdd']    


# Eddington Ratio
df_master['mdot_ratio'] = mdot_gs / MEdd #AKA Eddington Ratio
mdot_ratio = df_master['mdot_ratio']


#Isotropic Luminosities (3 methods described above)
df_master['XLsph'] = abs(2.2E39 * (m/10) * (mdot_ratio/10)**2 * (1 + np.log(mdot_ratio)))
df_master['XLsph2'] = 2.2E39 * (m/10) * (mdot_ratio/10)**2
df_master['LXtot'] = np.where(mdot_ratio > 1, LEdd * (1 + np.log(mdot_ratio)), LEdd * mdot_ratio)
#df_master['LXtot'] = LEdd * (1 + np.log(mdot_ratio))


# Beaming factor
# b = 1 for m_dot < 8.5
# b = 73/mdot**2 >= 8.5
# b = 3.2E-3 for >150
df_master['b'] = np.where(mdot_ratio >= 8.5, 73/mdot_ratio**2, 1)
df_master['b'] = np.where(mdot_ratio >= 150, 3.2E-3, df_master['b'])


# Observed Luminosity
df_master['Lx'] = df_master['LXtot']/df_master['b']

#Eddington Ratios for beamed and unbeamed.
df_master['ratio'] = df_master['LXtot'] / df_master['LEdd']
df_master['ratio_beamed'] = df_master['Lx'] / df_master['LEdd']


# Opening angles
df_master['theta'] = 2 * np.arccos(1-df_master['b']) #full opening angle in rad
df_master['theta_deg'] = df_master['theta'] * 180/np.pi #deg
df_master['theta_half_deg'] = df_master['theta_deg'] / 2 #Half opening angle



# Zeta is the opening angle of the wind
# we will also set a floor of zeta = 2
df_master['zeta'] = np.tan((pi/2) - np.arccos(1 - (73/(df_master['mdot_ratio']**2))))
df_master['zeta'] = np.where(df_master['zeta'] <= 2, 2, df_master['zeta'])


# Constants used for calculating r_out
# epsilon_wind is given by L_wind / L_tot and is normally between 0.25 and 0.95
epsilon_wind = 0.25 #Normally between 0.25 to 0.95 (check latex)
beta = 1.4 #Velocity of the wind, distinct from the beta used in ulxlc

'''
For now we will use:
    r_in = 6 for a = 0.001 (NS)
    r_in = 1.25 for a =  0.998 (BH)
'''

# General Relativity stuff
df_master['r_schw'] = 2 * G * m * Msol / c**2   #Schwarzschild radius (cm)
df_master['r_isco_nospin'] = (6 * G * m * Msol) / c**2 #ISCO (nospin) (cm)
df_master['r_isco'] = np.where(m < 2.5, 6, 1.25)    #Units of R_g (i think)
df_master['r_sph'] = df_master['r_isco'] * df_master['mdot_ratio']
df_master['r_out'] = 3 * epsilon_wind / (beta * df_master['zeta']) * df_master['mdot_ratio']**3/2 * df_master['r_isco']
df_master['P_wind'] = np.where(m < 2.5,
         G * m * Msol * pi / (3 * c**3 * a_ns) * df_master['r_out']**3 * ((1 - (df_master['r_isco']/df_master['r_out']))**3)/(np.log(df_master['r_out']/df_master['r_isco'])),
         G * m * Msol * pi / (3 * c**3 * a_bh) * df_master['r_out']**3 * ((1 - (df_master['r_isco']/df_master['r_out']))**3)/(np.log(df_master['r_out']/df_master['r_isco'])))
df_master['f_wind'] = 1 / df_master['P_wind']
df_master = df_master.sort_values('Z')
df_master = df_master.sort_values('tage')
df_master = df_master.reset_index()
df_master = df_master.drop('index', axis=1)
df_master.to_csv('../dataframe.csv')

# =============================================================================
# Plotting Functions
# =============================================================================
def Plot_Mass_Lx(df):
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

def Plot_Evolution(df):
    plt.figure()
    pivot = pd.pivot_table(df, index = ['Z','tage'], aggfunc='count', columns='is_bh')
    
if __name__ == '__main__':
    # df_master = df_master[df_master['Lx'] > 1E39]
    # df_master = df_master[df_master['b'] < 1]
    # df_master.to_csv('../dataframe.csv')
    
    
    # Counting pivot tables
    pd.pivot_table(df_master, index = ['Z','tage'], aggfunc='count')
    pd.pivot_table(df_master, index = ['Z','tage'], aggfunc='count', columns='is_bh')
    Plot_Mass_Lx(df_master)
    
