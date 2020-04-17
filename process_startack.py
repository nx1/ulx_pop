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

G_SI = 6.67408E-11 #N.(m^2)/(kg)^2 
G = 6.674E-8 #(cm)^3 g^-1 s^-2
pi = np.pi

#GR constants
epsilon_wind = 0.25 #Normally between 0.25 to 0.95 (check latex)
beta = 1.4 #Velocity of the wind, distinct from the beta used in ulxlc
NS_SPIN = 0.001
BH_SPIN = 0.998

'''
FILE NAMES:
===========
    eg: 'Z_0.002_tage_20.dat'
    'Z_0.002' = Metallicity = 0.002
    Metallicity range: 0.0002 - 0.02 
    'tage_20' = Age of system (time since ZAMS in Myr).
    Age range: 10 - 10000 


COLUMN NAMES:
=============
    ['mdot', 'm', 'Z', 'tage', 'is_bh', 'mdot_gs', 'LEdd', 'MEdd',
       'mdot_ratio', 'XLsph', 'XLsph2', 'LXtot', 'b', 'Lx', 'ratio',
       'ratio_beamed', 'theta', 'theta_deg', 'theta_half_deg', 'zeta',
       'r_schw', 'r_isco_nospin', 'r_isco', 'r_sph', 'r_out', 'P_wind',
       'f_wind'] 


COLUMN DESCRIPTONS:
===================
    mdot:
        mass transfer rate (in solar masses per Myr; to the disk, not to the accretor, 
        which is lower due to the mass loss in disk wind),
        units: M_sol / Myr
    
    m:
        mass of a compact object (in solar masses; we assumed that these with a mass 
        lower than 2.5 Msun are NSs, whereas the heavier ones are BHs).
        
        units: M_sol
    
    idum, iidd:
        identificators of the system in the simulations. They have no physical meaning,
        but will allow to find more info about some potentially odd, or interesting system.
    
    Z:
        Metallicity | range: 0.0002 - 0.02 
        units: Solar metallicity
        
    tage:
        Age of system time since Zero age main sequence (ZAMS)
        units: Myr
        
        
    Mass Accretion rates:
    =====================
        mdot_gs:
            see mdot
            units: grams / second
            
        mdot_ratio:
            Eddington Ratio
            units: None (Eddington)
            
        MEdd:
            Mass Accretion rate for eddington.
            defined by the equation LEdd / n c^2 where n = 1/12 for NS and BH
            * Greg used 0.2 for NS
            (See Grzegorz Wiktorowicz 2017 Appendix)
            units: None (Eddington)
            
    
    Luminosities:
    =============
        LEdd:
            Eddington Luminosity in erg/s from (1.2E38 * m)
            units: erg / s
        XLsph:
            A.R King method with log contribution 
            units: erg / s
        XLsph2:
            A.R King Method without log contribution
            units: erg / s
        LXtot:
            Grzegorz Wiktorowicz Method
            units: erg / s
        Lx:
            Beamed Luminosity from Lxtot/b
            In the 2018 paper by Greg they propose the following:
                Lx = Ledd * (1 + ln(mdot))  for mdot > 1
                Lx = Ledd * mdot            for mdot < 1
            units: erg / s
        ratio:
            Ratio of LXtot/Ledd
            units: None
        ratio_beamed:
            Ratio of Lx/Ledd
            units: None
            
            
    Cone Geometry:
    ==============
    b:
        Beaming parameter set by the following:
            b = 1 for m_dot < 8.5
            b = 73/mdot**2 >= 8.5
            b = 3.2E-3 for >150
        units: None
            
    theta:
        cone opening angle:
            b = 1 - cos(theta/2)
        units: radians
        
    theta_deg:
        cone opening angle 
        units: degrees
        
    theta_half_deg:
        cone half opening angle
        units: degrees
        


For now we will use:


    General Relativity:
    ===================
        zeta:
            tan of the opening angle of the wind, floor set at zeta=2
            zeta = tan[pi/2 - arccos(1 - 73 / mdot^2)]
            units: None
                
        r_g:
            Gravitational radius
            G*m / c^2
            units: cm
        
        a*:
            System spin, currently black holes are set to maximum spin
            while neutron stars are set to an average spin.
                a_bh = 0.998
                a_ns = 0.001
            units: None
            
        r_schw:
            Schwartzchild radius
                2*G*m / c^2
            units: R_g
                
        r_isco_nospin:
            Innermost stable circular orbit for non spinning object (r_in)
                6*G*m / c^2
                3*r_schw
            units: R_g
                
        r_isco:
            Innermost stable circular orbit (r_in)
                Set to be 6 for both NS and BH (no spin)
            units: R_g
                
        r_sph:
            spherization radius
                r_isco * mdot
            
            * alternatively from King 2009, we could use:
                R_sph = 27 * mdot_ratio * r_schw / 4
                
            units: R_g
                
        r_out:
            3 * epsilon_wind / (beta * zeta' * mdot_ratio^3/2 * r_isco
            units: R_g (from r_isco)
            where:
                beta = 1.4 | wind velocity
                units:(?)
                
                epsilon_wind = 0.25
                given by L_wind / L_tot and is normally between 0.25 and 0.95
                units: None

        P_wind_at_rsph:
            Precession period at r_sph
            (G * m * pi) / (3 * c^3 * a) * r_out^3 * [ (1 - (r_in/r_out)^3) / ( ln( r_out / r_in ) )]
            where all radii are expressed in units of r_g see Middleton 2017, 2019
            units: seconds
            
        P_envelope:
            precession period of the optically thick envolope
            (G * m * pi) / (3 * c^3 * a) * r_sph^3 * [ (1 - (r_in/r_out)^3) / ( ln( r_out / r_in ) )] * (r_out/r_sph)
            where all radii are expressed in units of r_g see Middleton 2017, 2019
            units: seconds
'''

df_master = pd.read_csv('../../data/processed/startrack_concat.csv')
df_master.index = df_master['Unnamed: 0']
df_master = df_master.drop(['Unnamed: 0'], axis=1)
df_master['is_bh'] = np.where(df_master['m'] < 2.5, 0 , 1) #Add type column

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
# where eta = 1/12 for bh and ns, since we are considering mass transfer to the disc and not the compact object
# This is distinctly different to what was done by Wiktorowicz who used 1/12 for BH and 0.2 for NS.
df_master['MEdd'] =LEdd / (1/12 * c**2)
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


df_master['zeta'] = np.tan((pi/2) - np.arccos(1 - (73/(df_master['mdot_ratio']**2))))
df_master['zeta'] = np.where(df_master['zeta'] <= 2, 2, df_master['zeta'])


# General Relativity stuff
df_master['R_g'] = (G * m*Msol) / c**2 #gravitational radii
df_master['a*'] = np.where(m<2.5, NS_SPIN, BH_SPIN)

df_master['r_schw'] = ((2 * G * m * Msol) / c**2) / df_master['R_g']
df_master['r_isco_nospin'] = ((6 * G * m * Msol) / c**2) / df_master['R_g']
df_master['r_isco'] = 6
# df_master['r_isco'] = np.where(m < 2.5, 6, 1.25)  #We were previously using r_isco = 1.25 for BH and 6 for NS
df_master['r_sph'] = df_master['r_isco'] * df_master['mdot_ratio']
df_master['r_out'] = 3 * epsilon_wind / (beta * df_master['zeta']) * df_master['mdot_ratio']**3/2 * df_master['r_isco']

df_master['P_inflow_at_rsph'] = (G * m * Msol * pi) / (3 * c**3 * df_master['a*']) * df_master['r_sph']**3 * ((1 - (df_master['r_isco']/df_master['r_sph']))**3)/(np.log(df_master['r_sph']/df_master['r_isco']))
df_master['P_envelope'] = (G * m * Msol * pi) / (3 * c**3 * df_master['a*']) * df_master['r_sph']**3 * ((1 - (df_master['r_isco']/df_master['r_sph']))**3)/(np.log(df_master['r_sph']/df_master['r_isco'])) * (df_master['r_out']/df_master['r_sph'])**2
df_master['P_wind'] = (G * m * Msol * pi) / (3 * c**3 * df_master['a*']) * df_master['r_out']**3 * ((1 - (df_master['r_isco']/df_master['r_out']))**3)/(np.log(df_master['r_out']/df_master['r_isco']))

df_master['P_inflow_days'] = df_master['P_inflow_at_rsph'] / (24*60*60)
df_master['P_envelope_days'] = df_master['P_envelope'] / (24*60*60)
df_master['P_wind_days'] = df_master['P_wind'] / (24*60*60)



if __name__ == '__main__':
    # df_master = df_master[df_master['Lx'] > 1E39]
    # df_master = df_master[df_master['b'] < 1]
    
    # Counting pivot tables
    pd.pivot_table(df_master, index = ['Z','tage'], aggfunc='count')
    pd.pivot_table(df_master, index = ['Z','tage'], aggfunc='count', columns='is_bh')
    
    
    # =========================================================================
    # UNCOMMENT THIS FOR CHANGES TO DATAFRAME
    # =========================================================================
    
    df_master.to_csv('../../data/processed/all_systems_df.csv')
