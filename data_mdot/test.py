#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Feb  4 09:41:07 2019

@author: nk7g14
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import glob


df_dict={}  #Dictionary for storing dataframes for m_dot


#CGS UNITS
eta = 0.1   #Efficiency factor for Eddington Mass
c = 3E10     #Speed of Light in cm/s
Msol = 1.989E33 #Mass of sun in g
Myr = 31557600 * 1E6 #1 Myr in Seconds
yr = 31557600 #1 yr in Seconds


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
files = glob.glob('*.dat')

for filename in files:
    print(filename)
    df_dict[filename]=pd.read_csv(filename, delim_whitespace=True, header=None,
           names=['mdot', 'm', 'idum', 'iidd'])

    cut_string = filename.split(sep='_', maxsplit=4)
    df_dict[filename]['Z'] = float(cut_string[1])
    df_dict[filename]['tage'] = float(cut_string[3][:-4])


df_master=pd.concat(df_dict, ignore_index=True)
m = df_master['m']
mdot = df_master['mdot']

df_master['mdot_gs'] = mdot * (Msol/Myr)   #Mass accretion rate in g/s
mdot_gs = df_master['mdot_gs']        

df_master['LEdd'] = 1.2E38 * m #ERG/S            
LEdd = df_master['LEdd']       
df_master['MEdd'] = np.where(m <= 2.5, LEdd / (0.2 * c**2), LEdd / (1/12 * c**2))   
#Ledd = eta * M_dot_eddington * c^2 where eta = 0.2 for NS and 1/12 for bh    
MEdd = df_master['MEdd']            
df_master['mdot_ratio'] = mdot_gs / MEdd #AKA Eddington Ratio
mdot_ratio = df_master['mdot_ratio'] 

df_master['XLsph'] = abs(2.2E39 * (m/10) * (mdot_ratio/10)**2 * (1 + np.log(mdot_ratio)))            
df_master['XLsph2'] = 2.2E39 * (m/10) * (mdot_ratio/10)**2
df_master['LXtot'] = np.where(mdot_ratio > 1, LEdd * (1 + np.log(mdot_ratio)), LEdd * mdot_ratio)
#df_master['LXtot'] = LEdd * (1 + np.log(mdot_ratio))


df_master['b'] = np.where(mdot_ratio >= 8.5, 73/mdot_ratio**2, 1)
df_master['b'] = np.where(mdot_ratio >= 150, 3.2E-3, df_master['b'])

df_master['Lx'] = df_master['LXtot']/df_master['b']

df_master['ratio'] = df_master['LXtot'] / df_master['LEdd']
df_master['ratio_beamed'] = df_master['Lx'] / df_master['LEdd']

df_master['theta'] = 2 * np.arccos(1-df_master['b']) #full opening angle in rad
df_master['theta_deg'] = df_master['theta'] * 180/np.pi #deg
df_master['theta_half_deg'] = df_master['theta_deg'] / 2 #Half opening angle



#==========================================================================#
tage_dict={}
Z_dict={}
df_master = df_master.sort_values('Z')
df_master = df_master.sort_values('tage')
#df_master = df_master[df_master['Lx'] > 1E39]
#df_master = df_master[df_master['b'] < 1]

fig, axarr = plt.subplots(3, 10)


for tage, i in zip(df_master.tage.unique(), range(len(df_master.tage.unique()))):
    for z, j in zip(df_master.Z.unique(), range(len(df_master.Z.unique()))):
        tage_dict[tage] = df_master[df_master['tage'] == tage]
        Z_dict[z] = tage_dict[tage][df_master['Z'] == z]
        is_ns = Z_dict[z]['m'] <= 2.5
        df_ns = Z_dict[z][is_ns]
        df_bh = Z_dict[z][~is_ns]
        
        #mass vs luminosity
        axarr[j, i].scatter(df_ns['m'], df_ns['Lx'], s=1.0, color='b')
        axarr[j, i].scatter(df_bh['m'], df_bh['Lx'], s=1.0, color='black')
        axarr[j, i].axhline(y=1E39, color='r')
        
        #Mass Distribution Histograms
#        axarr[j, i].hist(df_ns['m'])
#        axarr[j, i].hist(df_ns['m'])
        
        #axarr[j, i].text(5, 0.9E39, str(len(df_bh) + len(df_ns)), color='purple')
        
        axarr[0, i].set_title(tage)
        axarr[j, 0].set_ylabel(z)
        axarr[j, i].set_xscale('log')
        axarr[j, i].set_yscale('log')
        axarr[j, i].tick_params(axis='both', which = 'both', labelbottom=False, labelleft=False)

plt.show()




'''
#Add type column:
df_master['type'] = np.where(df_master['m'] < 2.5, 'NS' , 'BH')

#PIVOT TABLE TO COUNT
pd.pivot_table(df_master, index = ['Z','tage'], aggfunc='count')
pd.pivot_table(df_master, index = ['Z','tage'], aggfunc='count', columns='type')


#Mass Distribution Histrogram;
plt.hist(df_bh['m'], bins = 10)
plt.hist(df_ns['m'], bins = 10)

plt.xlabel('m ($M_{\odot}$)')
#plt.ylabel('m_dot ($M_{\odot} Myr^{-1}$)')
#plt.ylabel('b/x')
plt.ylabel('Luminosity $(erg s^{-1})$')

plt.xscale('log')
plt.yscale('log')

is_ulx = df_master['Lx'] > 1E39
df_ULXs = df_master[is_ulx]

is_ns = df_master['m'] <= 2.5
df_ns = df_master[is_ns]
df_bh = df_master[~is_ns]

plt.scatter(df_ns['m'],df_ns['Lx'], s=1.0, color='blue', label='Neutron Star = '+ str(len(df_ns)))
plt.scatter(df_bh['m'],df_bh['Lx'], s=1.0, color='black', label='Black Hole = '+ str(len(df_bh)))

is_ns_ulx = df_ns['Lx'] > 1E39
is_bh_ulx = df_bh['Lx'] > 1E39
df_ns_ulx = df_ns[is_ns_ulx]
df_bh_ulx = df_bh[is_bh_ulx]

plt.title('#NS: ' + str(len(df_ns_ulx)) + ' #BH: ' + str(len(df_bh_ulx)) + ' #T: ' + str(len(df_master)))
plt.axhline(y=1E39)
plt.legend()
plt.show()
'''