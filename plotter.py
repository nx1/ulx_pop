#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Feb 14 15:38:42 2019

@author: nk7g14

This script runs through the output curves generated from the xspec model 
ulxlc and performs alive/dead time analysis for each of them.

Lightcurves are stored in .txt files
Each lightcurve has an associated .info file stored in /info/
.
The filename format goes as:
    simulation_number - dincl - inclination.txt

where dincl is the precession angle

1) Load in all lightcurves into df_dict
2) Load in all info files into info_dict

For each simulation_number dincl combo there is a 0 inclination system and
a random inclination system.

3) Find the 0 inclination system and set the scaling factor, c, and the limiting
value, N_lim.

c = Lx / max(curves[key]['Flux']) #Scaling factor
N_lim = 1E39 / c                  #Limiting value
"""
import pandas as pd
import matplotlib.pyplot as plt
import glob
import numpy as np
import time

def LoadSystems():
    df_master = pd.read_csv('dataframe.csv')
    df = df_master[df_master['Lx'] < 1E39]
    df = df_master[df_master['b'] < 1]
    df = df.reset_index()
    df.columns
    df = df.drop(columns=['index', 'Unnamed: 0'])
    return df

def LoadCurves():
    #Importing Files
    curve_files = glob.glob('./curves/*.txt')
    print('Loading curve files...')
    df_dict={}      #Dictionary for storing dataframes for lightcurves
    for filename in curve_files:
        df_dict[filename[9:-4]] = pd.read_csv(filename, delimiter=' ',
               header=None, names=['Time', 'Time_Err', 'Flux'], skiprows=3)
    print('Loading curve files...DONE!')
    return df_dict

def FindCurves(df, num, dincl):
    '''
    For a given simulation number, this function returns all associated
    lightcurves
    '''
    CurveDict = {}
    for i in df.keys():
        split_key = i.split('-')
        sim_num = split_key[0]
        sim_dincl = split_key[1]
        sim_inclination = split_key[2]
#        print('key:', i)
#        print('num: {} sim_num: {}'.format(num,sim_num))
#        print('dincl {} sim_dincl {}'.format(dincl, sim_dincl))
#        print('inclination', sim_inclination)
        if str(num) == str(sim_num) and str(sim_dincl) == str(dincl):
            #print('found')
            CurveDict[i] = df[i]
    return CurveDict


def GetSimulationInfo(key):
    split_key = key.split('-')
    
    sim_number = int(split_key[0])
    dincl = split_key[1]
    inclination = split_key[2]
    
    row = df.loc[sim_number] 
    row['dincl'] = dincl
    row['inclination'] = inclination
    return row

def GetLx(num):
        return df.loc[num]['Lx']

def AliveTime(df_dict, key, limit):
    '''
    Calculates for a given curve the amount of time above and below a given
    limit
    '''
    
    if max(df_dict[key]['Flux']) < limit:
        Alive = 0
        Dead = df_dict[key]['Time'].iloc[-1]
        return Alive, Dead
    if min(df_dict[key]['Flux']) > limit:
        Alive = df_dict[key]['Time'].iloc[-1]
        Dead = 0
        return Alive, Dead
    
    curve = df_dict[key]
    curve['Above'] = curve['Flux'].map(lambda x: 1 if x >= limit else -1)
    
    df2 = curve[curve['Above'].diff() != 0]
    df2['time_diff'] = df2['Time'].diff()
    
    if len(df2) == 1:   #If the Curve never goes above the limit
        Alive = 0.0
        Dead = curve['Time'].iloc[-1]
    elif df2['Above'][0] == -1: #If the curve starts below the limit
        df_above = df2[df2['Above'] == 1]
        df_below = df2[df2['Above'] == -1]
        Dead =  np.sum(df_above['time_diff']) + df_above['Time'].iloc[0]
        Alive = np.sum(df_below['time_diff'])
    elif df2['Above'][0] == 1: #If the curve starts above the limit
        df_above = df2[df2['Above'] == 1]
        df_below = df2[df2['Above'] == -1]
        Dead =  np.sum(df_above['time_diff']) 
        Alive = np.sum(df_below['time_diff']) + df_below['Time'].iloc[0]
        
    return Alive, Dead

    
def PlotCurve(key):
    splitkey  = key.split('-')
    sim_num = splitkey[0]
    dincl =  splitkey[1]
    inclination =  splitkey[2]
    
    curves = FindCurves(df_dict, sim_num, dincl)
    
    N_lim = Normalise(curves)
    
    fig, ax = plt.subplots()
    
    for i in curves.keys():
        ax.plot(curves[i]['Time'],curves[i]['Flux'], label = i)
        
    ax.axhline(y=N_lim, color='r', linestyle='-', label = 'limit')
    ax.legend()
    return plt.show()
    

def Normalise(curves):
    '''
    Takes two curves and Normalises them based on inclination
    '''
    N_lim = 0
    for key in curves:
        splitkey = key.split('-')
        sim_num = int(splitkey[0])
        inclination = splitkey[-1]
        Lx = GetLx(sim_num)
        print('-----------------')
        print('Curve:', key)
        print('Lx:', Lx)
        if inclination == '0':
            
            max_flux = max(curves[key]['Flux']) #Maximum Flux
            c = Lx / max_flux                   #Scaling factor
            N_lim = 1E39 / c                    #Limiting value
            print('Curve Maximum:', max_flux)
            print('Normalisation factor c:', c)
            print('Resulting N_lim:', N_lim)
        else:
            pass
    return N_lim
            

try:
    df_dict
except:
    df_dict = LoadCurves()

df = LoadSystems()

#Paramaters we wish to test
dincl_list = [5.0, 10.0, 15.0, 20.0, 25.0, 30.0, 35.0, 40.0, 45.0] 
z_list = [0.02, 0.002, 0.0002]
tage_list = [10,20,50,100,200,500,1000,2000,5000,10000]

results_dict = {}

N_lim = 0
for i in range(len(df)):
    t1 = time.time()
    Lx = GetLx(i)
    for dincl in dincl_list:
        curves = FindCurves(df_dict, i, dincl)
        for key in curves:
            inclination = key.split('-')[-1]
            if inclination == '0':
                max_flux = max(curves[key]['Flux']) #Maximum Flux
                c = Lx / max_flux                   #Scaling factor
                N_lim = 1E39 / c                    #Limiting value
                
#                print('max flux:', max_flux)
#                print('min flux:', min(curves[key]['Flux']))
#                print('c:', c)
#                print('N_lim:', N_lim)
            else:
                pass
            
            Alive, Dead = AliveTime(df_dict, key, N_lim)
            results_dict[key] = Alive, Dead
    timeleft = round(((time.time() - t1) * 4320 - i)/3600, 2)
    print(i,'/', len(df), '| eta:', timeleft, 'hours')
    




df_a = pd.DataFrame.from_dict(results_dict, orient='index')


incls = []
dincls_list = []
for i, row in df_a.iterrows():
    split_key = i.split('-')
    inclination = split_key[-1]
    dincl = split_key[1]
    incls.append(float(inclination))
    dincls_list.append(float(dincl))

df_a['inclination'] = incls
df_a['dincl'] = dincls_list

df_a.columns = ['alive', 'dead', 'inclination', 'dincl']
df_a['ratio'] = df_a['alive']/(df_a['dead']+df_a['alive'])
df_a = df_a[df_a['inclination']!=0]
df_a_nonzero = df_a[df_a['alive']!=0]
df_a_nonzero = df_a_nonzero[df_a_nonzero['dead']!=0]

# =============================================================================
# ALL SYSTEMS HISTOGRAM
# =============================================================================
plt.title('Alive time distribution for {} sources'.format(len(df_a)))
plt.hist(df_a['alive'], bins=30)
plt.show()

#plt.savefig('figures/alive_hist_all.eps', format='eps', dpi=1000)
#plt.savefig('figures/alive_hist_all.png', format='png', dpi=1000)
#
# =============================================================================
# ALL NONZERO SYSTEMS HISTOGRAM
# =============================================================================
plt.title('Alive time distribution for {} sources'.format(len(df_a_nonzero)))
plt.hist(df_a_nonzero['alive'], bins=30)
plt.show()

#plt.savefig('figures/alive_hist_nonzero.eps', format='eps', dpi=1000)
#plt.savefig('figures/alive_hist_nonzero.png', format='png', dpi=1000)

# =============================================================================
# ALL SYSTEMS ALIVE TIME VS DINCL (PRECESSION ANGLE)
# =============================================================================
plt.scatter(df_a['dincl'], df_a['alive'])
#plt.savefig('figures/dincl_alive_all.eps', format='eps', dpi=1000)
#plt.savefig('figures/dincl_alive_all.png', format='png', dpi=1000)
# =============================================================================
# ALL NON-ZERO SYSTEMS ALIVE TIME VS DINCL (PRECESSION ANGLE)
# =============================================================================
df_a_nonzero_5 = df_a_nonzero[df_a['dincl']==5]
df_a_nonzero_10 = df_a_nonzero[df_a['dincl']==10]
df_a_nonzero_15 = df_a_nonzero[df_a['dincl']==15]
df_a_nonzero_20 = df_a_nonzero[df_a['dincl']==20]
df_a_nonzero_25 = df_a_nonzero[df_a['dincl']==25]
df_a_nonzero_30 = df_a_nonzero[df_a['dincl']==30]
df_a_nonzero_35 = df_a_nonzero[df_a['dincl']==35]
df_a_nonzero_40 = df_a_nonzero[df_a['dincl']==40]
df_a_nonzero_45 = df_a_nonzero[df_a['dincl']==45]

biglist = [df_a_nonzero_5 ,
df_a_nonzero_10,
df_a_nonzero_15,
df_a_nonzero_20,
df_a_nonzero_25,
df_a_nonzero_30,
df_a_nonzero_35,
df_a_nonzero_40,
df_a_nonzero_45]

print('{} {} {} {}'.format('dincl', 'mean', 'std', '# systems'))
for i in biglist:
    print('{} {} {} {}'.format(
            np.mean( i['dincl']),
            round(np.mean(i['ratio']), 3),
            round(np.std(i['ratio']), 3),
            len(i['ratio']),
            ))





for i in dincl_list:
    cut = df_a_nonzero[df_a['dincl']==i]
    plt.scatter(cut['dincl'], cut['alive'], label=i)

plt.xlabel('dincl')
plt.ylabel('alive/dead')
#plt.savefig('figures/dincl_alive_nonzero.eps', format='eps', dpi=1000)
#plt.savefig('figures/dincl_alive_nonzero.png', format='png', dpi=1000)

# =============================================================================
# ALL NON-ZERO SYSTEMS ALIVE TIME VS DINCL (PRECESSION ANGLE) MEANS
# =============================================================================
for i in dincl_list:
    cut = df_a_nonzero[df_a['dincl']==i]
    plt.scatter(np.mean(cut['dincl']), np.mean(cut['alive']), label=i)


plt.xlabel('dincl')
plt.ylabel('alive/dead')
plt.legend()
#plt.savefig('figures/dincl_alive_nonzero_mean.eps', format='eps', dpi=1000)
#plt.savefig('figures/dincl_alive_nonzero_mean.png', format='png', dpi=1000)

# =============================================================================
# SPECIFIC LIGHTCURVES and LIMIT
# =============================================================================
PlotCurve('0-10.0-0')

for i in np.arange(30,40):
    print(str(i)+'-10.0-0')
    PlotCurve(str(i)+'-10.0-0')

'''
Things you can plot:
    FOR DIFFERENT Z AND TAGE:
        Alive Time vs dincl
        Dead Time vs dincl
        Ratio vs dincl
        beaming vs alive time
        
'''
