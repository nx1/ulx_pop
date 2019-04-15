#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Feb 14 15:38:42 2019

@author: nk7g14

This script runs through the output curves generated from the xspec model 
ulxlc and performs alive/dead time analysis for each of them.

Lightcurves are stored in .txt files
Each lightcurve has an associated .info file stored in /info/

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

def FindCurves(df, index, dincl):
    '''
    Finds the two curves for a given simulation number and dincl
    corresponding two different inclinations, 0 and random.
    
    df -- dictionary containing all light curves
    '''
    CurveDict = {}
    for i in df:
        num = i[:-4].split('-')[0]
        incl = i[:-4].split('-')[1]
        if num == str(index) and incl == str(dincl):
            CurveDict[i] = df[i]
    return CurveDict


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

    

#Importing Files
info_files = glob.glob('./curves/info/*.info')
txt_files = glob.glob('./curves/*.txt')
               
print('Loading info files...')
try:
    info_dict #Saves time by checking if the dict is already loaded in
except:
    info_dict={}    #Dictonary for storing info associated info files
    for filename in info_files:
        F = open(filename, 'r')
        info_dict[filename[7:]] = F.read()
        F.close()
print('Loading info files...DONE!')

print('Loading curve files...')
try:
    df_dict     #Saves time by checking if the dict is already loaded in
except:
    df_dict={}      #Dictionary for storing dataframes for lightcurves
    for filename in txt_files:
        df_dict[filename] = pd.read_csv(filename, delimiter=' ',
               header=None, names=['Time', 'Time_Err', 'Flux'], skiprows=3)
print('Loading curve files...DONE!')


df_master = pd.read_csv('dataframe.csv')

df = df_master[df_master['Lx'] < 1E39]
df = df_master[df_master['b'] < 1]
df = df.reset_index()
df.columns
df = df.drop(columns=['index', 'Unnamed: 0'])

    
#Paramaters we wish to test
dincl_list = [5.0, 10.0, 15.0, 20.0, 25.0, 30.0, 35.0, 40.0, 45.0] 
z_list = [0.02, 0.002, 0.0002]
tage_list = [10,20,50,100,200,500,1000,2000,5000,10000]

df_alive = {}   #Dictionary for storing Alive/Dead times for each lightcurve
#fig, axarr = plt.subplots(3, 10)

N_lim = 0
for i in range(len(df_dict)):
    print('==========',i,'/',len(df_dict),'==========')
    for dincl in dincl_list:    #Loop over all precession angles
        curves = FindCurves(df_dict, i, dincl)
        print('hello')
        for key in curves:  #Loop through the found curves
            print('key')
            inclination = key[:-4].split('-')[2]  #Split key name for incl
            print('inclination:', inclination)
            if inclination == '0':
                max_flux = max(curves[key]['Flux']) #Maximum Flux
                c = Lx / max_flux                   #Scaling factor
                N_lim = 1E39 / c                    #Limiting value
                
                print('max flux:', max_flux)
                print('min flux:', min(curves[key]['Flux']))
                print('c:', c)
                print('N_lim:', N_lim)
            else:
                pass
            
#            print('c:', c)
#            print('N_lim:', N_lim)
            
            Alive, Dead = AliveTime(df_dict, key, N_lim)
#            print('Alive:', Alive, 'Dead:', Dead, 'total:', Alive+Dead)
#            print('=====================')
#            plt.scatter(b, np.divide(Alive, Dead),
#                        color=color_dict[dincl], marker='.')
            '''
            if np.divide(Alive, Dead) > 0.5:
                axarr[z_list.index(z), tage_list.index(tage)].scatter(dincl,
                      np.divide(Alive, Dead), marker='.')
                '''
            try:
                percent = Alive/(Alive+Dead)
                df_alive[key] = [Alive, Dead, percent,b,z,tage,Lx,dincl,inclination]
            except:
                df_alive[key] = [Alive, Dead, 'n/a',b,z,tage,Lx,dincl,inclination]
                
            
            
            #print(df_alive[key])


df_a = pd.DataFrame.from_dict(df_alive, orient='index')
df_a.columns = ['alive', 'dead', 'percent','b','z','tage','Lx','dincl','inclination']
df_a = df_a[df_a['inclination']!=0]













def PlotCurve(key):
    time = df_dict[key]['Time']
    flux = df_dict[key]['Flux']
    return plt.plot(time,flux)

def GetSimulationInfo(key):
    split_key = key.split('-')
    sim_number = int(split_key[0].split('/')[-1])
    dincl = split_key[1]
    inclination = split_key[2].split('.')[0]
    row = df.loc[sim_number] 
    row['dincl'] = dincl
    row['inclination'] = inclination
    return row

    

'''
Things you can plot:
    FOR DIFFERENT Z AND TAGE:
        Alive Time vs dincl
        Dead Time vs dincl
        Ratio vs dincl
        beaming vs alive time
        
'''
