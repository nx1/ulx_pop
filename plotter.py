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
from multiprocessing import Pool

# =============================================================================
# Functions
# =============================================================================

def LoadSystems():
    '''
    Loads all systems from STARTRACK from csv files
    '''
    df_master = pd.read_csv('dataframe.csv')
    df = df_master[df_master['Lx'] < 1E39]
    df = df_master[df_master['b'] < 1]
    df = df.reset_index()
    df = df.drop(columns=['index', 'Unnamed: 0'])
    return df


def LoadCurves(folder):
    '''
    Loads all simulation curves outputted from xspec model ulxlc
    '''
    #Importing Files
    curve_files = glob.glob('./{}/*.txt'.format(folder))
    print('Loading curve files...')
    df_dict = {}      #Dictionary for storing dataframes for lightcurves
    for filename in curve_files:
        split1 = filename.split('/')
        split2 = split1[-1][:-4]
        df_dict[split2] = pd.read_csv(filename, delimiter=' ',
               header=None, names=['Time', 'Time_Err', 'Flux'], skiprows=3)
    print('Loading curve files...DONE!')
    return df_dict


def FindCurves(df, num, dincl):
    '''
    Returns the lightcurves corresponding to a
    given simulation number and dincl
    '''
    CurveDict = {}
    for i in df.keys():
        split_key = i.split('-')
        sim_num = split_key[0]
        sim_dincl = split_key[1]
        #sim_inclination = split_key[2]
#        print('key:', i)
#        print('num: {} sim_num: {}'.format(num,sim_num))
#        print('dincl {} sim_dincl {}'.format(dincl, sim_dincl))
#        print('inclination', sim_inclination)
        if str(num) == str(sim_num) and str(sim_dincl) == str(dincl):
            #print('found')
            CurveDict[i] = df[i]
    return CurveDict


def GetSimulationInfo(key):
    '''
    Returns the row of simulation info for a given simulation key
    '''
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
        Dead = np.sum(df_above['time_diff']) + df_above['Time'].iloc[0]
        Alive = np.sum(df_below['time_diff'])
    elif df2['Above'][0] == 1: #If the curve starts above the limit
        df_above = df2[df2['Above'] == 1]
        df_below = df2[df2['Above'] == -1]
        Dead = np.sum(df_above['time_diff'])
        Alive = np.sum(df_below['time_diff']) + df_below['Time'].iloc[0]

    return Alive, Dead


def PlotCurve(key):
    #Plots a specific curve based on key
    splitkey = key.split('-')
    sim_num = splitkey[0]
    dincl = splitkey[1]
    inclination = splitkey[2]

    curves = FindCurves(df_dict, sim_num, dincl)

    N_lim = Normalise(curves)

    fig, ax = plt.subplots()

    for i in curves.keys():
        ax.plot(curves[i]['Time'], curves[i]['Flux'], label=i)

    ax.axhline(y=N_lim, color='r', linestyle='-', label='limit = '+str(N_lim))
    ax.legend()
    return plt.show()


def Normalise(curves):
    '''
    inputs :
        curves = dictionary of two curves
    Takes two curves and Normalises them based on inclination
    '''
    N_lim = 0
    for key in curves:
        splitkey = key.split('-')
        sim_num = int(splitkey[0])
        inclination = splitkey[-1]
        Lx = GetLx(sim_num)
#        print('-----------------')
#        print('Curve:', key)
#        print('Lx:', Lx)
        if inclination == '0':

            max_flux = max(curves[key]['Flux']) #Maximum Flux
            c = Lx / max_flux                   #Scaling factor
            N_lim = 1E39 / c                    #Limiting value
#            print('Curve Maximum:', max_flux)
#            print('Normalisation factor c:', c)
#            print('Resulting N_lim:', N_lim)
        else:
            pass
    return N_lim


def LookAtULX(df_dict, key):
    '''
    This method will attempt to calculate the alive/dead probability
    in a different manner.
    Instead of the previous method of calculating the ratio of time
    above and below the curve, this method will randomly select a time
    on the time series to observe the curve.
    We may either do this in an instantanous fashion or use a time interval
    that could correspond to E-ROSITA exposure time. the issue with this, is
    that the the time output from the system is in abitrary time units and so
    constraints on the period of Lense-Thirring precession would be required in
    order to normalise it.

        Pick random value in time range
    check to see if it's below or above N_lim
    return for each system the number of times looked vs the number of times
    alive
    '''
    split_key = key.split('-')
    sim_num = split_key[0]
    dincl = split_key[1]

    curves = FindCurves(df_dict, sim_num, dincl) #Two curves
    N_lim = Normalise(curves) #Find normalization limit
    for k in curves:
        inclination = k.split('-')[-1]
        if inclination == '0':
            pass
        else:
#            print('Found system of inclination:', inclination)
            time = curves[k]['Time']
            flux = curves[k]['Flux']
            rand = np.random.randint(low=0, high=len(time))
#            print('selected random time bin:', rand)
#            print('timebin corresponds to: t =', time[rand])
#            print('Associated flux is: f = ', flux[rand])
            if flux[rand] > N_lim:
                alive = 1
#                print('Flux:', flux[rand], 'N_lim:', N_lim)
#                print('Alive!')
            else:
#                print('Flux:', flux[rand], 'N_lim:', N_lim)
#                print('Dead :(')
                alive = 0
    return alive

def ResultsDictToPandas(r_dict):
    incls = []
    dincls_list = []
    folder_list = []
    df_a = pd.DataFrame.from_dict(results_dict, orient='index')

    for i, row in df_a.iterrows():
        split1 = i.split(':')
        folder_num = split1[0]
        sim_things = split1[1]

        split_key = sim_things.split('-')
        inclination = split_key[-1]
        dincl = split_key[1]
        incls.append(float(inclination))
        dincls_list.append(float(dincl))
        folder_list.append(int(folder_num))


    df_a['inclination'] = incls
    df_a['dincl'] = dincls_list
    df_a.columns = ['alive', 'dead', 'inclination', 'dincl']
    df_a['ratio'] = df_a['alive']/(df_a['dead']+df_a['alive'])
    df_a['folder_num'] = folder_list 
    return df_a

###############################################################################
###############################################################################
#########################====MAIN CODE====#####################################
###############################################################################
###############################################################################


# =============================================================================
# Loading curves & Systems df
# =============================================================================
try:
    df_dict
except:
    df_dict = LoadCurves('curves0')
df = LoadSystems()

# =============================================================================
# Calculate Alive/Dead time for systems
# =============================================================================
dincl_list = [5.0, 10.0, 15.0, 20.0, 25.0, 30.0, 35.0, 40.0, 45.0]
z_list = [0.02, 0.002, 0.0002]
tage_list = [10, 20, 50, 100, 200, 500, 1000, 2000, 5000, 10000]
results_dict = {}

for i in range(1):
    folder = 'curves{}'.format(i)   #Folder String
    df_dict = LoadCurves(folder)    #Load Curves from folder to dict
    #Paramaters we wish to test

    for sim_num in range(len(df)):  #Go through all systems
        Lx = GetLx(sim_num)         #Get Lx for system
        for dincl in dincl_list:    #Run through all precession angles
            curves = FindCurves(df_dict, sim_num, dincl) #Two curves
            N_lim = Normalise(curves) #Find normalization limit
            for key in curves:
                Alive, Dead = AliveTime(df_dict, key, N_lim)
                results_dict[str(i)+':'+key] = Alive, Dead
        print(sim_num, '/', len(df))

df_a = ResultsDictToPandas(results_dict)
df_a = df_a[df_a['inclination'] != 0]
df_a_nonzero = df_a[df_a['alive'] != 0]
df_a_nonzero = df_a_nonzero[df_a_nonzero['dead'] != 0]

pivot = pd.pivot_table(df_a_nonzero, index=['dincl'],
                       aggfunc=('count', 'mean', 'std'))

for index, mean in zip(pivot.index, pivot[('ratio', 'mean')]):
    plt.scatter(index, mean, label=index)

plt.legend()


# =============================================================================
# dincl vs ratio all for each simulation
# =============================================================================
plt.title('28 runs, 151 systems, ratio vs dincl')
plt.xlabel('dincl')
plt.ylabel('ratio')
for i in range(28):
    chunk = pivot[('ratio', 'mean')][i]
    plt.plot(chunk, label=str(i))

#dincl vs ratio for all simulations averaged
pivot2 = pd.pivot_table(df_a_nonzero, index =['dincl'], aggfunc=('count', 'mean', 'std'))
plt.plot(pivot2.index, pivot2['ratio', 'mean'],color='black', linewidth=5.0)


# =============================================================================
# Calculate explicit alive/dead for each system
# =============================================================================
results_dict = {}
for sim_num in range(len(df)):
    Lx = GetLx(sim_num) #Get Lx for system
    for dincl in dincl_list:    #Run through all precession angles
        curves = FindCurves(df_dict, sim_num, dincl) #Two curves
        N_lim = Normalise(curves) #Find normalization limit
        for key in curves:
            Alive, Dead = AliveTime(df_dict, key, N_lim)
            results_dict[key] = Alive, Dead
    print(sim_num,'/',len(df))



# =============================================================================
# Look at ULX method
# =============================================================================
alive_dict_looking = {}
observations = 1000 #Number of observations per system

for key, i in zip(df_dict.keys(),range(len(df_dict))):
    print(i, '/', len(df_dict))
    #print(key)
    alive_sum = 0
    for i in range(observations):
        alive = LookAtULX(df_dict, key)
        alive_sum+=alive
    alive_dict_looking[key] = [alive_sum, observations, alive_sum/observations]

# =============================================================================
# MULTIPROCESSING
# =============================================================================
'''
p = Pool(6)

def multiprocess(key):
    observations = 10
    alive_list = [0]*observations
    alive_sum = 0
    for i in range(observations):
        alive = LookAtULX(df_dict, key)
        alive_sum+=alive
    return alive_sum


for key in df_dict.keys():
    print(key)
    maped = p.map(multiprocess, key)
    print('Alive_sum:', list(maped))
'''
#Splitting by dincl
df_a_looking = pd.DataFrame.from_dict(alive_dict_looking, orient='index')
df_a_looking.columns = ['alive', 'observed', 'ratio']

incls = []
dincls_list = []
for i, row in df_a_looking.iterrows():
    split_key = i.split('-')
    inclination = split_key[-1]
    dincl = split_key[1]
    incls.append(float(inclination))
    dincls_list.append(float(dincl))

df_a_looking['inclination'] = incls
df_a_looking['dincl'] = dincls_list

df_a_looking = df_a_looking[df_a_looking['inclination'] != 0]
df_a_looking = df_a_looking[df_a_looking['alive'] != 0]
df_a_looking = df_a_looking[df_a_looking['alive'] != observations]


pivot = pd.pivot_table(df_a_looking, index = ['dincl'], aggfunc=('count', 'mean', 'std'))
plt.scatter(pivot.index, pivot['ratio']['mean'])
#Histogram of dincl
dincl_plot = 15.0
plt.hist(df_a_looking[df_a_looking['dincl']==dincl_plot]['ratio'], label=dincl_plot)
plt.legend()


# =============================================================================
# Dictionary to pandas and analysis    
# =============================================================================
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

piv1 = pd.pivot_table(df_a, index = ['dincl'], aggfunc=('count', 'mean', 'std'))
piv2 = pd.pivot_table(df_a_nonzero, index = ['dincl'], aggfunc=('count', 'mean', 'std'))

plt.xlabel('dincl $ \Delta i $')
plt.ylabel('Alive/Dead Ratio')
plt.plot(piv1.index, piv1[('ratio', 'mean')], label = 'all i=0 sources')
plt.plot(piv2.index, piv2[('ratio', 'mean')], label = 'all nonzero i=0 sources')
plt.legend()

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
    cut = df_a_nonzero[df_a['dincl'] == i]
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
    plt.scatter(np.mean(cut['dincl']), np.mean(cut['ratio']), label=i)


plt.xlabel('dincl')
plt.ylabel('alive/dead')
plt.legend()
#plt.savefig('figures/dincl_alive_nonzero_mean.eps', format='eps', dpi=1000)
#plt.savefig('figures/dincl_alive_nonzero_mean.png', format='png', dpi=1000)

# =============================================================================
# SPECIFIC LIGHTCURVES and LIMIT
# =============================================================================
PlotCurve('0-10.0-0')

'''
for i in np.arange(30,40):
    print(str(i)+'-10.0-0')
    PlotCurve(str(i)+'-10.0-0')
'''

'''
Things you can plot:
    FOR DIFFERENT Z AND TAGE:
        Alive Time vs dincl
        Dead Time vs dincl
        Ratio vs dincl
        beaming vs alive time
        
'''
