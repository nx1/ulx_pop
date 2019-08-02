#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
@author: nk7g14

This script is used for running the XSPEC model 'ulxlc' for a variety of
different model parameters.
"""

import pandas as pd
import subprocess
import numpy as np
import os
import glob
from multiprocessing import Pool
import shutil
from tqdm import tqdm

def CreateOutputFileName(i, parameters):
    params_rounded = np.round(parameters, 2)

    output_file_name = str(i) + '-' + str(params_rounded[4]) + '-' + str(params_rounded[3])
    return output_file_name


def MakeXCM(filename, parameters, index, save_folder):
    '''
    Creates xcm script file for use in xspec
        filename: name of xcm file with extension
        parameters: list of model parameters, listed below
    NOTE THESE COUNT UP FROM 1
    newpar 1  --> period
    newpar 2  --> phase
    newpar 3  --> theta
    newpar 4  --> incl
    newpar 5  --> dincl
    newpar 6  --> beta
    newpar 7  --> dopulse
    newpar 8  --> norm
    '''
    F = open(filename, 'w')

    params_rounded = np.round(parameters, 2)

    output_file_name = CreateOutputFileName(index, parameters)
    savedir = str(save_folder) + '/' + output_file_name
    cwd = os.getcwd()

    string = '''lmod ulxlc {}
model ulxlc & /*
newpar 1  {}
newpar 2  {}
newpar 3  {}
newpar 4  {}
newpar 5  {}
newpar 6  {}
newpar 7  {}
newpar 8  {}
cpd /null
plot model
setplot command wdata {}.txt
plot
exit''' .format(cwd, *params_rounded, savedir)

    F.write(string)
    F.close()
    # print('Made XCM file to : {:<10}'.format(filename))
    header = '{:<10} {:<10} {:<10} {:<10} {:<10} {:<10} {:<10} {:<10}'
    # print(header.format('period', 'phase', 'theta', 'incl', 'dincl', 'beta', 'dopulse', 'norm'))
    # print(header.format(*params_rounded))


def RunXCM(XCMfile):
    '''
    XCMfile: XCM file including extension
    '''
    devnull = open(os.devnull, 'w')
    subprocess.call(['xspec - {}'.format(XCMfile)], shell=True, stdout=devnull)
    # subprocess.call(['xspec - {}'.format(XCMfile)], shell=True)

def DeleteAllXCMFiles():
    cwd = os.getcwd()
    os.chdir(cwd)
    files = glob.glob('*.xcm')
    if files == []:
        print('No XCM files to delete')
    else:
        for file in files:
            os.remove(file)


def DeleteAlltxtFiles():
    cwd = os.getcwd()
    os.chdir(cwd)
    files = glob.glob('*.txt')
    if files == []:
        print('No txt files to delete')
    else:
        for file in files:
            os.remove(file)


def ReadLightCurve(filename):
    df = pd.read_csv(filename, delimiter=' ',
                     header=None, names=['Time', 'Time_Err', 'Flux'], skiprows=3)
    return df


def AliveTime(lc, limit):
    '''
    Calculates for a given curve the amount of time above and below a given
    limit
    '''

    if max(lc['Flux']) < limit:
        Alive = 0
        Dead = lc['Time'].iloc[-1]
        return Alive, Dead
    if min(lc['Flux']) > limit:
        Alive = lc['Time'].iloc[-1]
        Dead = 0
        return Alive, Dead

    curve = lc
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

def FilterNSBH(df):
    df_bh = df[df['is_bh'] == 1]
    df_ns = df[df['is_bh'] == 0]
    return df_bh, df_ns
    

def ChooseSystem(BH_NS, df_bh, df_ns):
    draw = np.random.random()
    if draw > BH_NS:
        choice = np.random.choice(df_ns.index)
    else:
        choice = np.random.choice(df_bh.index)
    return choice


def isAlwaysVisible(df, index):
    if (df['b'][index] >= 1) or (df['theta_half_deg'][index] > 45):
        return True
    else:
        return False

def simulate(dincl):
    
    parameters = [period, phase, theta, incl, dincl, beta, dopulse, norm]
    
    # print('Making XCM file')
    XCM_file_name = str(round(dincl,2)) + 'xspec.xcm'
    MakeXCM(XCM_file_name, parameters, c, simulation_number)
    
    # print('Calling xspec')
    RunXCM(XCM_file_name)
    os.remove(XCM_file_name)

# =============================================================================
# Main Code
# =============================================================================

if __name__ == '__main__':
    df_master = pd.read_csv('../dataframe.csv') #36420 SYSTEMS
    df = df_master[df_master['Lx'] > 1E39]  #Only ULX -     992 ULXs
    # df = df[df['b'] < 1]                    #Only Beamed -  227 Beamed ULXs
    # df = df[df['theta_half_deg'] < 45]      #thetha < 45 -  151 Beamed ULXs with half opening angles < 45
    df = df.reset_index()
    df = df.drop(columns=['index', 'Unnamed: 0'])
    
    number_of_simulations = 500
    
    period = 10.0
    phase = 0.0
    theta = None
    incl = None
    dincl = None
    beta = 0.2
    dopulse = 0
    norm = 1.0
    
    df_bh, df_ns = FilterNSBH(df)
    
    os.makedirs('./curves', exist_ok=True)
    df_bh, df_ns = FilterNSBH(df)       
    for BH_NS in np.arange(0.01, 0.1, 0.01):
        os.makedirs('./curves/{}'.format(BH_NS), exist_ok=True)
        pbar = tqdm(range(number_of_simulations))
        for simulation_number in pbar:
            os.makedirs('{}'.format(simulation_number), exist_ok=True)
            c = ChooseSystem(BH_NS, df_bh, df_ns)
            # print('c:', c)
            if isAlwaysVisible(df, c):
                pbar.set_description('%s Always Visible!' % c)
            else:
                pbar.set_description('%s transient!' % c)
                dincls = np.linspace(1.0, 45, 500)
                theta = df['theta_half_deg'][c]
                for isRandom in range(2):
                    if isRandom:
                        incl = np.random.uniform(0,90)
                    else:
                        incl = 0
                    pool = Pool()
                    pool.map(simulate, dincls)
                    pool.close()
            shutil.move('./{}'.format(simulation_number), './curves/{}'.format(BH_NS))
    


'''
# =============================================================================
# Single Run Test
# =============================================================================
t0 = time.time()
#It is possible to test any parameter by using it in a loop

period = 10.0
phase = 0.0
theta = 30
incl = 10.0
dincl = 20
beta = 0.4
dopulse = 0
norm = 1.0

i=5
limit = 5

plt.figure()



alive_deads = []

# incl = np.random.randint(0, high=90)
par = [period, phase, theta, incl, dincl, beta, dopulse, norm]

MakeXCM('xspec.xcm', par, i)
RunXCM('xspec.xcm')
output_file_name = CreateOutputFileName(i, par) + '.txt'

lc = ReadLightCurve(output_file_name)

alive, dead = AliveTime(lc, limit)
ratio = np.divide(alive,dead+alive)
alive_deads.append(ratio)
plt.plot(lc['Time'], lc['Flux'], label=(r'$di$ = '+str(dincl)+' : '+str(np.round(ratio, 2))))

DeleteAlltxtFiles()
time_taken = time.time()-t0
print('Time Taken:', time_taken)
print('Predicted time for 1000 iterations:', (time_taken*1000/60), 'minutes') 

plt.legend()
plt.axhline(y=limit, c='r')
plt.title("""period, phase, theta, incl, dincl, beta, dopulse, norm
          {}""".format(par))
'''
