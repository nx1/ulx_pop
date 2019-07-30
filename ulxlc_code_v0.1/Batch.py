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
import matplotlib.pyplot as plt
import time

def CreateOutputFileName(parameters):
    params_rounded = np.round(parameters, 2)

    output_file_name = str(i) + '-' + str(params_rounded[4]) + '-' + str(params_rounded[3])
    return output_file_name


def MakeXCM(filename, parameters):
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

    output_file_name = CreateOutputFileName(parameters)


    string = '''lmod ulxlc /home/nk7g14/Desktop/gitbox/ulx_pop/ulxlc_code_v0.1
model ulxlc & /*
newpar 1  {}
newpar 2  {}
newpar 3  {}
newpar 4  {}
newpar 5  {}
newpar 6  {}
newpar 7  {}
newpar 8  {}
cpd /xw
plot model
setplot command wdata {}.txt
plot
exit''' .format(*params_rounded, output_file_name)

    F.write(string)
    F.close()
    print('Made XCM file to : {:<10}'.format(filename))
    header = '{:<10} {:<10} {:<10} {:<10} {:<10} {:<10} {:<10} {:<10}'
    print(header.format('period', 'phase', 'theta', 'incl', 'dincl', 'beta', 'dopulse', 'norm'))
    print(header.format(*params_rounded))


def RunXCM(XCMfile):
    '''
    XCMfile: XCM file including extension
    '''
    subprocess.call(['xspec - {}'.format(XCMfile)], shell=True)


def DeleteAllXCMFiles():
    os.chdir('/home/nk7g14/Desktop/gitbox/ulx_pop/ulxlc_code_v0.1')
    files = glob.glob('*.xcm')
    if files == []:
        print('No XCM files to delete')
    else:
        for file in files:
            os.remove(file)


def DeleteAlltxtFiles():
    os.chdir('/home/nk7g14/Desktop/gitbox/ulx_pop/ulxlc_code_v0.1')
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

# =============================================================================
# Main Code
# =============================================================================
'''
if __name__ == '__main__':
    df_master = pd.read_csv('../dataframe.csv') #36420 SYSTEMS
    df = df_master[df_master['Lx'] > 1E39]  #Only ULX -     992 ULXs
    df = df[df['b'] < 1]                    #Only Beamed -  227 Beamed ULXs
    df = df[df['theta_half_deg'] < 45]      #thetha < 45 -  151 Beamed ULXs with half opening angles < 45
    df = df.reset_index()
    df = df.drop(columns=['index', 'Unnamed: 0'])
    
    
    period = 10.0
    phase = 0.0
    theta = None
    incl = None
    dincl = None
    beta = 0.2
    dopulse = 0
    norm = 1.0
    
    
    os.makedirs('./curves', exist_ok=True)
    # dincls = [5.0, 10.0, 15.0, 20.0, 25.0, 30.0, 35.0, 40.0, 45.0]
    dincls = np.random.uniform(low=1.0, high=45, size=1000)
    for i in range(147,152):
        theta = df['theta_half_deg'][i]
        for isRandom in range(2):
            for dincl in dincls:
                if isRandom:
                    incl = np.random.uniform(0,90)
                else:
                    incl = 0
                parameters = [period, phase, theta, incl, dincl, beta, dopulse, norm]
    
                print('Making XCM file')
    
                XCM_file_name = 'xspec.xcm'
                MakeXCM(XCM_file_name, parameters)
    
                print('Calling xspec')
                RunXCM(XCM_file_name)
        
        os.makedirs('./curves/{}'.format(i))
        curve_files = glob.glob('*.txt')
        for file in curve_files:
            os.rename(file, './curves/{}/{}'.format(i, file))
    
           
    
    
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

MakeXCM('xspec.xcm', par)
RunXCM('xspec.xcm')
output_file_name = CreateOutputFileName(par) + '.txt'

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
