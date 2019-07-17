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

#df_master contains the whole dataset from the STARTRACK population synthesis
#code

df_master = pd.read_csv('../dataframe.csv')
df = df_master[df_master['Lx'] > 1E39]  #Only ULX
df = df[df['b'] < 1]                    #Only Beamed
df = df[df['theta_half_deg'] < 45]      #Only half opening angle of > 45
df = df.reset_index()
df = df.drop(columns=['index', 'Unnamed: 0'])

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

    params_rounded = np.round(parameters)
    output_file_name = params_rounded

    output_file_name = str(i) + '-' + str(params_rounded[4]) + '-' + str(params_rounded[3])


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


'''
period = 10.0
phase = 0.0
theta = None
incl = None
dincl = None
beta = 0.2
dopulse = 0
norm = 1.0

dincls = [5.0, 10.0, 15.0, 20.0, 25.0, 30.0, 35.0, 40.0, 45.0]

for i in range(len(df)):
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
'''

# =============================================================================
# Single Run Test
# =============================================================================
period = 10.0
phase = 0.0
theta = 10
incl = 5.0
dincl = 10
beta = 0.2
dopulse = 0
norm = 1.0

i=1

par = [period, phase, theta, incl, dincl, beta, dopulse, norm]

MakeXCM('xspec.xcm', par)
RunXCM('xspec.xcm')

lc = ReadLightCurve('1-10.0-5.0.txt')
plt.step(lc['Time'], lc['Flux'])

