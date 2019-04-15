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


#df_master contains the whole dataset from the STARTRACK population synthesis
#code
df_master = pd.read_csv('dataframe.csv')

df = df_master[df_master['Lx'] < 1E39]
df = df_master[df_master['b'] < 1]
df = df.reset_index()



def MakeString(df, index, dincl):
    '''
    Creates long string of column names for appending to output file.
    
    Inputs:
        df = dataframe
        index = numbered index of item in question
        dincl = current dincl used
    '''
    longstring = ''
    for col in df.columns:
        longstring += col + ':' + str(df[col][index]) + '\n'
    longstring += 'dincl:' + str(dincl)
    return longstring

def MakeXCM(theta, incl, dincl, filename):
    '''
    Creates xcm script file for use in xspec
    '''
    F = open('xspec.xcm','w') 
    
    string='''lmod ulxlc /home/nk7g14/Desktop/gitbox/ulx_pop/ulxlc_code_v0.1
model ulxlc & /*
newpar 1  1.0
newpar 2  0.0
newpar 3  %s
newpar 4  %s
newpar 5  %s
newpar 6  0.5
newpar 7  0
newpar 8  1.0
cpd /xw
plot model
setplot command wdata %s.txt
plot
exit''' %(round(theta,2), round(incl,2), dincl, filename)
    
    F.write(string)
    F.close()
    
    
def MakeInfo(string, filename):
    F = open('ulxlc_code_v0.1/info/' + filename, 'w+')
    F.write(string)
    F.close()


dincls = [5.0, 10.0, 15.0, 20.0, 25.0, 30.0, 35.0, 40.0, 45.0]


for i in range(len(df)):
    
    for dincl in dincls:        
        incl = 0
        longstr = MakeString(df, i, dincl)
        fileName = str(i) + '-' + str(round(dincl,2)) + '-' + str(round(incl,2))
        
        print('Making XCM file')
        MakeXCM(df['theta_half_deg'][i], incl, dincl, fileName)
        
        print('Making info file')
        MakeInfo(longstr, fileName + '.info')
        
        print('Calling xspec')
        subprocess.call(['xspec - xspec.xcm'], shell=True)
        
    for dincl in dincls:        
        incl = np.random.uniform(0,90)
        longstr = MakeString(df, i, dincl)
        fileName = str(i) + '-' + str(round(dincl,2)) + '-' + str(round(incl,2))
        
        print('Making XCM file')
        MakeXCM(df['theta_half_deg'][i], incl, dincl, fileName)
        
        print('Making info file')
        MakeInfo(longstr, fileName + '.info')
        
        print('Calling xspec')
        subprocess.call(['xspec - xspec.xcm'], shell=True)
    