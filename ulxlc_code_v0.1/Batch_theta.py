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

df = df_master[df_master['Lx'] > 1E39]
df = df[df['b'] < 1]
df = df.reset_index()
df.columns
df = df.drop(columns=['index', 'Unnamed: 0'])

def MakeXCM(theta, incl, dincl, filename):
    '''
    Creates xcm script file for use in xspec
    newpar 1  --> period
    newpar 2  --> phase
    newpar 3  --> theta
    newpar 4  --> incl
    newpar 5  --> dincl
    newpar 6  --> beta
    newpar 7  --> dopulse
    newpar 8  --> norm
    '''
    F = open('xspec.xcm','w') 
    
    string='''lmod ulxlc /home/nk7g14/Desktop/gitbox/ulx_pop/ulxlc_code_v0.1
model ulxlc & /*
newpar 1  10.0
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

for i in range(45):
    MakeXCM(i,0,15,i)
    subprocess.call(['xspec - xspec.xcm'], shell=True)