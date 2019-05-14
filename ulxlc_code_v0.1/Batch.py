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
from multiprocessing import Pool
import os
import glob

#df_master contains the whole dataset from the STARTRACK population synthesis
#code
df_master = pd.read_csv('dataframe.csv')

df = df_master[df_master['Lx'] > 1E39]
df = df[df['b'] < 1]
df = df[df['theta_half_deg'] < 45]
df = df.reset_index()
df.columns
df = df.drop(columns=['index', 'Unnamed: 0'])

def MakeXCM(theta, incl, dincl, filename, i):
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
    F = open('xspec{}.xcm'.format(i),'w') 
    
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
    
    
'''
dincls = [5.0, 10.0, 15.0, 20.0, 25.0, 30.0, 35.0, 40.0, 45.0]

for i in range(len(df)):
#    print(incl, dincl, fileName)
    for dincl in dincls:
        incl = 0
        fileName = str(i) + '-' + str(round(dincl,2)) + '-' + str(round(incl,2))
        
        print('Making XCM file')
        MakeXCM(df['theta_half_deg'][i], incl, dincl, fileName)
        
        print('Calling xspec')
        subprocess.call(['xspec - xspec.xcm'], shell=True)
        
    for dincl in dincls:
        incl = np.random.uniform(0,90)
        fileName = str(i) + '-' + str(round(dincl,2)) + '-' + str(round(incl,2))
        
        print('Making XCM file')
        MakeXCM(df['theta_half_deg'][i], incl, dincl, fileName)
        
        print('Calling xspec')
        subprocess.call(['xspec - xspec.xcm'], shell=True)
     '''   
      
     
#MultiProcessing Function
def mp(i):
    
    dincls = [5.0, 10.0, 15.0, 20.0, 25.0, 30.0, 35.0, 40.0, 45.0]

#    print(incl, dincl, fileName)
    for dincl in dincls:
        incl = 0
        fileName = str(i) + '-' + str(round(dincl,2)) + '-' + str(round(incl,2))
        
        print('Making XCM file')
        MakeXCM(df['theta_half_deg'][i], incl, dincl, fileName, i)
        
        print('Calling xspec')
        subprocess.call(['xspec - xspec{}.xcm'.format(i)], shell=True)
        
    for dincl in dincls:
        incl = np.random.uniform(0,90)
        fileName = str(i) + '-' + str(round(dincl,2)) + '-' + str(round(incl,2))
        
        print('Making XCM file')
        MakeXCM(df['theta_half_deg'][i], incl, dincl, fileName, i)
        
        print('Calling xspec')
        subprocess.call(['xspec - xspec{}.xcm'.format(i)], shell=True)

# =============================================================================
# Main code
# =============================================================================
    
p = Pool(6)


for i in np.arange(32,36):
    p.map(mp, range(len(df)))
    xcm_files = glob.glob('*.xcm')  
    txt_files = glob.glob('*.txt')
    
    for filePath in xcm_files:
        try:
            os.remove(filePath)
        except:
            print("Error while deleting file : ", filePath)
            
    os.mkdir('curves{}'.format(i))
    for filePath in txt_files:    
        os.rename(filePath, './curves{}/{}'.format(i, filePath))  