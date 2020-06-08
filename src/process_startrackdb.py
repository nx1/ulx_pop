# -*- coding: utf-8 -*-
"""
Created on Fri May 29 10:27:46 2020

@author: norma


https://universeathome.pl/universe/bhdb.php

Process large .dat files from startrack database

The script currently runs through the standard model
Z = 0.02, 0.002, and 0.0002 files, there is a slight malformatting issue
with the extracted 0.002 dat file where the last column contains spaces which
are used as delimiters within the file.

Each record for one time-step of a bound binary has a format:
    
t dt Ma Mb Ka Kb a e Ra Rb La Lb aspina aspinb mt mttype Lxmt Lx Vsmx Vsmy Vsmz Mzamsa Mzamsb a0 e0 idum iidd evroute

    t is the age of the system (i.e. time since ZAMS) [Myr]
    dt is the length of the actual time-step [Myr]
    a is the semi-major axis of the system [Rsun]
    e is the eccentricity of the system
    Rx is the radius of the star x [Rsun]
    Lx is the luminosity of the star x [Lsun]
    aspinx is a spin of the star x (non-zero only if x is a black hole)
    mt is 1/0 if the mass transfer is on/off
    mttype is the timescale of mass transfer (1 - nuclear, 4 - thermal, 5 - mass transfer from white dwarf)
    Lxmt is a luminosity of an accretion disk during RLOF [erg/s] (non-zero only if mass transfer onto a compact accretor is present)
    Lx is a luminosity coming from wind accretion [erg/s] (non-zero only if a donor has a significant wind and accretor is a compact object)
    Vsmx, Vsmy, Vsmz is a 3D centre of mass velocity [Rsun/day]
    Mzamsa/Mzamsb is the mass on ZAMS for star A/B [Msun]
    a0 is the separation on ZAMS [Rsun]
    e0 is the eccentricity on ZAMS
    idum, iidd - identifiers of a system (i.e. all lines with identical idum and iidd refer to the same binary)
    evroute is a symbolical description of binary evolution (see below)

"""
import tarfile
import glob

import pandas as pd
import populations as pop
import numpy as np
import matplotlib.pyplot as plt
import corner
from scipy.optimize import curve_fit
import os


def fix_csv(infile, outfile):
    """Replaces the first 27 spaces in each line with commas so the file can be
    read properly as a csv"""
    
    with open(infile, 'r') as fin, open(outfile, 'w') as fout:
        for line in fin:
            fout.write(line.replace(' ', ',', 27))
            




def my_filter(df):    
    ns_a = df['Ka'] == 13 # ns
    bh_a = df['Ka'] == 14 # bh
    
    ns_b = df['Kb'] == 13
    bh_b = df['Kb'] == 14
    
    mt = df['mt'] == 1  # mass transfer on
    t_nuc = df['mttype'] == 1   # Nuclear timescale mass transfer
    t_therm = df['mttype'] == 4 # Thermal timescale mass transfer
    
    df = df[mt]
    return df
    
    
def read_dat(dat_file, delimiter=' ', use_dask=False,):
    if use_dask==True:
        return NotImplemented
        
    if use_dask == False:
        saved = pd.DataFrame()
        colnames = ['t', 'dt', 'Ma', 'Mb', 'Ka', 'Kb', 'a', 'e', 'Ra', 'Rb', 'La', 'Lb', 'aspina', 'aspinb', 'mt', 'mttype', 'Lxmt', 'Lx', 'Vsmx', 'Vsmy', 'Vsmz', 'Mzamsa', 'Mzamsb', 'a0', 'e0', 'idum', 'iidd', 'evroute']
        
        for chunk in pd.read_csv(dat_file,
                                 chunksize=chunksize,
                                 index_col=False,
                                 sep=delimiter,
                                 names=colnames,
                                 dtype={'Ka' : 'uint32',
                                        'Kb' : 'uint32',
                                        'mt' : 'uint32',
                                        'mttype': 'uint32',
                                        'evroute': 'category'}):
            
            # chunk = chunk[(chunk['idum'] == idum) & (chunk['iidd'] == iidd)]
            # chunk = chunk[(chunk['iidd'] == iidd)]
            # chunk = chunk.loc[chunk['iidd'].isin(st['iidd'])]
            # chunk = chunk.loc[chunk['iidd'].isin(st['iidd']) & (chunk['idum'].isin(st['idum']))]    
            # chunk = chunk[chunk['Mb'] < 2.5]
            
            chunk = my_filter(chunk)
            print(chunk)
            if len(chunk) > 0:
                saved = saved.append(chunk)
                print(len(saved))
                
        
        return saved

def extract_tar(tar_file, dat_file):
    
    dat_exists = os.path.isfile(dat_file)
    
    if dat_exists:
        print('tbz file already extracted.')
        os.chdir(home_dir)
        return True
    else:
        tf = tarfile.open(tar_file)
        os.chdir(extract_dir)
        tf.extractall()
        os.chdir(home_dir)
    

def get_dataframe(tbz_file, dat_file, csv_file):
        tbz_exists = os.path.isfile(tbz_file)    
        dat_exists = os.path.isfile(dat_file)
        csv_exists = os.path.isfile(csv_file)
        
        if not tbz_exists:
            raise FileNotFoundError('No tbz file found, download from https://universeathome.pl/universe/bhdb.php')
        
        
        if csv_exists:
            df = pd.read_csv(csv_file)
            return df
        
        if dat_exists:
            if 'fix' in dat_file:
                df = read_dat(dat_file, delimiter=',')
            else:
                df = read_dat(dat_file)
            df.to_csv(csv_file)
            os.remove(dat_file)
            return df
        
        else:
            extract_tar(tbz_file, dat_file)
            df = read_dat(dat_file)
            df.to_csv(csv_file)
            os.remove(dat_file)
            return df
        
        
        


    
# =============================================================================
# Read .dat files:
# =============================================================================







if __name__ == "__main__":
    home_dir = os.getcwd()
    extract_dir = 'D:/st/'
    os.makedirs('D:/st', exist_ok=True)
    
    
    
    chunksize = 10 ** 6
    
    
    
    
    # tbz_files = ['../data/external/startrack/std_binary.dat.tbz',
    #              '../data/external/startrack/midZ_binary.dat.tbz',
    #              '../data/external/startrack/lowZ_binary.dat.tbz']
    
    tbz_files = ['../data/external/startrack/std_binary.dat.tbz',
                 extract_dir + 'midZ_binary.dat.tbz',
                 extract_dir + 'lowZ_binary.dat.tbz']
        
    dat_files = [extract_dir+'bhdb180327_ZZ_0.02_Sal_-2.3_SS_1_BHSPIN_1_kick_6_data1.dat',
                extract_dir+'fix_bhdb180327_ZZ_0.002_Sal_-2.3_SS_1_BHSPIN_1_kick_6_data1.dat',
                extract_dir+'bhdb180327_ZZ_0.0002_Sal_-2.3_SS_1_BHSPIN_1_kick_6_data1.dat']
    
    csv_files = ['../data/processed/startrackdb/bhdb180327_ZZ_0.02_Sal_-2.3_SS_1_BHSPIN_1_kick_6_data1_mt_1.csv',
                  '../data/processed/startrackdb/bhdb180327_ZZ_0.002_Sal_-2.3_SS_1_BHSPIN_1_kick_6_data1_mt_1.csv',
                  '../data/processed/startrackdb/bhdb180327_ZZ_0.0002_Sal_-2.3_SS_1_BHSPIN_1_kick_6_data1_mt_1.csv']

    files = {'tbz': tbz_files,
             'dat' : dat_files,
             'csv': csv_files}
    
    N_files = len(tbz_files)
    
    
    
    # =============================================================================
    #     
    # =============================================================================
    
    
    dataframes = {}
    
    for i in range(N_files):
        df = get_dataframe(files['tbz'][i], files['dat'][i], files['csv'][i])
        dataframes[files['csv'][i]] = df
        print(dataframes)
