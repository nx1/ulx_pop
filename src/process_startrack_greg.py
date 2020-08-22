# -*- coding: utf-8 -*-
"""
Created on Mon Jun 29 12:51:09 2020

@author: norma

Sorry that it took so long. You can download the data files
(raw output from PS code) from
https://universeathome.pl/universe/pub/z*_data1.dat where * is the fractional
part of the Metallicity, viz. 02, 002,or 0002.

every line follows the following format (splitted for better readability):
            s.t, s.dt, s.A.M, s.B.M, s.A.K, s.B.K, s.a, s.e,
            s.A.R, s.B.R, s.A.L, s.B.L, s.A.aspin, s.B.aspin,
            s.mt, s.mttype, s.Lxmt, s.Lx, 
            s.A.dMmt,s.A.dMwind,s.B.dMmt,s.B.dMwind,
            s.Vsm[0], s.Vsm[1], s.Vsm[2],
            s.A.Mzams, s.B.Mzams, s.a0, s.e0,
            idum_run, iidd_old, s.evroute

A/B is a star A/B. K is the evolutionary type as described here
https://arxiv.org/pdf/astro-ph/0511811.pdf (Sec. 2.2).
dMmt is what you were looking for, i.e. mass loss from RLOF star [Msun/Myr].
Lxmt is calculated in a different way than in the previous datafile and it
should give similar values as your method before including beaming.
The files can be a bit unwieldy,but you can cut out all the information that
you need. 

Should you have any problem (e.g. files too big to download, unclear format)
I can help you to handle this. 

If you want to discuss your former results in more details,
we can do it also (if my quick response wasn't descriptive enough)

Cheers,
Greg
"""

import os
from pathlib import Path
from urllib.request import urlopen
from shutil import copyfileobj

import pandas as pd
from auxil import download_file


def fix_csv_file(path, fixed_path, number_of_commas):
    """
    Fix .dat file to csv file by inserting a set number of commas instead
    of spaces at each line.

    Parameters
    ----------
    path : path
        Input file
    fixed_path : path
        Output file.        
    number_of_commas : int
        Commas to insert.
    """
    if os.path.isfile(fixed_path):
        print(f'{fixed_path} already exists, not fixing.')
    else:
        with open(path, 'r') as fin, open(fixed_path, 'w') as fout:
            for line in fin:
                fout.write(line.replace(' ', ',', number_of_commas))

def process_fixed_csv_file(fixed_path, processed_csv_path):
    """
    Read, extract and save rows from fixed csv file.
    
    Parameters
    ----------
    fixed_path : path
        Input fixed csv filepath.
    processed_csv_path : path
        Output csv filepath.
    """
    if os.path.isfile(processed_csv_path):
        print(f'{processed_csv_path} already exists, not processing')
    else:
        df = pd.DataFrame()
        for chunk in pd.read_csv(fixed_path, chunksize=chunksize, names=cols):
            chunk = chunk[chunk['mt'] == 1] # Save rows with mass transfer on
            df = df.append(chunk)
            print(len(df))
        df.to_csv(processed_csv_path)

def run():
    for i in range(number_of_downloads):
        print(f'Downloading {data_urls[i]}')
        download_file(data_urls[i], save_paths[i])
        print(f'Fixing file {save_paths[i]}')
        fix_csv_file(save_paths[i], fixed_paths[i], number_of_commas)
        print(f'Reading file {fixed_paths[i]}')
        process_fixed_csv_file(fixed_paths[i], processed_csvs_paths[i])

if __name__ == "__main__":
    EXTERNAL_DATA_PATH = Path('D:/startrack/')
    INTERIM_DATA_PATH = Path('../data/interim/startrack/')
    
    data_urls = ['https://universeathome.pl/universe/pub/z02_data1.dat',
                 'https://universeathome.pl/universe/pub/z002_data1.dat',
                 'https://universeathome.pl/universe/pub/z0002_data1.dat']
    
    filenames            = ['z02_data1.dat', 'z002_data1.dat', 'z0002_data1.dat']
    save_paths           = [EXTERNAL_DATA_PATH / fn for fn in filenames]
    fixed_paths          = [EXTERNAL_DATA_PATH / (fn[:-4]+'_fixed.csv') for fn in filenames]
    processed_csvs_paths = [INTERIM_DATA_PATH / (fn[:-4]+'.csv') for fn in filenames]
    
    number_of_downloads = len(data_urls)
    
    cols = ['t', 'dt', 'M_a', 'M_b', 'K_a', 'K_b', 'a', 'e', 'R_a', 'R_b', 'L_a', 'L_b', 'spin_a', 'spin_b', 'mt', 'mttype', 'Lxmt', 'Lx',  'dMmt_a', 'dMwind_a', 'dMmt_b', 'dMwind_b', 'Vsm0', 'Vsm1', 'Vsm2', 'Mzams_a', 'Mzams_b', 'a0', 'e0', 'idum_run', 'iidd_old', 'evroute']
    number_of_commas = len(cols)-1
    
    chunksize = 1000000
    
    run()




        
    