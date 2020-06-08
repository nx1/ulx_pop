# -*- coding: utf-8 -*-
"""
Created on Tue May  5 16:09:06 2020

@author: norma
"""


import glob
from pathlib import Path
import numpy as np
import pandas as pd
from tqdm import tqdm
import sqlite3

from create_erass_curves import create_parent_population
from curvefunc import load_curve_file

def erass_sample_curve(curve, curve_period, N_lim, number_of_repeats):
    """Samples a lightcurve of a given period as would be done by eRASS
    (8 cycles of 6 months)
    parameters:
    -----------
    curve: lightcurve dataframe
    curve_period: period of the lightcurve
    sampling interval: how often to sample the lightcurve after selecting
    a random position somwhere in the first period.
    number_of_repeats: number of MC iterations
    """
    def truth_table_processor(truth):
        """Runs through the 8*N is ulx truth table
        to create a new dataframe with each column corresponding to an eRASS
        cycle and weather or not the source was observed as a transient ulx or not.
        thanks to based James for the help
        """
        observed_as_transient = []
        for row in truth:
            transient=8*[True]
            i=0 # Not observed as transient
            while row[i] == row[0] and i!=8:
                transient[i]=False
                i+=1 #observed as transient
            observed_as_transient.append(transient)
        
        df = pd.DataFrame(observed_as_transient)
        return df.mean()
    
    time_arr = np.array(curve['Time'])
    flux_arr = np.array(curve['Flux'])
    
    truth = []
    for n in range(number_of_repeats):
        start_time = np.random.uniform(0, curve_period)
        sample_times = np.array([(start_time + i*30*6)%curve_period for i in range(0,9)])
        sample_indexes = [np.argmin(np.abs(time_arr - t)) for t in sample_times]
        
        fluxes = flux_arr[sample_indexes]
        
        truth.append(list(np.where(fluxes > N_lim, True, False)))
    sample_results = truth_table_processor(truth)
    return sample_results

def filename_from_row(row):
    filename = str(row['theta']) + '-' + str(row['inclination']) + '-' + str(row['dincl'])
    return filename

def main():
    parent_population = create_parent_population()
    
    path = Path('./eRASS_sims/')
    # files = glob.glob()
    
    conn = sqlite3.Connection('erass.db')
    c = conn.cursor()
    c.execute("""CREATE TABLE IF NOT EXISTS sampResults (
    	curve_id int,
    	eRASS1 REAL,
        eRASS2 REAL,
        eRASS3 REAL,
        eRASS4 REAL,
        eRASS5 REAL,
        eRASS6 REAL,
        eRASS7 REAL,
        eRASS8 REAL
        )""")
    
    for index, row in tqdm(parent_population.iterrows()):
        # print(row)
        try:
            filename = filename_from_row(row)
            curve = load_curve_file('./eRASS_sims/' + filename +'.txt')
            N_lim = row['N_lim']
            res = erass_sample_curve(curve, 50, N_lim, 10000)
            #filename = f.stem
            c.execute("""INSERT INTO sampResults(
                curve_id, eRASS1, eRASS2, eRASS3, eRASS4, eRASS5, eRASS6, eRASS7, eRASS8)
                VALUES(?,?,?,?,?,?,?,?,?)""", (index, *list(res)))
            conn.commit()
        except:
            pass

if __name__ == "__main__":
    main()