# -*- coding: utf-8 -*-
"""
Created on Mon May  4 11:34:57 2020

@author: norma
"""

from os import path
import numpy as np
import pandas as pd
from tqdm import tqdm

from auxil import load_systems_dataframe
from batchfunc import run_ulxlc
from curvefunc import load_curve_file, scale_light_curve_period


def create_parent_population():    
    systems_df = load_systems_dataframe(True,True,True)
    curve_classifications = pd.read_csv('../data/processed/curve_classifications.csv')

    transient_curves = curve_classifications[curve_classifications['classification'] == 'transient']
    transient_curves['is_bh'] = transient_curves['system_id'].map(systems_df['is_bh'])
    transient_curves['P_wind_days'] = transient_curves['system_id'].map(systems_df['P_wind_days'])
    transient_curves['a*'] = transient_curves['system_id'].map(systems_df['a*'])
    transient_curves['Z'] = transient_curves['system_id'].map(systems_df['Z'])
    
    
    transient_curves = transient_curves[transient_curves['P_wind_days'] < 4*365]
    return transient_curves



def get_curve(row):
    def get_ulxlc_parameters(row):
        """Create parameters dictionary from df_a row"""
        parameters = {'period': 50,
                        'phase': 0,
                        'theta': row['theta'],
                        'inclination': row['inclination'],
                        'dincl': row['dincl'],
                        'beta': 0.2,
                        'dopulse': 1,
                        'norm': 1}
        return parameters

    parameters = get_ulxlc_parameters(row)
    filename = str(row['theta']) + '-' + str(row['inclination']) +'-' + str(row['dincl'])
    xcm_file = 'eRASS_sims/' + filename + '.xcm'
    lc_file = 'eRASS_sims/' + filename + '.txt'
    
    if not path.exists(lc_file):
        run_ulxlc(xcm_file, parameters, lc_file)


def main():
    parent_pop = create_parent_population()
    for index, row in tqdm(parent_pop.iterrows()):
        get_curve(row)

if __name__ == "__main__":
    main()
        
