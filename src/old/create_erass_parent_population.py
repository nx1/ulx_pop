# -*- coding: utf-8 -*-
"""
Created on Mon May 11 22:27:06 2020

@author: norma
"""

from os import path
import numpy as np
import pandas as pd
from tqdm import tqdm


import populations
from ulxlc import run_ulxlc
from curvefunc import load_curve_file, scale_light_curve_period


def main():    
    systems_df = populations.ulx_beamed_l_45()
    curve_classifications = pd.read_csv('../data/processed/curve_classifications.csv')

    transient_curves = curve_classifications[curve_classifications['classification'] == 'transient']
    transient_curves['is_bh'] = transient_curves['system_id'].map(systems_df['is_bh'])
    transient_curves['P_wind_days'] = transient_curves['system_id'].map(systems_df['P_wind_days'])
    transient_curves['a*'] = transient_curves['system_id'].map(systems_df['a*'])
    transient_curves['Z'] = transient_curves['system_id'].map(systems_df['Z'])
    
    transient_curves = transient_curves[transient_curves['P_wind_days'] < 4*365]
    
    return transient_curves

if __name__ == "__main__":
    parent_population = main()