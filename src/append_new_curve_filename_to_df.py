#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Feb 13 13:17:50 2020

@author: nk7g14
"""

import glob
import pandas as pd
from io import StringIO
from tqdm import tqdm
from process_new_curves import create_simulation_info_dict
import curvefunc

curve_files = glob.glob('./new_curves/*.txt')


csv_files = glob.glob('../src/new_curve_results/*.csv')
res = {}
for file in csv_files:
    res[file] = pd.read_csv(file)
df_new = pd.concat(res.values())
df_new = df_new.drop(['Unnamed: 0'], axis=1)
df_new['ratio'] = df_new['alive'] / (df_new['alive']+df_new['dead'])
#df_new.loc[np.isnan(df_new['ratio']),'ratio']=0



for curve_file in tqdm(curve_files):
    with open(curve_file, 'r') as file:
        data = file.read()
        simulation_info_list = data.splitlines()[-9:]
        simulation_info = create_simulation_info_dict(simulation_info_list)
        curve = curvefunc.load_curve_file_skip_footer(StringIO(data))
    #print(simulation_info)
    
    
def find_row(df_new, simulation_info):
    