#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
@author: Norman Khan

curve_classifier.py

Runs through output lightcurves from ULXLC and classifies them as alive/dead/transient

Results are saved in the following csv
/data/processed/curve_classications.csv

and contains the following columns:
    'system_id', 'theta', 'dincl', 'inclination',
    'lc_min', 'lc_max', 'N_lim', 'classification'
"""

from tqdm import tqdm
from pathlib import Path
import pickle
import pandas as pd
import numpy as np
import itertools

from auxil import load_systems_dataframe
from curvefunc import load_curve_file

# Load strongly beamed systems
systems_df = load_systems_dataframe(True, True, True)

# values gridsearched over
systems = systems_df.index
dincls = np.arange(0, 46, 1)
inclinations = np.arange(0, 91, 1)

# Create dictionary of system_ids : opening angles
system_id_theta_dict = systems_df['theta_half_deg'].round(2).to_dict()

# Create dictionary of system_ids : Luminosities
system_id_Lx_dict = systems_df['Lx'].to_dict()



# Load 0 inclination normalisation lookup table.
f = open('../data/processed/N_lim_dict.pickle','rb')
N_lim_dict = pickle.load(f)
f.close()

# Path where lightcurves are saved.
p = Path('../data/interim/curves/MC_curves_eta_0.08_ns/gridsearch/')

#Get the minimum and maximum fluxes of each lightcurve
lc_min_max = {}
for path in tqdm(p.glob('*.txt')):
    filename = path.stem
    curve = load_curve_file(path)
    lc_min_max[filename] = (curve['Flux'].min(), curve['Flux'].max())



def classify(N_lim, lc_min, lc_max):
    if lc_min > N_lim:
        classification = 'alive'
    elif lc_max < N_lim:
        classification = 'dead'
    elif (N_lim < lc_max) and (N_lim > lc_min):
        classification = 'transient'
    return classification

results = []
for system_id, dincl, i in tqdm(itertools.product(systems, dincls, inclinations)):
    theta = system_id_theta_dict[system_id]
    key = str(theta) + '-' + str(dincl) + '-' +str(i)
    lc_min, lc_max = lc_min_max[key]
    N_lim = N_lim_dict[str(system_id) + '-' + str(theta) + '-' + str(dincl)]
    classification = classify(N_lim, lc_min, lc_max)
    results.append([system_id, theta, dincl, i, lc_min, lc_max, N_lim, classification])
    
res_df = pd.DataFrame(results, columns=['system_id', 'theta', 'dincl', 'inclination', 'lc_min', 'lc_max', 'N_lim', 'classification'])
res_df.to_csv('../data/processed/curve_classications.csv', index=False)
