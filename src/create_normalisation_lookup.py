#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
@author: Norman Khan

create_normalisation_lookup.py

Run through 0 inclination lightcurves and calculate normalisation ULX
limits at 1E39 ergs equivilent, store the output in a .pickle file in 

data/processed/N_lim_dict.pickle

final dict is as follows N_lim_dict = {system_id-dincl-inclination: N_lim}

For a given:
    system_id (Lx)
    precessional angle
    inclinationation
    you get back a paticular value of N_lim
"""


from pathlib import Path
import pickle
from tqdm import tqdm
import numpy as np

from auxil import load_systems_dataframe
from curvefunc import load_curve_file


def calc_N_lim(c0_max, Lx):
    """
    Calculate where the 1E39 equivilent flux would be on a given lightcurve
    
    Parameters
    ----------
    c0_max : float
        Maximum lightcurve value at 0 inclination.
    Lx : float
        Beamed Luminosity.

    Returns
    -------
    N_lim : float
        Normalised 1E39 erg/s limit
    """
    c = Lx / lc_max #Curve Normalisation constant
    N_lim = 1E39 / c
    return N_lim


systems_df = load_systems_dataframe(True, True, True)
Lx_dict = systems_df['Lx'].to_dict()

# Round to match curve filenames
thetas_rounded = systems_df['theta_half_deg'].round(2)
id_theta_dict = thetas_rounded.to_dict()


p = Path('../data/interim/curves/MC_curves_eta_0.08_ns/zero_inclination_curves/')

lc_min_max = {}

for path in tqdm(p.glob('*.txt')):
    filename = path.stem
    curve = load_curve_file(path)
    lc_min_max[filename] = (curve['Flux'].min(), curve['Flux'].max())

N_lim_dict = {}
dincls = np.arange(0, 46, 1)
for system_id, theta in id_theta_dict.items():
    for d in dincls:
        key = str(theta)+'-'+str(d)+'-'+'0'
        
        Lx = Lx_dict[system_id]
        lc_max = lc_min_max[key][1]
        N_lim = calc_N_lim(lc_max, Lx)
        N_lim_dict[str(system_id) + '-' + str(theta) + '-' + str(d)] = N_lim


f = open('../data/processed/N_lim_dict.pickle','wb')
pickle.dump(N_lim_dict,f)
f.close()