# -*- coding: utf-8 -*-
"""
Created on Wed Dec  9 15:47:42 2020

@author: norma
"""
import numpy as np
from tqdm import tqdm


from ulxlc import ULXLC
import populations
from constants import params_default

df = populations.startrack_v2_mt_1_all(nrows=10000)
pop = populations.Population(df)



ulxlc = ULXLC()
ulxlc.set_params(*params_default)


N_mc = 10000

for j in range(N_mc):
    idx = pop.sample_ulxs(0.5, size=500)
    df_sampled = pop.df.loc[idx]
    thetas = df_sampled['theta_half_deg'].values
    for theta in tqdm(thetas):
    # for i, r in df_sampled.iterrows():
    #     theta = r['theta_half_deg']
        if theta > 45:
            continue
        else:
            ulxlc.set_theta(theta)
            
            for inclination in np.arange(0,91,5):
                ulxlc.set_inclination(inclination)
                for dincl in np.arange(0,46,5):
                    ulxlc.set_dincl(dincl)
                    ulxlc.ulxlc_model()
                    #print(max(ulxlc.lc_flux))
                    #print(r)
            
        
