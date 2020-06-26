# -*- coding: utf-8 -*-
"""
Created on Wed Apr 15 09:27:43 2020

@author: norma
"""

import pandas as pd
import numpy as np

import populations

def table(df_systems):
    piv1 = pd.pivot_table(df_systems, columns=['is_bh'], index=['Z', 'tage'], aggfunc='count')
    piv1 = piv1[piv1.columns[0:2]]
    piv1 = piv1.fillna(0)

    for Z in piv1.index.unique(level=0):
        Z_total = piv1.loc[Z].sum()

    n_ns = piv1[piv1.columns[0]]
    n_bh = piv1[piv1.columns[1]]
    
    tot = n_ns+n_bh
    piv1['%_NS'] = round(n_ns/tot*100, 2) 
    piv1['%_BH'] = round(n_bh/tot*100, 2)
    piv1['Total'] = tot
    
    sum_total = piv1.sum()
    n_ns_tot = sum_total[0]
    n_bh_tot = sum_total[1]
    
    sum_total['%_NS'] = round(n_ns_tot/(n_ns_tot + n_bh_tot)*100, 2) 
    sum_total['%_BH'] = round(n_bh_tot/(n_ns_tot + n_bh_tot)*100, 2) 
    sum_total = sum_total.rename(('SUM', ''))
    
    piv1 = piv1.append(sum_total)
    
    piv1[piv1.columns[0]] = piv1[piv1.columns[0]].astype('int32')
    piv1[piv1.columns[1]] = piv1[piv1.columns[1]].astype('int32')
    piv1[piv1.columns[4]] = piv1[piv1.columns[4]].astype('int32')
    return piv1


if __name__ == "__main__":
    import pdb; pdb.set_trace()
    all_populations = populations.all()    
    for name, df in all_populations.items():
        print(name)
        print(table(df).to_latex())
