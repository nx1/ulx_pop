#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Feb 21 13:00:57 2020

@author: nk7g14
"""

from pathlib import Path
import glob
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

if __name__ == "__main":
    folders = glob.glob('../data/interim/eROSITA_sampling_results/*')
    
    
    results = pd.DataFrame()
    
    for folder in folders:
        settings_file = folder + '/settings.info'
        results_file = folder + '/results.csv'
        systems_file = folder + '/df_systems.csv'
        sampled_systems_file = folder + '/sampled_systems.csv'
        
        
        with open(settings_file,"r") as f:
            lines = f.readlines()
            Z = lines[-1].split(':')[1]
            uuid = lines[0].split(':')[1]
            sampling_iterations = lines[4].split(':')[1].split(' ')[1]
    
        systems = pd.read_csv(systems_file)
    
        res = pd.read_csv(results_file)
        res = res.transpose()
        res = res.drop(['Unnamed: 0'])
        
        
    #    plt.figure()
    #    plt.title(uuid)
    #    for i, row in res.iterrows():
    #        plt.plot(row.index, row.values)
        if len(res) > 200 and int(sampling_iterations) == 100000:
            print('================')
            print(uuid)
            print(f'Z: {Z}')
            print(f'sampling_iterations: {sampling_iterations}')
            print(f'number of results: {len(res)}')
            print(f'number of systems: {len(systems)}')
            print(f'number of BH: {sum(systems.is_bh)}')
            print(f'number of NS: {len(systems) - sum(systems.is_bh)}')
            
    #        print(np.mean(res))
            #plt.yscale('log')        
    #        plt.errorbar(x=np.std(res).index, y=np.mean(res), yerr=np.std(res), label=Z)
            results[Z] = np.mean(res)
            
            
        """
        plt.figure()
        plt.violinplot(dataset = [res[0],
                          res[1],
                          res[2],
                          res[3],
                          res[4],
                          res[5],
                          res[6],
                          res[7]])
        """
    
    plt.figure(figsize=(3.5,3))
    plt.xlabel('eRASS cycle #')
    plt.ylabel('Probability observed as transient ULX')
    
    plt.plot(results.index, results[' all'], label='Z = all')
    plt.plot(results.index, results[' 0.02'], label='Z = 0.02')
    plt.plot(results.index, results[' 0.002'], label='Z = 0.002')
    plt.plot(results.index, results[' 0.0002'], label='Z = 0.0002')
    plt.tight_layout()
    plt.legend()
    plt.savefig('../reports/figures/eRASS_prob.png')
    plt.savefig('../reports/figures/eRASS_prob.eps')
    plt.savefig('../reports/figures/eRASS_prob.pdf')