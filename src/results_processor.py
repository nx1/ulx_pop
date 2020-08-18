# -*- coding: utf-8 -*-
"""
Created on Wed Aug 12 12:24:13 2020

@author: norma
"""
import glob
from pathlib import Path
import re
import sqlite3
import numpy as np
import pandas as pd

import populations


import matplotlib.pyplot as plt



class ulx:
    def __init__(self):
        self.is_ulx = None
        self.is_transient = None
        self.P_cycle_1_ulx = np.random.random()
        # self.P_transient = [None, 0.33, 0.46, 0.58, 0.69, 0.78, 0.87, 0.99]
        self.P_transient = 0.1*np.random.random(size=8)
        self.P_transient[0] = None
        self.transient_cycle = None
        
    def is_ulx_first_cycle(self):
        self.is_ulx = np.random.random() < self.P_cycle_1_ulx
        return self.is_ulx

    def observe(self, cycle):
        if cycle == 0:
            return [self.is_ulx_first_cycle(), self.is_transient]
        else:
            self.is_transient = np.random.random() < self.P_transient[cycle]
            if self.is_transient:
                self.is_ulx = not self.is_ulx
                if self.transient_cycle == None:
                    self.transient_cycle = cycle
            return [self.is_ulx, self.is_transient]


class Sim_bh_sample:
    def __init__(self, Population):
        self.pop = Population
    
    def ADT_bh_ratio_sampler(self):
        def create_keys(selected_systems, selected_dincls, selected_inclinations):
            a = np.core.defchararray.add(selected_systems.astype(str), '-')
            b = np.core.defchararray.add(selected_dincls.astype(str), '-')
            keys = np.core.defchararray.add(a,b)
            keys = np.core.defchararray.add(keys,selected_inclinations.astype(str))
            return keys
        
        def create_classification_dict():
            # Classication dict
            keys = (self.df_classifications['system_row_id'].astype(str) +
                    '-' + self.df_classifications['system_dincl'].astype(str) +
                    '-' + self.df_classifications['system_inclination'].astype(str)).values
            classications = self.df_classifications['lc_classification'].values
            class_dict = dict(zip(keys,classications))
            return class_dict
        
        class_dict = create_classification_dict()
        
        df_systems = self.df_pop.copy()
        bh_systems = df_systems[df_systems['is_bh'] == 1].index
        ns_systems = df_systems[df_systems['is_bh'] == 0].index
        
        res = []
 
        number_of_repeats = 100
        for i in range(number_of_repeats):
            print(i)
            for bh_ratio in np.arange(0,1.05,0.05):
                ns_ratio = 1 - bh_ratio
        
                bh_weights = [bh_ratio/len(bh_systems)]*len(bh_systems)
                ns_weights = [ns_ratio/len(ns_systems)]*len(ns_systems)
        
                selected_systems = np.random.choice([*bh_systems, *ns_systems], size=500, p=[*bh_weights, *ns_weights])
                selected_dincls = np.random.randint(0,46, size=500)
                selected_inclinations = np.random.randint(0,91, size=500)
                selected_keys = create_keys(selected_systems, selected_dincls, selected_inclinations)
        
                res_classications = [class_dict.get(key) for key in selected_keys]
        
                # None systems correspond to opening angles > 45 and are considered alive
                N_alive = res_classications.count(2) + res_classications.count(None)
                N_transient = res_classications.count(1)
                N_dead = res_classications.count(0)
        
                res.append([bh_ratio, N_alive, N_transient, N_dead])
        
        df_res = pd.DataFrame(res, columns=['bh_ratio', 'alive', 'transient','dead'])
        self.df_bh_ratio = df_res
        
        #savepath = '../data/interim/bh_ns_sampling_results/'
        #filename = f'Z = {Z}.csv'
        
        #df_res.to_csv(savepath+filename, index=False)
        # df_res.to_csv(savepath+filename, mode='a', header=False, index=False)
    
    
    
    
        

class Sim_eRASS:
    def __init__(self, systems):
        self.systems = systems
        self.size = len(self.systems)

        self.N_new_systems 				   = np.array([0,0,0,0,0,0,0,0])
        self.N_old_system_become_transient = np.array([0,0,0,0,0,0,0,0])
    
    def calc_secondary_quantities(self):
        self.N_delta_obs_ulx = self.N_new_systems - self.N_old_system_become_transient
        self.N_observed_ulxs = np.cumsum(self.N_delta_obs_ulx)
        self.N_transients = self.N_new_systems[1:] + self.N_old_system_become_transient[1:]
        self.N_transients = np.insert(self.N_transients, 0, 0)
        self.N_cum_transients = np.cumsum(self.N_transients)
        self.N_total_systems = np.cumsum(self.N_new_systems)

    def run(self):
        for cycle in range(8):
            new_systems = 0
            old_system_become_transient = 0

            for s in self.systems:
                is_ulx, is_transient = s.observe(cycle)
                if is_ulx and is_transient == None:
                    new_systems += 1
                    continue
                elif is_ulx and is_transient and s.transient_cycle == cycle:
                    new_systems += 1
                    continue
                elif not is_ulx and is_transient and s.transient_cycle == cycle:
                    old_system_become_transient += 1
                    continue
                
            self.N_new_systems[cycle] = new_systems
            if cycle>0:
                self.N_old_system_become_transient[cycle] = old_system_become_transient
        self.calc_secondary_quantities()

    def results(self):
        res = pd.DataFrame()
        res['N_new_systems'] = self.N_new_systems
        res['N_old_system_become_transient'] = self.N_old_system_become_transient
        res['N_observed_ulxs'] = self.N_observed_ulxs
        res['N_delta_obs_ulx'] = self.N_delta_obs_ulx
        res['N_transients'] = self.N_transients
        res['N_cum_transients'] = self.N_cum_transients
        res['N_total_systems'] = self.N_total_systems
        self.res = res
        return res



class ResultsProcessor:
    def __init__(self):
        self.get_sim_file_paths()

    def get_sim_file_paths(self):
        self.sim_file_paths = glob.glob('ulxlc/ulxlc_code_v0.1/sim*')
        if len(self.sim_file_paths) == 0:
            raise FileNotFoundError('No simulation files found')
    
    def set_active_db(self, path):
        self.db = path
        
    def set_active_sim_info(self, path):
        self.info_path = path
        
    def set_parent_population(self, df_pop):
        self.df_pop = df_pop

    def read_sim_info(self):
        with open(self.info_path, 'r') as f:
            info = f.readlines()
        return info
    
    def load_table(self, table_name):
        conn = sqlite3.connect(self.db)
        df = pd.read_sql_query(f"SELECT * from {table_name}", conn)
        conn.close()
        return df

    def ADT_load_table(self):
        self.df_classifications = self.load_table('ADT')

    def ADT_map_rows(self):
        df_systems = self.df_pop.copy()
    
        #transient_curves = self.df_classifications[self.df_classifications['lc_classification'] == 1]
        self.df_classifications['is_bh'] = self.df_classifications['system_row_id'].map(df_systems['is_bh'])
        self.df_classifications['P_wind_days'] = self.df_classifications['system_row_id'].map(df_systems['P_wind_days'])
        self.df_classifications['a*'] = self.df_classifications['system_row_id'].map(df_systems['a*'])
        self.df_classifications['Z'] = self.df_classifications['system_row_id'].map(df_systems['Z'])
        
        #self.transient_curves = transient_curves[transient_curves['P_wind_days'] < 4*365]
        
        
    def ADT_calc_intermediate(self):
        df_res = self.df_classifications.copy()
        
        self.df_d = df_res[df_res['lc_classification']==0]
        self.df_t = df_res[df_res['lc_classification']==1]
        self.df_a = df_res[df_res['lc_classification']==2]
        
        systems_per_dincl_bin = df_res['system_dincl'].value_counts().unique()[0]
        systems_per_i_bin     = df_res['system_inclination'].value_counts().unique()[0]
        
        # Number of systems for each dincl
        self.a_dincl_N = self.df_a['system_dincl'].value_counts()
        self.t_dincl_N = self.df_t['system_dincl'].value_counts()
        self.d_dincl_N = self.df_d['system_dincl'].value_counts()
        
        # Number of systems for each inclination
        self.a_i_N = self.df_a['system_inclination'].value_counts().sort_index()
        self.t_i_N = self.df_t['system_inclination'].value_counts().sort_index()
        self.d_i_N = self.df_d['system_inclination'].value_counts().sort_index()
        
        # Percentages
        self.a_dincl_percent = (self.a_dincl_N/systems_per_dincl_bin * 100)
        self.t_dincl_percent = (self.t_dincl_N/systems_per_dincl_bin * 100)
        self.d_dincl_percent = (self.d_dincl_N/systems_per_dincl_bin * 100)
        
        self.a_i_percent = (self.a_i_N/systems_per_i_bin * 100).sort_index()
        self.t_i_percent = (self.t_i_N/systems_per_i_bin * 100)
        self.d_i_percent = (self.d_i_N/systems_per_i_bin * 100).sort_index()


    def ADT_MC_table(self):
        df_res = self.df_classifications.copy()
        
        piv1 = pd.pivot_table(df_res, columns=['is_bh'], index=['Z', 'lc_classification'], aggfunc='count')
        piv1 = piv1.sort_values(by='Z', ascending=False)
        piv1 = piv1[piv1.columns[0:2]]
        
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
        
        self.MC_table = piv1
        return self.MC_table
    
    def ADT_plot_distributions(self, percent=False, savefig=False):
        fig, ax = plt.subplots(ncols=2, nrows=1, figsize=(8,3))
        plt.tight_layout()
        ax[0].set_xlim(0,max(self.df_d['system_dincl']))
        ax[0].set_xlabel(r'Precessional angle $\Delta i$')
        
        ax[0].minorticks_on()
        
        
        ax[1].set_xlim(0,max(self.df_d['system_inclination']))
        ax[1].set_xlabel(r'Inclination $i$')
        ax[1].minorticks_on()
        
        
        
        if percent == False:
            ax[0].set_ylabel(r'Number of systems')
            ax[1].set_ylabel(r'Number of systems')
            
            # dincl vs Number systems
            self.a_dincl_N.plot(label='Alive', linestyle='-', color='black', linewidth=0.8, ax=ax[0])
            self.t_dincl_N.plot(label='Transient', linestyle='--', color='black', linewidth=0.8, ax=ax[0])
            self.d_dincl_N.plot(label='Dead', linestyle='dotted', color='black', linewidth=0.8, ax=ax[0])
            
            # inclination vs Number systems
            self.a_i_N.plot(label='Alive', linestyle='-', color='black', linewidth=0.8, ax=ax[1])
            self.t_i_N.plot(label='Transient', linestyle='--', color='black', linewidth=0.8, ax=ax[1])
            self.d_i_N.sort_index().plot(label='Dead', linestyle='dotted', color='black', linewidth=0.8, ax=ax[1])
            
            ax[1].set_ylabel(r'Number of systems')
            
            if savefig:
                plt.savefig('../reports/figures/dincl_i_classifications.png', dpi=1000)
                plt.savefig('../reports/figures/dincl_i_classifications.eps')
            
        else:
            ax[0].set_ylabel(r'% of systems')
            ax[1].set_ylabel(r'% of systems')
            ax[0].set_ylim(0,100)
            ax[1].set_ylim(0,100)
            
            # dincl vs percentages
            self.a_dincl_percent.plot(label='Alive', linestyle='-', color='black', linewidth=0.8, ax=ax[0])
            self.t_dincl_percent.plot(label='Transient', linestyle='--', color='black', linewidth=0.8, ax=ax[0])
            self.d_dincl_percent.plot(label='Dead', linestyle='dotted', color='black', linewidth=0.8, ax=ax[0])
            
            # inclination vs percentages
            self.a_i_percent.plot(label='Alive', linestyle='-', color='black', linewidth=0.8, ax=ax[1])
            self.t_i_percent.plot(label='Transient', linestyle='--', color='black', linewidth=0.8, ax=ax[1])
            self.d_i_percent.plot(label='Dead', linestyle='dotted', color='black', linewidth=0.8, ax=ax[1])
            
            
            if savefig:
                plt.savefig('../reports/figures/dincl_i_classifications_percent.png', dpi=1000)
                plt.savefig('../reports/figures/dincl_i_classifications_percent.eps')
        ax[0].legend()
        ax[1].legend()
        plt.show()


    def ADT_bh_ratio_sampler(self):
        def create_keys(selected_systems, selected_dincls, selected_inclinations):
            a = np.core.defchararray.add(selected_systems.astype(str), '-')
            b = np.core.defchararray.add(selected_dincls.astype(str), '-')
            keys = np.core.defchararray.add(a,b)
            keys = np.core.defchararray.add(keys,selected_inclinations.astype(str))
            return keys
        
        def create_classification_dict():
            # Classication dict
            keys = (self.df_classifications['system_row_id'].astype(str) +
                    '-' + self.df_classifications['system_dincl'].astype(str) +
                    '-' + self.df_classifications['system_inclination'].astype(str)).values
            classications = self.df_classifications['lc_classification'].values
            class_dict = dict(zip(keys,classications))
            return class_dict
        
        class_dict = create_classification_dict()
        
        df_systems = self.df_pop.copy()
        bh_systems = df_systems[df_systems['is_bh'] == 1].index
        ns_systems = df_systems[df_systems['is_bh'] == 0].index
        
        res = []
 
        number_of_repeats = 100
        for i in range(number_of_repeats):
            print(i)
            for bh_ratio in np.arange(0,1.05,0.05):
                ns_ratio = 1 - bh_ratio
        
                bh_weights = [bh_ratio/len(bh_systems)]*len(bh_systems)
                ns_weights = [ns_ratio/len(ns_systems)]*len(ns_systems)
        
                selected_systems = np.random.choice([*bh_systems, *ns_systems], size=500, p=[*bh_weights, *ns_weights])
                selected_dincls = np.random.randint(0,46, size=500)
                selected_inclinations = np.random.randint(0,91, size=500)
                selected_keys = create_keys(selected_systems, selected_dincls, selected_inclinations)
        
                res_classications = [class_dict.get(key) for key in selected_keys]
        
                # None systems correspond to opening angles > 45 and are considered alive
                N_alive = res_classications.count(2) + res_classications.count(None)
                N_transient = res_classications.count(1)
                N_dead = res_classications.count(0)
        
                res.append([bh_ratio, N_alive, N_transient, N_dead])
        
        df_res = pd.DataFrame(res, columns=['bh_ratio', 'alive', 'transient','dead'])
        self.df_bh_ratio = df_res
        
        #savepath = '../data/interim/bh_ns_sampling_results/'
        #filename = f'Z = {Z}.csv'
        
        #df_res.to_csv(savepath+filename, index=False)
        # df_res.to_csv(savepath+filename, mode='a', header=False, index=False)
            
        
    def ADT_classification_mean_std(results):
        """Calculate for each bh_ratio the mean percentage and 1 sigma error
        for alive, transient AND DEAD systems"""
        res=[]
        states = ['alive', 'transient', 'dead']
        for bh in results['bh_ratio'].unique():
            for state in states:
                cut = results[results['bh_ratio']==bh]
                std_dev = np.std(cut[state])
                mean = np.mean(cut[state])
                res.append([bh, state, mean/500*100, std_dev/500*100])
                
        df_res = pd.DataFrame(res)
        df_res.columns=['bh_ratio', 'classification', 'mean', 'error']
        return df_res
    
    def process_AT_results(results):
        """Calculate for each bh_ratio the mean percentage and 1 sigma error
        for ONLY Alive & Transient results
        """
        res=[]
        states = ['alive', 'transient']
        for bh in results['bh_ratio'].unique():
            for state in states:
                cut = results[results['bh_ratio'] == bh]
                state_avg_number_of_systems = np.mean(cut[state])
                alive_and_dead_avg_number_of_systems = np.mean(cut['alive'] + cut['transient'])
                
                state_std_dev = np.std(cut[state])
    
                res.append([bh,
                            state,
                            state_avg_number_of_systems/alive_and_dead_avg_number_of_systems*100,
                            state_std_dev/alive_and_dead_avg_number_of_systems*100])
        df_res = pd.DataFrame(res)
        df_res.columns=['bh_ratio', 'classification', 'mean', 'error']
        return df_res
    
    def plot_ADT(df_res, Z, axarr, arr_row, arr_col):
        alive_fmt = '-'
        dead_fmt =  '--'
        transient_fmt = ':'
    
        df_res_alive = df_res[df_res['classification'] == 'alive']
        df_res_dead = df_res[df_res['classification'] == 'dead']
        df_res_trans = df_res[df_res['classification'] == 'transient']
        
        axarr[arr_row, arr_col].text(x=5, y=5, s=f'Z = {Z}')
        axarr[arr_row, arr_col].set_xlabel('$\%_{BH}$')
        axarr[arr_row, arr_col].set_ylabel('% of systems')
        axarr[arr_row, arr_col].set_xlim(0,100)
        axarr[arr_row, arr_col].set_ylim(0,100)
        
        axarr[arr_row, arr_col].plot(100*df_res_alive['bh_ratio'], df_res_alive['mean'],
             label='Alive', color='black',
             linestyle=alive_fmt, linewidth=0.8)
        
        axarr[arr_row, arr_col].plot(100*df_res_dead['bh_ratio'], df_res_dead['mean'],
             label='Dead', color='black',
             linestyle=dead_fmt, linewidth=0.8)
        
        axarr[arr_row, arr_col].plot(100*df_res_trans['bh_ratio'], df_res_trans['mean'],
             label='Transient', color='black',
             linestyle=transient_fmt, linewidth=0.8)
        
        
    #    axarr[arr_row, arr_col].errorbar(df_res_alive['BH_ratio'], df_res_alive['mean'],
    #         yerr=df_res_alive['error'], capsize=1.0, color='black',
    #         fmt='none', elinewidth=0.25)
    #    
    #    axarr[arr_row, arr_col].errorbar(df_res_dead['BH_ratio'], df_res_dead['mean'],
    #         yerr=df_res_dead['error'], capsize=1.0, color='black',
    #         fmt='none', elinewidth=0.25)
    #        
    #    axarr[arr_row, arr_col].errorbar(df_res_trans['BH_ratio'], df_res_trans['mean'],
    #         yerr=df_res_trans['error'], capsize=1.0, color='black',
    #         fmt='none', elinewidth=0.25)
        
        axarr[arr_row, arr_col].legend(loc='right', prop={'size': 8})
    #    plt.savefig('../reports/figures/ADT_BHNS_Z_{}.png'.format(Z), format='png', dpi=1000)
    #    plt.savefig('../reports/figures/ADT_BHNS_Z_{}.eps'.format(Z), format='eps')
    #    plt.savefig('../reports/figures/ADT_BHNS_Z_{}.pdf'.format(Z), format='pdf')
    
    def plot_AT(df_res, Z, axarr, arr_row, arr_col):
        alive_color = 'black'
        transient_color = 'grey'
        
        #Code for these figures are in /reports/investigations.ipynb
        PERCENT_ALIVE_EARNSHAW = 0.8271604938271605 * 100
        PERCENT_ALIVE_EARNSHAW_ERROR = 0.12256472421344072 * 100
        
        PERCENT_ALIVE_EARNSHAW_UPPER = PERCENT_ALIVE_EARNSHAW + PERCENT_ALIVE_EARNSHAW_ERROR
        PERCENT_ALIVE_EARNSHAW_LOWER = PERCENT_ALIVE_EARNSHAW - PERCENT_ALIVE_EARNSHAW_ERROR
        
        PERCENT_TRANS_EARNSHAW = 0.1728395061728395 * 100
        PERCENT_TRANS_EARNSHAW_ERROR = 0.03744750536124969 * 100
        PERCENT_TRANS_EARNSHAW_UPPER = PERCENT_TRANS_EARNSHAW + PERCENT_TRANS_EARNSHAW_ERROR
        PERCENT_TRANS_EARNSHAW_LOWER = PERCENT_TRANS_EARNSHAW - PERCENT_TRANS_EARNSHAW_ERROR
        
        #Filtering dfs by state
        df_res_alive = df_res[df_res['classification'] == 'alive']
        df_res_trans = df_res[df_res['classification'] == 'transient']
        
        #Curve Fitting
        def func(x, m, c):
            return m*x + c
        
        def func_inv(y, m, c):
            return (y - c) / m
    
        res_alive = curve_fit(func, df_res_alive['bh_ratio'], df_res_alive['mean'], sigma=df_res_alive['error'])
        res_trans = curve_fit(func, df_res_trans['bh_ratio'], df_res_trans['mean'], sigma=df_res_trans['error'])
        
        m_alive = res_alive[0][0]
        c_alive = res_alive[0][1]
        
        m_trans = res_trans[0][0]
        c_trans = res_trans[0][1]
        
        x = np.arange(0, 1, 0.001)
        y_alive = [func(i, m_alive, c_alive) for i in x]
        y_trans = [func(i, m_trans, c_trans) for i in x]
        
        # Plotting
        axarr[arr_row, arr_col].set_xlim(0,100)
        axarr[arr_row, arr_col].set_ylim(0,100)
        
        
        # Plot errorbars
        axarr[arr_row, arr_col].errorbar(100*df_res_alive['bh_ratio'], df_res_alive['mean'],
             yerr=df_res_alive['error'], label='Alive', capsize=1.0, color=alive_color,
             fmt='x', elinewidth=0.8, markersize=3)
        
        axarr[arr_row, arr_col].errorbar(100*df_res_trans['bh_ratio'], df_res_trans['mean'],
             yerr=df_res_trans['error'], label='Transient', capsize=1.0, color=transient_color,
             fmt='x', elinewidth=0.8, markersize=3)
        

    
    # Plot Line of best fits
    # axarr[arr_row, arr_col].plot(100*x, y_alive, color='black', linewidth=0.8)
    # axarr[arr_row, arr_col].plot(100*x, y_trans, color='grey', linewidth=0.8)
    
    # Plot Labels and annotations
    axarr[arr_row, arr_col].set_ylim(0, 100)
    
    axarr[arr_row, arr_col].set_xlabel('$\%_{BH}$')
    axarr[arr_row, arr_col].set_ylabel('% of systems')

    axarr[arr_row, arr_col].text(x=5, y=5, s=f'Z = {Z}', fontsize=10)
    # axarr[arr_row, arr_col].text(x=0, y=10, s=f'alive bf: y = {m_alive:.2f}x + {c_alive:.2f}', fontsize=6)
    # axarr[arr_row, arr_col].text(x=0, y=15, s=f'trans bf: y = {m_trans:.2f}x + {c_trans:.2f}', fontsize=6)
    
#    axarr[arr_row, arr_col].text(x=0, y=PERCENT_ALIVE_EARNSHAW+2,
#         s=f'Earnshaw % alive: {PERCENT_ALIVE_EARNSHAW:.2f}+-{PERCENT_ALIVE_EARNSHAW_ERROR:.2f}',
#         fontsize=6)
#    axarr[arr_row, arr_col].text(x=0, y=PERCENT_TRANS_EARNSHAW+2,
#         s=f'Earnshaw % trans: {PERCENT_TRANS_EARNSHAW:.2f}+-{PERCENT_TRANS_EARNSHAW_ERROR:.2f}',
#         fontsize=6)
    
    
    # X values interp
    alive_interp_x = func_inv(PERCENT_ALIVE_EARNSHAW, m_alive, c_alive)
    alive_upper_interp_x = func_inv(PERCENT_ALIVE_EARNSHAW_UPPER, m_alive, c_alive)
    alive_lower_interp_x = func_inv(PERCENT_ALIVE_EARNSHAW_LOWER, m_alive, c_alive)
    
    trans_interp_x = func_inv(PERCENT_TRANS_EARNSHAW, m_trans, c_trans)
    trans_upper_interp_x = func_inv(PERCENT_TRANS_EARNSHAW_UPPER, m_trans, c_trans)
    trans_lower_interp_x = func_inv(PERCENT_TRANS_EARNSHAW_LOWER, m_trans, c_trans)
    
    print(f'Z = {Z}')
    print(f'%BH from alive interp upper {alive_upper_interp_x}')
    print(f'%BH from alive interp mean {alive_interp_x}')
    print(f'%BH from alive interp lower {alive_lower_interp_x}')
    print('----')
    print(f'%BH from transient interp upper {trans_upper_interp_x}')
    print(f'%BH from transient interp mean {trans_interp_x}')
    print(f'%BH from transient interp lower {trans_lower_interp_x}')
    

    
#    axarr[arr_row, arr_col].hlines(PERCENT_ALIVE_EARNSHAW, xmin=0, xmax=alive_interp_x,
#         color='r', linewidth=0.8, linestyle='--')
#    axarr[arr_row, arr_col].hlines(PERCENT_ALIVE_EARNSHAW_UPPER, xmin=0, xmax=100*alive_upper_interp_x,
#         color='b', linewidth=0.8, linestyle='--')
#    axarr[arr_row, arr_col].hlines(PERCENT_ALIVE_EARNSHAW_LOWER, xmin=0, xmax=100*alive_lower_interp_x,
#         color='b', linewidth=0.8, linestyle='--', color='b', linewidth=0.8)    
#    axarr[arr_row, arr_col].axhline(PERCENT_ALIVE_EARNSHAW_UPPER, color='b', linewidth=0.6,
#                                    linestyle='--', color='b')
    
    axarr[arr_row, arr_col].axhline(PERCENT_ALIVE_EARNSHAW_UPPER, color='g', linewidth=0.6, label=r'Alive $1\sigma$')
    axarr[arr_row, arr_col].axhline(PERCENT_ALIVE_EARNSHAW_LOWER, color='g', linewidth=0.6)
    axarr[arr_row, arr_col].axhline(PERCENT_TRANS_EARNSHAW_UPPER, color='purple', linewidth=0.6, label=r'Transient $1\sigma$')
    axarr[arr_row, arr_col].axhline(PERCENT_TRANS_EARNSHAW_LOWER, color='purple', linewidth=0.6)
    
    
#    axarr[arr_row, arr_col].hlines(PERCENT_TRANS_EARNSHAW, xmin=0, xmax=trans_interp_x,
#         color='r', linewidth=0.8, linestyle='--')
#    axarr[arr_row, arr_col].hlines(PERCENT_TRANS_EARNSHAW_UPPER, xmin=0, xmax=trans_upper_interp_x,
#         color='b', linewidth=0.8, linestyle='--')
#    axarr[arr_row, arr_col].hlines(PERCENT_TRANS_EARNSHAW_LOWER, xmin=0, xmax=trans_lower_interp_x,
#         color='b', linewidth=0.8, linestyle='--')
    
    """
    axarr[arr_row, arr_col].vlines(100*alive_upper_interp_x, ymin=0, ymax=PERCENT_ALIVE_EARNSHAW_UPPER,
         color='b', linewidth=0.8, linestyle='--')
    axarr[arr_row, arr_col].vlines(100*alive_lower_interp_x, ymin=0, ymax=PERCENT_ALIVE_EARNSHAW_LOWER,
         color='b', linewidth=0.8, linestyle='--')
    
    axarr[arr_row, arr_col].text(x=100*alive_upper_interp_x, y=2,
         s=f'{alive_upper_interp_x:.2f}', c='red')
    axarr[arr_row, arr_col].text(x=100*alive_lower_interp_x, y=2,
         s=f'{alive_lower_interp_x:.2f}', c='red')
    """
    axarr[arr_row, arr_col].legend(loc='right', prop={'size': 8})
        
    def plot_all_ADT():
        import matplotlib
        matplotlib.rcParams['mathtext.fontset'] = 'stix'
        matplotlib.rcParams['font.family'] = 'STIXGeneral'
        fig, axarr = plt.subplots(2,2, figsize=(6,6))
        plt.gcf().subplots_adjust(bottom=0.15)
        
        plot_ADT(df_adt_Z_002, '0.02', axarr, 0, 0)
        plot_ADT(df_adt_Z_0002, '0.002', axarr, 0, 1)
        plot_ADT(df_adt_Z_00002, '0.0002', axarr, 1, 0)
        plot_ADT(df_adt_Z_all, 'all', axarr, 1, 1)
        plt.tight_layout()
        # plt.savefig('../reports/figures/ADT_BHNS_array.png', format='png', dpi=1000)
        # plt.savefig('../reports/figures/ADT_BHNS_array.eps', format='eps')
        # plt.savefig('../reports/figures/ADT_BHNS_array.pdf', format='pdf')
        
    def plot_all_AT():
        import matplotlib
        matplotlib.rcParams['mathtext.fontset'] = 'stix'
        matplotlib.rcParams['font.family'] = 'STIXGeneral'
        fig, axarr = plt.subplots(2,2, figsize=(6,6))
        plt.subplots_adjust(wspace=0.295)
        
        plot_AT(df_at_Z_002, '0.02', axarr, 0, 0)
        plot_AT(df_at_Z_0002, '0.002', axarr, 0, 1)
        plot_AT(df_at_Z_00002, '0.0002', axarr, 1, 0)
        plot_AT(df_at_Z_all, 'all', axarr, 1, 1)
        plt.savefig('../reports/figures/AT_BHNS_array.png', format='png', dpi=1000)
        plt.savefig('../reports/figures/AT_BHNS_array.eps', format='eps')
        plt.savefig('../reports/figures/AT_BHNS_array.pdf', format='pdf')
        
        
        
        
        
        
        
        
    def ERASS_load_table(self):
        self.df_erass = self.load_table('TRANSIENT')
        

        
        
    def ERASS_plot_param(self, param, subs):
        fig , axes  = plt.subplots(2, 4, figsize=(10,8))
        for i, ax in enumerate(axes.flat):
            ax.set_title(f'eRASS_{i}_{param}')
            if param == 'N_delta_obs_ulx':
                ax.hist(subs[i][param], bins=np.arange(-250,250,1))
            else:
                ax.hist(subs[i][param], bins=np.arange(0,500,1))
            ax.set_xlabel('N')
    
    def ERASS_plot_all_params(self):
        subs = []
        for cycle in range(8):
            sub = self.erass_mc_results[self.erass_mc_results.index == cycle]
            subs.append(sub)
        
        params = ['N_new_systems', 'N_old_system_become_transient', 'N_observed_ulxs',
               'N_delta_obs_ulx', 'N_transients', 'N_cum_transients', 'N_total_systems']
        for p in params:
            self.ERASS_plot_param(p, subs)
    
                

if __name__ == "__main__":
    df = populations.startrack_v2_mt_1_all(nrows=10000)
    pop = populations.Population(df, 'pop_10k')
    df_pop = pop.df_ulx_P_wind_l_4_years
    
    rp = ResultsProcessor()
    rp.get_sim_file_paths()
    print(rp.sim_file_paths)

    txt_path = 'ulxlc/ulxlc_code_v0.1/sim_1597740387_info.txt'
    db_path = 'ulxlc/ulxlc_code_v0.1/sim_1597740387.db'
    
    rp.set_active_db(db_path)
    rp.set_active_sim_info(txt_path)
    rp.set_parent_population(df_pop)
    
    rp.ADT_load_table()
    rp.ADT_calc_intermediate()
    rp.ADT_map_rows()
    rp.ADT_MC_table()
    rp.ADT_plot_distributions()
    rp.ADT_plot_distributions(percent=True)
    rp.ADT_bh_ratio_sampler()
    
    rp.ERASS_load_table()
    rp.ERASS_MC_sim(N=1000)
    rp.ERASS_plot_all_params()
    

