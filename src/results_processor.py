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
import os
from subprocess import run

from uuid import uuid4

class ulx:
    def __init__(self):
        self.is_ulx = None
        self.is_transient = None
        self.P_cycle_1_ulx = None
        self.P_transient = [None, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
        self.transient_cycle = None
    
    @classmethod
    def persistent_alive_system(cls):
        ulx = cls()
        ulx.P_transient = [None, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
        ulx.P_cycle_1_ulx = 1.0
        return ulx
    
    @classmethod
    def persistent_dead_system(cls):
        ulx = cls()
        ulx.P_transient = [None, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
        ulx.P_cycle_1_ulx = 0.0
        return ulx
    
    @classmethod
    def from_df_classifications_row(cls, Series, period):
        ulx = cls()
        
        ulx.P_cycle_1_ulx = Series.erass_1_ulx_prob
        ulx.P_transient = [None, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
        if period == 'P_wind':
            ulx.P_transient[1] = Series.erass_2_P_wind_transient_prob
            ulx.P_transient[2] = Series.erass_3_P_wind_transient_prob
            ulx.P_transient[3] = Series.erass_4_P_wind_transient_prob
            ulx.P_transient[4] = Series.erass_5_P_wind_transient_prob
            ulx.P_transient[5] = Series.erass_6_P_wind_transient_prob
            ulx.P_transient[6] = Series.erass_7_P_wind_transient_prob
            ulx.P_transient[7] = Series.erass_8_P_wind_transient_prob
        elif period == 'P_sup':
            ulx.P_transient[1] = Series.erass_2_P_sup_transient_prob
            ulx.P_transient[2] = Series.erass_3_P_sup_transient_prob
            ulx.P_transient[3] = Series.erass_4_P_sup_transient_prob
            ulx.P_transient[4] = Series.erass_5_P_sup_transient_prob
            ulx.P_transient[5] = Series.erass_6_P_sup_transient_prob
            ulx.P_transient[6] = Series.erass_7_P_sup_transient_prob
            ulx.P_transient[7] = Series.erass_8_P_sup_transient_prob
        return ulx
    

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
        self.N_transients_cum = np.cumsum(self.N_transients)
        self.N_total_systems = np.cumsum(self.N_new_systems)
        self.N_persistent_ulx_systems = self.N_new_systems[0] - np.cumsum(self.N_old_system_become_transient)

    @classmethod
    def from_run_id(cls, db, run_id, period):
        sql1 = f"""SELECT * FROM ERASS_MC_SAMPLED_SYSTEMS
                  WHERE run_id='{run_id}'"""
        sql2 = f"""SELECT * FROM CLASSIFICATIONS
          WHERE run_id='{run_id}'"""
        sql3 = f"""SELECT * FROM TRANSIENT
          WHERE run_id='{run_id}'"""
          

        conn = sqlite3.connect(db)
        df_sampled_systems = pd.read_sql_query(sql1, conn)
        df_classifications = pd.read_sql_query(sql2, conn)
        df_transient = pd.read_sql_query(sql3, conn)
        conn.close()
        
        class_count = df_classifications['lc_classification'].value_counts()
        
        size = len(df_sampled_systems)
        N_systems_not_simulated = size - len(df_classifications) # These systems are persistent
        try:
            N_alive_systems = class_count[1]
        except KeyError:
            N_alive_systems = 0
        try:
            N_dead_systems = class_count[0]
        except KeyError:
            N_dead_systems = 0
        N_transient_not_sampled = N_alive_systems - len(df_transient) # P_wind & P_sup > 4 year systems
        N_persistent_systems = N_systems_not_simulated + N_alive_systems + N_transient_not_sampled
        

        persistent_systems = [ulx.persistent_alive_system() for i in range(N_persistent_systems)]
        dead_systems = [ulx.persistent_dead_system() for i in range(N_dead_systems)]
        transient_systems = [ulx.from_df_classifications_row(row, period) for i, row in df_transient.iterrows()]
        
        systems = persistent_systems + dead_systems + transient_systems
        Sim_eRASS = cls(systems)
        return Sim_eRASS
        

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

    def collect_results(self):
        res = pd.DataFrame()
        res['erass_cycle'] = np.arange(1, 9)
        res['N_new_systems'] = self.N_new_systems
        res['N_old_system_become_transient'] = self.N_old_system_become_transient
        res['N_observed_ulxs'] = self.N_observed_ulxs
        res['N_delta_obs_ulx'] = self.N_delta_obs_ulx
        res['N_transients'] = self.N_transients
        res['N_transients_cum'] = self.N_transients_cum
        res['N_total_systems'] = self.N_total_systems
        res['N_persistent_ulx_systems'] = self.N_persistent_ulx_systems
        self.res = res
        return res
    
    def sim_write_to_sql(self):
        pass


class ResultsProcessor:
    def __init__(self):
        pass
    
    def set_active_db(self, path):
        self.db = path
        
    def set_parent_population(self, Population):
        self.pop = Population

    def table_create_erass_mc_info(self):
        conn = sqlite3.connect(self.db)
        sql = """CREATE TABLE IF NOT EXISTS ERASS_MC_INFO(
                run_id CHAR,
                size INT,
                bh_ratio REAL,
                dincl_cutoff INT,
                Z CHAR,
                erass_system_period_cutoff INT);"""
        conn.execute(sql)
        conn.close()
    
    def table_create_erass_mc_results(self):
        conn = sqlite3.connect(self.db)
        sql = """CREATE TABLE IF NOT EXISTS ERASS_MC_RESULTS(
                erass_cycle INT,
                N_new_systems INT,
                N_old_system_become_transient INT,
                N_observed_ulxs INT,
                N_delta_obs_ulx INT,
                N_transients INT,
                N_transients_cum INT,
                N_total_systems INT,
                N_persistent_ulx_systems INT,
                period CHAR,
                run_id CHAR);"""
        conn.execute(sql)
        conn.close()
    
    def table_create_erass_mc_sampled_systems(self):
        conn = sqlite3.connect(self.db)
        sql = """CREATE TABLE IF NOT EXISTS ERASS_MC_SAMPLED_SYSTEMS(
                row_id INT,
                theta_half_deg REAL,
                Lx1 REAL,
                P_wind_days REAL,
                P_sup_days REAL,
                inclination INT,
                dincl INT,
                idum_run INT,
                iidd_old INT,
                run_id CHAR);"""
        conn.execute(sql)
        conn.close()
    
    def table_load(self, table_name):
        conn = sqlite3.connect(self.db)
        df = pd.read_sql_query(f"SELECT * from {table_name}", conn)
        conn.close()
        return df
    
    def table_load_transient(self):
        self.df_transient = self.table_load('TRANSIENT')

    def table_load_classifications(self):
        self.df_classifications = self.table_load('CLASSIFICATIONS')
        
    def table_load_erass_mc_info(self):
        self.df_erass_mc_info = self.table_load('ERASS_MC_INFO')
        
    def table_load_erass_mc_results(self):
        self.df_erass_mc_results = self.table_load('ERASS_MC_RESULTS')
        
    def table_load_erass_mc_sampled_systems(self):
        self.df_erass_mc_sampled_systems = self.table_load('ERASS_MC_SAMPLED_SYSTEMS')
    
    def table_load_all(self):
        self.table_load_transient()
        self.table_load_classifications()
        self.table_load_erass_mc_info()
        self.table_load_erass_mc_results()
        self.table_load_erass_mc_sampled_systems()
        
    def table_map_erass_results_info(self):
        info = self.df_erass_mc_info.set_index('run_id')
        self.df_erass_mc_results['bh_ratio'] = self.df_erass_mc_results['run_id'].map(info['bh_ratio'])    
        self.df_erass_mc_results['size'] = self.df_erass_mc_results['run_id'].map(info['size'])    
        self.df_erass_mc_results['dincl_cutoff'] = self.df_erass_mc_results['run_id'].map(info['dincl_cutoff'])
        self.df_erass_mc_results['Z'] = self.df_erass_mc_results['run_id'].map(info['Z'])   
        

    def table_map_classifications_systems(self):
        self.df_classifications['is_bh'] = self.df_classifications['system_row_id'].map(self.pop.df['is_bh'])
        self.df_classifications['P_wind_days'] = self.df_classifications['system_row_id'].map(self.pop.df['P_wind_days'])
        self.df_classifications['a*'] = self.df_classifications['system_row_id'].map(self.pop.df['a*'])
        self.df_classifications['Z'] = self.df_classifications['system_row_id'].map(self.pop.df['Z'])
        
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


    def ADT_MC_bh_ratio_sampler(self, N=100, size=500):

        def create_classification_dict():
            df_classifications_reindexed = self.df_classifications.set_index(['system_row_id', 'system_dincl', 'system_inclination'])
            class_dict = dict(df_classifications_reindexed['lc_classification'])
            return class_dict
        
        try:
            self.class_dict
        except:
            self.class_dict = create_classification_dict()
        
        res = []
 
        for i in range(N):
            print(i)
            for bh_ratio in np.arange(0, 1.05, 0.05):
                
                selected_systems = self.pop.sample_ulxs(bh_ratio, size=size)
                selected_dincls = np.random.randint(0,46, size=size)
                selected_inclinations = np.random.randint(0,91, size=size)
                selected_keys = list(zip(selected_systems, selected_dincls, selected_inclinations))

                res_classications = [self.class_dict.get(key) for key in selected_keys]

                # None systems correspond to opening angles > 45 and are considered alive
                N_alive = res_classications.count(2) + res_classications.count(None)
                N_transient = res_classications.count(1)
                N_dead = res_classications.count(0)
        
                res.append([bh_ratio, N_alive, N_transient, N_dead, i])
        
        df_res = pd.DataFrame(res, columns=['bh_ratio', 'alive', 'transient', 'dead', 'iter'])
        self.df_bh_ratio = df_res
        
        #savepath = '../data/interim/bh_ns_sampling_results/'
        #filename = f'Z = {Z}.csv'
        
        #df_res.to_csv(savepath+filename, index=False)
        # df_res.to_csv(savepath+filename, mode='a', header=False, index=False)


    def ERASS_MC_run(self, N=10000, size=500, dincl_cutoff=46, Z='all', erass_system_period_cutoff=1460):
        print('Running erass MC')
        self.table_create_erass_mc_results()
        self.table_create_erass_mc_info()
        self.table_create_erass_mc_sampled_systems()
        
        # Get columns needed for simulation
        sim_cols = ['theta_half_deg', 'Lx1', 'P_wind_days', 'P_sup_days', 'idum_run', 'iidd_old']
        
        # Delete large dataframes to save memory
        del(self.pop.df)
        
        # Filter population df to only have columns we are interested in to save memory
        self.pop.df_ulx = self.pop.df_ulx[sim_cols]
        
        
        df_info = pd.DataFrame()
        for i in range(N):
            for bh_ratio in [0.0, 0.25, 0.5, 0.75, 1.0]:
                run_id = str(uuid4())
                
                df_info = pd.DataFrame()
                df_info['run_id'] = [run_id]
                df_info['size'] = [size]
                df_info['bh_ratio'] = [bh_ratio]
                df_info['dincl_cutoff'] = [dincl_cutoff]
                df_info['Z'] = [Z]
                df_info['erass_system_period_cutoff'] = [erass_system_period_cutoff]
                
                
                sampled_indexs = self.pop.sample_ulxs(bh_ratio, size=size)
                
                selected_inclinations = np.random.randint(0, 91, size=size)
                selected_dincls = np.random.randint(0, dincl_cutoff, size=size)
                
                df_sampled = self.pop.df_ulx.loc[sampled_indexs]
                df_sampled['Lx1'] = df_sampled['Lx1'] / 1e39
                df_sampled['inclination'] = selected_inclinations
                df_sampled['dincl'] = selected_dincls
                df_sampled['row_id'] = sampled_indexs
                df_sampled['run_id'] = run_id

                conn = sqlite3.connect(self.db)
                df_sampled.to_sql('ERASS_MC_SAMPLED_SYSTEMS', conn, if_exists='append', index=False)
                df_info.to_sql('ERASS_MC_INFO', conn, if_exists='append', index=False)
                conn.close()

                run(["ulxlc/ulxlc_code_v0.1/a.out", str(run_id), str(size), str(erass_system_period_cutoff)])
                
                for period in ['P_wind', 'P_sup']:
                    sim = Sim_eRASS.from_run_id(self.db, run_id, period)
                    sim.run()
                    res = sim.collect_results()
                    res['period'] = period
                    res['run_id'] = run_id
                    conn = sqlite3.connect(self.db)
                    res.to_sql('ERASS_MC_RESULTS', conn, if_exists='append', index=False)
                    conn.close()
    
    def calc_erass_mc_results_stats(self):
        res = self.df_erass_mc_results
        df_stats = res.groupby(['dincl_cutoff', 'Z', 'bh_ratio', 'period', 'erass_cycle']).agg(['min', 'mean', 'max', 'std', 'count'])
        self.df_erass_mc_stats = df_stats
        return df_stats

    def plot_erass_mc_results_hist(self, period='P_wind', by='cycle'):
        params = ['N_new_systems',
                  'N_old_system_become_transient',
                  'N_observed_ulxs',
                  'N_delta_obs_ulx',
                  'N_transients',
                  'N_transients_cum',
                  'N_total_systems']

        unique_bhs = self.df_erass_mc_results['bh_ratio'].unique()
        sub_p = self.df_erass_mc_results[self.df_erass_mc_results['period'] == period]
        
        if by == 'cycle':
            for p in params:
                print(f'Plotting {p}')
                fig, axes  = plt.subplots(2, 4, figsize=(10,8))
                for i, ax in enumerate(axes.flat):
                    cycle = i+1
                    ax.set_title(f'eRASS_{cycle}_{p}')
                    sub = sub_p[sub_p['erass_cycle'] == cycle]
                    for b in unique_bhs:
                        sub1 = sub[sub['bh_ratio'] == b]
                        ax.hist(sub1[p], bins=np.arange(sub1[p].min(), sub1[p].max()+5, 1), label=f'%bh={b}', edgecolor='black', histtype='stepfilled', alpha=0.5)
                        ax.set_xlabel('N')
                        ax.legend()

        if by == 'bh_ratio':
            for p in params:
                print(f'Plotting {p}')
                fig, axes  = plt.subplots(len(unique_bhs), 2, figsize=(10,8))
                for i, b in enumerate(unique_bhs):
                    sub = self.df_erass_mc_results[self.df_erass_mc_results['bh_ratio'] == b]
                    for c in range(1, 9):
                        sub1 = sub[sub['erass_cycle'] == c]
                        sub_wind = sub1[sub1['period'] == 'P_wind']
                        sub_sup = sub1[sub1['period'] == 'P_sup']
                        
                        axes[i][0].hist(sub_wind[p], bins=np.arange(sub_wind[p].min(), sub_wind[p].max()+1, 1), label=f'%cycle={c}', edgecolor='black', histtype='stepfilled', alpha=0.5)
                        axes[i][1].hist(sub_sup[p], bins=np.arange(sub_sup[p].min(), sub_sup[p].max()+1, 1), label=f'%cycle={c}', edgecolor='black', histtype='stepfilled', alpha=0.5)
                        
                        axes[i][0].set_title(f'{p} | P_wind')
                        axes[i][0].set_xlabel('N')
                        axes[i][0].legend()
                        
                        axes[i][1].set_title(f'{p} | P_sup')
                        axes[i][1].set_xlabel('N')
                        axes[i][1].legend()


if __name__ == "__main__":
    
    db_path = 'ulxlc.db'
    
    rp = ResultsProcessor()
    rp.set_active_db(db_path)
    
    # Insert population
    df = populations.startrack_v2_mt_1_all()
    pop = populations.Population(df, 'pop_10k')
    rp.set_parent_population(pop)
    

    # rp.plot_erass_mc_results_hist()
        

    # rp.table_load_classifications()
    # rp.ADT_calc_intermediate()
    # rp.table_map_classifications_systems()
    # rp.ADT_MC_table()
    # rp.ADT_plot_distributions()
    # rp.ADT_plot_distributions(percent=True)
    # rp.ADT_MC_bh_ratio_sampler(N=10)
    
    rp.ERASS_MC_run(erass_system_period_cutoff=5840)
    
    # rp.table_load_all()
    # rp.table_map_erass_results_info()
    
    # rp.plot_erass_mc_results_hist()
    # rp.plot_erass_mc_results_hist(period='P_sup')
    # rp.plot_erass_mc_results_hist(period='P_sup', by='bh_ratio')
    
    
    
    

