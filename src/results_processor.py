# -*- coding: utf-8 -*-
"""
Created on Wed Aug 12 12:24:13 2020

@author: norma
"""
import sqlite3
import numpy as np
import pandas as pd

import populations

import matplotlib.pyplot as plt
from matplotlib.transforms import Affine2D
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
        self.plot_set_latex_font()
        
        # matplotlib colors
        self.color_dead  = 'C3'
        self.color_trans = 'C1'
        self.color_alive = 'C2'
        
        self.linestyle_dead = 'dotted'
        self.linestyle_trans = '--'
        self.linestyle_alive = '-'
        
        
        #Code for these figures are in ../reports/investigations.ipynb
        self.PERCENT_ALIVE_EARNSHAW = 0.8271604938271605 * 100
        self.PERCENT_ALIVE_EARNSHAW_ERROR = 0.12256472421344072 * 100
        
        self.PERCENT_ALIVE_EARNSHAW_UPPER = self.PERCENT_ALIVE_EARNSHAW + self.PERCENT_ALIVE_EARNSHAW_ERROR
        self.PERCENT_ALIVE_EARNSHAW_LOWER = self.PERCENT_ALIVE_EARNSHAW - self.PERCENT_ALIVE_EARNSHAW_ERROR
        
        self.PERCENT_TRANS_EARNSHAW = 0.1728395061728395 * 100
        self.PERCENT_TRANS_EARNSHAW_ERROR = 0.03744750536124969 * 100
        self.PERCENT_TRANS_EARNSHAW_UPPER = self.PERCENT_TRANS_EARNSHAW + self.PERCENT_TRANS_EARNSHAW_ERROR
        self.PERCENT_TRANS_EARNSHAW_LOWER = self.PERCENT_TRANS_EARNSHAW - self.PERCENT_TRANS_EARNSHAW_ERROR
        
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
                erass_system_period_cutoff INT,
                erass_lmxrb_duty_cycle REAL);"""
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
                lmxrb INT,
                inclination INT,
                dincl INT,
                idum_run INT,
                iidd_old INT,
                run_id CHAR);"""
        conn.execute(sql)
        sql = """CREATE INDEX IF NOT EXISTS idx_sampled_systems_run_id
                 ON ERASS_MC_SAMPLED_SYSTEMS (run_id);"""
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
    
    def table_load_all(self, load_sampled=False):
        self.table_load_transient()
        self.table_load_classifications()
        self.table_load_erass_mc_info()
        self.table_load_erass_mc_results()
        if load_sampled == True:
            self.table_load_erass_mc_sampled_systems()
        
        
    def table_erass_mc_results_map_info(self):
        try:
            self.df_erass_mc_info
        except:
            self.table_load_erass_mc_info()
        try:
            self.df_erass_mc_results
        except:
            self.table_load_erass_mc_results()

        info = self.df_erass_mc_info.set_index('run_id')
        self.df_erass_mc_results['bh_ratio'] = self.df_erass_mc_results['run_id'].map(info['bh_ratio'])    
        self.df_erass_mc_results['size'] = self.df_erass_mc_results['run_id'].map(info['size'])    
        self.df_erass_mc_results['dincl_cutoff'] = self.df_erass_mc_results['run_id'].map(info['dincl_cutoff'])
        self.df_erass_mc_results['Z'] = self.df_erass_mc_results['run_id'].map(info['Z'])  
        self.df_erass_mc_results['erass_system_period_cutoff'] = self.df_erass_mc_results['run_id'].map(info['erass_system_period_cutoff'])
        

    def table_classifications_map_systems(self):
        try:
            self.df_classifications
        except:
            self.table_load_classifications()
            
        self.df_classifications['is_bh'] = self.df_classifications['system_row_id'].map(self.pop.df['is_bh'])
        self.df_classifications['P_wind_days'] = self.df_classifications['system_row_id'].map(self.pop.df['P_wind_days'])
        self.df_classifications['a*'] = self.df_classifications['system_row_id'].map(self.pop.df['a*'])
        self.df_classifications['Z'] = self.df_classifications['system_row_id'].map(self.pop.df['Z'])
        self.df_classifications['lmxrb'] = self.df_classifications['system_row_id'].map(self.pop.df['lmxrb'])
        
        
    def table_classifications_map_info(self):
        info = self.df_erass_mc_info.set_index('run_id')
        self.df_classifications['bh_ratio'] = self.df_erass_mc_results['run_id'].map(info['bh_ratio'])    
        self.df_classifications['size'] = self.df_erass_mc_results['run_id'].map(info['size'])    
        self.df_classifications['dincl_cutoff'] = self.df_erass_mc_results['run_id'].map(info['dincl_cutoff'])
        self.df_classifications['Z'] = self.df_erass_mc_results['run_id'].map(info['Z'])


    def table_classifications_pivot(self, margins=True, split_Z=True):
        try:
            self.df_classifications['Z']
        except KeyError:
            self.table_classifications_map_systems()
        if split_Z==True:
            piv = pd.pivot_table(self.df_classifications,
                                 columns=['is_bh'],
                                 index=['Z', 'lc_classification'],
                                 aggfunc='count',
                                 margins=margins,
                                 margins_name='total').run_id


        if split_Z==False:
            piv = pd.pivot_table(self.df_classifications,
                           columns=['is_bh'],
                           index=['lc_classification'],
                           aggfunc='count',
                           margins=margins,
                           margins_name='total').run_id
            
        if margins==True:
            piv['%_NS'] = piv[0]/piv['total']*100
            piv['%_BH'] = piv[1]/piv['total']*100

        self.df_pivot_lc_classifcations = piv
        return self.df_pivot_lc_classifcations
        
    def table_classifications_split_by_metallicity(self):
        try:
            self.df_classifications['Z']
        except:
            self.table_classifications_map_systems()
            
        self.df_c_02 = self.df_classifications[self.df_classifications['Z']==0.02]
        self.df_c_002 = self.df_classifications[self.df_classifications['Z']==0.002]
        self.df_c_0002 = self.df_classifications[self.df_classifications['Z']==0.0002]
    
    def table_classifications_calc_intermediate(self):
        self.df_d = self.df_classifications[self.df_classifications['lc_classification']==0]
        self.df_t = self.df_classifications[self.df_classifications['lc_classification']==1]
        self.df_a = self.df_classifications[self.df_classifications['lc_classification']==2]
        
        systems_per_dincl_bin = self.df_classifications['system_dincl'].value_counts().unique()[0]
        systems_per_i_bin     = self.df_classifications['system_inclination'].value_counts().unique()[0]
        
        # Number of systems for each dincl
        self.a_dincl_N = self.df_a['system_dincl'].value_counts().sort_index()
        self.t_dincl_N = self.df_t['system_dincl'].value_counts().sort_index()
        self.d_dincl_N = self.df_d['system_dincl'].value_counts().sort_index()
        
        # Number of systems for each inclination
        self.a_i_N = self.df_a['system_inclination'].value_counts().sort_index()
        self.t_i_N = self.df_t['system_inclination'].value_counts().sort_index()
        self.d_i_N = self.df_d['system_inclination'].value_counts().sort_index()
        
        # Percentages
        self.a_dincl_percent = (self.a_dincl_N/systems_per_dincl_bin * 100)
        self.t_dincl_percent = (self.t_dincl_N/systems_per_dincl_bin * 100)
        self.d_dincl_percent = (self.d_dincl_N/systems_per_dincl_bin * 100)
        
        self.a_i_percent = (self.a_i_N/systems_per_i_bin * 100).sort_index()
        self.t_i_percent = (self.t_i_N/systems_per_i_bin * 100).sort_index()
        self.d_i_percent = (self.d_i_N/systems_per_i_bin * 100).sort_index()


    def MC_bh_ratio_classifications_sampler(self, N=10000, size=500, dincl_cutoff=46, Z='all'):
        def create_classification_dict():
            df_classifications_reindexed = self.df_classifications.set_index(['system_row_id', 'system_dincl', 'system_inclination'])
            class_dict = dict(df_classifications_reindexed['lc_classification'])
            return class_dict
        
        try:
            self.class_dict
        except:
            self.class_dict = create_classification_dict()
            
        if Z!='all':
            self.pop.df_ulx[self.pop.df_ulx['Z'] == Z]
        
        res = []
 
        for i in range(N):
            print(i)
            for bh_ratio in np.arange(0, 1.05, 0.05):
                
                selected_systems = self.pop.sample_ulxs(bh_ratio, size=size)
                selected_dincls = np.random.randint(0, dincl_cutoff, size=size)
                selected_inclinations = np.random.randint(0, 91, size=size)
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

    
    def MC_ERASS_simulation_run(self,
                                size=500,
                                dincl_cutoff=46,
                                Z='all',
                                bh_ratio=0.5,
                                erass_system_period_cutoff=1460,
                                erass_lmxrb_duty_cycle=0.1):
        """
        eRASS Monte-Carlo Simulation.

        Parameters
        ----------
        size : int, optional
            Sampled population size.
            The default is 500.
        dincl_cutoff : int, optional
            Maximum precessional angle.
            The default is 46.
        Z : float
            Metallicity to use, one of 'all', 0.02, 0.002, 0.0002.
            The default is 'all'.
        bh_ratio : float
            Ratio of black holes to neutron stars.
            The default is 0.5.
        erass_system_period_cutoff : int, optional
            Maximum period in days after which systems will be considered as persistent.
            The default is 1460.
        erass_lmxrb_duty_cycle : float, optional
            Duty cycle to use for LMXRB systems.
            The default is 0.1.

        Returns
        -------
        None.

        """
        print('Running erass MC...')
        
        # Create SQLite tables
        self.table_create_erass_mc_results()
        self.table_create_erass_mc_info()
        self.table_create_erass_mc_sampled_systems()
        
        # Get columns needed for simulation
        sim_cols = ['theta_half_deg', 'Lx1', 'P_wind_days', 'P_sup_days', 'idum_run', 'iidd_old', 'lmxrb']
        
        # Delete large dataframes to save memory
        # del(self.pop.df)
        
        # Check if we have the correct Z population
        if self.pop.df_ulx_Z_subset != float(Z):
            self.pop.filter_df_ulx_by_Z(Z)
            self.pop.calc_bh_ns_ulx_sub_populations()

        # Filter population df to only have columns we are interested in to save memory
        self.pop.df_ulx = self.pop.df_ulx[sim_cols]
        
        # Create run info dataframe for storing id, size etc..
        df_info = pd.DataFrame()

        # Random 
        run_id = str(uuid4())
        
        df_info = pd.DataFrame()
        df_info['run_id'] = [run_id]
        df_info['size'] = [size]
        df_info['bh_ratio'] = [bh_ratio]
        df_info['dincl_cutoff'] = [dincl_cutoff]
        df_info['Z'] = [Z]
        df_info['erass_system_period_cutoff'] = [erass_system_period_cutoff]
        df_info['erass_lmxrb_duty_cycle'] = [erass_lmxrb_duty_cycle]
        
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

        run(["ulxlc/ulxlc_code_v0.1/a.out", str(run_id), str(size), str(erass_system_period_cutoff), str(erass_lmxrb_duty_cycle)])

        for period in ['P_wind', 'P_sup']:
            sim = Sim_eRASS.from_run_id(self.db, run_id, period)
            sim.run()
            res = sim.collect_results()
            res['period'] = period
            res['run_id'] = run_id
            conn = sqlite3.connect(self.db)
            res.to_sql('ERASS_MC_RESULTS', conn, if_exists='append', index=False)
            conn.close()


    def MC_calc_N_persistent(self, size=500):
        """Calculate the number of systems that had half opening angles > 45
        and thus were not put through the lightcurve simulation"""
        self.df_N_persistent = size - self.df_classifications['run_id'].value_counts()

    
    def MC_get_run_ids(self, group_precession_cuttoff=True):
        try:
            self.df_erass_mc_info
        except AttributeError:
            self.table_load_erass_mc_info()
        
        info = self.df_erass_mc_info.copy()
        info = info.set_index('run_id')
        if group_precession_cuttoff:
            self.dict_MC_run_ids = info.groupby(by=['size', 'bh_ratio', 'dincl_cutoff', 'Z']).groups
        else:
            self.dict_MC_run_ids = info.groupby(by=['size', 'bh_ratio', 'dincl_cutoff', 'Z', 'erass_system_period_cutoff']).groups
    
    def MC_get_run_counts(self):
        count_dict = {}
        for k, v in self.dict_MC_run_ids.items():
            count_dict[k] = len(v)
        return count_dict

    def MC_get_classifications(self, key):
        """
        Parameters
        ----------
        key : tuple from key in self.dict_MC_run_ids
        """
        try:
            self.dict_MC_run_ids
        except:
            self.MC_get_run_ids()
        
        ids = self.dict_MC_run_ids[key]
        df_classifications = self.df_classifications.set_index('run_id')
        # See https://pandas.pydata.org/pandas-docs/stable/user_guide/indexing.html#deprecate-loc-reindex-listlike
        df_classifications = df_classifications.loc[df_classifications.index.intersection(ids).unique()]
        return df_classifications
    
    def MC_get_erass_mc_results(self, key):
        """
        Parameters
        ----------
        key : tuple from key in self.dict_MC_run_ids
        """
        try:
            self.dict_MC_run_ids
        except:
            self.MC_get_run_ids(group_precession_cuttoff=False)
        
        ids = self.dict_MC_run_ids[key]
        df_erass_mc_results = self.df_erass_mc_results.set_index('run_id')
        df_erass_mc_results = df_erass_mc_results.loc[df_erass_mc_results.index.intersection(ids).unique()]
        return df_erass_mc_results
    
    def MC_ADT_result_stats(self):
        try:
            self.df_classifications['Z']
        except:
            self.table_classifications_map_info()
        self.df_classifications_stats = self.df_classifications.groupby(['dincl_cutoff', 'Z', 'bh_ratio']).agg(['min', 'mean', 'max', 'std', 'count'])
        
    def MC_ERASS_result_stats(self):
        try:
            self.df_erass_mc_results['dincl_cutoff']
        except:
            self.table_erass_mc_results_map_info()
            
        self.df_erass_mc_stats = self.df_erass_mc_results.groupby(['dincl_cutoff', 'Z', 'bh_ratio', 'period', 'erass_cycle']).agg(['min', 'mean', 'max', 'std', 'count'])


    def MC_calc_classification_counts(self, key):
        try:
            self.df_N_persistent
        except AttributeError:
            self.MC_calc_N_persistent()
            
        if len(key) != 4:
            print('key length !=4, probably seperated by precessional angle cut.')
            
        df_c = self.MC_get_classifications(key)
        
        df_d = df_c[df_c['lc_classification']==0]
        df_t = df_c[df_c['lc_classification']==1]
        df_a = df_c[df_c['lc_classification']==2]
        
        df_N_d = df_d.index.value_counts().rename('N_dead')
        df_N_t = df_t.index.value_counts().rename('N_trans')
        df_N_a = df_a.index.value_counts().rename('N_alive')
        df_N_p = self.df_N_persistent.rename('N_persistent')

        df_classification_counts = pd.concat([df_N_d, df_N_t, df_N_a, df_N_p], axis=1).dropna()
        
        # Add in those systems with opening angles > 45 that weren't simulated
        df_classification_counts['N_alive_persistent'] = df_classification_counts['N_alive'] + df_classification_counts['N_persistent']
        
        df_classification_counts['frac_alive_visible'] = 100 * df_classification_counts['N_alive_persistent'] /  (df_classification_counts['N_alive_persistent'] + df_classification_counts['N_trans'])
        df_classification_counts['frac_trans_visible'] = 100 * df_classification_counts['N_trans'] / (df_classification_counts['N_alive_persistent'] + df_classification_counts['N_trans'])
        return df_classification_counts
    
    
    def MC_calc_all_classifications_count_stats(self):
        res = pd.DataFrame()
        for k in self.dict_MC_run_ids.keys():
            print(k)
            df = self.MC_calc_classification_counts(k)
            df['size'] = k[0]
            df['bh_ratio'] = k[1]
            df['dincl_cutoff'] = k[2]
            df['Z'] = k[3]
            res = res.append(df)
        gb = res.groupby(['Z', 'dincl_cutoff', 'bh_ratio']).agg(['min', 'mean', 'max', 'std', 'count'])
        gb = gb.drop('size', axis=1)
        self.df_classifications_counts_stats = gb


    def plot_set_latex_font(self):
        import matplotlib
        matplotlib.rcParams['mathtext.fontset'] = 'stix'
        matplotlib.rcParams['font.family'] = 'STIXGeneral'    
    
    def plot_classifications_hist(self, Z, dincl_cut, frac_visible=False, save=False):
        
        fig, ax = plt.subplots(5,1, sharex=True, sharey=True, figsize=(6, 4.5))
        ax[0].set_title(fr'Z = {Z} | $\Delta i_{{max}} =$ {dincl_cut}$^{{\circ}}$')

        
        for i, bh in enumerate([0.0, 0.25, 0.5, 0.75, 1.0]):
            print(i, bh)
            key = (500, bh, dincl_cut, Z)
            df_classification_counts = self.MC_calc_classification_counts(key)
            if frac_visible==False:
                df_classification_counts['N_dead'].hist(bins=np.arange(0,500,5), label='Dead', alpha=0.8, edgecolor='black', linestyle=self.linestyle_dead, histtype='step', ax=ax[i], grid=False)
                df_classification_counts['N_trans'].hist(bins=np.arange(0,500,5), label='Transient', alpha=0.8, edgecolor='black', linestyle=self.linestyle_trans, histtype='step', ax=ax[i], grid=False)
                df_classification_counts['N_alive_persistent'].hist(bins=np.arange(0,500,5), label='Alive', alpha=0.8, edgecolor='black', linestyle=self.linestyle_alive, histtype='step', ax=ax[i], grid=False)
                ax[-1].set_xlabel('Number of observed ULXs')
            else:
                ax[i].axvspan(self.PERCENT_ALIVE_EARNSHAW_LOWER, self.PERCENT_ALIVE_EARNSHAW_UPPER, alpha=0.3, color='grey')
                ax[i].axvspan(self.PERCENT_TRANS_EARNSHAW_LOWER, self.PERCENT_TRANS_EARNSHAW_UPPER, alpha=0.7, color='grey')
                
                df_classification_counts['frac_alive_visible'].hist(bins=np.arange(0,100,1), label='% Alive', alpha=0.8, edgecolor='black', linestyle=self.linestyle_alive, histtype='step', ax=ax[i], grid=False)
                df_classification_counts['frac_trans_visible'].hist(bins=np.arange(0,100,1), label='% Transient', alpha=0.8, edgecolor='black', linestyle=self.linestyle_trans, histtype='step', ax=ax[i], grid=False)
                ax[-1].set_xlabel('% of observed ULXs')
                
            ax[i].text(x=350, y=ax[i].get_ylim()[1]*(1/2), s=r'$\%_{{BH}} = $ {}'.format((bh)*100), fontsize=9)    
            ax[i].legend(fontsize=7, loc='right')
        
        if frac_visible:
            ax[0].text(x=(self.PERCENT_ALIVE_EARNSHAW_LOWER)+0.25, y=ax[0].get_ylim()[1]+5, s='Alive $1 \sigma$ XMM', fontsize=7)
            ax[0].text(x=(self.PERCENT_TRANS_EARNSHAW_LOWER)+0.25, y=ax[0].get_ylim()[1]+5, s='Transient $1 \sigma$ XMM', fontsize=7)                


        
        plt.tight_layout()
        plt.subplots_adjust(hspace=0.0)
        
        if save:
            if frac_visible==False:
                plt.savefig(f'../reports/figures/ADT_{Z}_{dincl_cut}.png', dpi=500)
                plt.savefig(f'../reports/figures/ADT_{Z}_{dincl_cut}.eps')
                plt.savefig(f'../reports/figures/ADT_{Z}_{dincl_cut}.pdf')
            else:
                plt.savefig(f'../reports/figures/ADT_frac_{Z}_{dincl_cut}.png', dpi=500)
                plt.savefig(f'../reports/figures/ADT_frac_{Z}_{dincl_cut}.eps')
                plt.savefig(f'../reports/figures/ADT_frac_{Z}_{dincl_cut}.pdf')
        
        
        
    def plot_classifications_dincl_i(self, percent=False, savefig=False):
        try:
            self.df_d
        except:
            self.table_classifications_calc_intermediate()
        
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
            self.a_dincl_N.plot(label='Alive', linestyle=self.linestyle_alive, color='black', linewidth=0.8, ax=ax[0])
            self.t_dincl_N.plot(label='Transient', linestyle=self.linestyle_trans, color='black', linewidth=0.8, ax=ax[0])
            self.d_dincl_N.plot(label='Dead', linestyle=self.linestyle_dead, color='black', linewidth=0.8, ax=ax[0])
            
            # inclination vs Number systems
            self.a_i_N.plot(label='Alive', linestyle=self.linestyle_alive, color='black', linewidth=0.8, ax=ax[1])
            self.t_i_N.plot(label='Transient', linestyle=self.linestyle_trans, color='black', linewidth=0.8, ax=ax[1])
            self.d_i_N.sort_index().plot(label='Dead', linestyle=self.linestyle_dead, color='black', linewidth=0.8, ax=ax[1])
            
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
            self.a_dincl_percent.plot(label='Alive', linestyle=self.linestyle_alive, color='black', linewidth=0.8, ax=ax[0])
            self.t_dincl_percent.plot(label='Transient', linestyle=self.linestyle_trans, color='black', linewidth=0.8, ax=ax[0])
            self.d_dincl_percent.plot(label='Dead', linestyle=self.linestyle_dead, color='black', linewidth=0.8, ax=ax[0])
            
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


    def plot_erass_mc_results_hist(self, key, by='cycle', save=False):
        try:
            self.df_erass_mc_results['Z']
        except:
            self.table_erass_mc_results_map_info()
        
        params = ['N_new_systems',
                  'N_old_system_become_transient',
                  'N_observed_ulxs',
                  'N_delta_obs_ulx',
                  'N_transients',
                  'N_transients_cum',
                  'N_total_systems']
        
        res = self.MC_get_erass_mc_results(key)
        
        if by == 'cycle':
            for p in params:
                print(f'Plotting {p}')
                fig, axes  = plt.subplots(2, 4, figsize=(12,10))
                for i, ax in enumerate(axes.flat):
                    cycle = i+1
                    ax.set_title(f'eRASS_{cycle}_{p}')
                    sub = res[res['erass_cycle'] == cycle]
                    for b in unique_bhs:
                        sub1 = sub[sub['bh_ratio'] == b]
                        ax.hist(sub1[p], bins=np.arange(sub1[p].min(), sub1[p].max()+5, 1), label=f'%bh={b}', edgecolor='black', histtype='stepfilled', alpha=0.5)
                        ax.set_xlabel('N')
                        ax.legend()

        if by == 'bh_ratio':
            for p in params:
                print(f'Plotting {p}')
                fig, axes  = plt.subplots(len(unique_bhs), 2, figsize=(12,10), sharex=True)
                for i, b in enumerate(unique_bhs):
                    sub = res[res['bh_ratio'] == b]
                    for c in range(1, 9):
                        sub1 = sub[sub['erass_cycle'] == c]
                        sub_wind = sub1[sub1['period'] == 'P_wind']
                        sub_sup = sub1[sub1['period'] == 'P_sup']
                        
                        axes[i][0].tick_params(axis='both', labelsize=8, labelbottom=True)
                        axes[i][1].tick_params(axis='both', labelsize=8, labelbottom=True)
                        
                        axes[i][0].hist(sub_wind[p], bins=np.arange(sub_wind[p].min(), sub_wind[p].max()+1, 1), label=f'cycle={c}', edgecolor='black', histtype='stepfilled', alpha=0.5)
                        axes[i][1].hist(sub_sup[p], bins=np.arange(sub_sup[p].min(), sub_sup[p].max()+1, 1), label=f'cycle={c}', edgecolor='black', histtype='stepfilled', alpha=0.5)
                        
                        axes[i][0].set_title(f'{p} | bh_ratio = {b} |P_wind')
                        # axes[i][0].set_xlabel('N')
                        axes[i][0].legend(fontsize=7)
                        
                        axes[i][1].set_title(f'{p} | bh_ratio = {b} |P_sup')
                        # axes[i][1].set_xlabel('N')
                        axes[i][1].legend(fontsize=7)
                plt.tight_layout()
                if save==True:
                    plt.savefig(f'../reports/figures/e_{p}_{by}_{erass_system_period_cutoff}_{dincl_cutoff}_{Z}.png',dpi=500)
                    plt.savefig(f'../reports/figures/e_{p}_{by}_{erass_system_period_cutoff}_{dincl_cutoff}_{Z}.eps')
                    plt.savefig(f'../reports/figures/e_{p}_{by}_{erass_system_period_cutoff}_{dincl_cutoff}_{Z}.pdf')

    def plot_bar_classifications_Z(self):
        piv = self.table_classifications_pivot(margins=False)
        piv.plot(kind='bar')
        
    def plot_erass_transients(self, erass_system_period_cutoff=1460, Z='all', save=False):
    
        cycles = ['1', '2', '3', '4', '5', '6', '7', '8']
        dincls = [46, 20]
        bh_percents = [0, 0.25, 0.5, 0.75, 1.0]
        
        nbars = 8
        spacing = 0.1
        linewidth = 1.0
        
        fig, ax = plt.subplots(ncols=2, nrows=len(dincls), figsize=(6,6), sharey='row')
        
        ax[-1][0].set_xlabel('eRASS cycle')
        ax[-1][1].set_xlabel('eRASS cycle')
        
        clist = ["#ff595e", "#ffca3a", "#8ac926", "#1982c4", "#6a4c93"]
        
        for i, dincl in enumerate(dincls):
            for j, bh in enumerate(bh_percents):
                key = (500, bh, dincl, Z, erass_system_period_cutoff)
                df = self.MC_get_erass_mc_results(key)
                
                sub_wind = df[df['period'] == 'P_wind']
                sub_sup = df[df['period'] == 'P_sup']
                
                
                stats_wind = sub_wind.groupby(['erass_cycle']).agg(['mean', 'std'])['N_transients']
                stats_sup = sub_sup.groupby(['erass_cycle']).agg(['mean', 'std'])['N_transients']
                
                bh_label = f'$\%_{{BH}}$ = {bh*100}'            
                
                trans1 = Affine2D().translate(-spacing*(nbars)/4 + spacing*j, 0.0) + ax[0][i].transData
                trans2 = Affine2D().translate(-spacing*(nbars)/4 + spacing*j, 0.0) + ax[1][i].transData
        
                ax[0][i].errorbar(cycles, stats_wind['mean'], yerr=stats_wind['std'], linestyle="none",
                                  linewidth=linewidth, capsize=1.0, transform=trans1, label=bh_label, c=clist[j])
                
                ax[1][i].errorbar(cycles, stats_sup['mean'], yerr=stats_sup['std'], linestyle="none",
                      linewidth=linewidth, capsize=1.0, transform=trans2, label=bh_label, c=clist[j])
        
    
                ax[0][i].set_title(f'$\Delta i_{{max}} = {dincl}^{{\circ}}$ | $P_{{wind}}$ | Z = {Z}')
                ax[1][i].set_title(f'$\Delta i_{{max}} = {dincl}^{{\circ}}$ | $P_{{sup}}$ | Z = {Z}')
        
                
                ax[i][0].grid(axis='y')
                ax[i][1].grid(axis='y')
                #ax[1][i].grid(axis='y')
                #ax[0][i].tick_params(labelrotation=90)
                #ax[1][i].tick_params(axis='x', labelrotation=90)
                
        ax[0][0].set_ylabel('N Transients')
        ax[1][0].set_ylabel('N Transients')
        
        ax[0][0].set_ylim(0,35)
        ax[1][0].set_ylim(0,95)
        
        ax[0][0].legend()
        #ax[0][0].legend()
        
        plt.subplots_adjust(wspace=0)
        
        plt.tight_layout()
        if save:
            plt.savefig(f'../reports/figures/erass_N_transients_Z={Z}.eps', bbox_inches='tight')
            plt.savefig(f'../reports/figures/erass_N_transients_Z={Z}.png', bbox_inches='tight', dpi=1000)



if __name__ == "__main__":
    # Load population
    df = populations.startrack_v2_mt_1_all()
    pop = populations.Population(df)
    
    # Create sim object    
    rp = ResultsProcessor()
    
    # Set Database
    db_path = 'ulxlc.db'
    rp.set_active_db(db_path)
    
    # Insert population
    rp.set_parent_population(pop)

    # # Run Simulations
    # for Z in ['0.02', '0.002', '0.0002']:
    #     pop = populations.Population(df)
    #     rp.set_parent_population(pop)
    #     for i in range(10):
    #         for dincl in [21, 46]:
    #             for dc in [0.1, 0.2, 0.3, 1.0]:
    #                 for bh_ratio in [0.0, 0.25, 0.5, 0.75, 1.0]:
    #                     print(i, bh_ratio, dincl, Z, dc)
    #                     # import pdb; pdb.set_trace()
    #                     rp.MC_ERASS_simulation_run(size=500, dincl_cutoff=dincl, Z=Z, erass_lmxrb_duty_cycle=dc, erass_system_period_cutoff=999999)


    # rp.MC_bh_ratio_classifications_sampler(N=10)

    # Load tables
    # rp.table_load_all()
    rp.table_load_classifications()
    # rp.table_load_transient()
    rp.table_load_erass_mc_info()
    # rp.table_load_erass_mc_results()
    
    # Process tables
    # rp.table_erass_mc_results_map_info()
    # rp.table_classifications_map_systems()
    # rp.table_classifications_pivot()
    
    # rp.table_classifications_calc_intermediate()
    # rp.table_classifications_map_systems()
    # rp.table_classifications_pivot()
    
    # MC related functions
    # rp.MC_get_run_ids(group_precession_cuttoff=True)
    # rp.MC_get_run_counts()

    # Plot results
    
    # rp.plot_set_latex_font()
    # rp.plot_classifications_hist(Z='all', dincl_cut=20, frac_visible=True, save=True)
    # rp.plot_erass_transients()
    
    
    # for key in rp.dict_MC_run_ids.keys():
    #     print(key)

    
    # rp.plot_erass_mc_results_hist(key)
    # rp.plot_classifications_hist('0.002', 46, frac_visible=True, save=False)
    # rp.plot_classifications_dincl_i()
    # rp.plot_classifications_dincl_i(percent=True)   
    # rp.plot_classifications_hist('all', 20, frac_visible=True, save=True)

    # rp.plot_erass_mc_results_hist(period='P_sup',
    #                               by='bh_ratio',
    #                               erass_system_period_cutoff=1460,
    #                               dincl_cutoff=46,
    #                               Z='0.002',
    #                               save=True)

    # rp.MC_calc_N_persistent()
    