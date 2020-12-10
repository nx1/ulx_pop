# -*- coding: utf-8 -*-
"""
Created on Wed Aug 12 12:24:13 2020

@author: norma
"""
import logging
from subprocess import run
from uuid import uuid4

import sqlite3
import numpy as np
import pandas as pd
import matplotlib
import matplotlib.pyplot as plt
from matplotlib.transforms import Affine2D
from scipy.stats import norm
from scipy.optimize import curve_fit

import populations

class ResultsProcessor:
    def __init__(self):
        pass
    
    def set_active_db(self, path):
        logging.debug('Setting active db to %s', path)
        self.db = path
        
    def set_parent_population(self, Population):
        logging.debug('Setting parent population')
        self.pop = Population

    def table_create_erass_mc_info(self):
        logging.debug('Creating table ERASS_MC_INFO')
        conn = sqlite3.connect(self.db)
        sql = """CREATE TABLE IF NOT EXISTS ERASS_MC_INFO(
                run_id CHAR,
                size INT,
                bh_ratio REAL,
                dincl_cutoff INT,
                Z CHAR,
                erass_system_period_cutoff INT);"""
        conn.execute(sql)
        sql = """CREATE INDEX IF NOT EXISTS idx_mc_info_run_id
         ON ERASS_MC_INFO (run_id);"""
        conn.execute(sql)
        conn.close()

    def table_create_classifications(self):
        logging.debug('Creating table CLASSIFICATIONS')
        conn = sqlite3.connect(self.db)
        sql = """CREATE TABLE IF NOT EXISTS CLASSIFICATIONS(
                  system_row_id INT,
                  system_theta REAL,
                  system_dincl INT,
                  system_inclination INT,
                  lc_min_flux REAL,
                  lc_max_flux REAL,
                  lc_ulx_lim REAL,
                  lc_classification INT,
                  run_id CHAR);"""
        conn.execute(sql)
        sql = """CREATE INDEX IF NOT EXISTS idx_classifications_run_id
         ON CLASSIFICATIONS (run_id);"""
        conn.execute(sql)
        conn.close()
        
    def table_create_transient(self):
        logging.debug('Creating table TRANSIENT')
        conn = sqlite3.connect(self.db)
        sql = """CREATE TABLE IF NOT EXISTS TRANSIENT( 
                  system_row_id INT, 
                  system_dincl INT, 
                  system_inclination INT, 
                  erass_1_ulx_prob REAL, 
                  erass_2_P_wind_transient_prob REAL, 
                  erass_3_P_wind_transient_prob REAL, 
                  erass_4_P_wind_transient_prob REAL, 
                  erass_5_P_wind_transient_prob REAL, 
                  erass_6_P_wind_transient_prob REAL, 
                  erass_7_P_wind_transient_prob REAL, 
                  erass_8_P_wind_transient_prob REAL, 
                  erass_2_P_sup_transient_prob REAL, 
                  erass_3_P_sup_transient_prob REAL, 
                  erass_4_P_sup_transient_prob REAL, 
                  erass_5_P_sup_transient_prob REAL, 
                  erass_6_P_sup_transient_prob REAL, 
                  erass_7_P_sup_transient_prob REAL, 
                  erass_8_P_sup_transient_prob REAL, 
                  erass_P_wind_persistent_prob REAL, 
                  erass_P_sup_persistent_prob REAL, 
                  run_id CHAR);"""
        conn.execute(sql)
        sql = """CREATE INDEX IF NOT EXISTS idx_transient_run_id
         ON TRANSIENT (run_id);"""
        conn.execute(sql)
    
    def table_create_erass_mc_results(self):
        logging.debug('Creating table ERASS_MC_RESULTS')
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
        sql = """CREATE INDEX IF NOT EXISTS idx_mc_rusults_run_id
         ON ERASS_MC_RESULTS (run_id);"""
        conn.execute(sql)
        conn.close()
    
    def table_create_erass_mc_sampled_systems(self):
        logging.debug('Creating table ERASS_MC_SAMPLED_SYSTEMS')
        conn = sqlite3.connect(self.db)
        sql = """CREATE TABLE IF NOT EXISTS ERASS_MC_SAMPLED_SYSTEMS(
                system_row_id INT,
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
        
    def table_create_erass_evolution(self):
        conn = sqlite3.connect(self.db)
        sql = """CREATE TABLE IF NOT EXISTS ERASS_EVOLUTION(
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
                duty_cycle REAL,
                run_id CHAR);"""
        conn.execute(sql)
        sql = """CREATE INDEX IF NOT EXISTS idx_erass_evolution_run_id
         ON ERASS_EVOLUTION (run_id);"""
        conn.execute(sql)
        conn.close()
    
    def table_load(self, table_name):
        """
        Load all rows from SQLite table.
        """
        logging.debug('Loading table %s', table_name)
        conn = sqlite3.connect(self.db)
        df = pd.read_sql_query(f"SELECT * from {table_name}", conn)
        conn.close()
        return df
    
    def table_load_erass_mc_info(self):
        self.df_erass_mc_info = self.table_load('ERASS_MC_INFO')
    
    def table_load_erass_mc_sampled_systems(self):
        self.df_erass_mc_sampled_systems = self.table_load('ERASS_MC_SAMPLED_SYSTEMS')
    
    def table_load_classifications(self):
        self.df_classifications = self.table_load('CLASSIFICATIONS')
    
    def table_load_transient(self):
        self.df_transient = self.table_load('TRANSIENT')
        
    def table_load_erass_mc_results(self):
        self.df_erass_evolution = self.table_load('ERASS_MC_RESULTS')
        
    def table_load_erass_evolution(self):
        self.df_erass_evolution = self.table_load('ERASS_EVOLUTION')
        
    def table_delete_erass_evolution_non_1_duty_cycle(self):
        """Delete all rows in erass_evolution that do not have a duty cycle of 1."""
        logging.debug('Deleteing all non duty=1.0 rows from ERASS_EVOLUTION')
        conn = sqlite3.connect(self.db)
        sql = """DELETE FROM ERASS_EVOLUTION
                 WHERE duty_cycle != 1.0;"""
        conn.execute(sql)
        conn.close()

    def table_erass_mc_results_map_info(self):
        logging.debug('Mapping ERASS_MC_INFO to ERASS_MC_RESULTS')
        try:
            self.df_erass_mc_info
        except:
            self.table_load_erass_mc_info()
        try:
            self.df_erass_evolution
        except:
            self.table_load_erass_evolution()

        info = self.df_erass_mc_info.set_index('run_id')
        self.df_erass_evolution['bh_ratio'] = self.df_erass_evolution['run_id'].map(info['bh_ratio'])    
        self.df_erass_evolution['size'] = self.df_erass_evolution['run_id'].map(info['size'])    
        self.df_erass_evolution['dincl_cutoff'] = self.df_erass_evolution['run_id'].map(info['dincl_cutoff'])
        self.df_erass_evolution['Z'] = self.df_erass_evolution['run_id'].map(info['Z'])  
        self.df_erass_evolution['erass_system_period_cutoff'] = self.df_erass_evolution['run_id'].map(info['erass_system_period_cutoff'])
        

    def table_classifications_map_systems(self):
        logging.debug('Mapping systems to CLASSIFICATIONS')
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
        logging.debug('Mapping ERASS_MC_INFO to CLASSIFICATIONS')
        info = self.df_erass_mc_info.set_index('run_id')
        self.df_classifications['bh_ratio'] = self.df_erass_evolution['run_id'].map(info['bh_ratio'])    
        self.df_classifications['size'] = self.df_erass_evolution['run_id'].map(info['size'])    
        self.df_classifications['dincl_cutoff'] = self.df_erass_evolution['run_id'].map(info['dincl_cutoff'])
        self.df_classifications['Z'] = self.df_erass_evolution['run_id'].map(info['Z'])


    def table_classifications_pivot(self, margins=True, split_Z=True, split_dincl=True):
        logging.debug('Creating classifications pivot table')
        try:
            self.df_classifications['Z']
        except KeyError:
            self.table_classifications_map_systems()
        
        if split_dincl and split_Z:
            self.table_classifications_map_info()
            piv = pd.pivot_table(self.df_classifications,
                     columns=['is_bh'],
                     index=['Z', 'dincl_cutoff', 'lc_classification'],
                     aggfunc='count',
                     margins=margins,
                     margins_name='total').run_id
        
        if split_Z:
            piv = pd.pivot_table(self.df_classifications,
                                 columns=['is_bh'],
                                 index=['Z', 'lc_classification'],
                                 aggfunc='count',
                                 margins=margins,
                                 margins_name='total').run_id

        if not split_Z:
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
        logging.debug('Splitting CLASSIFICATIONS table by metallicity')
        try:
            self.df_classifications['Z']
        except:
            self.table_classifications_map_systems()
            
        self.df_c_02 = self.df_classifications[self.df_classifications['Z']==0.02]
        self.df_c_002 = self.df_classifications[self.df_classifications['Z']==0.002]
        self.df_c_0002 = self.df_classifications[self.df_classifications['Z']==0.0002]
    
    def table_classifications_calc_intermediate(self):
        logging.debug('calculating intermediate CLASSIFICATIONS tables')
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
        
        
    
    def sim_erass_run(self,
                      size=500,
                      dincl_cutoff=46,
                      Z='all',
                      bh_ratio=0.5,
                      erass_system_period_cutoff=1460):
        """
        eRASS Monte-Carlo Simulation
        ----------------------------
        Samples a given parent population to create ULX lightcurves
        to evaulate the number of alive/dead/transient lightcurves.
        
        as well as samples each transient lightcurve to obtain the probability
        of the source being transient for a given eRASS cycle.

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
        Returns
        -------
        None.

        """
        logging.debug('Running eRASS MC')
        
        # Create SQLite tables
        self.table_create_erass_mc_info()
        self.table_create_classifications()
        self.table_create_transient()
        self.table_create_erass_mc_sampled_systems()
        self.table_create_erass_mc_results()
        
        # Get columns needed for simulation
        logging.debug('Getting simulation columns')
        sim_cols = ['theta_half_deg', 'Lx1', 'P_wind_days', 'P_sup_days', 'idum_run', 'iidd_old', 'lmxrb']
        
        # Delete large dataframes to save memory
        # del(self.pop.df)
        
        logging.debug('Checking population Z')
        self.pop.filter_df_ulx_by_Z(Z)

        # Filter population df to only have columns we are interested in to save memory
        logging.debug('Filtering out unused columns')
        self.pop.df_ulx = self.pop.df_ulx[sim_cols]
        
        # Create run info dataframe for storing id, size etc..
        logging.debug('Creating info table')
        df_info = pd.DataFrame()

        # Random run id
        run_id = str(uuid4())
        logging.debug('run_id: %s', run_id)

        # Create info table        
        df_info = pd.DataFrame()
        df_info['run_id'] = [run_id]
        df_info['size'] = [size]
        df_info['bh_ratio'] = [bh_ratio]
        df_info['dincl_cutoff'] = [dincl_cutoff]
        df_info['Z'] = [Z]
        df_info['erass_system_period_cutoff'] = [erass_system_period_cutoff]
        
        # Sample Systems
        logging.debug('Sampling systems, inclinations and dincl')
        sampled_indexs = self.pop.sample_ulxs(bh_ratio, size=size)
        selected_inclinations = np.random.randint(0, 91, size=size)
        selected_dincls = np.random.randint(0, dincl_cutoff, size=size)
        
        # Get sampled systems from populations df
        logging.debug('Retrieving sampled systems')
        df_sampled = self.pop.df_ulx.loc[sampled_indexs]
        df_sampled['Lx1'] = df_sampled['Lx1'] / 1e39
        df_sampled['inclination'] = selected_inclinations
        df_sampled['dincl'] = selected_dincls
        df_sampled['system_row_id'] = sampled_indexs
        df_sampled['run_id'] = run_id
        
        # Write sampled systems to sql
        logging.debug('Writing sampled systems')
        conn = sqlite3.connect(self.db)
        df_sampled.to_sql('ERASS_MC_SAMPLED_SYSTEMS', conn, if_exists='append', index=False)
        df_info.to_sql('ERASS_MC_INFO', conn, if_exists='append', index=False)
        conn.close()
        
        # Call ULXLC
        logging.debug('Calling ulxlc')
        run(["ulxlc/ulxlc_code_v0.1/a.out", str(run_id), str(size), str(erass_system_period_cutoff)])


    def sim_classifications_sampler_with_duty_cycle(self, duty_cycle):
        """
        Run through classifications and apply a duty cycle to calculate
        'lc_classification_new'.
        
        Parameters
        ----------
        duty_cycle : float
            between 0.0 and 1.0.
        """
        self.table_load_classifications()
        self.table_classifications_map_systems()
         
        logging.debug('Recalculating classifcations using duty cycle=%s', duty_cycle)
        
        self.df_classifications['duty_cycle'] = np.where(self.df_classifications['lmxrb']==1, duty_cycle, 1.0)
        self.df_classifications[['lc_classification', 'duty_cycle']]
        self.df_classifications['rand'] = np.random.random(size=len(self.df_classifications))
        self.df_classifications['rand>dc'] =  self.df_classifications['rand'] > self.df_classifications['duty_cycle']#System is off
        
        # rolled duty cycle = dead and if alive or transient then make dead
        cond = (self.df_classifications['rand>dc']==True) & ( (self.df_classifications['lc_classification']==1) | (self.df_classifications['lc_classification']==2))
        self.df_classifications['lc_classification_new'] = np.where(cond, 0, self.df_classifications['lc_classification'])

    def sim_bh_ratio_classifications_sampler(self, N=10000, size=500, dincl_cutoff=46, Z='all'):
        """
        Re-sample lightcurve classifications.

        Parameters
        ----------
        N : int
            Number of repeats. The default is 10000.
        size : int
            Number of clasifications per run. The default is 500.
        dincl_cutoff : TYPE, optional
            DESCRIPTION. The default is 46.
        Z : TYPE, optional
            DESCRIPTION. The default is 'all'.

        Returns
        -------
        class_dict : TYPE
            DESCRIPTION.

        """
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


    def calc_N_persistent(self, size=500):
        """Calculate the number of systems that had half opening angles > 45
        and thus were not put through the lightcurve simulation"""
        self.df_N_persistent = size - self.df_classifications['run_id'].value_counts()


    def get_run_ids(self, group_period_cutoff=True):
        """
        Get simulation parameters run_id dictionary.

        Parameters
        ----------
        group_period_cutoff : bool
            Group results together irrespective of period cutoff.
            The default is True.

        Returns
        -------
        dict_MC_run_ids : dict
            dictionary with keys as tuple of parameters and entries as
            associated run_ids.
        """
        try:
            self.df_erass_mc_info
        except AttributeError:
            self.table_load_erass_mc_info()
        
        info = self.df_erass_mc_info.copy()
        info = info.set_index('run_id')
        
        sim_cols = list(info.columns)
        
        if group_period_cutoff:
            sim_cols.remove('erass_system_period_cutoff')
            
        self.dict_MC_run_ids = info.groupby(by=sim_cols).groups
        return self.dict_MC_run_ids
    
    def get_run_counts(self):
        count_dict = {}
        for k, v in self.dict_MC_run_ids.items():
            count_dict[k] = len(v)
        return count_dict

    def get_classifications(self, key):
        """
        Parameters
        ----------
        key : tuple from key in self.dict_MC_run_ids
        """
        try:
            self.dict_MC_run_ids
        except:
            self.get_run_ids()
        
        ids = self.dict_MC_run_ids[key]
        df_classifications = self.df_classifications.set_index('run_id')
        # See https://pandas.pydata.org/pandas-docs/stable/user_guide/indexing.html#deprecate-loc-reindex-listlike
        df_classifications = df_classifications.loc[df_classifications.index.intersection(ids).unique()]
        return df_classifications
    
    def get_erass_evolution_from_key(self, key):
        """
        Parameters
        ----------
        key : tuple from key in self.dict_MC_run_ids
        """
        try:
            self.dict_MC_run_ids
        except:
            self.get_run_ids(group_period_cutoff=False)
        
        ids = self.dict_MC_run_ids[key]
        df_erass_evolution = self.df_erass_evolution.set_index('run_id')
        df_erass_evolution = df_erass_evolution.loc[df_erass_evolution.index.intersection(ids).unique()]
        return df_erass_evolution
    


    def calc_classification_counts(self, key):
        """
        Count the number of alive, dead and persistent sources for a given
        key.

        Parameters
        ----------
        key : tuple
            key from dict_MC_run_ids.

        Returns
        -------
        df_classification_counts : TYPE
            DESCRIPTION.

        """
        try:
            self.df_N_persistent
        except AttributeError:
            self.calc_N_persistent()
            
        df_c = self.get_classifications(key)
        
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
    
    
    def calc_bh_ratio_classfication_correlation_constant(self):
        """
        Calculate the result of fitting a line with y=mx+c
        to bh_ratio vs number in each classification

        Returns
        -------
        pd.DataFrame
        """
        def line(x, m, c):
            return m*x+c
        
        res = []
        for Z in ['0.02', '0.002', '0.0002']:
            for dincl in [46, 21]:
                sub = self.df_classifications_counts_stats.loc[Z,dincl,:]
                for classification in ['N_alive_persistent',
                          'N_trans',
                          'N_dead']:
    
                    sub2 = sub[classification]
                    xdata = [0.0, 0.25, 0.5, 0.75, 1.0]
                    ydata  = sub2['mean']
                    ystd  = sub2['std']
                    
                    popt, pcov = curve_fit(line, xdata, ydata, sigma=ystd)
                    perr = np.sqrt(np.diag(pcov))
                    
                    m = popt[0]
                    c = popt[1]
                    
                    m_std = perr[0]
                    c_std = perr[1]
                    
                    res.append((Z, dincl, classification, m, m_std))
                    
        df_res = pd.DataFrame(res, columns=['Z', 'dincl', 'classification', 'm', 'm_std'])
        df_res = df_res.set_index(['Z', 'dincl', 'classification'])
        return df_res
    
    def calc_all_classifications_count_stats(self):
        res = pd.DataFrame()
        for k in self.dict_MC_run_ids.keys():
            print(k)
            df = self.calc_classification_counts(k)
            df['size'] = k[0]
            df['bh_ratio'] = k[1]
            df['dincl_cutoff'] = k[2]
            df['Z'] = k[3]
            res = res.append(df)
        gb = res.groupby(['Z', 'dincl_cutoff', 'bh_ratio']).agg(['min', 'mean', 'max', 'std', 'count'])
        gb = gb.drop('size', axis=1)
        self.df_classifications_counts_stats = gb
        return self.df_classifications_counts_stats
        
        
    def ADT_result_stats(self):
        try:
            self.df_classifications['Z']
        except:
            self.table_classifications_map_info()
        self.df_classifications_stats = self.df_classifications.groupby(['dincl_cutoff', 'Z', 'bh_ratio']).agg(['min', 'mean', 'max', 'std', 'count'])
        return self.df_classifications_stats
        
    def ERASS_result_stats(self):
        try:
            self.df_erass_evolution['dincl_cutoff']
        except:
            self.table_erass_mc_results_map_info()
            
        self.df_erass_mc_stats = self.df_erass_evolution.groupby(['dincl_cutoff', 'Z', 'duty_cycle', 'bh_ratio', 'period', 'erass_cycle']).agg(['min', 'mean', 'max', 'std', 'count'])
        return self.df_erass_mc_stats


class Ulx:
    """
    This class is used to represent a single ULX, and used for
    simulating the number of detected ULXs over the course of eRASS
    via sampling the cycle dependent transient probabilities obtained via
    MC of lightcurves produced by ULXLC.
    """
    def __init__(self):
        self.is_ulx = None        # Is the source currently a ULX?
        self.is_transient = None  # Is the source currently transient?
        self.P_cycle_1_ulx = None # Probability of source being ULX on cycle 1
        self.P_transient = [None, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0] # Transient probabilities for each erass cycle
        self.transient_cycle = None # Which cycle was the source identifed as transient?
        self.is_lmxrb = None        # Is the source classified as a LMXRB?
    
    @classmethod
    def persistent_alive_system(cls):
        """
        Create ULX system that is persistently a ULX
        over the course of eRASS."""
        # logging.debug('Initializing persistent alive ULX')
        Ulx = cls()
        Ulx.P_cycle_1_ulx = 1.0
        return Ulx
    
    @classmethod
    def persistent_dead_system(cls):
        # logging.debug('Initializing persistent dead ULX')
        Ulx = cls()
        Ulx.P_cycle_1_ulx = 0.0
        return Ulx
    
    @classmethod
    def from_table_transient_row(cls, Series, period):
        """
        Create ULX from a row from transient probability row.

        Parameters
        ----------
        Series : pd.Series
            Row from CLASSFICATIONS table
        period : string
            Which period to use, either 'P_wind' or 'P_sup'
        """
        # logging.debug('initializing ULX from row, period=%s', period)
        
        Ulx = cls()

        Ulx.P_cycle_1_ulx = Series.erass_1_ulx_prob

        if period == 'P_wind':
            Ulx.P_transient[1] = Series.erass_2_P_wind_transient_prob
            Ulx.P_transient[2] = Series.erass_3_P_wind_transient_prob
            Ulx.P_transient[3] = Series.erass_4_P_wind_transient_prob
            Ulx.P_transient[4] = Series.erass_5_P_wind_transient_prob
            Ulx.P_transient[5] = Series.erass_6_P_wind_transient_prob
            Ulx.P_transient[6] = Series.erass_7_P_wind_transient_prob
            Ulx.P_transient[7] = Series.erass_8_P_wind_transient_prob
        elif period == 'P_sup':
            Ulx.P_transient[1] = Series.erass_2_P_sup_transient_prob
            Ulx.P_transient[2] = Series.erass_3_P_sup_transient_prob
            Ulx.P_transient[3] = Series.erass_4_P_sup_transient_prob
            Ulx.P_transient[4] = Series.erass_5_P_sup_transient_prob
            Ulx.P_transient[5] = Series.erass_6_P_sup_transient_prob
            Ulx.P_transient[6] = Series.erass_7_P_sup_transient_prob
            Ulx.P_transient[7] = Series.erass_8_P_sup_transient_prob
        return Ulx
    
    def observe(self, cycle):
        """
        Simulate observing ULX at a given eRASS cycle.
            
        Parameters
        ----------
        cycle : int
            eRASS cycle to observe the source at

        Returns
        -------
        result : list
            [self.is_ulx, self.is_transient]
            at the given cycle, was the source observed as transient, and was it a ulx or not?
        """
        logging.debug('Simulating observing ULX on eRASS cycle: %s', cycle)

        # If first cycle, then check to see if the source was a ULX or not.
        if cycle == 0:
            r = np.random.random()
            self.is_ulx = r < self.P_cycle_1_ulx
            logging.debug('cycle: %s, roll: %s, P_cycle_1_ulx: %s, is_ulx: %s', cycle, r, self.P_cycle_1_ulx, self.is_ulx)
        else:
            r = np.random.random()
            self.is_transient = r < self.P_transient[cycle]
            logging.debug('roll: %s, P_transient[%s]: %s, is_transient: %s', r, cycle, self.P_transient[cycle], self.is_transient)
            if self.is_transient:
                logging.debug('ULX identified as transient, flipping self.is_ulx')
                self.is_ulx = not self.is_ulx
                if self.transient_cycle == None:
                    self.transient_cycle = cycle
                    logging.debug('First transient cycle found at cycle: %s', self.transient_cycle)

        return [self.is_ulx, self.is_transient]



class ErassTransientSampler:
    """
    This class is used to sample the transient probabilities outputted from /ulxlc/ulxlcmod.c
    to determine how the observational nature of ULXs change over the course of eRASS.
    """
    def __init__(self):
        logging.debug('Initializing eRASS transient sampler')
        
        self.db = None              # SQL Database path
        self.run_id = None          # run_id to retrieve
        self.period = None          # Which period to use for simulation 'P_sup' or 'P_wind'
        self.duty_cycle = None
        self.systems = None         # Systems is provided as a list containing Ulxs [Ulx(), Ulx(), ...]

        # These are the two directly observable quantities we are concerning outselves with.
        self.N_new_systems 				   = np.array([0,0,0,0,0,0,0,0])
        self.N_old_system_become_transient = np.array([0,0,0,0,0,0,0,0])

    def set_db(self, db):
        self.db = db

    def set_period(self, period):
        """Either 'P_wind or P_sup'"""
        self.period = period
        
    def set_run_id(self, run_id):
        self.run_id = run_id

    def set_duty_cycle(self, duty_cycle):
        self.duty_cycle = duty_cycle

    def set_active_population(self, pop):
        self.pop = pop
    
    def load_table(self, table_name):
        logging.debug('Getting run_id: %s from table: %s', self.run_id, table_name)
        sql = f"""SELECT * FROM {table_name}
                  WHERE run_id='{self.run_id}'"""
        conn = sqlite3.connect(self.db)
        df = pd.read_sql_query(sql, conn)
        conn.close()
        return df
    
    def table_load(self, table_name):
        logging.debug('Loading table %s', table_name)
        conn = sqlite3.connect(self.db)
        df = pd.read_sql_query(f"SELECT * from {table_name}", conn)
        conn.close()
        return df
    
    
    def populate_systems(self):
        """
        Populate systems list from a given run_id and period.

        Parameters
        ----------
        run_id : str
            run_id to retrieve
            
        Returns
        -------
        systems : list
            list of Ulx() objects
        """
        logging.debug('Getting systems with db=%s, run_id=%s, period=%s', self.db, self.run_id, self.period)
       
        # Retrieve rows from SQL database
        df_sampled_systems = self.load_table('ERASS_MC_SAMPLED_SYSTEMS')
        df_classifications = self.load_table('CLASSIFICATIONS')
        df_transient = self.load_table('TRANSIENT')
        
        
        # Calculate number of systems in each classification
        class_count = df_classifications['lc_classification'].value_counts()
        logging.debug('class_count: %s', class_count)
        
        
        # On the off-chance that we got none of a paticular classification 
        if 0 not in class_count.index:
            logging.debug('No dead systems found, setting N_class_dead=0')
            N_class_dead = 0
        else:
            N_class_dead = class_count[0]
        
        if 1 not in class_count.index:
            logging.debug('No transient systems found, setting N_class_trans=0')
            N_class_trans = 0
        else:
            N_class_trans = class_count[0]
            
        if 2 not in class_count.index:
            logging.debug('No alive systems found, setting N_class_alive=0')
            N_class_alive = 0
        else:
            N_class_alive = class_count[2]
                
        N_samp = len(df_sampled_systems)                   # Total sampled systems
        N_class = len(df_classifications)                  # Systems with opening angles < 45
        N_class_not_simulated = N_samp - N_class           # Systems with opening angles > 45 (These systems are persistent)
        N_trans = len(df_transient)                        # Systems sampled for eRASS
        N_trans_not_simulated = N_class - N_trans          # Systems with period > erass_systems_period_cuttoff 
        N_alive = N_class_not_simulated + N_class_alive    # Total number of persistent systems
        
        logging.debug('N_samp = %s', N_samp)
        logging.debug('N_class = %s', N_class)
        logging.debug('N_class_not_simulated = %s', N_class_not_simulated)
        logging.debug('N_trans = %s', N_trans)
        logging.debug('N_trans_not_simulated = %s', N_trans_not_simulated)
        logging.debug('N_alive = %s', N_alive)

        # Check for LMXRB systems
        if self.duty_cycle != None or self.duty_cycle != 1.0:
            logging.debug('duty cyle=%s | not None or not 1.0, looking for lmxrb systems', self.duty_cycle)
            
            df_lmxrb = df_sampled_systems.loc[df_sampled_systems['lmxrb'] == 1]
            N_samp_lmxrb = len(df_lmxrb)
        
        
            if N_samp_lmxrb > 1:
                logging.debug('lmxrb systems found N_samp_lmxrb=%s', N_samp_lmxrb)
                
                logging.debug('mapping lmxrb flag to df_transient')
                df_transient.loc[:, 'lmxrb'] = df_transient['system_row_id'].map(self.pop.df['lmxrb'])
        
                # get lmxrb transient systems
                df_transient_lmxrb = df_transient.loc[df_transient['lmxrb'] == 1].copy()
                df_transient_non_lmxrb = df_transient.loc[df_transient['lmxrb'] == 0].copy()
                
                N_trans_lmxrb = len(df_transient_lmxrb)
                
                logging.debug('N_trans_lmxrb=%s',N_trans_lmxrb)
                
                logging.debug('Sampling lmxrb transients via duty cycle')
                
                # Sample lmxrb transients using duty cycle
                df_transient_lmxrb['duty_cycle'] = self.duty_cycle
                df_transient_lmxrb['rand'] = np.random.random(size=N_trans_lmxrb)
                df_transient_lmxrb['rand<duty_cycle'] = df_transient_lmxrb['rand'] < df_transient_lmxrb['duty_cycle']
                df_transient_lmxrb_active = df_transient_lmxrb[df_transient_lmxrb['rand<duty_cycle']==True]
                
                # Also get the number of transients thrown out (these are dead)
                N_trans_lmxrb_active = len(df_transient_lmxrb_active)
                N_trans_lmxrb_dead = N_trans_lmxrb - N_trans_lmxrb_active
                
                logging.debug('N_trans_lmxrb_active=%s', N_trans_lmxrb_active)
                logging.debug('N_trans_lmxrb_dead=%s', N_trans_lmxrb_dead)


        systems_alive = [Ulx.persistent_alive_system() for i in range(N_alive)]
        systems_dead = [Ulx.persistent_dead_system() for i in range(N_class_dead)]
        
        
        if self.duty_cycle != None or self.duty_cycle!=1.0:
            systems_transient = [Ulx.from_table_transient_row(row, self.period) for i, row in df_transient_non_lmxrb.iterrows()]
            systems_transient_lmxrb = [Ulx.from_table_transient_row(row, self.period) for i, row in df_transient_lmxrb_active.iterrows()]
            systems_dead_lmxrb = [Ulx.persistent_dead_system() for i in range(N_trans_lmxrb_dead)]
            self.systems = systems_alive + systems_dead + systems_transient + systems_transient_lmxrb + systems_dead_lmxrb
        else:
            systems_transient = [Ulx.from_table_transient_row(row, self.period) for i, row in df_transient.iterrows()]
            self.systems = systems_alive + systems_dead + systems_transient
        return self.systems
     
        
    def calc_secondary_quantities(self):
        logging.debug('Calculating secondary observable quantities')
        self.N_delta_obs_ulx = self.N_new_systems - self.N_old_system_become_transient
        self.N_observed_ulxs = np.cumsum(self.N_delta_obs_ulx)
        self.N_transients = self.N_new_systems[1:] + self.N_old_system_become_transient[1:]
        self.N_transients = np.insert(self.N_transients, 0, 0)
        self.N_transients_cum = np.cumsum(self.N_transients)
        self.N_total_systems = np.cumsum(self.N_new_systems)
        self.N_persistent_ulx_systems = self.N_new_systems[0] - np.cumsum(self.N_old_system_become_transient)

     
    def run(self):
        """
        Run eRASS Transient probability sampler.
        """
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
        res['period'] = self.period
        res['duty_cycle'] = self.duty_cycle
        res['run_id'] = self.run_id
        self.res = res
        return res
    
    def write_results_to_erass_evolution(self):
        logging.debug('Writing results to ERASS_MC_RESULTS')
        conn = sqlite3.connect(self.db)
        self.res.to_sql('ERASS_EVOLUTION', conn, if_exists='append', index=False)
        conn.close()
        



class Plotter:
    def __init__(self):
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
        
        self.set_latex_font()
        self.set_savefig(False)
    
    def set_latex_font(self):
        matplotlib.rcParams['mathtext.fontset'] = 'stix'
        matplotlib.rcParams['font.family'] = 'STIXGeneral'
    
    def set_results_processor(self, ResultsProcessor):
        self.rp = ResultsProcessor
        
    def set_savefig(self, savefig=False):
        """Turn saving figures on or off"""
        self.savefig = savefig
    
    
    def bar_classifications_Z(self):
        piv = self.rp.table_classifications_pivot(margins=False)
        piv.plot(kind='bar')
        
    
    def classifications_dincl_i(self, percent=False):
        try:
            self.rp.df_d
        except:
            self.rp.table_classifications_calc_intermediate()
        
        fig, ax = plt.subplots(ncols=2, nrows=1, figsize=(8,3))
        plt.tight_layout()
        ax[0].set_xlim(0,max(self.rp.df_d['system_dincl']))
        ax[0].set_xlabel(r'Precessional angle $\Delta i$')
        
        ax[0].minorticks_on()
        
        
        ax[1].set_xlim(0,max(self.rp.df_d['system_inclination']))
        ax[1].set_xlabel(r'Inclination $i$')
        ax[1].minorticks_on()
        
        
        if percent == False:
            ax[0].set_ylabel(r'Number of systems')
            ax[1].set_ylabel(r'Number of systems')
            
            # dincl vs Number systems
            self.rp.a_dincl_N.plot(label='Alive', linestyle=self.linestyle_alive, color='black', linewidth=0.8, ax=ax[0])
            self.rp.t_dincl_N.plot(label='Transient', linestyle=self.linestyle_trans, color='black', linewidth=0.8, ax=ax[0])
            self.rp.d_dincl_N.plot(label='Dead', linestyle=self.linestyle_dead, color='black', linewidth=0.8, ax=ax[0])
            
            # inclination vs Number systems
            self.rp.a_i_N.plot(label='Alive', linestyle=self.linestyle_alive, color='black', linewidth=0.8, ax=ax[1])
            self.rp.t_i_N.plot(label='Transient', linestyle=self.linestyle_trans, color='black', linewidth=0.8, ax=ax[1])
            self.rp.d_i_N.sort_index().plot(label='Dead', linestyle=self.linestyle_dead, color='black', linewidth=0.8, ax=ax[1])
            
            ax[1].set_ylabel(r'Number of systems')
            
            if self.savefig:
                plt.savefig('../reports/figures/dincl_i_classifications.png', dpi=1000)
                plt.savefig('../reports/figures/dincl_i_classifications.eps')
            
        else:
            ax[0].set_ylabel(r'% of systems')
            ax[1].set_ylabel(r'% of systems')
            ax[0].set_ylim(0,100)
            ax[1].set_ylim(0,100)
            
            # dincl vs percentages
            self.rp.a_dincl_percent.plot(label='Alive', linestyle=self.linestyle_alive, color='black', linewidth=0.8, ax=ax[0])
            self.rp.t_dincl_percent.plot(label='Transient', linestyle=self.linestyle_trans, color='black', linewidth=0.8, ax=ax[0])
            self.rp.d_dincl_percent.plot(label='Dead', linestyle=self.linestyle_dead, color='black', linewidth=0.8, ax=ax[0])
            
            # inclination vs percentages
            self.rp.a_i_percent.plot(label='Alive', linestyle='-', color='black', linewidth=0.8, ax=ax[1])
            self.rp.t_i_percent.plot(label='Transient', linestyle='--', color='black', linewidth=0.8, ax=ax[1])
            self.rp.d_i_percent.plot(label='Dead', linestyle='dotted', color='black', linewidth=0.8, ax=ax[1])
            
            
            if self.savefig:
                plt.savefig('../reports/figures/dincl_i_classifications_percent.png', dpi=1000)
                plt.savefig('../reports/figures/dincl_i_classifications_percent.eps')
        ax[0].legend()
        ax[1].legend()
        plt.show()
        
    
    
    def erass_mc_results_hist(self, key, by='cycle'):
        try:
            self.rp.df_erass_evolution['Z']
        except:
            self.rp.table_erass_mc_results_map_info()
        
        params = ['N_new_systems',
                  'N_old_system_become_transient',
                  'N_observed_ulxs',
                  'N_delta_obs_ulx',
                  'N_transients',
                  'N_transients_cum',
                  'N_total_systems']
        
        res = self.rp.get_erass_evolution_from_key(key)
        
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
                if self.savefig==True:
                    plt.savefig(f'../reports/figures/e_{p}_{by}_{erass_system_period_cutoff}_{dincl_cutoff}_{Z}.png',dpi=500)
                    plt.savefig(f'../reports/figures/e_{p}_{by}_{erass_system_period_cutoff}_{dincl_cutoff}_{Z}.eps')
                    plt.savefig(f'../reports/figures/e_{p}_{by}_{erass_system_period_cutoff}_{dincl_cutoff}_{Z}.pdf')


    def classifications_hist(self, key, duty_cycle=1.0, frac_visible=False, normdist=False):
        """
        Plot classifications histogram for a given simulation key.

        Parameters
        ----------
        key : tuple
            Simulation key eg (500, 0.5, 21, '0.0002')
        frac_visible : bool
            Plot the % of visible systems i.e alive and transient.
            The default is False.
        normdist : bool
            Plot fitted normal histograms
        Returns
        -------
        None.
        """
        def make_norm(df, col):
            """Create scipy.stats.norm object"""
            n = norm(loc=df_classification_counts[col].mean(),
                         scale=df_classification_counts[col].std())
            return n
        
        # Confidence intervals ranges for plotting normdists
        lim_lower = 0.00001
        lim_upper = 0.99999
        
        fig, ax = plt.subplots(5,1, sharex=True, sharey=True, figsize=(6, 4.5))
        
        Z = key[3]
        dincl_cut = key[2]
        
        ax[0].set_title(fr'Z = {Z} | $\Delta i_{{max}} ={dincl_cut-1}^{{\circ}}$ | d = {duty_cycle}')

        for i, bh in enumerate([0.0, 0.25, 0.5, 0.75, 1.0]):
            print(i, bh)
            
            key = list(key)
            key[1] = bh
            key = tuple(key)
            
            df_classification_counts = self.rp.calc_classification_counts(key)
            
        
 
            if frac_visible==False:
                xmax=500
                ax[-1].set_xlabel('Number of sources')
                
                if normdist==False:
                    xbins = np.arange(0, xmax, 5)
                    
                    df_classification_counts['N_dead'].hist(bins=xbins, label='Dead', alpha=0.8, edgecolor='black', linestyle=self.linestyle_dead, histtype='step', ax=ax[i], grid=False)
                    df_classification_counts['N_trans'].hist(bins=xbins, label='Transient', alpha=0.8, edgecolor='black', linestyle=self.linestyle_trans, histtype='step', ax=ax[i], grid=False)
                    df_classification_counts['N_alive_persistent'].hist(bins=xbins, label='Alive', alpha=0.8, edgecolor='black', linestyle=self.linestyle_alive, histtype='step', ax=ax[i], grid=False)
                    
                if normdist==True:
                    norm_dead = make_norm(df_classification_counts, 'N_dead')
                    norm_trans = make_norm(df_classification_counts, 'N_trans')
                    norm_alive = make_norm(df_classification_counts, 'N_alive_persistent')
                    
                    x_d = np.linspace(norm_dead.ppf(lim_lower), norm_dead.ppf(lim_upper), 100)
                    x_t = np.linspace(norm_trans.ppf(lim_lower), norm_trans.ppf(lim_upper), 100)
                    x_a = np.linspace(norm_alive.ppf(lim_lower), norm_alive.ppf(lim_upper), 100)
                    
                    ax[i].plot(x_d, norm_dead.pdf(x_d), linestyle=self.linestyle_dead, color='black', label='Dead')
                    ax[i].plot(x_t, norm_trans.pdf(x_t), linestyle=self.linestyle_trans, color='black', label='Transient')
                    ax[i].plot(x_a, norm_alive.pdf(x_a), linestyle=self.linestyle_alive, color='black', label='Alive')
                    
            else:
                xmax=100
                ax[-1].set_xlabel('% of observed ULXs')
                
                ax[i].axvspan(self.PERCENT_ALIVE_EARNSHAW_LOWER, self.PERCENT_ALIVE_EARNSHAW_UPPER, alpha=0.3, color='grey')
                ax[i].axvspan(self.PERCENT_TRANS_EARNSHAW_LOWER, self.PERCENT_TRANS_EARNSHAW_UPPER, alpha=0.7, color='grey')
                
                if normdist==False:
                    xbins = np.arange(0, xmax, 1)
                    df_classification_counts['frac_alive_visible'].hist(bins=xbins, label='% Alive', alpha=0.8, edgecolor='black', linestyle=self.linestyle_alive, histtype='step', ax=ax[i], grid=False, density=True)
                    df_classification_counts['frac_trans_visible'].hist(bins=xbins, label='% Transient', alpha=0.8, edgecolor='black', linestyle=self.linestyle_trans, histtype='step', ax=ax[i], grid=False, density=True)
                    
                    
                if normdist==True:
                    norm_alive = make_norm(df_classification_counts, 'frac_alive_visible')
                    norm_trans = make_norm(df_classification_counts, 'frac_trans_visible')
                    
                    x_a = np.linspace(norm_alive.ppf(lim_lower), norm_alive.ppf(lim_upper), 100)
                    x_t = np.linspace(norm_trans.ppf(lim_lower), norm_trans.ppf(lim_upper), 100)
                    
                    ax[i].plot(x_a, norm_alive.pdf(x_a), linestyle=self.linestyle_alive, color='black', label='% Alive')
                    ax[i].plot(x_t, norm_trans.pdf(x_t), linestyle=self.linestyle_trans, color='black', label='% Transient')
            
            ax[i].set_xlim(0, xmax)
            ax[i].set_ylim(0, 1.1*ax[0].get_ylim()[1])
            ax[i].text(x=ax[i].get_xlim()[1]*(1/2), y=ax[i].get_ylim()[1]*(1/2), s=r'$\%_{{BH}} = $ {}'.format((bh)*100), fontsize=9)    
            ax[i].legend(fontsize=7, loc='right')
            ax[i].tick_params(axis='y',
                              which='both',
                              left=False,
                              labelleft=False)
        
        if frac_visible:
            ymax = ax[0].get_ylim()[1]
            ax[0].text(x=(self.PERCENT_ALIVE_EARNSHAW_LOWER)+0.25, y=ymax+0.05*ymax, s='Alive $1 \sigma$ XMM', fontsize=7)
            ax[0].text(x=(self.PERCENT_TRANS_EARNSHAW_LOWER)+0.25, y=ymax+0.05*ymax, s='Transient $1 \sigma$ XMM', fontsize=7)


        plt.tight_layout()
        plt.subplots_adjust(hspace=0.0)
        
        if self.savefig:
            if frac_visible==False:
                plt.savefig(f'../reports/figures/ADT_{Z}_{dincl_cut}_{duty_cycle}_normdist_{normdist}.png', dpi=500)
                plt.savefig(f'../reports/figures/ADT_{Z}_{dincl_cut}_{duty_cycle}_normdist_{normdist}.eps')
                plt.savefig(f'../reports/figures/ADT_{Z}_{dincl_cut}_{duty_cycle}_normdist_{normdist}.pdf')
            else:
                plt.savefig(f'../reports/figures/ADT_frac_{Z}_{dincl_cut}_{duty_cycle}_normdist_{normdist}.png', dpi=500)
                plt.savefig(f'../reports/figures/ADT_frac_{Z}_{dincl_cut}_{duty_cycle}_normdist_{normdist}.eps')
                plt.savefig(f'../reports/figures/ADT_frac_{Z}_{dincl_cut}_{duty_cycle}_normdist_{normdist}.pdf')
        
        
    def erass_transients(self, erass_system_period_cutoff=999999, Z='all', duty_cycle=0.1, percent_transient=False):
        cycles = ['1', '2', '3', '4', '5', '6', '7', '8']
        dincls = self.rp.df_erass_mc_info['dincl_cutoff'].unique()
        bh_percents = np.sort(self.rp.df_erass_mc_info['bh_ratio'].unique())
        
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
                df = self.rp.get_erass_evolution_from_key(key)
                df = df[df['duty_cycle']==duty_cycle]


                df['percent_obs'] = (df['N_total_systems']/500) * 100
                df['percent_trans_to_obs'] = df['N_transients_cum'] / df['N_total_systems'] * 100
            
                sub_wind = df[df['period'] == 'P_wind']
                sub_sup = df[df['period'] == 'P_sup']
            
                agg_wind = sub_wind.groupby(['erass_cycle']).agg(['mean', 'std'])
                agg_sup = sub_sup.groupby(['erass_cycle']).agg(['mean', 'std'])
                
                if not percent_transient:
                    y_wind = agg_wind['N_transients']['mean']
                    y_wind_std = agg_wind['N_transients']['std']
                    y_sup = agg_sup['N_transients']['mean']
                    y_sup_std = agg_sup['N_transients']['std']
                
                if percent_transient:
                    percent_obs_wind = agg_wind['percent_obs']
                    percent_obs_sup = agg_sup['percent_obs']

                    y_wind = agg_wind['percent_trans_to_obs']['mean']
                    y_wind_std = agg_wind['percent_trans_to_obs']['std']
                    y_sup = agg_sup['percent_trans_to_obs']['mean']
                    y_sup_std = agg_sup['percent_trans_to_obs']['std']
                
                bh_label = f'$\%_{{BH}}$ = {bh*100}' 
                
                trans1 = Affine2D().translate(-spacing*(nbars)/4 + spacing*j, 0.0) + ax[0][i].transData
                trans2 = Affine2D().translate(-spacing*(nbars)/4 + spacing*j, 0.0) + ax[1][i].transData
        
                ax[0][i].errorbar(cycles, y_wind, yerr=y_wind_std, linestyle="none",
                                  linewidth=linewidth, capsize=1.0, transform=trans1, label=bh_label, c=clist[j])
                
                ax[1][i].errorbar(cycles, y_sup, yerr=y_sup_std, linestyle="none",
                      linewidth=linewidth, capsize=1.0, transform=trans2, label=bh_label, c=clist[j])
                
                # ax[0][i].errorbar(cycles, percent_obs_wind['mean'], yerr=percent_obs_wind['std'], linestyle="none",
                #                   linewidth=linewidth, capsize=1.0, transform=trans1, c=clist[j])
                
                # ax[1][i].errorbar(cycles, percent_obs_sup['mean'], yerr=percent_obs_sup['std'], linestyle="none",
                #                   linewidth=linewidth, capsize=1.0,  transform=trans2, c=clist[j])
        
    
                ax[0][i].set_title(f'$P_{{wind}}$ | $\Delta i_{{max}} = {dincl-1}^{{\circ}}$ | Z = {Z} | d = {duty_cycle}')
                ax[1][i].set_title(f'$P_{{sup}}$ | $\Delta i_{{max}} = {dincl-1}^{{\circ}}$ | Z = {Z} | d = {duty_cycle}')
        
                
                ax[i][0].grid(axis='y')
                ax[i][1].grid(axis='y')
                #ax[1][i].grid(axis='y')
                #ax[0][i].tick_params(labelrotation=90)
                #ax[1][i].tick_params(axis='x', labelrotation=90)
        
        if not percent_transient:
            ax[0][0].set_ylabel('N Transients')
            ax[1][0].set_ylabel('N Transients')
            
            ax[0][0].set_ylim(0,35)
            ax[1][0].set_ylim(0,95)
            
        if percent_transient:
            ax[0][0].set_ylabel('% Transient')
            ax[1][0].set_ylabel('% Transient')
            # ax[0][0].set_ylim(0,100)
            # ax[1][0].set_ylim(0,100)
        
        
        ax[0][0].legend()
        #ax[0][0].legend()
        
        plt.subplots_adjust(wspace=0)
        
        plt.tight_layout()
        if self.savefig:
            if not percent_transient:
                plt.savefig(f'../reports/figures/erass_N_transients_Z={Z}_d={duty_cycle}.eps', bbox_inches='tight')
                plt.savefig(f'../reports/figures/erass_N_transients_Z={Z}_d={duty_cycle}.png', bbox_inches='tight', dpi=1000)
            if percent_transient:
                plt.savefig(f'../reports/figures/erass_perc_transients_Z={Z}_d={duty_cycle}.eps', bbox_inches='tight')
                plt.savefig(f'../reports/figures/erass_perc_transients_Z={Z}_d={duty_cycle}.png', bbox_inches='tight', dpi=1000)



if __name__ == "__main__":
    # logging.basicConfig(level=logging.DEBUG)
    # logging.basicConfig(level=0)
    
    
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

    # rp.sim_bh_ratio_classifications_sampler(N=10)
    # rp.sim_classifications_sampler_with_duty_cycle(duty_cycle=0.1)

    
    # =============================================================================
    # Table related functions
    # =============================================================================
    # Load tables
    rp.table_load_erass_mc_info()
    rp.table_load_classifications()
    rp.table_load_erass_evolution()

    
    # rp.table_load_transient()
    # rp.table_load_erass_mc_results()    

    # rp.table_erass_mc_results_map_info()
    # rp.table_classifications_map_systems()
    # rp.table_classifications_pivot()
    
    # rp.table_classifications_calc_intermediate()
    # rp.table_classifications_map_systems()
    
    # tab = rp.table_classifications_pivot()
    # with open('../reports/table_mc_class_beamed.txt', mode='w+') as f:
    #     f.write(tab.to_latex(float_format="%.2f"))
    

    # =============================================================================
    # Run ULXLC over simulation parameter grid
    # =============================================================================
    
    # for Z in ['0.002', '0.0002', '0.02']:
    #     pop = populations.Population(df)
    #     rp.set_parent_population(pop)
    #     for i in range(300):
    #         for dincl in [21, 46]:
    #             for bh_ratio in [0.0, 0.25, 0.75, 1.0]:
    #                 print(i, bh_ratio, dincl, Z)
    #                 # import pdb; pdb.set_trace()
    #                 rp.sim_erass_run(size=500,
    #                                  dincl_cutoff=dincl,
    #                                  Z=Z,
    #                                  bh_ratio=bh_ratio,
    #                                  erass_system_period_cutoff=999999)
    
    # =============================================================================
    # Simulate transient probabilities
    # =============================================================================
    # rp.table_load_erass_evolution()
    # ids_unique_info = rp.df_erass_mc_info.run_id.unique()
    # ids_unique_erass_evolution = rp.df_erass_evolution.run_id.unique()
    # ids_to_simulate = np.setdiff1d(ids_unique_info, ids_unique_erass_evolution)
    
    # from tqdm import tqdm
    
    # for run_id in tqdm(ids_to_simulate):
    #     for duty_cycle in [1.0, 0.1, 0.2, 0.3]:
    #         for period in ['P_wind', 'P_sup']:
    #             ets = ErassTransientSampler()
                
    #             ets.set_period(period)
    #             ets.set_db('ulxlc.db')
    #             ets.set_duty_cycle(duty_cycle)
    #             ets.set_active_population(pop)
    #             ets.set_run_id(run_id)
    #             systems = ets.populate_systems()
                
    #             ets.run()
                
    #             res = ets.collect_results()
                
    #             ets.write_results_to_erass_evolution()
        
    # =============================================================================
    # Run over Classifications with duty cycle
    # =============================================================================
    
    # rp.sim_classifications_sampler_with_duty_cycle(0.1)
    # rp.df_classifications['lc_classification'] = rp.df_classifications['lc_classification_new']
    
    # from tqdm import tqdm
    
    # rp.get_run_ids(True)
    
    # to_latex=True
    
    # df2 = pd.DataFrame()
    # for d in [1.0, 0.3, 0.2, 0.1]:
    #     print(d)
    #     if d != 1.0:
    #         rp.table_load_classifications()
    #         rp.sim_classifications_sampler_with_duty_cycle(d)
    #         rp.df_classifications['lc_classification'] = rp.df_classifications['lc_classification_new']
            
        
    #     for k in tqdm(rp.dict_MC_run_ids.keys()):
    #         res = rp.calc_classification_counts(k)
    #         df = res.agg(['mean', 'std'])
    #         df['duty_cycle'] = d
    #         df['bh_ratio'] = k[1]
    #         df['dincl_cutoff'] = k[2]
    #         df['Z'] = k[3]
    #         df2 = df2.append(df)
        
    
    # if to_latex:
    #     df2 = df2.drop(['N_alive', 'N_persistent'], axis=1)
        
        
    
    # df2['quantity'] = df2.index
    # cols = ['dincl_cutoff', 'duty_cycle', 'Z', 'bh_ratio']
    # df2 = df2.pivot_table(index=cols, columns='quantity')
    # if to_latex:
    #     with open('../reports/table_classifications.txt', mode='w+') as f:
    #         f.write(df2.to_latex(float_format="%.2f"))

    
    # =============================================================================
    # Print Result statistics
    # =============================================================================
    

    
    # rp.get_run_ids(False)
    # rp.get_run_counts()

    # df_classifications_stats = rp.ADT_result_stats()
    # print('Classification stats:')
    # print(df_classifications_stats)
    
    
    # rp.get_run_ids()
    
    # res = pd.DataFrame()

    # rp.calc_all_classifications_count_stats()
    # bh_ratio_cc = rp.calc_bh_ratio_classfication_correlation_constant()
    
    # df_erass_mc_stats = rp.ERASS_result_stats()
    # print('Erass MC stats:')
    # print(df_erass_mc_stats)
    

    # =============================================================================
    # Plotting Functions
    # =============================================================================
    
    # Plot eRASS number of transients
    
    p = Plotter()
    p.set_results_processor(rp)
    p.set_savefig(True)
    rp.get_run_ids(False)
    
    
    for Z in ['0.02', '0.002', '0.0002']:
        for dc in [1.0, 0.3, 0.2, 0.1]:
            
            # Plot num transients
            # p.erass_transients(erass_system_period_cutoff=999999, Z=Z, duty_cycle=dc, percent_transient=False)
            
            # Plot percent transients
            p.erass_transients(erass_system_period_cutoff=999999, Z=Z, duty_cycle=dc, percent_transient=True)

    # p.classifications_dincl_i(percent=False)
    
    
    # Plot all classification hists with duty cycle
    # rp.sim_classifications_sampler_with_duty_cycle()
    
    
    # Plot ADT Histograms
    
    # for d in [1.0, 0.3, 0.2, 0.1]:
    #     print(d)
    #     if d != 1.0:
    #         rp.table_load_classifications()
    #         rp.sim_classifications_sampler_with_duty_cycle(d)
    #         rp.df_classifications['lc_classification'] = rp.df_classifications['lc_classification_new']
            
    #     for Z in ['0.02', '0.002', '0.0002']:
    #         for dincl in [46, 21]:
    #             key = (500, 0.0, dincl, Z, 999999)
    #             p.classifications_hist(key, duty_cycle=d, frac_visible=False, normdist=True)
    #             p.classifications_hist(key, duty_cycle=d, frac_visible=True, normdist=True)
