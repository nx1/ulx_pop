# -*- coding: utf-8 -*-
"""
Created on Tue May 19 12:16:45 2020

@author: norma

See:
https://universeathome.pl/universe/bhdb.php

The second dataset is enormous and not fully provided here, instead a subset
is provided of all the systems accross all metallicities that are undergoing
mass transfer 'mt = 1', this provides 9,372,542 rows to work with.

4,624,270 of these rows are ULXs (88,748 binaries)
4,161,301 of these rows have opening angles < 45 deg ()

We sample from the set of unique binaries and weight them by the total
amount of time that they are a ULX.

"""
import logging
from pathlib import Path

import numpy as np
import pandas as pd

import matplotlib.pyplot as plt
from array import array

from constants import G, c, Myr, R_sol, M_sol, L_sol, sigma, NS_SPIN, BH_SPIN, NS_ISCO, BH_ISCO, epsilon, beta, eta


def startrack():
    """
    Old Population.
    """
    systems_df_path = Path('../data/processed/all_systems_df.csv')
    df = pd.read_csv(systems_df_path)
    df = df.drop(['Unnamed: 0'], axis=1)
    return df

def ulx():
    """
    Old Population.
    """
    df = startrack()
    df = df[df['Lx'] > 1E39]
    return df

def ulx_beamed():
    """
    Old Population.
    """
    df = ulx()
    df = df[df['b'] < 1]
    return df


def ulx_beamed_l_45():
    """
    Old Population.
    """
    df = ulx_beamed()
    df = df[df['theta_half_deg'] < 45]
    return df

def ulx_beamed_l_45_P_wind_l_4_years():
    df = ulx_beamed_l_45()
    df = df[df['P_wind_days'] < 4*365]
    return df

def startrack_v2_mt_1_z02(**kwargs):
    systems_df_path = Path('../data/interim/startrack/z02_data1.csv')
    df = pd.read_csv(systems_df_path, **kwargs)
    return df

def startrack_v2_mt_1_z002(**kwargs):
    systems_df_path = Path('../data/interim/startrack/z002_data1.csv')
    df = pd.read_csv(systems_df_path, **kwargs)
    return df

def startrack_v2_mt_1_z0002(**kwargs):
    systems_df_path = Path('../data/interim/startrack/z0002_data1.csv')
    df = pd.read_csv(systems_df_path, **kwargs)
    return df

def startrack_v2_mt_1_all(**kwargs):
    systems_df_path = Path('../data/interim/startrack/data_mt=1.csv')
    df = pd.read_csv(systems_df_path, **kwargs)
    return df

def startrack_v2_mt_1_test_subset(**kwargs):
    systems_df_path = Path('../data/interim/startrack/data_mt=1_test_subset.csv')
    df = pd.read_csv(systems_df_path, **kwargs)
    return df

    
def set_latex_font():
    import matplotlib
    matplotlib.rcParams['mathtext.fontset'] = 'stix'
    matplotlib.rcParams['font.family'] = 'STIXGeneral'

class Population:
    """
    This class contains functions for describing, plotting, sampling and outputting populations.
    """
    def __init__(self, df):
        self.df = df
        self.calc_columns()
        self.calc_sub_populations()
        
        # XLF settings
        self.bin_min = 38
        self.bin_max = 44
        self.bin_width = 0.25
        self.bins = np.arange(self.bin_min, self.bin_max, self.bin_width)
        self.nbins = len(self.bins)
        self.bin_centers = 0.5 * (self.bins[:-1] + self.bins[1:])
        
        self.df_ulx_Z_subset = 'all'
        
    def calc_columns(self):
        logging.debug('Calculating population columns')
        
        self.df['is_bh'] = np.where(self.df['K_a'] == 14, 1, 0)
        
        # Make mass accretion rates positive
        self.df['dMmt_a'] = abs(self.df['dMmt_a'])
        self.df['dMmt_b'] = abs(self.df['dMmt_b'])
    
        #Convert to grams per second
        self.df['mdot_gs_a'] = self.df['dMmt_a'] * (M_sol/Myr)
        self.df['mdot_gs_b'] = self.df['dMmt_b'] * (M_sol/Myr)
    
        #Calculate Eddington luminosity and Eddington accretion rate
        self.df['LEdd'] = 1.3E38 * self.df['M_a']
        self.df['mdot_Edd'] = self.df['LEdd'] / (eta * c**2)
    
        self.df['mdot_ratio'] = self.df['mdot_gs_b'] / self.df['mdot_Edd']
    
        self.df['Lx_iso'] = np.where(self.df['mdot_ratio'] > 1, self.df['LEdd'] * (1 + np.log(self.df['mdot_ratio'])), self.df['LEdd'] * self.df['mdot_ratio'])
        self.df['log_Lx_iso'] = np.log10(self.df['Lx_iso'])
        #self.df['Lx_iso'] = self.df['LEdd'] * (1 + np.log(self.df['mdot_ratio']))
    
        # Beaming factor
        # b = 1 for m_dot < 8.5
        # b = 73/mdot**2 >= 8.5
        # b = 3.2E-3 for >150
        self.df['b'] = np.where(self.df['mdot_ratio'] >= 8.5, 73/self.df['mdot_ratio']**2, 1)
        self.df['b'] = np.where(self.df['mdot_ratio'] >= 150, 3.2E-3, self.df['b'])
    
        # Observed Luminosity
        self.df['Lx1'] = self.df['Lx_iso']/self.df['b']
        self.df['log_Lx1'] = np.log10(self.df['Lx1'])
    
        self.df['theta'] = 2 * np.arccos(1-self.df['b'])     # Full opening angle in rad
        self.df['theta_deg'] = self.df['theta'] * 180/np.pi  # degrees
        self.df['theta_half_deg'] = self.df['theta_deg'] / 2 # Half opening angle (theta / 2)
    
        self.df['zeta'] = np.tan((np.pi/2) - np.arccos(1 - (73/(self.df['mdot_ratio']**2))))
        self.df['zeta'] = np.where(self.df['zeta'] <= 2, 2, self.df['zeta'])
    
        # General Relativity stuff
        self.df['R_g']           = (G * self.df['M_a']*M_sol) / c**2 # Gravitational radii
        self.df['a*']            = np.where(self.df['M_a']<2.5, NS_SPIN, BH_SPIN)
    
        self.df['r_schw']        = ((2 * G * self.df['M_a'] * M_sol) / c**2) / self.df['R_g']
        self.df['r_isco_nospin'] = ((6 * G * self.df['M_a'] * M_sol) / c**2) / self.df['R_g']
        self.df['r_isco']        = np.where(self.df['M_a']<2.5, NS_ISCO, BH_ISCO)
    
        self.df['r_sph'] = self.df['r_isco'] * self.df['mdot_ratio']
        self.df['r_out'] = 3 * epsilon / (beta * self.df['zeta']) * self.df['mdot_ratio']**3/2 * 6 # The 6 R_g makes the calculation correct for eta = 1/12
    
        self.df['P_inflow_at_rsph'] = (G * self.df['M_a'] * M_sol * np.pi) / (3 * c**3 * self.df['a*']) * self.df['r_sph']**3 * ((1 - (self.df['r_isco']/self.df['r_sph']))**3)/(np.log(self.df['r_sph']/self.df['r_isco']))
        self.df['P_envelope']       = (G * self.df['M_a'] * M_sol * np.pi) / (3 * c**3 * self.df['a*']) * self.df['r_sph']**3 * ((1 - (self.df['r_isco']/self.df['r_sph']))**3)/(np.log(self.df['r_sph']/self.df['r_isco'])) * (self.df['r_out']/self.df['r_sph'])**2
        self.df['P_wind']           = (G * self.df['M_a'] * M_sol * np.pi) / (3 * c**3 * self.df['a*']) * self.df['r_out']**3 * ((1 - (self.df['r_isco']/self.df['r_out']))**3)/(np.log(self.df['r_out']/self.df['r_isco']))
    
        self.df['P_inflow_days']   = self.df['P_inflow_at_rsph'] / (24*60*60)
        self.df['P_envelope_days'] = self.df['P_envelope'] / (24*60*60)
        self.df['P_wind_days']     = self.df['P_wind'] / (24*60*60)
    
        self.df['P_orb'] = 2 * np.pi * np.sqrt((R_sol*self.df['a'])**3 / (G*M_sol*(self.df['M_b'] + self.df['M_a'])) )
        self.df['P_sup'] = 22.1 * self.df['P_orb']
        self.df['P_sup_err'] = self.df['P_sup'] * np.sqrt((0.1/22.1)**2)
        # In days
        self.df['P_orb_days'] = self.df['P_orb'] / (60*60*24)
        self.df['P_sup_days'] = self.df['P_sup'] / (60*60*24)
        self.df['P_sup_err_days'] = self.df['P_sup_err'] / (60*60*24)
        
        # Effective temp star a from Stefan-Boltzmann Law
        self.df['T_eff_b'] = ((self.df['L_b']*L_sol) / (4*np.pi*(self.df['R_b']*R_sol)**2 * sigma))**(1/4)
        
        # LMXRB nuclear mt, and dominated by disc accretion
        self.df['lmxrb'] = np.where((self.df['mttype'] == 1) &
                                    (self.df['dMmt_b'] > self.df['dMwind_b']) &
                                    (self.df['T_eff_b'] < 7000) &
                                    (self.df['M_b'] < 5), 1, 0)
        
    def calc_sub_populations(self):
        logging.debug('Calculating subpopulation')
        self.df_ulx = self.df[self.df['Lx1'] > 1e39]
        self.df_ulx = self.df_ulx[self.df_ulx['mttype'] != 5]   # Exclude WD mass transfer systems
        #self.df_ulx_opening_angle_le_45 = self.df_ulx[self.df_ulx['theta_half_deg'] <= 45]
        #self.df_ulx_P_wind_l_4_years = self.df_ulx_opening_angle_le_45[self.df_ulx_opening_angle_le_45['P_wind_days'] < 365*4]
        #self.df_ulx_P_sup_l_4_years = self.df_ulx_opening_angle_le_45[self.df_ulx_opening_angle_le_45['P_sup_days'] < 365*4]
        
    def filter_df_ulx_by_Z(self, Z):
        if self.df_ulx_Z_subset != Z:
            logging.debug('Population not filtered by Z')
            logging.debug('Filtering population by Z')
            self.df_ulx = self.df_ulx[self.df_ulx['Z'] == float(Z)]
            self.df_ulx_Z_subset = Z
            self.pop.calc_bh_ns_ulx_sub_populations()
        else:
            logging.debug('Population alread filtered by Z')
            logging.debug('Not filtering')
        
    
    def split_ns_bh(self, df):
        logging.debug('Splitting df into ns and bh')
        df_ns = df[df['K_a'] == 13]
        df_bh = df[df['K_a'] == 14]
        return df_ns, df_bh
    
    def calc_bh_ns_ulx_sub_populations(self):
        logging.debug('Calculating df_ulx_bh and df_ulx_ns')
        self.df_ulx_ns, self.df_ulx_bh = self.split_ns_bh(self.df_ulx)
         
    def gb_sys(self, df):
        """Get pandas GroupBy object for unique binaries"""
        logging.debug('Getting system groupby')
        gb = df.groupby(['idum_run', 'iidd_old'])
        return gb    
        
    
    def get_system_ids(self, df):
        """Get all unique idum_run, iidd_old pairs in a given df.
        returns a list of tuples."""
        logging.debug('Getting system ids')
        gb = self.gb_sys(df)
        ids = list(gb.groups)
        return ids
    
    def get_system(self, idum_run, iidd_old):
        logging.debug('Getting system idum_run: %s iidd_old: %s', idum_run, iidd_old)
        gb = self.gb_sys(self.df)
        df_system = gb.get_group((idum_run, iidd_old))
        return df_system


    def calc_sampling_weights(self, df):
        """
        P_samp = system alive time / all system alive times.
        """
        logging.debug('Calculating sampling weights')
        gb = self.gb_sys(df)
        sys_time = gb[['dt']].sum() # system_on_time
        sys_time['P_samp'] = sys_time['dt'] / sys_time['dt'].sum()
        samp_weights = sys_time['P_samp']
        return samp_weights
    
    def calc_ulx_sampling_weights(self):
        logging.debug('Calculating BH and NS ULX sampling weights')
        try:
            self.df_ulx_bh
        except(AttributeError):
            self.calc_bh_ns_ulx_sub_populations()
        self.ulx_bh_samp_weights = self.calc_sampling_weights(self.df_ulx_bh)
        self.ulx_ns_samp_weights = self.calc_sampling_weights(self.df_ulx_ns)
    
    
    def calc_ulx_binary_dict(self):
        logging.debug('Calculating ulx binary dictionary')
        gb = self.gb_sys(self.df_ulx)
        self.binary_dict = gb.groups
        
    def calc_system_averages(self, df):
        logging.debug('Calculating system_averages')
        gb = self.gb_sys(df)
        df_sys_avg = gb.mean().reset_index()
        return df_sys_avg

    def sample_ulxs(self, bh_ratio, size=500):
        logging.debug('Sampling ulxs bh_ratio=%s size=%s', bh_ratio, size)
        try:
            self.ulx_bh_samp_weights
            self.binary_dict
        except AttributeError:
            self.calc_ulx_sampling_weights()
            self.calc_ulx_binary_dict()
        
        N_bh = int(size*bh_ratio)
        N_ns = size-N_bh
        
        sampled_bh = np.random.choice(self.ulx_bh_samp_weights.index, size=N_bh, p=self.ulx_bh_samp_weights.values)
        sampled_ns = np.random.choice(self.ulx_ns_samp_weights.index, size=N_ns, p=self.ulx_ns_samp_weights.values)
        
        sampled_ulxs_idum_iidd_pairs = np.concatenate((sampled_bh, sampled_ns))
        
        sampled_indexs = np.array([np.random.choice(self.binary_dict[sampled_ulxs_idum_iidd_pairs[i]]) for i in range(size)])
        return sampled_indexs
    

   
    
    def sample_visible(self, df, include_beaming=True, include_duty_cycle=False, lmxrb_duty_cycle=0.1):
        """
        Randomly sample based on observation probability b, and or duty cycle
        and return systems that were visible.
        """
        logging.debug('Sampling visible systems')
        logging.debug('N of input rows: %s', len(df))
        
        dfw = df.copy()
        dfw['bin'] = pd.cut(np.log10(dfw['Lx1']), bins=self.bins)
        if include_beaming:
            dfw['rand'] = np.random.random(size=len(dfw))
            dfw['visible_b'] = dfw['rand'] < dfw['b']
            
        if include_duty_cycle:
            dfw['duty_cycle'] = np.where(dfw['lmxrb']==1, lmxrb_duty_cycle, 1)
            dfw['rand2'] = np.random.random(size=len(dfw))
            dfw['visible_dc'] = dfw['rand2'] < dfw['duty_cycle']
            df_vis = dfw[(dfw['visible_b'] == True) & (dfw['visible_dc'] == True)]
        else:
            df_vis = dfw[(dfw['visible_b'] == True)]
            
        logging.debug('N of visible rows: %s', len(df_vis))
        return df_vis
    
    
    def XLF_mc_histogram(self,
                         df,
                         include_beaming=True,
                         include_duty_cycle=False,
                         lmxrb_duty_cycle=0.1,
                         mc_iters=10000):
        """
        Monte-Carlo the XLF in order to obtain errors due to beaming
        
        Returns
        -------
        hist_results : np.array()
            mc_iters x len(bins) array containing the numbers in each bin.
        """
        
        hist_results = np.ndarray((mc_iters, len(self.bins)-1))
        for i in range(mc_iters):
            df_vis = self.sample_visible(df, include_duty_cycle, lmxrb_duty_cycle)
            if i%100==0:
                print(f'MC XLF | duty_cycle={lmxrb_duty_cycle} include={include_duty_cycle} | iteration: {i}/{mc_iters}')
            h = np.histogram(np.log10(df_vis['Lx1']), bins=self.bins)
            hist_results[i] = h[0]
        return hist_results


    def XLF_plot(self,
                 df,
                 by='NSBH',
                 average_systems=False,
                 include_beamed=False,
                 include_duty_cycle=False,
                 save=False):
        
        fontsize = 10
        
        if average_systems:
            df = self.calc_system_averages(df)
            
        if include_beamed:
            fig, ax = plt.subplots(3, 1, sharex=True,  sharey=True, figsize=(6, 6))
            ax[2].set_xlabel(r'log $L_{x}$ $(\mathrm{erg \ s^{-1}})$', fontsize=fontsize)
            ax[2].axvline(39, c='r')
            ax[2].set_ylabel(r'N', rotation=0, fontsize=fontsize)
        else:
            fig, ax = plt.subplots(2, 1, sharex=True, sharey=True, figsize=(6, 4.5))
        
        ax[0].tick_params(axis='both', labelsize=fontsize)
        ax[0].tick_params(axis='both', labelsize=fontsize)
        ax[0].axvline(39, c='r')
        ax[1].axvline(39, c='r')
    
        ax[0].set_yscale('log')
        ax[1].set_yscale('log')
        ax[2].set_yscale('log')
        
        ax[0].set_xlabel(r'log $L_{iso}$ $(\mathrm{erg \ s^{-1}})$', fontsize=fontsize)
        ax[1].set_xlabel(r'log $L_{x}$ $(\mathrm{erg \ s^{-1}})$', fontsize=fontsize)
        ax[0].set_ylabel(r'N', rotation=0, fontsize=fontsize)
        ax[1].set_ylabel(r'N', rotation=0, fontsize=fontsize)
        

        if by=='NSBH':
            df_ns, df_bh = self.split_ns_bh(df)
            
            np.log10(df_ns['Lx_iso']).hist(ax=ax[0], bins=self.bins, label='NS | Isotropic Emission',  edgecolor='black', histtype='step', alpha=1.0, grid=False)
            np.log10(df_bh['Lx_iso']).hist(ax=ax[0], bins=self.bins, label='BH | Isotropic Emission',  edgecolor='grey', histtype='step', alpha=1.0, grid=False, linestyle='--')
        
            np.log10(df_ns['Lx1']).hist(ax=ax[1], bins=self.bins, label='NS | With Beaming',  edgecolor='black', histtype='step', alpha=1.0, grid=False)
            np.log10(df_bh['Lx1']).hist(ax=ax[1], bins=self.bins, label='BH | With Beaming',  edgecolor='grey', histtype='step', alpha=1.0, grid=False, linestyle='--')
            if include_beamed:
                ns_hist_results = self.XLF_mc_histogram(df_ns, include_duty_cycle)
                bh_hist_results = self.XLF_mc_histogram(df_bh, include_duty_cycle)
            
            df_vis_ns = self.sample_visible(df_ns)
            df_vis_bh = self.sample_visible(df_bh)
            
            mean_ns = np.mean(ns_hist_results, axis=0)
            std_ns  = np.std(ns_hist_results, axis=0)
            
            mean_bh = np.mean(bh_hist_results, axis=0)
            std_bh  = np.std(bh_hist_results, axis=0)
            
            np.log10(df_vis_ns['Lx1']).hist(ax=ax[2], bins=self.bins, label='NS | With Beaming | Visible only',  edgecolor='black', histtype='step', alpha=1.0, grid=False)
            np.log10(df_vis_bh['Lx1']).hist(ax=ax[2], bins=self.bins, label='BH | With Beaming | Visible only',  edgecolor='grey', histtype='step', alpha=1.0, grid=False, linestyle='--')
            
            ax[2].errorbar(self.bin_centers, mean_ns, yerr=std_ns, linestyle='', color='black')
            ax[2].errorbar(self.bin_centers, mean_bh, yerr=std_bh, linestyle='', color='grey')
            
            ax[2].legend(loc='upper left')
        

        elif by=='mttype':
            df_nuc = df[df['mttype'] == 1]
            df_therm = df[df['mttype'] == 4]
            # df_wd = df[df['mttype'] == 5]
            
            df_nuc_sys = self.calc_system_averages(df_nuc)
            df_therm_sys = self.calc_system_averages(df_therm)
            # df_wd_sys = self.calc_system_averages(df_wd)
        
            np.log10(df_nuc_sys['Lx_iso']).hist(ax=ax[0], bins=self.bins, label='Nuclear MT | Isotropic Emission',  edgecolor='black', histtype='stepfilled', fc='green', alpha=0.5, grid=False)
            np.log10(df_therm_sys['Lx_iso']).hist(ax=ax[0], bins=self.bins, label='Thermal MT | Isotropic Emission',  edgecolor='black', histtype='stepfilled', fc='red', alpha=0.5, grid=False)
            # np.log10(df_wd_sys['Lx_iso']).hist(ax=ax[0], bins=self.bins, label='WD MT | Isotropic Emission',  edgecolor='black', histtype='stepfilled', fc='white', alpha=0.5, grid=False)
        
            np.log10(df_nuc_sys['Lx1']).hist(ax=ax[1], bins=self.bins, label='Nuclear MT | With Beaming',  edgecolor='black', histtype='stepfilled', fc='green', alpha=0.5, grid=False)
            np.log10(df_therm_sys['Lx1']).hist(ax=ax[1], bins=self.bins, label='Thermal MT | With Beaming',  edgecolor='black', histtype='stepfilled', fc='red', alpha=0.5, grid=False)
            # np.log10(df_wd_sys['Lx1']).hist(ax=ax[1], bins=self.bins, label='WD MT | With Beaming',  edgecolor='black', histtype='stepfilled', fc='white', alpha=0.5, grid=False)
            
        elif by=='both':
            df_ns, df_bh = self.split_ns_bh(df)
            
            df_nuc_bh   = df_bh[df_bh['mttype'] == 1]
            df_therm_bh = df_bh[df_bh['mttype'] == 4]
            df_nuc_ns   = df_ns[df_ns['mttype'] == 1]
            df_therm_ns = df_ns[df_ns['mttype'] == 4]
            
            df_nuc_bh_sys = self.calc_system_averages(df_nuc_bh)
            df_therm_bh_sys = self.calc_system_averages(df_therm_bh)
            df_nuc_ns_sys = self.calc_system_averages(df_nuc_ns)
            df_therm_ns_sys = self.calc_system_averages(df_therm_ns)
        
            np.log10(df_nuc_bh_sys['Lx_iso']).hist(ax=ax[0], bins=self.bins, label='BH | Nuclear MT | Isotropic Emission',  edgecolor='black', histtype='stepfilled', fc='cyan', alpha=0.5, grid=False)
            np.log10(df_therm_bh_sys['Lx_iso']).hist(ax=ax[0], bins=self.bins, label='BH | Thermal MT | Isotropic Emission',  edgecolor='black', histtype='stepfilled', fc='purple', alpha=0.5, grid=False)
            np.log10(df_nuc_ns_sys['Lx_iso']).hist(ax=ax[0], bins=self.bins, label='NS | Nuclear MT | Isotropic Emission',  edgecolor='black', histtype='stepfilled', fc='green', alpha=0.5, grid=False)
            np.log10(df_therm_ns_sys['Lx_iso']).hist(ax=ax[0], bins=self.bins, label='NS | Thermal MT | Isotropic Emission',  edgecolor='black', histtype='stepfilled', fc='red', alpha=0.5, grid=False)
        
            np.log10(df_nuc_bh_sys['Lx1']).hist(ax=ax[1], bins=self.bins, label='BH | Nuclear MT | With Beaming',  edgecolor='black', histtype='stepfilled', fc='cyan', alpha=0.5, grid=False)
            np.log10(df_therm_bh_sys['Lx1']).hist(ax=ax[1], bins=self.bins, label='BH | Thermal MT | With Beaming',  edgecolor='black', histtype='stepfilled', fc='purple', alpha=0.5, grid=False)
            np.log10(df_nuc_ns_sys['Lx1']).hist(ax=ax[1], bins=self.bins, label='NS | Nuclear MT | With Beaming',  edgecolor='black', histtype='stepfilled', fc='green', alpha=0.5, grid=False)
            np.log10(df_therm_ns_sys['Lx1']).hist(ax=ax[1], bins=self.bins, label='NS | Thermal MT | With Beaming',  edgecolor='black', histtype='stepfilled', fc='red', alpha=0.5, grid=False)
        
        ax[0].legend(fontsize=fontsize, loc='upper right')
        ax[1].legend(fontsize=fontsize, loc='upper left')
        
        plt.tight_layout()
        if save==True:
            if by=='NSBH':
                plt.savefig('../reports/figures/XLF_by_BHNS.png', dpi=500)
                plt.savefig('../reports/figures/XLF_by_BHNS.eps')
                plt.savefig('../reports/figures/XLF_by_BHNS.pdf')
            elif by=='mttype':
                    plt.savefig('../reports/figures/XLF_by_MT.png', dpi=500)
                    plt.savefig('../reports/figures/XLF_by_MT.eps')
                    plt.savefig('../reports/figures/XLF_by_MT.pdf')
            elif by=='both':
                    plt.savefig('../reports/figures/XLF_by_MT_NSBH.png', dpi=500)
                    plt.savefig('../reports/figures/XLF_by_MT_NSBH.eps')
                    plt.savefig('../reports/figures/XLF_by_MT_NSBH.pdf')
                    
    def XLF_plot_by_time(self):
        times = [(0,100),(100,500), (500,1000), (1000,2000), (2000,8000)]
        for i in range(len(times)):
    
            low = times[i][0]
            high = times[i][1]
            print(i, high, low)
            cut = self.df[(pop.df['t']>low) & (pop.df['t']<high)]
            cut = self.calc_system_averages(cut)
            plt.yscale('log')
            plt.hist(np.log10(cut['Lx1']), bins=self.bins, label=f'{low} - {high} Myr', histtype='step')
            
        plt.legend()
    

    def pivot_binaries_count(self, df):
        piv = pd.pivot_table(df, columns=['is_bh'], index=['Z'], aggfunc='count', margins=True, margins_name='total')['theta_half_deg']
        return piv


    def plot_mass_accretion_rates(self, df, save=False):
        df_ns = df[df['K_a'] == 13]
        df_bh = df[df['K_a'] == 14]
        
        fig, ax = plt.subplots(3,1)
        ax[0].hist(np.log10(df_ns['dMmt_b']), label='NS', histtype='step', bins=50)
        ax[0].hist(np.log10(df_bh['dMmt_b']), label='BH', histtype='step', bins=50, linestyle='--')
        ax[0].legend()
        ax[0].set_xlabel(r'log10(dMmt_b) ($M_{\odot} \ Myr^{-1}$)')
        plt.tight_layout()
        
        ax[1].hist(np.log10(df_ns['mdot_gs_b']), label='NS', histtype='step', bins=50)
        ax[1].hist(np.log10(df_bh['mdot_gs_b']), label='BH', histtype='step', bins=50, linestyle='--')
        ax[1].legend()
        ax[1].set_xlabel(r'log10(mdot_gs_b) ($g \ s^{-1}$)')
    
        ax[2].hist(np.log10(df_ns['mdot_ratio']), label='NS', histtype='step', bins=50)
        ax[2].hist(np.log10(df_bh['mdot_ratio']), label='BH', histtype='step', bins=50, linestyle='--')
        ax[2].legend()
        ax[2].set_xlabel(r'log10($\dot{m}_{0}$) (Dimensionless)')
        plt.tight_layout()
    
        plt.tight_layout()
        if save:
            plt.savefig('../reports/figures/mass_accretion_rates_by_row.png', dpi=500)
            plt.savefig('../reports/figures/mass_accretion_rates_by_row.eps')
            plt.savefig('../reports/figures/mass_accretion_rates_by_row.pdf')

    def plot_system_luminosity_evolution(self, idum_run, iidd_old):
        df_system = self.get_system(idum_run, iidd_old)
        plt.figure()
        plt.title(f'System: {idum_run} {iidd_old}')
        plt.scatter(df_system['t'], np.log10(df_system['Lx1']))
        plt.axhline(39, c='r') 
        plt.xlabel('Time (Myr)')
        plt.ylabel('log10 (Lx)')
        plt.show()
    
    
    def export_to_binary_data(self, df):
        """
        Create binary datafiles used for importing into C code.
        """
        N = len(df)
        system_id = array('i', df.index)
        system_theta = array('d', df.theta_half_deg)
        system_Lx = array('d', df.Lx1 / 1e39)    # In units of x10^39 erg/s
        system_P_wind_days = array('d', df.P_wind_days)
        system_P_sup_days = array('d', df.P_sup_days)
        
        
        def save_to_binary_file(array, param_name):
            """
            The tofile() function takes one argument which is known as the file object
            The file object is the object that is created when you call open()

            """
            dtype = array.typecode
            N = len(array)
            f = open(f'data_{param_name}_N={N}_{dtype}.bin', 'wb+')
            array.tofile(f)
            f.close()
            
        save_to_binary_file(system_id, 'id')
        save_to_binary_file(system_theta, 'theta')
        save_to_binary_file(system_Lx, 'Lx')
        save_to_binary_file(system_P_wind_days, 'P_wind_days')
        save_to_binary_file(system_P_sup_days, 'P_sup_days')
    
    def export_df(self, df, identifier):
        df.to_csv(f'{identifier}.csv', index=False)
        
    
    def export_for_ulxlc(self, df, identifier):
        self.describe()
        self.export_population_info(df, identifier)
        self.export_df(df, identifier)
        self.export_to_binary_data(df)
        
    def describe(self, df, df_name):        
        df_bh = df[df['K_a'] == 14]
        df_ns = df[df['K_a'] == 13]
        
        df_nuc = df[df['mttype'] == 1]
        df_therm = df[df['mttype'] == 4]
        df_wd = df[df['mttype'] == 5]
        
        df_lmxrb = df[df['lmxrb']==1]
        df_beamed = df[df['b']>=1]
        df_unbeamed = df[df['b']<1]
        
        N_rows = len(df)
        N_rows_bh = len(df_bh)
        N_rows_ns = len(df_ns)
        N_rows_nuc = len(df_nuc)
        N_rows_therm = len(df_therm)
        N_rows_wd = len(df_wd)
        N_rows_lmxrb = len(df_lmxrb)
        N_rows_unbeamed = len(df_unbeamed)
        N_rows_beamed = len(df_beamed)
        
        
        N_binaries = len(df.groupby(['idum_run', 'iidd_old']))
        N_binaries_bh = len(df_bh.groupby(['idum_run', 'iidd_old']))
        N_binaries_ns = len(df_ns.groupby(['idum_run', 'iidd_old']))
        N_binaries_nuc = len(df_nuc.groupby(['idum_run', 'iidd_old']))
        N_binaries_therm = len(df_therm.groupby(['idum_run', 'iidd_old']))
        N_binaries_wd = len(df_wd.groupby(['idum_run', 'iidd_old']))
        N_binaries_lmxrb = len(df_lmxrb.groupby(['idum_run', 'iidd_old']))
        N_binaries_unbeamed = len(df_unbeamed.groupby(['idum_run', 'iidd_old']))
        N_binaries_beamed = len(df_beamed.groupby(['idum_run', 'iidd_old']))
        
        print(f'Information for {df_name}')
        print('==========================')
        print(f'N rows: {N_rows} ({N_rows/N_rows*100:0.2f}%)')
        print(f'N BH rows: {N_rows_bh} ({N_rows_bh/N_rows*100:0.2f}%)')
        print(f'N NS rows: {N_rows_ns} ({N_rows_ns/N_rows*100:0.2f}%)')
        print(f'N nuclear mt rows: {N_rows_nuc} ({N_rows_nuc/N_rows*100:0.2f}%)')
        print(f'N thermal mt rows: {N_rows_therm} ({N_rows_therm/N_rows*100:0.2f}%)')
        print(f'N WD mt rows: {N_rows_wd} ({N_rows_wd/N_rows*100:0.2f}%)')
        print(f'N lmxrb rows: {N_rows_lmxrb} ({N_rows_lmxrb/N_rows*100:0.2f}%)')
        print(f'N b>=1 rows: {N_rows_unbeamed} ({N_rows_unbeamed/N_rows*100:0.2f}%)')
        print(f'N b<1 rows: {N_rows_beamed} ({N_rows_beamed/N_rows*100:0.2f}%)')
    
        print('--------------------------')
        print(f'N binaries: {N_binaries} ({N_binaries/N_binaries*100:0.2f}%)')
        print(f'N BH binaries: {N_binaries_bh} ({N_binaries_bh/N_binaries*100:0.2f}%)')
        print(f'N NS binaries: {N_binaries_ns} ({N_binaries_ns/N_binaries*100:0.2f}%)')
        print(f'N nuclear mt binaries: {N_binaries_nuc} ({N_binaries_nuc/N_binaries*100:0.2f}%)')
        print(f'N thermal mt binaries: {N_binaries_therm} ({N_binaries_therm/N_binaries*100:0.2f}%)')
        print(f'N WD mt binaries: {N_binaries_wd} ({N_binaries_wd/N_binaries*100:0.2f}%)')
        print(f'N lmxrb binaries: {N_binaries_lmxrb} ({N_binaries_lmxrb/N_binaries*100:0.2f}%)')
        print(f'N b>=1 rows binaries: {N_binaries_unbeamed} ({N_binaries_unbeamed/N_binaries*100:0.2f}%)')
        print(f'N b<1 binaries: {N_binaries_beamed} ({N_binaries_beamed/N_binaries*100:0.2f}%)')
        print('--------------------------')

if __name__ == "__main__":
    logging.basicConfig(level=logging.DEBUG)
    
    df = startrack_v2_mt_1_all() #nrows=100000
    # df = startrack_v2_mt_1_test_subset()
    
    pop = Population(df)
    
    # df_sys = pop.calc_system_averages(pop.df)
    
    # Population Statistics
    # pop.describe(pop.df, 'df')
    # pop.describe(pop.df_ulx, 'df_ulx')
    
    # Plotting
    set_latex_font()
    
    # pop.XLF_plot(df=pop.df,
    #              by='NSBH',
    #              average_systems=True,
    #              include_beamed=True,
    #              include_duty_cycle=True,
    #              save=False)


    # =============================================================================
    # MONTE CARLO XLF
    # =============================================================================


    # def sample_l_iso_with_bh_ns(df):
    #     df_samp = df[['K_a', 'log_Lx_iso']].sample(n=5000, replace=True)
    #     L = df_samp['log_Lx_iso'].values
    #     K_a = df_samp['K_a'].values
    #     return np.array([L, K_a])

    # def sample_with_bh_ns(df):
    #     df_samp = df[['K_a', 'log_Lx1']].sample(n=5000, replace=True)
    #     L = df_samp['log_Lx1'].values
    #     K_a = df_samp['K_a'].values
    #     return np.array([L, K_a])
    
    # def sample_with_beaming_and_duty_cycle(df):
    #     """For use using GroupBy.apply(sample_with_beaming_and_duty_cycle)"""
    #     duty_cycle = 0.3
    #     df = df[['K_a', 'log_Lx1', 'lmxrb', 'b']]
    #     N_sys = len(df)
    #     df['dc'] = np.where(df['lmxrb']==1, duty_cycle, 1.0)
    #     df['p_obs'] = df['b'] * df['dc']    # Observation probability
    #     df['1-p_obs'] = 1 - df['p_obs']     # Non-Observation prob
        
    #     # print(df)
    #     df_2 = pd.DataFrame()
    #     df_2['p_obs'] = df['1-p_obs']
    #     df_2['log_Lx1'] = np.zeros(shape=N_sys) # L=0 for non obs
    #     df_2['K_a'] = df['K_a']
    #     # print(df_2)
        
    #     df_c = pd.concat([df, df_2], axis=0)
    #     df_c['p_obs'] = df_c['p_obs']/N_sys
    #     df_c['p_obs'] = df_c['p_obs'].abs()
    #     df_c['p_obs'] /= df_c['p_obs'].sum()
    #     # print(df_c)
        
    #     df_samp = df_c[['K_a', 'log_Lx1']].sample(n=5000, replace=True, weights=df_c['p_obs'])
    #     L = df_samp['log_Lx1'].values
    #     K_a = df_samp['K_a'].values
    #     return np.array([L, K_a])
        



    # gb = pop.gb_sys(pop.df)
    # samp = gb.apply(sample_with_beaming_and_duty_cycle)
    
    # def samp_2_means_std(samp):
    #     val = np.stack(samp.values, axis=0) 
        
    #     val_L = val[:,0,:]
    #     val_K = val[:,1,:]
        
    #     ns_mask = ~(val_K==13)
    #     bh_mask = ~(val_K==14)
        
    #     ns_L_ma = np.ma.masked_array(data=val_L, mask=ns_mask)
    #     bh_L_ma = np.ma.masked_array(data=val_L, mask=bh_mask)
        
    #     ns_L = ns_L_ma.filled(fill_value=0)
    #     bh_L = bh_L_ma.filled(fill_value=0)
        
    #     ns_L = ns_L.T
    #     bh_L = bh_L.T
        
    #     ns_hists = [np.histogram(ns_L[i], bins=pop.bins)[0] for i in range(len(ns_L))]
    #     bh_hists = [np.histogram(bh_L[i], bins=pop.bins)[0] for i in range(len(bh_L))]
        
    #     # ns_hists = np.where(ns_hists==0, np.nan, ns_hists)
    #     # bh_hists = np.where(bh_hists==0, np.nan, bh_hists)
        
    #     bh_means = np.nanmean(bh_hists, axis=0)
    #     bh_stds  = np.nanstd(bh_hists, axis=0)
        
    #     ns_means = np.nanmean(ns_hists, axis=0)
    #     ns_stds  = np.nanstd(ns_hists, axis=0)
    #     return bh_means, bh_stds, ns_means, ns_stds
    
    
    
    # from tqdm import tqdm
    
    # tqdm.pandas()
    
    # samp = gb.progress_apply(sample_with_bh_ns)
    # bh_means, bh_stds, ns_means, ns_stds = samp_2_means_std(samp)
    # del(samp)
    
    # samp = gb.progress_apply(sample_with_beaming_and_duty_cycle)
    # bh_means_obs, bh_stds_obs, ns_means_obs, ns_stds_obs = samp_2_means_std(samp)
    # del(samp)
    
    # samp = gb.progress_apply(sample_l_iso_with_bh_ns)
    # bh_means_iso, bh_stds_iso, ns_means_iso, ns_stds_iso = samp_2_means_std(samp)
    
    
    
    # fig, ax = plt.subplots(1,3, sharey=True, figsize=(10,4))
    
    # ax[0].errorbar(x=pop.bin_centers, y=bh_means_iso, yerr=bh_stds_iso, linestyle='-', capsize=1.0, label='BH | isotropic emission', c='black')
    # ax[0].errorbar(x=pop.bin_centers, y=ns_means_iso, yerr=ns_stds_iso, linestyle='--', capsize=1.0, label='NS | isotropic emission', c='grey')
    # ax[0].set_xlabel(r'log $L_{iso}$ $(\mathrm{erg \ s^{-1}})$', fontsize=10)
    # ax[0].axvline(39, c='r', linestyle='--')
    # ax[0].legend(loc='lower right', fontsize=10)
    
    
    # ax[1].errorbar(x=pop.bin_centers, y=bh_means, yerr=bh_stds, linestyle='-', capsize=1.0, label='BH | with beaming', c='black')
    # ax[1].errorbar(x=pop.bin_centers, y=ns_means, yerr=ns_stds, linestyle='--', capsize=1.0, label='NS | with beaming', c='grey')
    # ax[1].set_xlabel(r'log $L_{x}$ $(\mathrm{erg \ s^{-1}})$', fontsize=10)
    # ax[1].axvline(39, c='r', linestyle='--')
    # ax[1].legend(loc='lower right', fontsize=10)
    
    
    # ax[2].errorbar(x=pop.bin_centers, y=bh_means_obs, yerr=bh_stds_obs, linestyle='-', capsize=1.0, label='BH | observed | d = 0.3', c='black')
    # ax[2].errorbar(x=pop.bin_centers, y=ns_means_obs, yerr=ns_stds_obs, linestyle='--', capsize=1.0, label='NS | observed | d = 0.3', c='grey')
    # ax[2].set_xlabel(r'log $L_{x}$ $(\mathrm{erg \ s^{-1}})$', fontsize=10)
    # ax[2].axvline(39, c='r', linestyle='--')
    # ax[2].set_yscale('log')
    # ax[2].legend(loc='lower right', fontsize=10)
    # plt.tight_layout()
    # plt.subplots_adjust(wspace=0)
    
    # plt.savefig('../reports/figures/XLF_by_bh_ns_samp_5000.png', dpi=500)
    # plt.savefig('../reports/figures/XLF_by_bh_ns_samp_5000.eps')
    # plt.savefig('../reports/figures/XLF_by_bh_ns_samp_5000.pdf')
    
    
    
    


