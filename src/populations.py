# -*- coding: utf-8 -*-
"""
Created on Tue May 19 12:16:45 2020

@author: norma

This file contains two datasets ontained from startrack.
The first an old dataset containing #36k rows which we have ommited using as
the time sampling was not performed uniformly.

The second dataset is enormous and not fully provided here, instead a subset
is provided of all the systems accross all metallicities that are undergoing
mass transfer 'mt = 1', this provides 9,372,542 rows to work with.

4,624,270 of these rows are ULXs (88,748 binaries)
4,161,301 of these rows have opening angles < 45 deg ()

We sample from the set of unique binaries and weight them by the total
amount of time that they are a ULX.
    
# I feel like I have discovered forbidden mathematics
    print(sys_df['t'].diff().sum())
    print(sys_df['t'].max() - sys_df['t'].min())

"""
from pathlib import Path

import numpy as np
import pandas as pd

import matplotlib.pyplot as plt
from array import array

from constants import G, c, Myr, R_sol, M_sol, NS_SPIN, BH_SPIN, NS_ISCO, BH_ISCO, epsilon, beta, LMXRB_MASS, eta


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
    df = pd.read_csv(systems_df_path, index_col=0, **kwargs)
    return df
    
    
class Population:
    """
    This class contains functions for describing, plotting, sampling and outputting populations.
    """
    def __init__(self, df):
        self.df = df
        self.calc_columns()
        self.calc_sub_populations()
        
        # XLF settings
        self.bin_min = 34
        self.bin_max = 44
        self.bin_width = 0.1
        self.bins = np.arange(self.bin_min, self.bin_max, self.bin_width)
        self.bin_centers = 0.5 * (self.bins[:-1] + self.bins[1:])
        
        self.df_ulx_Z_subset = None
        
    def calc_columns(self):   
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
        #self.df['Lx_iso'] = self.df['LEdd'] * (1 + np.log(self.df['mdot_ratio']))
    
        # Beaming factor
        # b = 1 for m_dot < 8.5
        # b = 73/mdot**2 >= 8.5
        # b = 3.2E-3 for >150
        self.df['b'] = np.where(self.df['mdot_ratio'] >= 8.5, 73/self.df['mdot_ratio']**2, 1)
        self.df['b'] = np.where(self.df['mdot_ratio'] >= 150, 3.2E-3, self.df['b'])
    
        # Observed Luminosity
        self.df['Lx1'] = self.df['Lx_iso']/self.df['b']
    
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
        
        # LMXRB defined as m_b < 1.5, nuclear mt, and dominated by disc accretion.
        self.df['lmxrb'] = np.where((self.df['M_b'] < LMXRB_MASS) & (self.df['mttype'] == 1) & (self.df['dMmt_b'] > self.df['dMwind_b']), 1, 0)
        
    def calc_sub_populations(self):
        self.df_ulx = self.df[self.df['Lx1'] > 1e39]
        self.df_ulx = self.df_ulx[self.df_ulx['mttype'] != 5]
        #self.df_ulx_opening_angle_le_45 = self.df_ulx[self.df_ulx['theta_half_deg'] <= 45]
        #self.df_ulx_P_wind_l_4_years = self.df_ulx_opening_angle_le_45[self.df_ulx_opening_angle_le_45['P_wind_days'] < 365*4]
        #self.df_ulx_P_sup_l_4_years = self.df_ulx_opening_angle_le_45[self.df_ulx_opening_angle_le_45['P_sup_days'] < 365*4]
    
    def filter_df_ulx_by_Z(self, Z):
        self.df_ulx = self.df_ulx[self.df_ulx['Z'] == float(Z)]
        self.df_ulx_Z_subset = float(Z)
    
    def calc_bh_ns_ulx_sub_populations(self):
        self.df_ulx_bh = self.df_ulx[self.df_ulx['is_bh'] == 1]
        self.df_ulx_ns = self.df_ulx[self.df_ulx['is_bh'] == 0]
       
    
    def get_system(self, idum_run, iidd_old):
        df_system = self.df[self.df['idum_run']==idum_run]
        df_system = df_system[df_system['iidd_old']==iidd_old]
        return df_system
    
    def get_unique_system_counts(self, df):
        df_unique_systems_count = df.groupby(['idum_run', 'iidd_old']).size().reset_index().rename(columns={0:'count'})
        df_unique_systems_count = df_unique_systems_count.sort_values('count', ascending=False)
        return df_unique_systems_count


    def calc_sampling_weights(self, df):
        gb = df.groupby(['idum_run', 'iidd_old'])[['t']].agg(['min', 'max']) # Get min and max t for each system
        gb[('t', 'range')] = gb[('t', 'max')] - gb[('t', 'min')] #Calculate t_range (i.e system alive time)
        t_tot = gb[('t', 'range')].sum() #Calculate total time alive for all systems
        gb['P_samp'] = gb[('t', 'range')] / t_tot
        samp_weights = gb['P_samp']
        return samp_weights
    
    def calc_ulx_sampling_weights(self):
        try:
            self.df_ulx_bh
        except(AttributeError):
            self.calc_bh_ns_ulx_sub_populations()
        self.ulx_bh_samp_weights = self.calc_sampling_weights(self.df_ulx_bh)
        self.ulx_ns_samp_weights = self.calc_sampling_weights(self.df_ulx_ns)
    
    
    def calc_ulx_binary_dict(self):
        self.binary_dict = self.df_ulx.groupby(['idum_run', 'iidd_old']).groups
    
    def sample_ulxs(self, bh_ratio, size=500):
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
        
    
    def sample_visible(self, df, include_duty_cycle=False, lmxrb_duty_cycle=0.1):
        """
        Randomly sample based on observation probability b and return systems
        that were visible.

        """        
        dfw = df.copy()
        dfw['bin'] = pd.cut(np.log10(dfw['Lx1']), bins=self.bins)
        dfw['rand'] = np.random.random(size=len(dfw))
        dfw['visible_b'] = dfw['rand'] < dfw['b']
        
        if include_duty_cycle:
            dfw['duty_cycle'] = np.where(dfw['lmxrb']==1, lmxrb_duty_cycle, 1)
            dfw['rand2'] = np.random.random(size=len(dfw))
            dfw['visible_dc'] = dfw['rand2'] < dfw['duty_cycle']
            df_vis = dfw[(dfw['visible_b'] == True) & (dfw['visible_dc'] == True)]
            return df_vis
        else:
            df_vis = dfw[(dfw['visible_b'] == True)]
            return df_vis
    
    def XLF_mc_histogram(self, df, include_duty_cycle=False, lmxrb_duty_cycle=0.1, mc_iters=1000):
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
                print(f'MC XLF | iteration: {i}/{mc_iters}')
            h = np.histogram(np.log10(df_vis['Lx1']), bins=self.bins)
            hist_results[i] = h[0]
        return hist_results
    
    def set_latex_font(self):
        import matplotlib
        matplotlib.rcParams['mathtext.fontset'] = 'stix'
        matplotlib.rcParams['font.family'] = 'STIXGeneral'
        
    
    def XLF_plot(self, df, by='NSBH', save=False, include_beamed=False, include_duty_cycle=False):
        fontsize = 10
        
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
            df_ns = df[df['K_a'] == 13]
            df_bh = df[df['K_a'] == 14]
        
            df_ns_sys = df_ns.groupby(['idum_run', 'iidd_old']).mean().reset_index()
            df_bh_sys = df_bh.groupby(['idum_run', 'iidd_old']).mean().reset_index()
            
            np.log10(df_ns_sys['Lx_iso']).hist(ax=ax[0], bins=self.bins, label='NS | Isotropic Emission',  edgecolor='black', histtype='step', alpha=1.0, grid=False)
            np.log10(df_bh_sys['Lx_iso']).hist(ax=ax[0], bins=self.bins, label='BH | Isotropic Emission',  edgecolor='grey', histtype='step', alpha=1.0, grid=False, linestyle='--')
        
            np.log10(df_ns_sys['Lx1']).hist(ax=ax[1], bins=self.bins, label='NS | With Beaming',  edgecolor='black', histtype='step', alpha=1.0, grid=False)
            np.log10(df_bh_sys['Lx1']).hist(ax=ax[1], bins=self.bins, label='BH | With Beaming',  edgecolor='grey', histtype='step', alpha=1.0, grid=False, linestyle='--')
            if include_beamed:
                ns_hist_results = self.XLF_mc_histogram(df_ns_sys, include_duty_cycle)
                bh_hist_results = self.XLF_mc_histogram(df_bh_sys, include_duty_cycle)
            
                df_vis_ns = self.sample_visible(df_ns_sys)
                df_vis_bh = self.sample_visible(df_bh_sys)
                
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
            
            df_nuc_sys = df_nuc.groupby(['idum_run', 'iidd_old']).mean().reset_index()
            df_therm_sys = df_therm.groupby(['idum_run', 'iidd_old']).mean().reset_index()
            # df_wd_sys = df_wd.groupby(['idum_run', 'iidd_old']).mean().reset_index()
        
            np.log10(df_nuc_sys['Lx_iso']).hist(ax=ax[0], bins=self.bins, label='Nuclear MT | Isotropic Emission',  edgecolor='black', histtype='stepfilled', fc='green', alpha=0.5, grid=False)
            np.log10(df_therm_sys['Lx_iso']).hist(ax=ax[0], bins=self.bins, label='Thermal MT | Isotropic Emission',  edgecolor='black', histtype='stepfilled', fc='red', alpha=0.5, grid=False)
            # np.log10(df_wd_sys['Lx_iso']).hist(ax=ax[0], bins=self.bins, label='WD MT | Isotropic Emission',  edgecolor='black', histtype='stepfilled', fc='white', alpha=0.5, grid=False)
        
            np.log10(df_nuc_sys['Lx1']).hist(ax=ax[1], bins=self.bins, label='Nuclear MT | With Beaming',  edgecolor='black', histtype='stepfilled', fc='green', alpha=0.5, grid=False)
            np.log10(df_therm_sys['Lx1']).hist(ax=ax[1], bins=self.bins, label='Thermal MT | With Beaming',  edgecolor='black', histtype='stepfilled', fc='red', alpha=0.5, grid=False)
            # np.log10(df_wd_sys['Lx1']).hist(ax=ax[1], bins=self.bins, label='WD MT | With Beaming',  edgecolor='black', histtype='stepfilled', fc='white', alpha=0.5, grid=False)
            
        elif by=='both':
            df_ns = df[df['K_a'] == 13]
            df_bh = df[df['K_a'] == 14]
            
            df_nuc_bh   = df_bh[df_bh['mttype'] == 1]
            df_therm_bh = df_bh[df_bh['mttype'] == 4]
            df_nuc_ns   = df_ns[df_ns['mttype'] == 1]
            df_therm_ns = df_ns[df_ns['mttype'] == 4]
            
            df_nuc_bh_sys = df_nuc_bh.groupby(['idum_run', 'iidd_old']).mean().reset_index()
            df_therm_bh_sys = df_therm_bh.groupby(['idum_run', 'iidd_old']).mean().reset_index()
            df_nuc_ns_sys = df_nuc_ns.groupby(['idum_run', 'iidd_old']).mean().reset_index()
            df_therm_ns_sys = df_therm_ns.groupby(['idum_run', 'iidd_old']).mean().reset_index()
        
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
        plt.scatter(df_system['t'], np.log10(df_system['Lx1']))
        plt.xlabel('Time')
        plt.ylabel('log (Lx)')
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
    df = startrack_v2_mt_1_all() #nrows=100000
    pop = Population(df)
    df_sys = pop.df.groupby(['idum_run', 'iidd_old']).mean().reset_index()
    # pop.describe(pop.df, 'df')
    # pop.describe(pop.df_ulx, 'df_ulx') 
    pop.set_latex_font()
    # pop.XLF_plot(df_sys, include_beamed=True, include_duty_cycle=True, save=True)
    
    lmxrb = pop.df[pop.df['lmxrb'] == 1]
    # lmxrb = lmxrb.groupby(['idum_run', 'iidd_old']).mean().reset_index()

    plt.yscale('log')
    
    plt.hist(np.log10(lmxrb['Lx1']), bins=pop.bins, histtype='step', label='lmxrb')
    
    df_ns = lmxrb[lmxrb['K_a'] == 13]
    df_bh = lmxrb[lmxrb['K_a'] == 14]

    df_ns_sys = df_ns.groupby(['idum_run', 'iidd_old']).mean().reset_index()
    df_bh_sys = df_bh.groupby(['idum_run', 'iidd_old']).mean().reset_index()
    
    for dc in [0.1, 0.2, 0.3, 1.0]:
        
        vis = pop.sample_visible(lmxrb, include_duty_cycle=True, lmxrb_duty_cycle=dc)
        
        ns_hist_results = pop.XLF_mc_histogram(df_ns_sys, True, lmxrb_duty_cycle=dc)
        bh_hist_results = pop.XLF_mc_histogram(df_bh_sys, True, lmxrb_duty_cycle=dc)
    
        df_vis_ns = pop.sample_visible(df_ns_sys)
        df_vis_bh = pop.sample_visible(df_bh_sys)
        
        mean_ns = np.mean(ns_hist_results, axis=0)
        std_ns  = np.std(ns_hist_results, axis=0)
        
        mean_bh = np.mean(bh_hist_results, axis=0)
        std_bh  = np.std(bh_hist_results, axis=0)
        
        np.log10(df_vis_ns['Lx1']).hist(bins=pop.bins, label=f'NS | With Beaming | Visible only | d={dc}', histtype='step', alpha=1.0, grid=False)
        np.log10(df_vis_bh['Lx1']).hist(bins=pop.bins, label=f'BH | With Beaming | Visible only | d={dc}',histtype='step', alpha=1.0, grid=False, linestyle='--')
        
        plt.errorbar(pop.bin_centers, mean_ns, yerr=std_ns, linestyle='', capsize=1.0, label=f'NS | With Beaming | Visible only | d={dc}')
        plt.errorbar(pop.bin_centers, mean_bh, yerr=std_bh, linestyle='', capsize=1.0, label=f'BH | With Beaming | Visible only | d={dc}')
        
                
        # plt.hist(np.log10(vis['Lx1']), bins=pop.bins, histtype='step', label=f'd = {dc}')
        plt.legend()
        
        
    plt.xlabel(r'log $L_{x}$ $(\mathrm{erg \ s^{-1}})$')
    plt.legend()
    plt.savefig('../reports/figures/XLF_lmxrb_sys_inc_dc.png', dpi=500)
    plt.savefig('../reports/figures/XLF_lmxrb_sys_inc_dc.eps')
    plt.savefig('../reports/figures/XLF_lmxrb_sys_inc_dc.pdf')
    