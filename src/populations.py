# -*- coding: utf-8 -*-
"""
Created on Tue May 19 12:16:45 2020

@author: norma

This file contains two datasets ontained from startrack.
The first an old dataset containing #36k rows which we have ommited using as
the time sampling was not performed uniformly.

The second dataset is enormous and not fully provided here, instead a subset
is provided of all the systems accross all metallicities that are undergoing
mass trasnfer 'mt = 1', this provides 9,372,542 rows to work with.

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

from constants import G, c, Myr, R_sol, M_sol, NS_SPIN, BH_SPIN, NS_ISCO, BH_ISCO, epsilon, beta


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
    def __init__(self, df, name):
        self.df = df
        self.name = name
        self.calc_columns()
        self.calc_sub_populations()
        
    def calc_columns(self):
        self.df['is_bh'] = np.where(self.df['K_a'] == 13, 1, 0)
        
        # Make mass accretion rates positive
        self.df['dMmt_a'] = abs(self.df['dMmt_a'])
        self.df['dMmt_b'] = abs(self.df['dMmt_b'])
    
        #Convert to grams per second
        self.df['mdot_gs_a'] = self.df['dMmt_a'] * (M_sol/Myr)
        self.df['mdot_gs_b'] = self.df['dMmt_b'] * (M_sol/Myr)
    
    
        #Calculate Eddington luminosity and Eddington accretion rate
        self.df['LEdd'] = 1.3E38 * self.df['M_a']
        self.df['mdot_Edd'] = self.df['LEdd'] / (1/12 * c**2)
    
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
        self.df['r_isco']        = np.where(self.df['M_a'] < 2.5, NS_ISCO, BH_ISCO)
    
        self.df['r_sph'] = self.df['r_isco'] * self.df['mdot_ratio']
        self.df['r_out'] = 3 * epsilon / (beta * self.df['zeta']) * self.df['mdot_ratio']**3/2 * self.df['r_isco']
    
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

    
    def calc_sub_populations(self):
        self.df_ulx = self.df[self.df['Lx1'] > 1e39]
        self.df_ulx = self.df_ulx[self.df_ulx['mttype'] != 5]
        self.df_ulx_bh = self.df_ulx[self.df_ulx['is_bh'] == 1]
        self.df_ulx_ns = self.df_ulx[self.df_ulx['is_bh'] == 0]
        self.df_ulx_opening_angle_le_45 = self.df_ulx[self.df_ulx['theta_half_deg'] <= 45]
        self.df_ulx_P_wind_l_4_years = self.df_ulx_opening_angle_le_45[self.df_ulx_opening_angle_le_45['P_wind_days'] < 365*4]
        self.df_ulx_P_sup_l_4_years = self.df_ulx_opening_angle_le_45[self.df_ulx_opening_angle_le_45['P_sup_days'] < 365*4]
        
    
    def get_system(self, df, idum_run, iidd_old):
        df_system = df.copy()
        df_system = df[df['idum_run']==idum_run]
        df_system = df_system[df_system['iidd_old']==iidd_old]
        return df_system
    
    def get_unique_system_counts(self, df):
        df_unique_systems_count = df.groupby(['idum_run', 'iidd_old']).size().reset_index().rename(columns={0:'count'})
        df_unique_systems_count = df_unique_systems_count.sort_values('count', ascending=False)
        return df_unique_systems_count

    def calc_sampling_weights(self, df):
        gb = df.groupby(['idum_run', 'iidd_old'])[['t']].agg(['min', 'max']) # Get min and max t for each system
        gb[('t', 'range')] = gb[('t', 'max')] - gb[('t', 'min')] #Calculate t_range (i.e system alive time)
        t_tot = gb[('t', 'range')].sum()#Calculate total time alive for all systems
        gb['P_samp'] = gb[('t', 'range')] / t_tot
        samp_weights = gb['P_samp']
        return samp_weights
    
    def calc_ulx_sampling_weights(self):                
        self.ulx_bh_samp_weights = self.calc_sampling_weights(self.df_ulx_bh)
        self.ulx_ns_samp_weights = self.calc_sampling_weights(self.df_ulx_ns)
    
    def sample_ulxs(self, bh_ratio, size=500):
        try:
            self.ulx_bh_samp_weights
        except:
            self.calc_ulx_sampling_weights()
            
        N_bh = int(size*bh_ratio)
        N_ns = size-N_bh
        
        
        sampled_bh = np.random.choice(self.ulx_bh_samp_weights.index, size=N_bh , p=self.ulx_bh_samp_weights.values)
        sampled_ns = np.random.choice(self.ulx_ns_samp_weights.index, size=N_ns , p=self.ulx_ns_samp_weights.values)
        
        sampled_ulxs = np.concatenate((sampled_bh, sampled_ns))
        return sampled_ulxs
        
        
    
    def plot_XLF_by_NSBH(self, df):
        df_ns = df[df['K_a'] == 13]
        df_bh = df[df['K_a'] == 14]
    
        fig, ax = plt.subplots(2,1, sharex=True, figsize=(10,8))
        np.log10(df_ns['Lx_iso']).hist(ax=ax[0], bins=100, label='NS | Isotropic Emission',  edgecolor='black', histtype='stepfilled', fc='blue', alpha=0.5)
        np.log10(df_bh['Lx_iso']).hist(ax=ax[0], bins=100, label='BH | Isotropic Emission',  edgecolor='black', histtype='stepfilled', fc='gray', alpha=0.8)
    
        np.log10(df_ns['Lx1']).hist(ax=ax[1], bins=100, label='NS | With Beaming',  edgecolor='black', histtype='stepfilled', fc='blue', alpha=0.5)
        np.log10(df_bh['Lx1']).hist(ax=ax[1], bins=100, label='BH | With Beaming',  edgecolor='black', histtype='stepfilled', fc='gray', alpha=0.8)
    
        ax[0].axvline(39, c='r')
        ax[1].axvline(39, c='r')
    
        ax[0].set_title('XLF')
        ax[0].set_yscale('log')
        ax[1].set_yscale('log')
        ax[1].set_xlabel(r'log Luminosity $(\mathrm{erg \ s^{-1}})$')
        ax[0].set_ylabel(r'N')
        ax[1].set_ylabel(r'N')
    
        ax[0].legend()
        ax[1].legend()
        # plt.savefig('../reports/figures/luminosity_distributions_by_type.png', dpi=500)
        # plt.savefig('../reports/figures/luminosity_distributions_by_type.eps')
        
    def plot_XLF_by_mttype(self, df):
        df_nuc = df[df['mttype'] == 1]
        df_therm = df[df['mttype'] == 4]
        df_wd = df[df['mttype'] == 5]
    
    
        fig, ax = plt.subplots(2,1, sharex=True, figsize=(10,8))
        np.log10(df_nuc['Lx_iso']).hist(ax=ax[0], bins=100, label='Nuclear MT | Isotropic Emission',  edgecolor='black', histtype='stepfilled', fc='green', alpha=0.5)
        np.log10(df_therm['Lx_iso']).hist(ax=ax[0], bins=100, label='Thermal MT | Isotropic Emission',  edgecolor='black', histtype='stepfilled', fc='red', alpha=0.5)
        np.log10(df_wd['Lx_iso']).hist(ax=ax[0], bins=100, label='Thermal MT | Isotropic Emission',  edgecolor='black', histtype='stepfilled', fc='white', alpha=0.5)
    
    
        np.log10(df_nuc['Lx1']).hist(ax=ax[1], bins=100, label='Nuclear MT | With Beaming',  edgecolor='black', histtype='stepfilled', fc='green', alpha=0.5)
        np.log10(df_therm['Lx1']).hist(ax=ax[1], bins=100, label='Thermal MT | With Beaming',  edgecolor='black', histtype='stepfilled', fc='red', alpha=0.5)
        np.log10(df_wd['Lx1']).hist(ax=ax[1], bins=100, label='Thermal MT | With Beaming',  edgecolor='black', histtype='stepfilled', fc='white', alpha=0.5)
    
        ax[0].axvline(39, c='r')
        ax[1].axvline(39, c='r')
    
    
        ax[0].set_yscale('log')
        ax[1].set_yscale('log')
        ax[1].set_xlabel(r'log Luminosity $(\mathrm{erg \ s^{-1}})$')
        ax[0].set_ylabel(r'N')
        ax[1].set_ylabel(r'N')
    
        ax[0].legend()
        ax[1].legend()
        #plt.savefig('../reports/figures/luminosity_distributions_by_MT.png', dpi=500)
        #plt.savefig('../reports/figures/luminosity_distributions_by_MT.eps')
    
    def plot_XLF_by_mttype_and_bh(self, df):
        df_ns = df[df['K_a'] == 13]
        df_bh = df[df['K_a'] == 14]
        df_nuc_bh   = df_bh[df_bh['mttype'] == 1]
        df_therm_bh = df_bh[df_bh['mttype'] == 4]
        df_nuc_ns   = df_ns[df_ns['mttype'] == 1]
        df_therm_ns = df_ns[df_ns['mttype'] == 4]
    
        fig, ax = plt.subplots(2,1, sharex=True, figsize=(10,8))
        np.log10(df_nuc_bh['Lx_iso']).hist(ax=ax[0], bins=100, label='BH | Nuclear MT | Isotropic Emission',  edgecolor='black', histtype='stepfilled', fc='cyan', alpha=0.5)
        np.log10(df_therm_bh['Lx_iso']).hist(ax=ax[0], bins=100, label='BH | Thermal MT | Isotropic Emission',  edgecolor='black', histtype='stepfilled', fc='purple', alpha=0.5)
        np.log10(df_nuc_ns['Lx_iso']).hist(ax=ax[0], bins=100, label='NS | Nuclear MT | Isotropic Emission',  edgecolor='black', histtype='stepfilled', fc='green', alpha=0.5)
        np.log10(df_therm_ns['Lx_iso']).hist(ax=ax[0], bins=100, label='NS | Thermal MT | Isotropic Emission',  edgecolor='black', histtype='stepfilled', fc='red', alpha=0.5)
    
    
    
        np.log10(df_nuc_bh['Lx1']).hist(ax=ax[1], bins=100, label='BH | Nuclear MT | With Beaming',  edgecolor='black', histtype='stepfilled', fc='cyan', alpha=0.5)
        np.log10(df_therm_bh['Lx1']).hist(ax=ax[1], bins=100, label='BH | Thermal MT | With Beaming',  edgecolor='black', histtype='stepfilled', fc='purple', alpha=0.5)
        np.log10(df_nuc_ns['Lx1']).hist(ax=ax[1], bins=100, label='NS | Nuclear MT | With Beaming',  edgecolor='black', histtype='stepfilled', fc='green', alpha=0.5)
        np.log10(df_therm_ns['Lx1']).hist(ax=ax[1], bins=100, label='NS | Thermal MT | With Beaming',  edgecolor='black', histtype='stepfilled', fc='red', alpha=0.5)
    
        ax[0].axvline(39, c='r')
        ax[1].axvline(39, c='r')
    
    
        ax[0].set_yscale('log')
        ax[1].set_yscale('log')
        ax[1].set_xlabel(r'log Luminosity $(\mathrm{erg \ s^{-1}})$')
        ax[0].set_ylabel(r'N')
        ax[1].set_ylabel(r'N')
    
        ax[0].legend(loc='upper left')
        ax[1].legend(loc='upper left')
        #plt.savefig('../reports/figures/luminosity_distributions.png', dpi=500)
        #plt.savefig('../reports/figures/luminosity_distributions.eps')
        
    def plot_system_luminosity_evolution(self, df, idum_run, iidd_old):
        df_system = self.get_system(df, idum_run, iidd_old)
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
        
    def export_population_info(self, df, identifier):
        """Population Info
           ---------------
           name : {self.name}
           identifier : {identifier}
           N_rows :
           N_rows_ulx
           N_rows_beamed :
           N_rows_opening_angle_l_45 :
           N_rows_P_wind_l_4_years :
        
           N_rows_bh :
           N_rows_ns :
        
           N_binaries :
           N_ulx_binaries :
           N_beamed_ulx_binaries :
           N_opening_angle_l_45_binaries :
           N_P_wind_l_4_years_binaries :
             
           unique_Z :
        """
    
    def export_df(self, df, identifier):
        df.to_csv(f'{identifier}.csv', index=False)
        
    
    def export_for_ulxlc(self, df, identifier):
        self.describe()
        self.export_population_info(df, identifier)
        self.export_df(df, identifier)
        self.export_to_binary_data(df)
        
    def describe(self):
        try:
            self.df_ulx
        except(AttributeError):
            self.calc_sub_populations()

        self.N_rows = len(self.df)
        self.N_rows_ulx = len(self.df_ulx)
        self.N_rows_opening_angle_le_45 = len(self.df_ulx_opening_angle_le_45)
        self.N_rows_P_wind_l_4_years = len(self.df_ulx_P_wind_l_4_years)
        self.N_rows_P_sup_l_4_years = len(self.df_ulx_P_sup_l_4_years)
        
        self.N_rows_bh = len(df[df['K_a'] == 14])
        self.N_rows_ns = len(df[df['K_a'] == 13])
        
        self.N_binaries = len(self.df.groupby(['idum_run', 'iidd_old']).size().reset_index().rename(columns={0:'count'}))
        self.N_beamed_ulx_binaries = len(self.df_ulx.groupby(['idum_run', 'iidd_old']).size().reset_index().rename(columns={0:'count'}))
        self.N_opening_angle_l_45_binaries = len(self.df_ulx_opening_angle_le_45.groupby(['idum_run', 'iidd_old']).size().reset_index().rename(columns={0:'count'}))
        self.N_P_wind_l_4_years_binaries = len(self.df_ulx_P_wind_l_4_years.groupby(['idum_run', 'iidd_old']).size().reset_index().rename(columns={0:'count'}))
        self.N_P_sup_l_4_years_binaries = len(self.df_ulx_P_sup_l_4_years.groupby(['idum_run', 'iidd_old']).size().reset_index().rename(columns={0:'count'}))
        
                
        self.unique_Z =self. df['Z'].unique()

    
if __name__ == "__main__":
    df = startrack_v2_mt_1_all(10000)
    
    pop = Population(df, 'startrack_mt_1')
    
    #pop.describe()
    #pop.calc_ulx_sampling_weights()
