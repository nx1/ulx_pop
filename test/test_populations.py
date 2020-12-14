# -*- coding: utf-8 -*-
"""
Created on Thu Dec 10 11:17:06 2020

@author: norma
"""



import sys
sys.path.insert(0, '../src')

import pytest
import populations

def test_load_population():
    populations.startrack_v2_mt_1_all(nrows=1000)
    populations.startrack_v2_mt_1_z02(nrows=1000)
    populations.startrack_v2_mt_1_z002(nrows=1000)
    populations.startrack_v2_mt_1_z0002(nrows=1000)


@pytest.fixture
def pop():
    df = populations.startrack_v2_mt_1_all(nrows=1000)
    pop = populations.Population(df) 
    return pop

def test_calc_sub_populations(pop):
    pop.calc_sub_populations()

def test_filter_df_ulx_by_Z(pop):
    pop.filter_df_ulx_by_Z(0.002)
    
def test_split_ns_bh(pop):
    pop.split_ns_bh(pop.df)

def test_calc_bh_ns_ulx_sub_populations(pop):
    pop.calc_bh_ns_ulx_sub_populations()
    
def test_gb_sys(pop):
    pop.gb_sys(pop.df)

def test_get_system_ids(pop):
    pop.get_system_ids(pop.df)
    
def test_get_system(pop):
    idum_run = -100000
    iidd_old = 11058
    pop.get_system(idum_run, iidd_old)
    
def test_calc_sampling_weights(pop):
    pop.calc_sampling_weights(pop.df)

def test_calc_ulx_sampling_weights(pop):
    pop.calc_ulx_sampling_weights()
    
def test_calc_ulx_binary_dict(pop):
    pop.calc_ulx_binary_dict()

def test_calc_system_averages(pop):
    pop.calc_system_averages(pop.df)

def test_sample_ulxs(pop):
    pop.sample_ulxs(0.5, size=500)
      
    
    
def test_pivot_binaries_count(pop):
    pop.pivot_binaries_count(pop.df)
    
def test_plot_mass_accretion_rates(pop):
    pop.plot_mass_accretion_rates(pop.df, save=False)

def test_plot_system_luminosity_evolution(pop):
    idum_run = -100000
    iidd_old = 11058
    pop.plot_system_luminosity_evolution(idum_run, iidd_old)
    
def test_describe(pop):
    pop.describe(pop.df, 'df')
