# -*- coding: utf-8 -*-
"""
Created on Thu Dec 10 11:17:06 2020

@author: norma
"""
import sys
sys.path.insert(0, '../src')

import numpy as np
import matplotlib.pyplot as plt
    
import pytest
import populations

def test_load_population():
    populations.startrack_v2_mt_1_all(nrows=1000)
    populations.startrack_v2_mt_1_z02(nrows=1000)
    populations.startrack_v2_mt_1_z002(nrows=1000)
    populations.startrack_v2_mt_1_z0002(nrows=1000)


@pytest.fixture
def pop():
    df = populations.startrack_v2_mt_1_all(nrows=100000)
    pop = populations.Population(df) 
    return pop

def test_calc_sub_populations(pop):
    pop.calc_sub_populations()

def test_filter_df_ulx_by_Z(pop):
    # TODO
    """
    b = pop()
    a = pop()

    a.filter_df_ulx_by_Z(0.002)
    b.Z_filter(Z, ) ?
    
    assert a.df_ulx == assert b.df_ulx 
    assert a.df == assert b.df
    
    assert a.df_ulx == b.df_ulx
    assert a.df == b.df 
    """
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

def test_calc_system_averages(pop):
    pop.calc_system_averages(pop.df)


def test_sample_systems(pop):
    np.random.seed(500)
    p_bh = 0.5
    N = 500
    sampled_indexs2 = pop.sample_systems(p_bh, size=N, subset='ulx') # Sample ULXs only
    sampled_indexs3 = pop.sample_systems(p_bh, size=N, subset='all') # Sample all
    assert len(sampled_indexs2)==N
    assert len(sampled_indexs3)==N
    
    with pytest.raises(KeyError):
        pop.sample_systems(0.0,subset='bad_subset')

    df_samp = pop.sample_systems(0.0,subset='all', return_df=True)   
    assert 14 not in df_samp['K_a'].value_counts()
    assert df_samp['Lx1'].min() < 1e39

    df_samp = pop.sample_systems(1.0,subset='all', return_df=True)
    assert 13 not in df_samp['K_a'].value_counts()
    assert df_samp['Lx1'].min() < 1e39

    df_samp = pop.sample_systems(0.0, return_df=True)
    assert 14 not in df_samp['K_a'].value_counts()
    assert df_samp['Lx1'].min() > 1e39

    df_samp = pop.sample_systems(1.0, return_df=True)
    assert 13 not in df_samp['K_a'].value_counts()
    assert df_samp['Lx1'].min() > 1e39



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

def test_arccos(pop):
    """
    
    https://www.wolframalpha.com/input/?i=y+%3D+tan%28+%28pi%2F2%29+-+arccos%281+-+%2873+%2F+x**2%29%29+%29
    https://www.wolframalpha.com/input/?i=tan%28+%28pi%2F2%29+-+arccos%281+-+%2873+%2F+x**2%29%29+%29+%3D+2
    https://www.wolframalpha.com/input/?i=y+%3D+arccos%281+-+%2873+%2F+x**2%29%29
    
    for zeta = 2 mdot_0 = sqrt(73*(5+2sqrt(5))) = 26.295739668527474
    
    arccos term domain: x >= +-sqrt(146) / 2 = +-6.041522986797286
    """
    #Test code removed, it's sorted now.
    pass
    

    
