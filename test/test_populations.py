# -*- coding: utf-8 -*-
"""
Created on Thu Dec 10 11:17:06 2020

@author: norma
"""
import sys
sys.path.insert(0, '../src')

import numpy as np

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
    # TODO
    """
    b = pop()
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
    
def test_calc_ulx_binary_dict(pop):
    pop.calc_ulx_binary_dict()

def test_calc_system_averages(pop):
    pop.calc_system_averages(pop.df)

def test_sample_ulxs(pop):
    p_bh = 0.5
    N = 500
    sampled_indexs = pop.sample_ulxs(p_bh, size=N)
    assert len(sampled_indexs)==N

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
    arccos() domain is from -1 to 1 acrcos(-1) = pi ; arccos(1) = 0; 
    therefore 73/mdot^2 should also be between 2 to 0 or it will give bad values
    
    if mdot_0 < sqrt(73/2) (6.04) then nan
    mdot_0 == 8.5 gives zeta = 0
    mdot_0 < 8.5 gives zeta < 0
    as mdot_0 --> inf arccos(arg) --> 0 zeta --> inf
    
    however we have a lower limit of zeta = 2 which corresponds to mdot 6.452
    
    b>=1
    
    if mdot <= 8.5
        set zeta 
    """
    
#     assert pop.df['b'].min() >= 0
#     assert pop.df['b'].max() <= 1.0
    
#     import pandas as pd
    
#     pop.df['cos_arg']  = 1 - (73/(pop.df['mdot_ratio']**2))
#     pop.df['tan_arg']  = (np.pi / 2) - np.arccos(pop.df['cos_arg'])
#     # self.df['zeta']    = np.tan((np.pi/2) - np.arccos(1 - (73/(self.df['mdot_ratio']**2))))
#     pop.df['zeta_new'] = np.tan(pop.df['tan_arg'])
    
#     df_test = pd.DataFrame()
    
#     df_test['mdot_ratio'] = np.random.normal(loc=0, scale=5, size=10000)
#     df_test['cos_arg']    = 1 - (73/(df_test['mdot_ratio']**2))
#     df_test['tan_arg']    = (np.pi / 2) - np.arccos(df_test['cos_arg'])
#     df_test['zeta_new']   = (np.pi / 2) - np.arccos(df_test['cos_arg'])
#     df_test['zeta']       = np.tan((np.pi/2) - np.arccos(1 - (73/(df_test['mdot_ratio']**2))))
    
    
#     import matplotlib.pyplot as plt
    
#     def plot_hists(cols):
#         for c in cols:
#             plt.figure()
#             plt.title(c)
#             plt.hist(df_test[c], bins=100)
            
#     plot_hists(['mdot_ratio', 'cos_arg', 'tan_arg', 'zeta_new', 'zeta'])

#     plt.show()
#     pd.testing.assert_series_equal(df_test['zeta'], df_test['zeta_new'])
    
#     import pandas as pd
    
#     pd.testing.assert_series_equal(pop.df['zeta'], pop.df['zeta_new'])

    

    
