# -*- coding: utf-8 -*-
"""
Created on Mon Dec 14 18:18:10 2020

@author: norma
"""
import sys
import os
import sqlite3
import pytest
sys.path.insert(0, '../src')
os.chdir('../src')

import populations
from results_processor import ResultsProcessor

@pytest.fixture
def pop():
    df = populations.startrack_v2_mt_1_all(nrows=1000)
    pop = populations.Population(df) 
    return pop



def test_set_active_db():
    rp = ResultsProcessor()
    db1 = 'xlf.db'
    db2 = 'ulxlc.db'
    db3 = 'old/ulxlc3.db'
    db4 = 'old_db/ulxlc2.db'
    db5 = 'old_db/ulxlc3.db'
    dbs = [db1, db2, db3, db4, db5]
    for db in dbs:
        rp.set_active_db(db)
        print(rp.db)
        conn = sqlite3.connect(rp.db)
        conn.close()


def test_set_parent_population(pop):
    rp = ResultsProcessor()
    rp.set_parent_population(pop)


# def test_ADT_result_stats():
# def test_ERASS_result_stats():
# def test_calc_N_persistent():
# def test_calc_all_classifications_count_stats():
# def test_calc_bh_ratio_classfication_correlation_constant():
# def test_calc_classification_counts():
# def test_get_classifications():
# def test_get_erass_evolution_from_key():
# def test_get_run_counts():
# def test_get_run_ids():
# def test_set_active_db():

# def test_sim_bh_ratio_classifications_sampler():
# def test_sim_classifications_sampler_with_duty_cycle():
# def test_sim_erass_run():
    
    
# def test_table_classifications_calc_intermediate():
# def test_table_classifications_map_info():
# def test_table_classifications_map_systems():
# def test_table_classifications_pivot():
# def test_table_classifications_split_by_metallicity():
# def test_table_create_classifications():
# def test_table_create_erass_evolution():
# def test_table_create_erass_mc_info():
# def test_table_create_erass_mc_results():
# def test_table_create_erass_mc_sampled_systems():
# def test_table_create_transient():
# def test_table_delete_erass_evolution_non_1_duty_cycle():
# def test_table_erass_mc_results_map_info():
# def test_table_load():
# def test_table_load_classifications():
# def test_table_load_erass_evolution():
# def test_table_load_erass_mc_info():
# def test_table_load_erass_mc_results():
# def test_table_load_erass_mc_sampled_systems():
# def test_table_load_transient():


        
