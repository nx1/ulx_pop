# -*- coding: utf-8 -*-
"""
Created on Tue May 19 12:31:13 2020

@author: norma

eRASS Simulation

"""
import numpy as np
import os
import time
import sqlite3

import populations
import auxil
import curvefunc
import batchfunc

def create_classification_dict():
    # Classication dict
    keys = (curve_classifications['system_id'].astype(str) +
            '-' + curve_classifications['dincl'].astype(str) +
            '-' + curve_classifications['inclination'].astype(str)).values
    classications = curve_classifications['classification'].values
    class_dict = dict(zip(keys, classications))
    return class_dict


def create_N_lim_dict():
    # Classication dict
    keys = (curve_classifications['system_id'].astype(str) +
            '-' + curve_classifications['dincl'].astype(str) +
            '-' + curve_classifications['inclination'].astype(str)).values
    N_lim = curve_classifications['N_lim'].values
    N_lim_dict = dict(zip(keys, N_lim))
    return N_lim_dict


def create_keys(selected_systems, selected_dincls, selected_inclinations,):
    a = np.core.defchararray.add(selected_systems.astype(str), '-')
    b = np.core.defchararray.add(selected_dincls.astype(str), '-')
    keys = np.core.defchararray.add(a,b)
    keys = np.core.defchararray.add(keys,selected_inclinations.astype(str))
    return keys


def get_erass_curve(theta, inclination, dincl):
    theta = round(theta, 2)
    inclination = int(inclination)
    dincl = int(dincl)
    
    
    key = f'{theta}-{inclination}-{dincl}'
    
    lc_path = erass_curve_folder + f'{key}.txt'
    isFile = os.path.isfile(lc_path)
    
    if isFile:        
        curve = curvefunc.load_curve_file(lc_path)
        return curve
    else:
        
        xcm_path = erass_curve_folder + f'{key}.xcm'
        parameters = {'period': 50,
                      'phase': 0,
                      'theta': theta,
                      'inclination': inclination,
                      'dincl': dincl,
                      'beta': 0.2,
                      'dopulse': 1,
                      'norm': 1}
        batchfunc.run_ulxlc(xcm_path, parameters, lc_path)
        curve = curvefunc.load_curve_file(lc_path)
        return curve


def is_ulx_on_first_cycle(curve, curve_period, P_wind, N_lim, sampling_interval=30*6):
    time_arr = np.array(curve['Time'])
    flux_arr = np.array(curve['Flux'])
    # Time scaling constant
    k = curve_period / P_wind
    
    t_start = np.random.uniform(P_wind)
    t_start_curve = k * t_start

    sample_index = np.argmin(np.abs(time_arr - t_start_curve))
    flux = flux_arr[sample_index]
    
    if flux > N_lim:
        return True
    else:
        return False


def sample_curve(curve, curve_period, P_wind, N_lim, sampling_interval=30*6):
    time_arr = np.array(curve['Time'])
    flux_arr = np.array(curve['Flux'])

    # Time scaling constant
    k = curve_period / P_wind
    
    # Sampling in days
    t_start = np.random.uniform(P_wind)
    t_interval = sampling_interval
    t_sample = [(t_start + i * t_interval) for i in range(8)]
    
    # Sampling points scaled for curve
    t_sample_curve = k * np.array(t_sample) % curve_period
    
    # Get corresponding indexs
    sample_indexes = [np.argmin(np.abs(time_arr - t)) for t in t_sample_curve]
    
    # Find corresponding fluxes
    fluxes = flux_arr[sample_indexes]
    return fluxes


def get_first_transient_cycle(alive_arr):
    for i, a in enumerate(alive_arr):
        if a != alive_arr[0]:
            return i
    raise ValueError('No transient cycle found, should not have got here')


def erass_classifier(is_ulx_erass):
    arr_sum = sum(is_ulx_erass)
    if  arr_sum == 0:
        return 'dead'
    if arr_sum == len(is_ulx_erass):
        return 'alive'
    else:
        return 'transient'

def run(sample_size=500, bh_ratio=0.5, dincl_cut=46):
    # eRASS Simulation Variables
    sample_size = sample_size
    bh_ratio = bh_ratio
    dincl_cut = dincl_cut
    
    # Treat sources that are over a wind period length as persistent sources
    P_wind_persistent_level = 4 * 365
    
    # Inputs
    key_class_dict = create_classification_dict()
    key_N_lim_dict = create_N_lim_dict()
    id_theta_dict = dict(df_ulx['theta_half_deg'])
    id_P_wind_dict = dict(df_ulx['P_wind_days'])
    
    
    selected_systems       = auxil.sample_by_bh_ratio(df_ulx, bh_ratio, sample_size)
    selected_dincls        = np.random.randint(0, dincl_cut, size=sample_size)
    selected_inclinations  = np.random.randint(0,91, size=sample_size)
    selected_keys          = create_keys(selected_systems, selected_dincls, selected_inclinations)
    selected_classications = [key_class_dict.get(key) for key in selected_keys]
    selected_N_lims        = [key_N_lim_dict.get(key) for key in selected_keys]
    selected_thetas        = [id_theta_dict.get(key) for key in selected_systems]
    selected_P_winds       = [id_P_wind_dict.get(key) for key in selected_systems]


    # Sample Result variables
    N_total     = sample_size
    N_alive     = selected_classications.count('alive') + selected_classications.count(None) #None --> no lightcurve --> alive
    N_transient = selected_classications.count('transient')
    N_dead      = selected_classications.count('dead')

    # Sample Check    
    assert N_alive + N_transient + N_dead == N_total

    # eRASS Results
    N_alive_persisitent = 0
    N_dead_persisitent  = 0
    N_transient_erass   = [0] * 8


    # for i in range(sample_size):
    #     erass.run()
    

    for i in range(sample_size):
        
        curve_classification = selected_classications[i]
        
        if curve_classification == None:
            N_alive_persisitent += 1
    
        if curve_classification == 'alive':
            N_alive_persisitent += 1
    
        if curve_classification == 'dead':
            N_dead_persisitent += 1
    
        if curve_classification == 'transient':
            
            # Get variables for calculation
            system_id   = selected_systems[i]
            dincl       = selected_dincls[i]
            inclination = selected_inclinations[i]
            theta       = selected_thetas[i]
            N_lim       = selected_N_lims[i]
            P_wind      = selected_P_winds[i]
            
            curve       = get_erass_curve(theta, inclination, dincl)

            # Long period systems will be treated as persistent
            if P_wind > P_wind_persistent_level:
                is_ulx = is_ulx_on_first_cycle(curve, 50, P_wind, N_lim, sampling_interval=30*6)
                if is_ulx:
                    N_alive_persisitent += 1
                else:
                    N_dead_persisitent += 1


            else:
                sampled_fluxes = sample_curve(curve, 50, P_wind, N_lim, sampling_interval=30*6)
                is_ulx_erass = list(np.where(sampled_fluxes > N_lim, True, False))
                
                erass_sample_classification = erass_classifier(is_ulx_erass)
                
                if erass_sample_classification == 'alive':
                    N_alive_persisitent += 1
                    
                if erass_sample_classification == 'dead':
                    N_dead_persisitent += 1
    
                if erass_sample_classification == 'transient':
                    first_transient_cycle = get_first_transient_cycle(is_ulx_erass)
                    N_transient_erass[first_transient_cycle] += 1
            
        
    N_transient_erass_cumulative = np.cumsum(N_transient_erass).tolist()
    
    assert N_alive_persisitent + N_dead_persisitent + N_transient_erass_cumulative[-1] == sample_size         
                    
    return N_alive_persisitent, N_dead_persisitent, N_transient_erass, N_transient_erass_cumulative



def sql_create_db():
    create_table_sql = """CREATE TABLE IF NOT EXISTS SAMPLE_RESULTS (
                                id integer PRIMARY KEY AUTOINCREMENT,
                                dincl integer,
                                bh_percent real,
                                N_alive_persisitent integer,
                                N_dead_persisitent integer,
                                N_transient_erass1 integer,
                                N_transient_erass2 integer,
                                N_transient_erass3 integer,
                                N_transient_erass4 integer,
                                N_transient_erass5 integer,
                                N_transient_erass6 integer,
                                N_transient_erass7 integer,
                                N_transient_erass8 integer,
                                N_t_cum1 integer,
                                N_t_cum2 integer,
                                N_t_cum3 integer,
                                N_t_cum4 integer,
                                N_t_cum5 integer,
                                N_t_cum6 integer,
                                N_t_cum7 integer,
                                N_t_cum8 integer
                                );"""
    
    conn = sqlite3.connect('erass.db')
    c = conn.cursor()
    c.execute(create_table_sql)
    conn.commit()
    conn.close()
    
    
def sql_save_result(unpacked_tuple):
    
    insert = """ INSERT INTO SAMPLE_RESULTS (
    dincl,
    bh_percent,
    N_alive_persisitent,
    N_dead_persisitent,
    N_transient_erass1,
    N_transient_erass2,
    N_transient_erass3,
    N_transient_erass4,
    N_transient_erass5,
    N_transient_erass6,
    N_transient_erass7,
    N_transient_erass8,
    N_t_cum1,
    N_t_cum2,
    N_t_cum3,
    N_t_cum4,
    N_t_cum5,
    N_t_cum6,
    N_t_cum7,
    N_t_cum8
    )
    VALUES(?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?)"""
    
    conn = sqlite3.connect('../../../src/erass.db')
    c = conn.cursor()
    print(unpacked_tuple)
    c.execute(insert, unpacked_tuple)
    conn.commit()
    conn.close()
    
    
if __name__ == "__main__":
    sql_create_db()
    df_ulx = populations.ulx()
    curve_classifications = auxil.load_curve_classifications()
    
    
    # This script struggles pretty hard when making it change directory.
    os.chdir('../data/interim/curves/')
    erass_curve_folder = 'eRASS_sims/'
    
    import itertools

    dincls = [46, 40, 30, 20, 10]
    bh_ratios = [0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0]
    repeats = np.arange(100)
    
    results = {}
    for i, variables in enumerate(itertools.product(dincls, bh_ratios, repeats)):
        print(i)
        dincl, bh, repeat = variables
        N_alive_persisitent, N_dead_persisitent, N_transient_erass, N_transient_erass_cumulative = run(sample_size=500, bh_ratio=bh, dincl_cut=dincl)
        results[i] =(dincl,
                     bh,
                     N_alive_persisitent,
                     N_dead_persisitent,
                     N_transient_erass,
                     N_transient_erass_cumulative)
        
        print(results[i])

        unpacked_tuple = (dincl,
                         bh,
                         N_alive_persisitent,
                         N_dead_persisitent,
                         *N_transient_erass,
                         *N_transient_erass_cumulative)
        
        sql_save_result(unpacked_tuple)




    
    
    
    
