# -*- coding: utf-8 -*-
"""
Created on Tue Apr 21 21:38:52 2020

@author: norma
"""

from os import path
import numpy as np
import pandas as pd
from tqdm import tqdm

from ulxlc import run_ulxlc
from curvefunc import load_curve_file, scale_light_curve_period
import create_erass_parent_population

class Population:
    def __init__(self):
        self.dataframe = None
        
    def split_bh_ns_systems(self):
        self.bh_systems = self.dataframe[self.dataframe['is_bh'] == 1].index
        self.ns_systems = self.dataframe[self.dataframe['is_bh'] == 0].index
        
    def sample(self, bh_ratio, n):
        ns_ratio = 1 - bh_ratio
        bh_weights = [bh_ratio/len(self.bh_systems)]*len(self.bh_systems)
        ns_weights = [ns_ratio/len(self.ns_systems)]*len(self.ns_systems)
        selected_systems = np.random.choice([*self.bh_systems, *self.ns_systems], size=n, p=[*bh_weights, *ns_weights])
        sampled_df = self.dataframe.loc[selected_systems]
        return sampled_df


def get_curve(row):
    
    def get_ulxlc_parameters(row):
        """Create parameters dictionary from df_a row"""
        parameters = {'period': 50,
                        'phase': 0,
                        'theta': row['theta'],
                        'inclination': row['inclination'],
                        'dincl': row['dincl'],
                        'beta': 0.2,
                        'dopulse': 1,
                        'norm': 1}
        return parameters

    P_wind = row['P_wind_days']
    parameters = get_ulxlc_parameters(row)
    filename = str(row['theta']) + '-' + str(row['inclination']) +'-' + str(row['dincl'])
    xcm_file = 'eRASS_sims/' + filename + '.xcm'
    lc_file = 'eRASS_sims/' + filename + '.txt'
    
    if not path.exists(lc_file):
        run_ulxlc(xcm_file, parameters, lc_file)

    curve = load_curve_file(lc_file)
    curve = scale_light_curve_period(curve, 50, P_wind)
    return curve


def erass_sample_curve(curve, curve_period, N_lim, number_of_repeats):
    """Samples a lightcurve of a given period as would be done by eRASS
    (8 cycles of 6 months)
    parameters:
    -----------
    curve: lightcurve dataframe
    curve_period: period of the lightcurve
    sampling interval: how often to sample the lightcurve after selecting
    a random position somwhere in the first period.
    number_of_repeats: number of MC iterations
    """
    def truth_table_processor(truth):
        """Runs through the 8*N is ulx truth table
        to create a new dataframe with each column corresponding to an eRASS
        cycle and weather or not the source was observed as a transient ulx or not.
        thanks to based James for the help
        """
        observed_as_transient = []
        for row in truth:
            transient=8*[True]
            i=0 # Not observed as transient
            while row[i] == row[0] and i!=8:
                transient[i]=False
                i+=1 #observed as transient
            observed_as_transient.append(transient)
        
        df = pd.DataFrame(observed_as_transient)
        return df.mean()
    
    time_arr = np.array(curve['Time'])
    flux_arr = np.array(curve['Flux'])
    
    truth = []
    for n in tqdm(range(number_of_repeats)):
        start_time = np.random.uniform(0, curve_period)
        sample_times = np.array([(start_time + i*30*6)%curve_period for i in range(0,9)])
        sample_indexes = [np.argmin(np.abs(time_arr - t)) for t in sample_times]
        
        fluxes = flux_arr[sample_indexes]
        
        truth.append(list(np.where(fluxes > N_lim, True, False)))
    sample_results = truth_table_processor(truth)
    return sample_results


def main():
    transient_curves = create_erass_parent_population.main()
    
    parent_population = Population()
    parent_population.dataframe = transient_curves
    parent_population.split_bh_ns_systems()
    
    parent_systems.dataframe.to_csv(f'../data/interim/erass_simulations/parent_population.csv')
    for i in range(1000):
        for bh_ratio in [0, 0.25, 0.5, 0.75, 1.0]:
            df_results = pd.DataFrame()
            sampled_systems = parent_systems.sample(bh_ratio, 500)
            for index, row in sampled_systems.iterrows():
                curve = get_curve(row)
                N_lim = row['N_lim']
                curve_period = row['P_wind_days']
                df_results[index] = erass_sample_curve(curve, curve_period, N_lim, number_of_repeats=10000)
                
            print(f'cycle: {i}')
            print('results saved to:')
            print(f'../data/interim/erass_simulations/')
            df_results.to_csv(f'../data/interim/erass_simulations/{i}_bh_r={bh_ratio}.csv')
            sampled_systems.to_csv(f'../data/interim/erass_simulations/{i}_bh_r={bh_ratio}_sampled_systems.csv')


if __name__ == "__main__":
    main()