"""
parameter_distibutions.py

Monte-Carlo simulation to obtain the input distribution of binary
system parameters used in simulations.
"""
import numpy as np
import matplotlib.pyplot as plt
import populations
from tqdm import tqdm

if __name__ == "__main__":
    print('Loading population...')
    df = populations.startrack_v2_mt_1_all() # 
    pop = populations.Population(df)

    pop.filter_non_bh_ns()
    pop.filter_non_thermal_or_nuclear_mt()
    pop.calc_sub_populations()
    
    N_samp = 500
    N_mc = 10000
    N_bins = 25
    bh_ratios = [0.0, 0.25, 0.5, 0.75, 1.00]
    
    params = ['mdot_ratio', 'Lx_iso', 'b', 'Lx1', 'theta_half_deg', 'zeta'
              'r_out', 'P_wind_days', 'P_sup_days']
    for param in params:
        print(param)
        plt.figure()
        for bh in bh_ratios:
        
            hists = np.ndarray((N_mc, N_bins))
            
            for i in tqdm(range(N_mc)):
                idx = pop.sample_systems(bh_ratio=bh, size=N_samp, subset='ulx')
                df_sampled = df.loc[idx]
                vals = df_sampled[param].values
                # print(param, bh, i)
                hists[i] = np.histogram(vals, bins=N_bins)[0]
            
            hist_bins = np.histogram(vals, bins=N_bins)[1]
            hist_mean = np.mean(hists, axis=0) 
            hist_std = np.std(hists, axis=0)

            plt.errorbar(hist_bins[:-1], hist_mean, yerr=hist_std, label=f'%BH = {bh}')
            # plt.hist(df_sampled['b'], bins=N_bins)
        
        plt.legend()
        plt.xlabel(param)
        plt.ylabel('N')
        plt.savefig(f'../reports/figures/paramdist/MC_{param}_N_mc={N_mc}')
        plt.yscale('log')
        plt.savefig(f'../reports/figures/paramdist/MC_{param}_N_mc={N_mc}_ylog')
