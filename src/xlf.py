"""
Panel 1 : Unbeamed emission split NS/BH
Panel 2 : Beamed Emission split NS/BH
Panel 3 : Beamed + duty_cycle (observed) split NS/BH
Panel 4 : Beamed + duty_cycle + precession (observed) split NS/BH
Panel 5 : Beamed + duty_cycle + prcession (observed) combined NS/BH different %_BH
"""


import ctypes
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from tqdm import tqdm

import populations
from ulxlc import ULXLC

from constants import set_latex_font, params_default


def calc_cum_hist(data, bins):
    hist = np.histogram(data, bins=bins)[0]
    hist_cum = np.cumsum(hist)
    hist_surv = sum(hist) - hist_cum
    return hist_surv

def sample(pop, df, p_BH, N):
    sampled_indexs = pop.sample_ulxs(p_BH,N)
    df_sampled = df.loc[sampled_indexs]
    
    # df_sampled['N'] = np.arange(0,N,1)
    # df_sampled['run_id'] = run_id
    
    
    df_sampled['r'] = np.random.random(size=N)
    # df_sampled['d_vis'] = df_sampled['r'] < df_sampled['d']
    df_sampled['vis'] = df_sampled['r'] < df_sampled['p_obs']
    df_sampled['log_Lx1_vis'] = np.where(df_sampled['vis']==True, df_sampled['log_Lx1'], 0)
    df_sampled['inclination'] = np.random.randint(0, 91, size=N)
    df_sampled['dincl'] = np.random.randint(0, 46, size=N)
    return df_sampled

def main(pop, N=500, N_mc=10000, save=False):    
    duty_cycle = 0.3
    p_BH = 0.5
    bh_ratios = [0, 0.25, 0.5, 0.75, 1.0]
    
    # XLF settings
    N = N #Population size
    N_mc = N_mc
    
    bin_min = 38
    bin_max = 43
    bin_width = 0.25
    bins = np.arange(bin_min, bin_max, bin_width)
    nbins = len(bins)
    bin_centers = 0.5 * (bins[:-1] + bins[1:])
    
    panels = ['L_iso',
              'Lx1',
              'Lx1_obs']
    
    
    # df = df.set_index(['idum_run', 'iidd_old'])
    cols = ['idum_run', 'iidd_old', 't', 'dt', 'K_a', 'mttype', 'log_Lx_iso', 'log_Lx1', 'Lx1', 'theta_half_deg', 'b', 'lmxrb']
    df = pop.df[cols]
    
    df = df[(df['K_a']==13) | (df['K_a']==14)] # Only include BH/NS systems
    df = df[(df['mttype'] == 1.0) | (df['mttype'] == 4.0)] # Remove systems that aren't thermal or nuclear mass transfer
    
    df['d'] = np.where(df['lmxrb']==1, duty_cycle, 1.0) # Duty cycle
    df['p_obs'] = df['b'] * df['d']     # Observation probability
    df['1-p_obs'] = 1 - df['p_obs']     # Non-Observation prob

    
    Lx = np.ndarray((N_mc, N))
    
    iso_ns_hists = np.ndarray((N_mc, nbins-1))
    iso_bh_hists = np.ndarray((N_mc, nbins-1))
    
    lx1_ns_hists = np.ndarray((N_mc, nbins-1))
    lx1_bh_hists = np.ndarray((N_mc, nbins-1))
    
    obs_b_d_ns_hists = np.ndarray((N_mc, nbins-1))
    obs_b_d_bh_hists = np.ndarray((N_mc, nbins-1))

    obs_prec_ns_hists = np.ndarray((N_mc, nbins-1))
    obs_prec_bh_hists = np.ndarray((N_mc, nbins-1))
    
    obs_prec_0 = np.ndarray((N_mc, nbins-1))
    obs_prec_025 = np.ndarray((N_mc, nbins-1))
    obs_prec_05 = np.ndarray((N_mc, nbins-1))
    obs_prec_075 = np.ndarray((N_mc, nbins-1))
    obs_prec_1 = np.ndarray((N_mc, nbins-1))

    
    ulxlc = ULXLC(100, 1.0)
    ulxlc.set_params(*params_default)
    
    for i in tqdm(range(N_mc)):
        # run_id = str(uuid4())
        df_sampled = sample(pop, df, p_BH, N)
            
    
        Lx[i] = df_sampled['log_Lx_iso'].values
        c_Lx_prec = (ctypes.c_double * N)()    
        c_Lx      = (ctypes.c_double * N)(*df_sampled['Lx1'].values)
        c_thetas  = (ctypes.c_double * N)(*df_sampled['theta_half_deg'].values)
        c_incls   = (ctypes.c_double * N)(*df_sampled['inclination'].values)
        c_dincls  = (ctypes.c_double * N)(*df_sampled['dincl'].values)
        
        ulxlc.xlf_calc_L_prec(c_Lx_prec, c_Lx, c_thetas, c_incls, c_dincls, N)
        
        Lx_prec = np.ctypeslib.as_array(c_Lx_prec)
        df_sampled['log_Lx1_prec'] = np.log10(Lx_prec)
        
        
        # df_sampled.loc[:,'L_prec'] = df_L_samp['LSAMP'].values
        # df_sampled['log_L_prec'] = np.log10(df_sampled['L_prec'])
        # df_sampled['log_L_prec_vis'] = np.where(df_sampled['d_vis']==True, df_sampled['log_L_prec'], 0)
        
        df_ns = df_sampled[df_sampled['K_a'] == 13]
        df_bh = df_sampled[df_sampled['K_a'] == 14]
        

        iso_ns_hists[i]      = calc_cum_hist(df_ns['log_Lx_iso'], bins)
        iso_bh_hists[i]      = calc_cum_hist(df_bh['log_Lx_iso'], bins)
        
        lx1_ns_hists[i]      = calc_cum_hist(df_ns['log_Lx1'], bins)
        lx1_bh_hists[i]      = calc_cum_hist(df_bh['log_Lx1'], bins)
        
        obs_b_d_ns_hists[i]  = calc_cum_hist(df_ns['log_Lx1_vis'], bins)
        obs_b_d_bh_hists[i]  = calc_cum_hist(df_bh['log_Lx1_vis'], bins)
        
        obs_prec_ns_hists[i] = calc_cum_hist(df_ns['log_Lx1_prec'], bins)
        obs_prec_bh_hists[i] = calc_cum_hist(df_bh['log_Lx1_prec'], bins)
        
        obs_prec_0[i]        = calc_cum_hist(df_bh['log_Lx1_prec'], bins)
        obs_prec_025[i]      = calc_cum_hist(df_bh['log_Lx1_prec'], bins)
        obs_prec_05[i]       = calc_cum_hist(df_bh['log_Lx1_prec'], bins)
        obs_prec_075[i]      = calc_cum_hist(df_bh['log_Lx1_prec'], bins)
        obs_prec_1[i]        = calc_cum_hist(df_bh['log_Lx1_prec'], bins)

    # =============================================================================
    # Calculate Means and stds
    # =============================================================================
    
    mean_bh_iso = np.mean(iso_bh_hists, axis=0) 
    std_bh_iso = np.std(iso_bh_hists, axis=0) 
    
    mean_ns_iso = np.mean(iso_ns_hists, axis=0) 
    std_ns_iso = np.std(iso_ns_hists, axis=0) 
    
    mean_bh_lx1 = np.mean(lx1_bh_hists, axis=0) 
    std_bh_lx1 = np.std(lx1_bh_hists, axis=0) 
    
    mean_ns_lx1 = np.mean(lx1_ns_hists, axis=0) 
    std_ns_lx1 = np.std(lx1_ns_hists, axis=0) 
    
    mean_bh_obs_b_d = np.mean(obs_b_d_bh_hists, axis=0) 
    std_bh_obs_b_d = np.std(obs_b_d_bh_hists, axis=0) 
    
    mean_ns_obs_b_d = np.mean(obs_b_d_ns_hists, axis=0) 
    std_ns_obs_b_d = np.std(obs_b_d_ns_hists, axis=0) 
    
    mean_obs_prec_bh = np.mean(obs_prec_bh_hists, axis=0)
    std_obs_prec_bh = np.std(obs_prec_bh_hists, axis=0)
    
    mean_obs_prec_ns = np.mean(obs_prec_ns_hists, axis=0)
    std_obs_prec_ns = np.std(obs_prec_ns_hists, axis=0)
    
    
    
    # mean_ns_obs_b_d_p = np.mean(obs_b_d_p_ns_hists, axis=0)
    # std_ns_obs_b_d_p = np.std(obs_b_d_p_ns_hists, axis=0)
    
    
    # =============================================================================
    # Plotting
    # =============================================================================
    
    def plot_area(ax, x, mean, std, **kwargs):
        ax.fill_between(x, mean-std, mean+std, **kwargs)
    
    nrows = 4
    w = 6.93
    h = w/nrows
    
    fig, ax = plt.subplots(1,nrows, sharey=True)
    fig.set_size_inches([w,h])
    plt.tight_layout()
    ax[0].set_ylabel(r'$N( > L)$')
    ax[0].set_xlabel(r'log $L_{\mathrm{iso}}$ $(\mathrm{erg \ s^{-1}})$', fontsize=10)
    ax[1].set_xlabel(r'log $L_{x}$ $(\mathrm{erg \ s^{-1}})$', fontsize=10)
    ax[2].set_xlabel(r'log $L_{x}$ $(\mathrm{erg \ s^{-1}})$', fontsize=10)
    ax[3].set_xlabel(r'log $L_{x}$ $(\mathrm{erg \ s^{-1}})$', fontsize=10)
    
    plt.subplots_adjust(top=0.99,
                        bottom=0.215,
                        left=0.088,
                        right=0.996,
                        hspace=0.0,
                        wspace=0.0)
    
    
    ax[0].errorbar(bins[:-1], mean_ns_iso, yerr=std_ns_iso, color='blue', capsize=1.0, linewidth=1.0, ds='steps-mid', label='NS | Isotropic Emission')
    ax[0].errorbar(bins[:-1], mean_bh_iso, yerr=std_bh_iso, color='black', capsize=1.0, linewidth=1.0, ds='steps-mid', label='BH | Isotropic Emission')
    
    ax[1].errorbar(bins[:-1], mean_ns_lx1, yerr=std_ns_lx1, color='blue', capsize=1.0, linewidth=1.0, ds='steps-mid', label='NS | Beamed Emission')
    ax[1].errorbar(bins[:-1], mean_bh_lx1, yerr=std_bh_lx1, color='black', capsize=1.0, linewidth=1.0, ds='steps-mid', label='BH | Beamed Emission')
    
    ax[2].errorbar(bins[:-1], mean_ns_obs_b_d, yerr=std_ns_obs_b_d, color='blue', capsize=1.0, linewidth=1.0, ds='steps-mid', label=r'NS | Observed $(b \times d)$')
    ax[2].errorbar(bins[:-1], mean_bh_obs_b_d, yerr=std_bh_obs_b_d, color='black', capsize=1.0, linewidth=1.0, ds='steps-mid', label=r'BH | Observed $(b \times d)$')

    ax[3].errorbar(bins[:-1], mean_obs_prec_ns, yerr=std_obs_prec_ns, color='blue', capsize=1.0, linewidth=1.0, ds='steps-mid', label=r'NS | Observed Precession')
    ax[3].errorbar(bins[:-1], mean_obs_prec_bh, yerr=std_obs_prec_bh, color='black', capsize=1.0, linewidth=1.0, ds='steps-mid', label=r'BH | Observed Precession')
        
    for a in ax:
        a.legend(loc='upper right', prop={'size': 7})
        a.set_yscale('log')
        a.set_ylim(1)

    # plot_area(ax[0], bin_centers, mean_ns_iso, mean_ns_iso, label='NS | Isotropic Emission', color='blue', alpha=0.8)
    # plot_area(ax[0], bin_centers, mean_bh_iso, std_bh_iso, label='BH | Isotropic Emission', color='black', alpha=0.8)
    
    # plot_area(ax[1], bin_centers, mean_ns_lx1, std_ns_lx1, label='NS | Beamed Emission', color='blue', alpha=0.8)
    # plot_area(ax[1], bin_centers, mean_bh_lx1, std_bh_lx1, label='BH | Beamed Emission', color='black', alpha=0.8)
    
    # plot_area(ax[2], bin_centers, mean_ns_obs_b_d, std_ns_obs_b_d, label=r'NS | Observed $(b \times d)$', color='blue', alpha=0.8)
    # plot_area(ax[2], bin_centers, mean_bh_obs_b_d, std_bh_obs_b_d, label=r'BH | Observed $(b \times d)$', color='black', alpha=0.8)
    
    # plot_area(ax[3], bin_centers, mean_ns_obs_b_d, std_ns_obs_b_d, label='NS | Observed b+d+p', color='blue', alpha=0.8)
    # plot_area(ax[3], bin_centers, mean_bh_obs_b_d, std_bh_obs_b_d, label='BH | Observed b+d+p', color='black', alpha=0.8)
    
    
    # ax[0].errorbar(bin_centers, mean_ns_iso, yerr=std_ns_iso, color='blue', linestyle='', capsize=1.0)
    # ax[0].errorbar(bin_centers, mean_bh_iso, yerr=std_bh_iso, color='black', linestyle='', capsize=1.0)
    
    # ax[1].errorbar(bin_centers, mean_ns_lx1, yerr=std_ns_lx1, color='blue', linestyle='', capsize=1.0)
    # ax[1].errorbar(bin_centers, mean_bh_lx1, yerr=std_bh_lx1, color='black', linestyle='', capsize=1.0)
    
    # ax[2].errorbar(bin_centers, mean_ns_obs_b_d, yerr=std_ns_obs_b_d, color='blue', linestyle='', capsize=1.0)
    # ax[2].errorbar(bin_centers, mean_bh_obs_b_d, yerr=std_bh_obs_b_d, color='black', linestyle='', capsize=1.0)
    if save:
        plt.savefig('../reports/figures/XLF_cumsum.png', dpi=1000)
        plt.savefig('../reports/figures/XLF_cumsum.eps')
        plt.savefig('../reports/figures/XLF_cumsum.pdf')
    
    
    from scipy.optimize import curve_fit
    
    def powlaw(L, A, gamma):
        return A * L**gamma
    
    
    # =============================================================================
    # Other
    # =============================================================================
    
    # print('removing NS <--> BH systems')
    # #Find & remove systems that didn't spend entire lifetime as NS or BH (161 binaries)
    # K_a_mean = gb['K_a'].mean()
    # idx = K_a_mean[(K_a_mean != 13.0) & (K_a_mean != 14.0)].index
    # df = df.drop(idx, axis=0)
    
    
    # print('Getting transitioning mttype systems')
    # # Systems that don't spend their entire time in one mass transfer state.
    # mttype_mean = gb['mttype'].mean()
    # idx = mttype_mean[(mttype_mean != 1.0) & (mttype_mean != 4.0)].index
    # variable = df.loc[idx]
    
    # gb = pop.gb_sys(variable)
    
    
    """
    for i, df in tqdm(gb):
        samp = df.sample(n=N_samples, replace=True)
        samp['dincl'] = np.random.randint(0, 46, size=N_samples)
        samp['inclination'] = np.random.randint(0,91, size=N_samples)
        samp['sample_id'] = np.arange(0, N_samples, 1)
        samp['system_row_id'] = samp.index
    
    """
    
    
if __name__=="__main__":
    set_latex_font()
    df = populations.startrack_v2_mt_1_all() #nrows=50000
    pop = populations.Population(df)
    main(pop, save=True)