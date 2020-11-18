"""
Panel 1 : Unbeamed emission split NS/BH
Panel 2 : Beamed Emission split NS/BH
Panel 3 : Beamed + duty_cycle (observed) split NS/BH
Panel 4 : Beamed + duty_cycle + precession (observed) split NS/BH
Panel 5 : Beamed + duty_cycle + prcession (observed) combined NS/BH different %_BH
"""
import numpy as np
import matplotlib.pyplot as plt
import sqlite3
import subprocess
from uuid import uuid4

from tqdm import tqdm

import populations
from constants import R_sol


def lmxrb(df):
    df['T_eff'] = df['R_a']*R_sol
    return df


db = 'xlf.db'

def table_create_xlf_sampled(db):
    conn = sqlite3.connect(db)
    sql = """CREATE TABLE IF NOT EXISTS XLF_SAMPLED(
             system_row_id INT,
             idum_run INT,
             iidd_old INT,
             K_a INT,
			 Lx_iso REAL,
             Lx1 REAL,
             theta_half_deg REAL,
             dincl INT,
             inclination INT,
             sample_id INT);"""
    conn.execute(sql)
    sql = """CREATE INDEX IF NOT EXISTS idx_xlf_sample_id
             ON XLF_SAMPLED (sample_id);"""
    conn.execute(sql)
    conn.close()


print('loading csv')
df = populations.startrack_v2_mt_1_all()

# cols = ['idum_run', 'iidd_old', 't', 'dt', 'K_a', 'mttype', 'log_Lx_iso', 'log_Lx1', 'theta_half_deg', 'lmxrb']

pop = populations.Population(df)

pop.df = pop.df[(pop.df['K_a']==13) | (pop.df['K_a']==14)] # Only include BH/NS systems

# df = pop.df[cols]

print('removing NS <--> BH systems')
#Find & remove systems that didn't spend entire lifetime as NS or BH (161 binaries)
gb = pop.gb_sys(df)
df = df.set_index(['idum_run', 'iidd_old'])
K_a_mean = gb['K_a'].mean()
idx = K_a_mean[(K_a_mean != 13.0) & (K_a_mean != 14.0)].index
df = df.drop(idx, axis=0)


print('Removing non mttype=1 or 4')
#Remove systems that aren't thermal or nuclear mass transfer
df = df[(df['mttype'] == 1.0) | (df['mttype'] == 4.0)]
gb = pop.gb_sys(df)

print('Getting transitioning mttype systems')
# Systems that don't spend their entire time in one mass transfer state.
mttype_mean = gb['mttype'].mean()
idx = mttype_mean[(mttype_mean != 1.0) & (mttype_mean != 4.0)].index

variable = df.loc[idx]

gb = pop.gb_sys(variable)


s=10
e=20
c=0
for i, df in gb:
    if c>s:
        print(f'SID: {i}')
        N = len(df)
        print(f'samples: {N}')
    
        mttype_uniq = df['mttype'].unique()
        K_a_uniq = df['K_a'].unique()
        K_b_uniq = df['K_b'].unique()
    
        t_min = df['t'].min()
        t_max = df['t'].max()
        
        dt_sum = df['dt'].sum()
        
        print(f'Time extent: {t_max - t_min:.2f} Myr')
        print(f'Time mt on: {dt_sum:.2f} Myr')
        
        print(f'Unique K_a: {K_a_uniq} ({len(K_a_uniq)})')
        print(f'Unique K_b: {K_b_uniq} ({len(K_b_uniq)})')
        print(f'Unique mttype: {mttype_uniq}')
        
        print('')
        
        # print(df)
        fig, ax = plt.subplots(5,1, figsize=(8,10))
        ax[0].plot(df['t'], df['mttype'], label='mttype 1=nuc, 4=thermal')
        ax[0].set_ylabel('mttype')
        
        ax[1].plot(df['t'], df['dMmt_a'], label='mdot_a')
        ax[1].set_ylabel(r'dMmt_a ($M_{\odot} Myr^{-1}$)')
        
        ax[2].plot(df['t'], df['dMmt_b'], label='mdot_b')
        ax[2].set_ylabel(r'dMmt_b ($M_{\odot} Myr^{-1}$)')
        
        ax[3].plot(df['t'], df['M_a'], label='M_a')
        ax[3].set_ylabel(r'$M_{\odot}$')
        
        ax[4].plot(df['t'], df['M_b'], label='M_b')
        ax[4].set_ylabel(r'$M_{\odot}$')
        
        
        
        for a in ax:
            a.legend()
        
    if c>e:
        break
        
    c+=1
    

# N_samples = 10000




# Distribution of system number of system samples
# c = gb.count()
# c['K_a'].hist(bins=np.arange(0,1000,1))



"""
table_create_xlf_sampled(db)
 
for i, df in tqdm(gb):
    samp = df.sample(n=N_samples, replace=True)
    samp['dincl'] = np.random.randint(0, 46, size=N_samples)
    samp['inclination'] = np.random.randint(0,91, size=N_samples)
    samp['sample_id'] = np.arange(0, N_samples, 1)
    samp['system_row_id'] = samp.index
    conn = sqlite3.connect(db)
    samp.to_sql('XLF_SAMPLED', conn, if_exists='append', index=False) 
    conn.close()    


class XLF:
    def __init__(self, Population):
        self.pop = Population
        self.percent_bhs = [0.00, 0.25, 0.5, 0.75, 1.00]
        self.max_dincls = [21, 46]

    def sample_luminosities():
        pass
"""