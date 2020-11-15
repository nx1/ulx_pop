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



df = populations.startrack_v2_mt_1_all()

cols = ['idum_run', 'iidd_old', 'K_a', 'Lx_iso', 'Lx1', 'theta_half_deg']

pop = populations.Population(df)
df = pop.df[cols]

gb = pop.gb_sys(df)

N_samples = 10000

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
