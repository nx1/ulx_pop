# -*- coding: utf-8 -*-
"""
Created on Tue Apr 28 11:34:19 2020

@author: norma

"""
import csv
import glob
import sqlite3

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from pathlib import Path


#CGS UNITS
c    = 3E10         # Speed of Light in cm/s
M_sol = 1.989E33     # Mass of sun in g
yr   = 31557600     # 1 yr in Seconds
Myr  = yr * 1E6     # 1 Myr in Seconds

G_SI = 6.67408E-11  # N.(m^2)/(kg)^2 
G    = 6.674E-8     # (cm)^3 g^-1 s^-2

#GR constants
epsilon = 0.25 # Normally between 0.25 to 0.95 (check latex)
beta    = 1.4  # Velocity of the wind, distinct from the beta used in ulxlc
BH_ISCO = 1.25
BH_SPIN = 0.998
NS_ISCO = 6
NS_SPIN = 0.01


class System:

    def __init__(self, system_id, Z, tage, mass, mdot_grams_sec):
        self.system_id = system_id
        self.Z    = Z
        self.tage = tage
        self.mass = mass
        self.mdot_grams_sec = mdot_grams_sec
        
        self.is_bh    = mass >= 2.5

        self.L_Edd    = 1.2E38 * mass
        self.eta      = 1 / 12
        self.mdot_Edd = self.L_Edd / (self.eta * c**2)
        self.mdot     = self.mdot_grams_sec / self.mdot_Edd
        self.Lx_iso   = self.calc_isotropic_luminosity()        
        self.b        = self.calc_beaming_factor()
        self.Lx       = self.Lx_iso / self.b
        
        self.theta          = 2 * np.arccos(1-self.b) 
        self.theta_half_deg = (self.theta * 180 / np.pi) / 2
        
        self.zeta = self.calc_zeta()
        self.spin = BH_SPIN if self.is_bh else NS_SPIN
        
        self.R_g    = (G * self.mass * M_sol) / c**2
        self.r_schw = 2 * self.R_g
        self.r_isco = BH_ISCO if self.is_bh else NS_ISCO
        self.r_in   = self.r_isco
        self.r_sph  = self.r_isco * self.mdot
        self.r_out  = 3 * epsilon / (beta * self.zeta) * self.mdot**3/2 * self.r_isco
        
        self.P_inflow_at_rsph = (G * self.mass * M_sol * np.pi) / (3 * c**3 * self.spin) * self.r_sph**3 * ((1 - (self.r_isco/self.r_sph))**3)/(np.log(self.r_sph/self.r_isco))
        self.P_envelope       = (G * self.mass * M_sol * np.pi) / (3 * c**3 * self.spin) * self.r_sph**3 * ((1 - (self.r_isco/self.r_sph))**3)/(np.log(self.r_sph/self.r_isco)) * (self.r_out/self.r_sph)**2
        self.P_wind           = (G * self.mass * M_sol * np.pi) / (3 * c**3 * self.spin) * self.r_out**3 * ((1 - (self.r_isco/self.r_out))**3)/(np.log(self.r_out/self.r_isco))
        
        if not self.is_bh:
            self.free_precession_period = self.calc_NS_free_precession()
    
    def __repr__(self):
        return repr(f'System: system_id: {self.system_id}, mass: {self.mass}, mdot: {self.mdot:2f}')
    
    def calc_beaming_factor(self):
        if self.mdot < 8.5:
            b = 1
        elif self.mdot <= 150:
            b = 73 / self.mdot**2
        else:
            b = 3.2E-2
        return b
    
    def calc_isotropic_luminosity(self):
        if self.mdot > 1:
            Lx_iso = self.L_Edd * (1 + np.log(self.mdot))
        else:
            Lx_iso = self.L_Edd * self.mdot
        return Lx_iso
    
    def calc_zeta(self):
        zeta = np.tan( np.pi / 2 - np.arccos(1 - (73 / (self.mdot**2))) )
        if zeta <= 2:
            zeta = 2
        return zeta
    
    def calc_NS_free_precession(self):
        # See Ravenhall(1994)
        R = 10
        NS_PERIOD = 1
        MOI_angle = 10 * (np.pi / 180)
            
        Lambda_R = (1 - 2.953 * self.mass / R)**-1
        delta_R = 1.74E-2 * R**2 / (self.mass * Lambda_R**1.5)
        Q_I = 1 - 0.62 * self.mass * delta_R / R
        fractional_moi = 9.06E-4 * ((R - delta_R)**6 * Q_I) / (self.mass**2 * R**2 * Lambda_R)
        fractional_moi = fractional_moi / 100
        period = NS_PERIOD * fractional_moi * 1/np.cos(MOI_angle)
        return period
    
    def create_table(self):
        conn = sqlite3.connect('ulx_pop.db')
        c = conn.cursor()
        c.execute("""CREATE TABLE IF NOT EXISTS systems(
            system_id INTEGER PRIMARY_KEY,
            Z REAL,
            tage INTEGER,
            mass REAL,
            mdot_grams_sec REAL,
            
            is_bh INTEGER,
            
            L_Edd REAL,
            mdot_Edd REAL,
            mdot REAL, 
            Lx_iso REAL, 
            b REAL, 
            Lx REAL, 
            
            theta REAL, 
            theta_half_deg REAL, 
            
            zeta REAL,
            spin REAL,
            
            R_g REAL,
            r_schw REAL,
            r_isco REAL,
            r_in REAL,
            r_sph REAL,
            r_out REAL,
            
            P_inflow_at_rsph REAL,
            P_envelope REAL,
            P_wind REAL
            )""")
        
            
    def sql_insert_to_table(cursor, table):
        cursor.execute("""INSERT INTO ? VALUES(
            )""", (table))
    
    def to_series(self):
        series = pd.Series(vars(self))
        return series
    


class Population:
    def __init__(self, systems):
        self.systems = systems
        self.size = len(self.systems)
        self.N_bh = sum([s.is_bh for s in self.systems])
        self.N_ns = self.size - self.N_bh
        self.bh_percent = self.N_bh / self.size
        self.ns_percent = self.N_ns / self.size
        self.parent = None
        self.Z = set([s.Z for s in systems])
        self.tage = set([s.tage for s in systems])
    
        
    def get_bh_systems(self):
        bh_systems = []
        for s in self.systems:
            if s.is_bh:
                bh_systems.append(s)
                
        return bh_systems
        
    def get_ns_systems(self):
        ns_systems = []
        for s in self.systems:
            if not s.is_bh:
                ns_systems.append(s)
        return ns_systems
        
    def get_ulx_systems(self):
        ulx_systems = []
        for s in self.systems:
            if s.Lx > 1E39:
                ulx_systems.append(s)
        return ulx_systems

    def get_beamed_systems(self):
        beamed_systems = []
        for s in ulx_systems:
            if s.b < 1 and s.Lx > 1E39:
                beamed_systems.append(s)
        return beamed_systems
        
    def get_beamed_l_45_systems(self):
        strongly_beamed_systems = []
        for s in ulx_systems:
            if s.theta_half_deg < 45 and s.b < 1 and s.Lx > 1E39:
                strongly_beamed_systems.append(s)
        return strongly_beamed_systems


    def get_systems_by_z(self, Z):
        Z_systems = []
        for s in self.systems:
            if s.Z == Z:
                Z_systems.append(s)
        return Z_systems
    
    def sample_with_bh_percent(self, n, bh_percent):
        ns_percent = 1 - bh_percent
        
        bh_systems = self.get_bh_systems()
        ns_systems = self.get_ns_systems()
        
        bh_weights = [bh_percent/len(bh_systems)]*len(bh_systems)
        ns_weights = [ns_percent/len(ns_systems)]*len(ns_systems)
        
        sampled_systems = np.random.choice([*bh_systems, *ns_systems],
                                           p=[*bh_weights, *ns_weights],
                                           size=n)
        return list(sampled_systems)
    
    def to_dataframe(self):
        system_vars = []
        for s in self.systems:
            system_vars.append(vars(s))
        df = pd.DataFrame(system_vars)
        return df
    
    def to_SQL(self, table_name):
        
        return table_name
    

class Curve:
    def __init__(self, path):
        self.path = Path(path)
        self.filename = self.path.stem
        
        self.parent_uuid = None
        
        # self.set_ulxlc_params()
        
        self.max  = None
        self.min  = None
        self.mean = None
        
        self.time = []
        self.flux = []
    
    def load_curve_data(self):
        curve = pd.read_csv(self.path, delimiter=' ', header=None,
                        names=['time', 'time_err', 'flux'], skiprows=3)
        self.time = curve['time'].to_numpy()
        self.flux = curve['flux'].to_numpy()
        
    
    def set_ulxlc_params(self, period, phase, theta, inclination, dincl, beta, dopulse, norm):
        self.period = period
        self.phase = phase
        self.theta = theta
        self.inclination = inclination
        self.dincl = dincl
        self.beta = beta
        self.dopulse = dopulse
        self.norm = norm
    
    def set_ulxlc_params_from_filename(self):
        theta, dincl, i = self.filename.split('-')
        self.set_ulxlc_params(self, 10, 0.0, float(theta),
                         float(i), float(dincl), 0.2, 0, 1.0)
    
    def get_ulxlc_parmeters(self):
        ulxlc_parameters = {}
        ulxlc_parameters['period'] = self.period
        ulxlc_parameters['phase'] = self.phase
        ulxlc_parameters['theta'] = self.theta
        ulxlc_parameters['inclination'] = self.inclination
        ulxlc_parameters['dincl'] = self.dincl
        ulxlc_parameters['beta'] = self.beta
        ulxlc_parameters['dopulse'] = self.dopulse
        ulxlc_parameters['norm'] = self.norm
        return ulxlc_parameters
    
    
    def plot(self):
        plt.figure()
        plt.plot(self.time, self.flux)
        plt.show()



def load_full_population():
    systems = []
    with open('../data/processed/startrack_concat.csv', 'r') as f:
        reader = csv.reader(f)
        next(reader)
        for row in reader:
            uuid = int(row[0])
            mdot = float(row[1]) * (M_sol/Myr)
            mass = float(row[2])
            Z = float(row[3])
            tage = float(row[4])
            systems.append(System(uuid, Z, tage, mass, mdot))
            
    pop_all = Population(systems)
    return pop_all

class eRASS_simulation:
    def __init__(self):
        self.parent_population = None
        self.number_of_sampled_systems = 500
        self.number_of_repeats = 100
        
        self.number_of_lightcurve_samples = 10000
        self.erass_cycle_duration_days = 6*30
        self.erass_number_of_cycles = 8
        
    def run():
        pass

    def get_curve(System, inclination, dincl):
    
        def get_ulxlc_parameters(System):
            """Create parameters dictionary from df_a row"""
            parameters = {'period': 50,
                            'phase': 0,
                            'theta': round(System.theta_half_deg, 2),
                            'inclination': inclination,
                            'dincl': dincl,
                            'beta': 0.2,
                            'dopulse': 1,
                            'norm': 1}
            return parameters

        P_wind_days = System.P_wind / (24*60*60)
        
        parameters = get_ulxlc_parameters(System)
        filename = str(row['theta']) + '-' + str(row['inclination']) +'-' + str(row['dincl'])
        
        xcm_file = 'eRASS_sims/' + filename + '.xcm'
        lc_file = 'eRASS_sims/' + filename + '.txt'
        
        if not path.exists(lc_file):
            run_ulxlc(xcm_file, parameters, lc_file)
    
        curve = load_curve_file(lc_file)
        curve = scale_light_curve_period(curve, 50, P_wind)
        return curve
    
    def scale_six_months(self):
        """Scale 6 months to the length of the curve"""
        six_month_time = curve.period * (self.erass_cycle_duration_days / P_wind_days)
        return six_month_time

    def run():
        create_lightcurves()
        calculate_erass_probabilities()


if __name__ == "__main__":
    df_master = pd.read_csv('../data/processed/startrack_concat.csv')
    path =  '../data/interim/curves/MC_curves_eta_0.08_ns/gridsearch\\11.95-18-84.txt'
    filename = Path(path).stem
    theta, dincl, inclination = filename.split('-')
    
    curve = Curve(path)
    import tqdm
    pop_all = load_full_population()
    ulx_systems = pop_all.get_ulx_systems()
    beamed = pop_all.get_beamed_systems()
    s_beamed = pop_all.get_beamed_l_45_systems()
    
    pop_ulx = Population(ulx_systems)
    df_ulx = pop_ulx.to_dataframe()
    df_all = pop_all.to_dataframe()
    
    curve_files = glob.glob('../data/interim/curves/MC_curves_eta_0.08_ns/gridsearch/*.txt')
    
    
    
    