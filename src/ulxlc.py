"""
ulxlc.py

This file allows for the running of ulxlc within python using a shared library

/ulxlc/ulxlc_code_v0.1/ulxlc_shared.c contains the file to compile.

It may be compiled using the following command:
    gcc -shared -Wl,-soname,ulxlc -o ulxlc_shared.so -fPIC ulxlcbase.c ulxlc_shared.c -L /home/x1/cfitsio/lib -lcfitsio -lm -lnsl

This create a file called ulxlc_shared.so that may be run using python's
ctypes module.

"""
import ctypes
import numpy as np
import pandas as pd

N_sys = 500     #Number of systems, must be the same as in C code.
N_double_arr = N_sys * ctypes.c_double
eight_int_arr = 8 * ctypes.c_int

class MC_input(ctypes.Structure):
    _fields_ = [("id", N_double_arr),
                ("theta", N_double_arr),
                ("inclination", N_double_arr),
                ("dincl", N_double_arr),
                ("Lx", N_double_arr),
                ("period", N_double_arr),
                ("phase", N_double_arr)]
    
    @classmethod
    def from_df(cls, df_sampled, dincl_max, period):
        N_sys = len(df_sampled)
        N_double_arr = N_sys * ctypes.c_double
        
        c_id     = (N_double_arr)(*df_sampled.index.values)
        c_theta  = (N_double_arr)(*df_sampled['theta_half_deg'].values)
        c_incl   = (N_double_arr)(*np.random.randint(low=0, high=91, size=N_sys))
        c_dincl  = (N_double_arr)(*np.random.randint(low=0, high=dincl_max, size=N_sys))
        c_Lx     = (N_double_arr)(*df_sampled['Lx1'].values)
        c_period = (N_double_arr)(*df_sampled[period].values)
        
        c_phase  = (N_double_arr)(*np.random.random(size=N_sys))
        
        return cls(c_id, c_theta, c_incl, c_dincl, c_Lx, c_period, c_phase)
        
        
    
class MC_output(ctypes.Structure):
    _fields_ = [("N_alive", ctypes.c_int),
                ("N_dead", ctypes.c_int),
                ("N_transient", ctypes.c_int),
                ("N_alive_unsimmed", ctypes.c_int),
                ("N_alive_tot", ctypes.c_int),
                ("N_ulx", eight_int_arr),
                ("N_not_ulx", eight_int_arr),
                ("N_new", eight_int_arr),
                ("N_dip", eight_int_arr),
                ("N_delta_ulx", eight_int_arr),
                ("N_transients", eight_int_arr)]
    
    @classmethod
    def initialize(cls):
        eight_int_arr = 8 * ctypes.c_int
        
        # Output params
        # Classification counts
        c_N_alive = 0
        c_N_dead = 0
        c_N_transient = 0
        c_N_alive_unsimmed = 0
        c_N_alive_tot = 0
        
        # eRASS evolution
        c_N_ulx        = (eight_int_arr)(*np.zeros(8, dtype=int))
        c_N_not_ulx    = (eight_int_arr)(*np.zeros(8, dtype=int))
        c_N_new        = (eight_int_arr)(*np.zeros(8, dtype=int))
        c_N_dip        = (eight_int_arr)(*np.zeros(8, dtype=int))
        c_N_delta_ulx  = (eight_int_arr)(*np.zeros(8, dtype=int))
        c_N_transients = (eight_int_arr)(*np.zeros(8, dtype=int))
        
        return cls(c_N_alive, c_N_dead, c_N_transient, c_N_alive_unsimmed,
                   c_N_alive_tot, c_N_ulx, c_N_not_ulx, c_N_new, c_N_dip,
                   c_N_delta_ulx, c_N_transients)
    
    def _collect_counts(self):
        res = {}
        res['N_alive'] = self.N_alive
        res['N_dead'] = self.N_dead
        res['N_transient'] = self.N_transient
        res['N_alive_unsimmed'] = self.N_alive_unsimmed
        res['N_alive_tot'] = self.N_alive_tot
        return res
    
    def _collect_erass_evolution(self):
        res = {}
        res['N_ulx']        = np.array(self.N_ulx)
        res['N_not_ulx']    = np.array(self.N_not_ulx)
        res['N_new']        = np.array(self.N_new)
        res['N_dip']        = np.array(self.N_dip)
        res['N_delta_ulx']  = np.array(self.N_delta_ulx)
        res['N_transients'] = np.array(self.N_transients)
        
        res['N_new_cum']        = np.cumsum(res['N_new'])
        res['N_dip_cum']        = np.cumsum(res['N_dip'])
        res['N_transients_cum'] = np.cumsum(res['N_transients'])
        return res
    
    def collect(self):
        self.res_counts = self._collect_counts()
        self.res_erass  = self._collect_erass_evolution()

class ULXLC:
    def __init__(self, lc_nt=5000, lc_timestep=0.01):
        self.path = './ulxlc/ulxlc_code_v0.1/ulxlc_shared.so'
        self.libc = ctypes.CDLL(self.path)
        
        self.params= (ctypes.c_double * 7)()
        self.lc_timestep = lc_timestep
        self.lc_nt = lc_nt
        self.lc_t = (ctypes.c_double*self.lc_nt)(*[self.lc_timestep*i for i in range(self.lc_nt)])
        self.lc_flux = (ctypes.c_double * self.lc_nt)()  

        #libc objects default to return ints, this forces the functions we want to return doubles.
        self.libc.max.restype = ctypes.c_double   
        self.libc.min.restype = ctypes.c_double   
        self.libc.lc_boost.restype = ctypes.c_double
        self.libc.calc_Lx_prec.restype= ctypes.c_double
    
    def set_params(self, period, phase, theta, incl, dincl, beta, dopulse):
        self.params= (ctypes.c_double * 7)(period, phase, theta, incl, dincl, beta, dopulse)
        
    def set_theta(self, theta):
        self.params[2] = theta
        
    def set_inclination(self, incl):
        self.params[3] = incl
        
    def set_dincl(self, dincl):
        self.params[4] = dincl
    
    def get_N_sys(self):
        N_sys = self.libc.get_N_sys()
        return N_sys
    
    def get_lc_max(self):
        return max(self.lc_flux)
    
    def get_lc_min(self):
        return min(self.lc_flux)
    
    def xlf_calc_L_prec(self, Lx_prec, lc_classification, Lxs, thetas, incls, dincls, N):
        # Random indexes to sample from the lc (since C rand() is awful.)
        lc_idx = (ctypes.c_int*N)(*np.random.randint(self.lc_nt, size=N))
        self.libc.xlf_calc_L_prec(self.lc_t, self.lc_nt, self.lc_flux, lc_idx, Lx_prec, lc_classification, Lxs, thetas, incls, dincls, N)
    
    def calc_lc_flux_scaling_constant(self, lc_zero_incl_max_flux, Lx):
        lc_flux_scaling_constant = Lx / lc_zero_incl_max_flux;
        return lc_flux_scaling_constant
    
    def calc_lc_ulx_lim(self, lc_flux_scaling_constant):
        lc_ulx_lim = 1e39 / lc_flux_scaling_constant
        return lc_ulx_lim

    def classify_curve(self, lc_ulx_lim, lc_max_flux, lc_min_flux):
        return self.libc.classify_curve(lc_ulx_lim, lc_max_flux, lc_min_flux)
    
    def grid_ulxlc_model(self, theta, Lx, N):
        self.libc.grid_ulxlc_model(theta, Lx, N, self.lc_t, self.lc_nt, self.lc_flux, self.params, 7)
    
    def ulxlc_model(self):
        self.libc.ulxlc_model(self.lc_t, self.lc_nt, self.lc_flux, self.params, 7)
        return self.lc_flux
    
    def sim(self, s_id, theta, incl, dincl, Lx, period, phase, N):
        """
        Parameters
        ----------
        s_id    : array
        theta : array
        incl : array
        dincl : array
        Lx : array
        P : array
        phase : array
        N : int
            array length
        """
        self.libc.sim(s_id, theta, incl, dincl, Lx, period, phase, N)
