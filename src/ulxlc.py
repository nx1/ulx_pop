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

class ULXLC:
    def __init__(self, lc_nt=5000, lc_timestep=0.01):
        self.path = './ulxlc/ulxlc_code_v0.1/ulxlc_shared.so'
        self.libc = ctypes.CDLL(self.path)
        
        self.params= (ctypes.c_double * 7)()
        self.lc_timestep = lc_timestep
        self.lc_nt = lc_nt
        self.lc_t = (ctypes.c_double*self.lc_nt)(*[self.lc_timestep*i for i in range(self.lc_nt)])
        self.lc_flux = (ctypes.c_double * self.lc_nt)()            
    
    def set_params(self, period, phase, theta, incl, dincl, beta, dopulse):
        self.params= (ctypes.c_double * 7)(period, phase, theta, incl, dincl, beta, dopulse)
        
    def set_theta(self, theta):
        self.params[2] = theta
        
    def set_inclination(self, incl):
        self.params[3] = incl
        
    def set_dincl(self, dincl):
        self.params[4] = dincl
    
    def get_lc_max(self):
        return max(self.lc_flux)
    
    def get_lc_min(self):
        return min(self.lc_flux)
    
    def calc_lc_flux_scaling_constant(self, lc_zero_incl_max_flux, Lx):
        lc_flux_scaling_constant = Lx / lc_zero_incl_max_flux;
        return lc_flux_scaling_constant
    
    def calc_lc_ulx_lim(self, lc_flux_scaling_constant):
        lc_ulx_lim = 1e39 / lc_flux_scaling_constant
        return lc_ulx_lim

    def classify_curve(self, lc_ulx_lim, lc_max_flux, lc_min_flux):
        return self.libc.classify_curve(lc_ulx_lim, lc_max_flux, lc_min_flux)
    
    def grid_ulxlc_model(self):
        self.libc.grid_ulxlc_model(self.theta, self.Lx, self.lc_t, self.lc_nt, self.lc_flux, self.params, 7)
    
    def ulxlc_model(self):
        self.libc.ulxlc_model(self.lc_t, self.lc_nt, self.lc_flux, self.params, 7)
        return self.lc_flux
    
    
        
    