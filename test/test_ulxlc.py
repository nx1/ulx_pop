#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
@author: nk7g14

test_batchfunc.py

Test file for batchfunc.py which contains files for running ulxlc

# ULXLC

# C Methods
# Dauser
ulxlc._set_check_par_model()
ulxlc._init_par_model()
ulxlc._interp_betTheta()
ulxlc._calc_orbit_ang()
ulxlc._get_flux_from_angle()
ulxlc._create_lc()
ulxlc._ulxlc_base()
ulxlc._free_ulxTable()

# Khan
ulxlc.max()
ulxlc.min()
ulxlc.get_first_transient_cycle_index()
ulxlc.callback()

# Python Methods
    ulxlc.get_lc_max()
    ulxlc.get_lc_min()
    calc_lc_flux_scaling_constant()
    calc_lc_ulx_lim()
    classify_curve()
    grid_ulxlc_model()
    ulxlc.ulxlc_model()
    ulxlc.set_params()
    ulxlc.classify_curve()
    
"""
import numpy as np
import sys
import pytest
import os
import time
import matplotlib.pyplot as plt
import ctypes
import subprocess

sys.path.insert(0, '../src')


from ulxlc import ULXLC, MC_input, MC_output
from constants import params_default, incl_default

os.chdir('../src')



@pytest.fixture
def ulxlc():
   ulxlc = ULXLC(100, 1.0)
   ulxlc.set_params(*params_default)
   return ulxlc


def test_init(ulxlc):
    assert ulxlc.lc_t[0] == 0.0
    assert ulxlc.lc_t[-1] == 99.00

def test_set_param(ulxlc):    
    ulxlc.set_theta(3.44)
    assert ulxlc.params[2] == 3.44
    
    ulxlc.set_inclination(76.3)
    assert ulxlc.params[3] == 76.3
    
    ulxlc.set_dincl(23.3)
    assert ulxlc.params[4] == 23.3
    

def test_get_lc_max(ulxlc):
    ulxlc.ulxlc_model()
    assert ulxlc.get_lc_max() == pytest.approx(35.26, 0.1)

def test_get_lc_min(ulxlc):
    ulxlc.ulxlc_model()
    assert ulxlc.get_lc_min() == pytest.approx(1.19, 0.1)
    
def test_calc_flux_scaling_constant(ulxlc):
    ulxlc.set_inclination(0)
    ulxlc.ulxlc_model()
    lc_zero_incl_max_flux = ulxlc.get_lc_max()
    
    ulxlc.set_inclination(incl_default)
    ulxlc.ulxlc_model()
    lc_fsc = ulxlc.calc_lc_flux_scaling_constant(lc_zero_incl_max_flux, 1.3e39)
    assert lc_fsc == pytest.approx(3.69e+37, 0.1e37)

def test_calc_lc_ulx_lim(ulxlc):
    ulxlc.ulxlc_model()
    lc_ulx_lim = ulxlc.calc_lc_ulx_lim(3.686509317961678e+37)
    assert lc_ulx_lim == pytest.approx(27.13, 0.1)
    
def test_classify_curve(ulxlc):
    ulxlc.ulxlc_model()
    lc_ulx_lim = ctypes.c_double(27.125931708017863)
    lc_max = ctypes.c_double(35.26444571367667)
    lc_min = ctypes.c_double(1.1942565810840877)
    
    assert ulxlc.classify_curve(lc_ulx_lim,lc_max,lc_min) == 1
    assert ulxlc.classify_curve(ctypes.c_double(0.25),lc_max,lc_min) == 2
    assert ulxlc.classify_curve(ctypes.c_double(37),lc_max,lc_min) == 0

def test_ulxlc_model(ulxlc):
    ulxlc.ulxlc_model()
    assert ulxlc.lc_flux[50] == pytest.approx(1.19, 0.1)
    ulxlc.set_inclination(30)
    ulxlc.ulxlc_model()
    assert ulxlc.lc_flux[50] == pytest.approx(0.15, 0.1)
    
def test_calc_lc_flux_scaling_constant(ulxlc):
    r"""
    C = Lx / max(lc_i=0)
    $C = L_{x} / \mathrm{max}(L_{\mathrm{curve, \ i=0}})$
    
    i vs C
    dincl vs C
    
    lc_i=0:
    lc_i=

    
    Example: Lx = 1.43e39, max(lc_i=43) < max(lc_i=0) = 12.2,
           C = 1.43e39 / 12.2 = 1.17e+38 (1.1721311475409836e+38)
           C = 2.37e40 / 0.54 = 4.39e40 (4.388888888888889e+40)
    
    C          = Lx         /  max(curve)
    [erg s^-1] = [erg s^-1] / [None]
    [10e39 * erg s^-1] = [10e39 * erg s^-1] / [None] <-- In old code
    
    None.
    """
    Lx = 2e39
    ulxlc.set_inclination(0)
    ulxlc.ulxlc_model()
    lc_zero_incl_max_flux = ulxlc.get_lc_max()
    lc_fsc = ulxlc.calc_lc_flux_scaling_constant(lc_zero_incl_max_flux, Lx)
    lc_ulx_lim = ulxlc.calc_lc_ulx_lim(lc_fsc)
    ulxlc.set_inclination(15)
    
    assert lc_fsc == pytest.approx(5.67e+37, 0.1e37)
    assert lc_ulx_lim == pytest.approx(17.63, 0.1)


def test_grid_ulxlc_model(ulxlc):
    N = 2
    theta = np.random.normal(20, size=N)
    theta = (ctypes.c_double * N)(*theta)
    
    Lx = np.random.normal(5e39, scale=1e38, size=N)
    Lx = (ctypes.c_double * N)(*Lx)
    
    ulxlc.libc.grid_ulxlc_model(theta, Lx, N, ulxlc.lc_t, ulxlc.lc_nt, ulxlc.lc_flux, ulxlc.params, 7)

def test_dincl_vs_lc(ulxlc):
    dincls = range(0,46)
    
    lc_max = []
    lc_min = []
    for dincl in dincls:
        ulxlc.set_dincl(dincl)
        ulxlc.ulxlc_model()
        lc_max.append(ulxlc.get_lc_max())
        lc_min.append(ulxlc.get_lc_min())

    # Plot 
    plt.figure()    
    plt.plot(dincls, lc_min)
    plt.title(f'params={list(ulxlc.params)}')
    plt.xlabel(r'Precessional angle ($\Delta i$)')
    plt.ylabel('lc_min')
    plt.savefig('../test/figures/dincl_vs_lc_min.png')
    plt.savefig('../test/figures/dincl_vs_lc_min.eps')
    
    plt.figure()
    plt.plot(dincls, lc_max)
    plt.title(f'params={list(ulxlc.params)}')
    plt.xlabel(r'Precessional angle ($\Delta i$)')
    plt.ylabel('lc_max')
    plt.savefig('../test/figures/dincl_vs_lc_max.png')
    plt.savefig('../test/figures/dincl_vs_lc_max.eps')

def test_incl_vs_lc(ulxlc):
    inclinations = range(0,91)
    
    lc_max = []
    lc_min = []
    for i in inclinations:
        ulxlc.set_inclination(i)
        ulxlc.ulxlc_model()
        lc_max.append(ulxlc.get_lc_max())
        lc_min.append(ulxlc.get_lc_min())

    # Plot 
    plt.figure()    
    plt.plot(inclinations, lc_min)
    plt.title(f'params={list(ulxlc.params)}')
    plt.xlabel('Inclination (i)')
    plt.ylabel('lc_min')
    plt.savefig('../test/figures/inclination_vs_lc_min.png')
    plt.savefig('../test/figures/inclination_vs_lc_min.eps')
    
    plt.figure()
    plt.plot(inclinations, lc_max)
    plt.title(f'params={list(ulxlc.params)}')
    plt.xlabel('Inclination (i)')
    plt.ylabel('lc_max')
    plt.savefig('../test/figures/inclination_vs_lc_max.png')
    plt.savefig('../test/figures/inclination_vs_lc_max.eps')
            
def test_dincl_vs_incl_vs_lc(ulxlc):
    dincls = np.arange(0,46)
    inclinations = np.arange(0,91)
    lc_max = np.ndarray((46,91))
    lc_min = np.ndarray((46,91))
    
    for i, dincl in enumerate(dincls):
        ulxlc.set_dincl(dincl)
        for j, incl in enumerate(inclinations):
            ulxlc.set_inclination(incl)
            ulxlc.ulxlc_model()
            lc_max[i][j] = ulxlc.get_lc_max()
            lc_min[i][j] = ulxlc.get_lc_min()
    
    plt.figure()
    plt.imshow(lc_max)
    plt.title(f'lc_max | params={list(ulxlc.params)}')
    plt.xlabel('Inclination (i)')
    plt.ylabel(r'Precessional angle ($\Delta i$)')
    plt.colorbar()
    plt.savefig('../test/figures/dincl_incl_vs_lc_max.png')
    
    
    plt.figure()    
    plt.imshow(lc_min)
    plt.title(f'lc_min | params={list(ulxlc.params)}')
    plt.xlabel('Inclination (i)')
    plt.ylabel(r'Precessional angle ($\Delta i$)')
    plt.colorbar()
    plt.savefig('../test/figures/dincl_incl_vs_lc_min.png')
    
def test_c_min_max(ulxlc):
    """Check if min and max functions work same in python and C."""
    ulxlc.ulxlc_model()
    assert ulxlc.get_lc_max() == ulxlc.libc.max(ulxlc.lc_flux, ulxlc.lc_nt)
    assert ulxlc.get_lc_min() == ulxlc.libc.min(ulxlc.lc_flux, ulxlc.lc_nt)

def test_calc_Lx_prec(ulxlc):
    lc_max_flux_zero_incl = ctypes.c_double(35)
    lc_flux = ctypes.c_double(21)
    Lx = ctypes.c_double(3e39)
    
    Lx_prec = ulxlc.libc.calc_Lx_prec(Lx, lc_max_flux_zero_incl, lc_flux)
    assert Lx_prec == 1.8e+39
    
def test_xlf_calc_L_prec(ulxlc):
    N = 5
    c_Lx_prec = (ctypes.c_double * N)()
    c_Lx      = (ctypes.c_double * N)(*np.random.normal(loc=2e39, scale=1e38, size=N))
    c_thetas  = (ctypes.c_double * N)(*np.random.randint(low=2.25e-01, high=46, size=N))
    c_incls   = (ctypes.c_double * N)(*np.random.randint(low=0, high=91, size=N))
    c_dincls  = (ctypes.c_double * N)(*np.random.randint(low=0, high=46, size=N))
    
    ulxlc.xlf_calc_L_prec(c_Lx_prec, c_Lx, c_thetas, c_incls, c_dincls, N)
    
    for n in range(N):
        assert c_Lx_prec[n] != 0.0

def test_sim(ulxlc):
    N_sys = ulxlc.get_N_sys()     #Number of systems, must be the same as in C code.
    N_double_arr = N_sys * ctypes.c_double
    
    # Input params
    c_id     = (N_double_arr)(*np.random.randint(low=0, high=9000000, size=N_sys))
    c_theta  = (N_double_arr)(*np.random.randint(low=7, high=90, size=N_sys))
    c_incl   = (N_double_arr)(*np.random.randint(low=0, high=91, size=N_sys))
    c_dincl  = (N_double_arr)(*np.random.randint(low=0, high=46, size=N_sys))
    c_Lx     = (N_double_arr)(*np.random.normal(loc=2e39, scale=1e38, size=N_sys))
    c_period = (N_double_arr)(*np.random.randint(low=1, high=1000, size=N_sys))
    c_phase  = (N_double_arr)(*np.random.random(size=N_sys))


    # Pass to params to structs
    inp = MC_input(c_id, c_theta, c_incl, c_dincl, c_Lx, c_period, c_phase)
    out = MC_output.initialize()
    
    ulxlc.libc.sim(ctypes.byref(inp), ctypes.byref(out))
    
    
    
    

# def test_calc_scaling_constant_grid():
#     # Grid values
#     incl_min = 0
#     incl_max = 91
#     incl_step = 1
    
#     dincl_min = 0
#     dincl_max = 46
#     dincl_step = 1
    
#     theta_min = 0
#     theta_max = 46
#     theta_step = 1
    
    # incls  = np.arange(0,91,1)
    # dincls = np.arange(0,46,1)
    # thetas = np.arange(0,46,1)
    
    # N_incls  = len(incls)
    # N_dincls = len(dincls)
    # N_thetas = len(thetas)
    
    # N_tot = N_incls*N_dincls*N_thetas
    
    # save_par = ['dincl', 'theta', 'incl', 'lc_min', 'lc_max', 'lc_boost']
    # N_save_par = len(save_par)
    
    # # arr = np.ndarray((N_tot, N_save_par))
    # # c_arr = ((ctypes.c_double * N_save_par) * N_tot)()
    
    
    # ulxlc.libc.lc_boost(ulxlc.lc_t, ulxlc.lc_nt, ulxlc.lc_flux, N_tot, N_save_par)
    # # assert c_arr = 

# def test_ulxlc_model_bad_param():
#     ulxlc = ULXLC()
#     ulxlc.set_params('a', phase, theta, incl, dincl, beta, dopulse)
#     ulxlc.set_params(-5, phase, theta, incl, dincl, beta, dopulse)
#     ulxlc.set_params(period, phase, 95, incl, dincl, beta, dopulse)
#     ulxlc.ulxlc_model()

                   
    

