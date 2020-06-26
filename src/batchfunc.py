#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
@author: nk7g14

This script is used for running the XSPEC model 'ulxlc' for a variety of
different model parameters.

It requires that one has properly set up local models and installed ulxlc
so that one may use the command "lmod ulxlc" in XSPEC.

This script for some reason does not work in spyder interactive window.
"""
import subprocess
import os

def make_xcm(xcm_filename, parameters, lightcurve_filename):
    '''
    Creates xcm script file for use in XSPEC.
    
    Parameters
    ----------
    xcm_filename : string
        Name of xcm file
    parameters : dictionary 
        ULXLC parameters
    lightcurve_filename : string
        Name of output lightcurve filename
    '''
    F = open(xcm_filename, 'w')


    string = f'''lmod ulxlc
model ulxlc & /*
newpar 1  {parameters['period']}
newpar 2  {parameters['phase']}
newpar 3  {parameters['theta']}
newpar 4  {parameters['inclination']}
newpar 5  {parameters['dincl']}
newpar 6  {parameters['beta']}
newpar 7  {parameters['dopulse']}
newpar 8  {parameters['norm']}
dummy 0 50 5000 lin
cpd /null
plot model
setplot command wdata {lightcurve_filename}
plot
exit'''

    F.write(string)
    F.close()


def run_xcm(xcm_filename):
    '''
    Run xcm file using XSPEC

    Parameters
    ----------
    xcm_filename : string
    '''
    devnull = open(os.devnull, 'w') # Use Null device so no printing to the terminal occurs, siginificantly faster
    subprocess.call([f'xspec - {xcm_filename}'], shell=True, stdout=devnull)


def remove_xcm(xcm_filename):
    try:
        os.remove(xcm_filename)
    except:
        print('xcm file not found')
        
    

def run_ulxlc(xcm_filename, parameters, lightcurve_filename):
    if not os.path.isfile(lightcurve_filename):
        make_xcm(xcm_filename, parameters, lightcurve_filename)
        run_xcm(xcm_filename)
        remove_xcm(xcm_filename)


