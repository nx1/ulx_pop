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
    
    xcm_filename: name of xcm file with extension
    parameters: dictionary of model parameters
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
cpd /null
plot model
setplot command wdata {lightcurve_filename}
plot
exit'''

    F.write(string)
    F.close()


def run_xcm(xcm_filename):
    '''
    XCMfile: XCM file including extension
    '''
    devnull = open(os.devnull, 'w')
    subprocess.call([f'xspec - {xcm_filename}'], shell=True, stdout=devnull)

def append_parameters_to_lightcurve_files(lightcurve_filename, system_id ,parameters):
    with open(lightcurve_filename, "a") as myfile:
        myfile.write(f"system_id:{system_id}\n")
        for k, v in parameters.items():
            myfile.write(f'{k}: {v}\n')

def run_ulxlc(xcm_filename, parameters, system_id, lightcurve_filename, append_to_file=False):
    make_xcm(xcm_filename, parameters, lightcurve_filename)
    run_xcm(xcm_filename)
    if append_to_file:
        append_parameters_to_lightcurve_files(lightcurve_filename, system_id ,parameters)
    os.remove(xcm_filename)
