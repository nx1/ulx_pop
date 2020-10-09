#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Apr 24 11:52:31 2019

@author: nk7g14
"""
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.animation as animation
import numpy as np

import os

import sys
sys.path.append('../../')

import curvefunc


def init():
    '''
    Cleans off the graph sometimes that's a lot more important than other times
    otherwise you might get a line stuck there -JL
    '''
    d_line.set_ydata([np.nan]*len(d_t))
    i_line.set_ydata([np.nan]*len(i_t))
    
    d_text.set_text('')
    i_text.set_text('')
    return d_line, 

def animate(i):
    d_filename = 'dincl_plot/'+str(i)+'.txt'
    i_filename = 'Inclination_plot/'+str(i)+'.txt'
    
    d_curve = curvefunc.load_curve_file(d_filename)
    i_curve = curvefunc.load_curve_file(i_filename)
    
    d_line.set_ydata(d_curve['Flux'])
    d_line.set_xdata(d_curve['Time'])
    
    i_line.set_ydata(i_curve['Flux'])
    i_line.set_xdata(i_curve['Time'])
    
    d_text.set_text(rf'Precessional angle $\Delta i = {{}}$'.format(i%len(inclination_range)))
    i_text.set_text(rf'Inclination $i = {{}}$'.format(i))
    
    return d_line, i_line, d_text, i_text





import matplotlib
import itertools

matplotlib.rcParams['mathtext.fontset'] = 'stix'
matplotlib.rcParams['font.family'] = 'STIXGeneral'

animation_fps = 5

dincl_range = np.arange(0,46,1)
inclination_range = np.arange(0,90,1)

grid_variables = [dincl_range, inclination_range]
gridsearch = itertools.product(*grid_variables)

# number_of_frames = len(gridsearch)

for i, j in gridsearch:
    print(i, j)
    
    
curves = []


for i in gridsearch:
    print(i)
    
    

frames = []


inclination = 30
dincl = 15


# Load Curves
dincl_file = f'dincl_plot/{inclination}.txt'
inclination_file = f'Inclination_plot/{dincl}.txt'

d_curve = curvefunc.load_curve_file(dincl_file)
i_curve = curvefunc.load_curve_file(inclination_file)

# Set up Axis
fig, ax = plt.subplots(2,1, figsize=(10,6), sharex=True)
d_ax, i_ax = ax[0], ax[1]
d_ax.set_ylim(0,40)
i_ax.set_ylim(0,40)



# Plot Labels

d_ax.set_xlabel('Time')

d_ax.set_ylabel(r'Luminosity')
i_ax.set_ylabel(r'Luminosity')

i_ax.set_xlabel('Time')




plt.style.use('dark_background')
plt.tight_layout()


# Animated variables
d_t = d_curve['Time']
d_flux = d_curve['Flux']

i_t = i_curve['Time']
i_flux = i_curve['Flux']


d_line, = d_ax.plot(d_t, d_flux, c='w')
i_line, = i_ax.plot(i_t, i_flux, c='w')

d_text = d_ax.text(0, 35, '', fontsize=15)
i_text = i_ax.text(0, 35, '', fontsize=15)



ani = animation.FuncAnimation(fig, animate, init_func = init, interval = 1000/animation_fps, blit=True, repeat=True, frames=45)
# Set up formatting for the movie files
plt.rcParams['animation.ffmpeg_path'] = 'ffmpeg-20200623-ce297b4-win64-static/ffmpeg-20200623-ce297b4-win64-static/bin/ffmpeg.exe'
FFwriter = animation.FFMpegWriter(fps=30)
ani.save('basic_animation.mp4', writer=FFwriter)



