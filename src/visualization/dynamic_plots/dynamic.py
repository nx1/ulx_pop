# -*- coding: utf-8 -*-

#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Apr 24 11:52:31 2019

@author: nk7g14
"""
from pathlib import Path

import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.animation as animation
import numpy as np


dincl_path = Path('dincl_plot/')
inclination_path = Path('inclination_plot/')
theta_path = Path('theta_plot/')



def read_curve(path):
    curve = pd.read_csv(path, delimiter=' ', header=None, names=['Time', 'Time_Err', 'Flux'], skiprows=3)
    time = curve['Time']
    flux = curve['Flux']
    return time, flux    

def load_curves(path):
    curves_list = list(path.rglob('*.txt'))
    curves = []
    for p in curves_list:
        curve = read_curve(p)

        stem = int(p.stem)
        print(stem)
        curves.append(curve)
        return curves





fig, ax = plt.subplots(3, 1)
ax[0].set_ylim(-3,40)
ax[0].set_xlim(0,50)
ax[1].set_ylim(-3,40)
ax[1].set_xlim(0,50)
ax[2].set_ylim(-3,40)
ax[2].set_xlim(0,50)

line1, = ax[0].plot([],[])
line2, = ax[1].plot([],[])
line3, = ax[2].plot([],[])


text = ax[0].text(-1.5, 1.0, '0')

def init():
    '''
    Cleans off the graph sometimes that's a lot more important than other times
    otherwise you might get a line stuck there -JL
    '''
    line1.set_ydata([])
    line2.set_ydata([])
    line3.set_ydata([])
    
    return line1, line2, line3, 

def animate(i):
    filename = 'dincl_plot/' + str(i)+'.txt'
    
    curve = pd.read_csv(filename, delimiter=' ',
              header=None, names=['Time', 'Time_Err', 'Flux'], skiprows=3)
    line1.set_ydata()
    line1.set_xdata()
    
    xlim=ax[0].get_xlim()
    ylim=ax[0].get_ylim()
    
    text.set_position((xlim[0]+(xlim[1]-xlim[0])/100,ylim[0]+(ylim[1]-ylim[0])/100))
    text.set_text('inclination:'+str(i))
    
    return line,  text


anim = animation.FuncAnimation(fig, animate, init_func = init, interval = 50, blit=True, repeat=True, frames=89)