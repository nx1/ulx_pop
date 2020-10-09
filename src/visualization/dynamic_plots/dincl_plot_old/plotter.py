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


fig, ax = plt.subplots()
ax.set_ylim(-1,50)

curve_init = pd.read_csv('0.txt', delimiter=' ',
              header=None, names=['Time', 'Time_Err', 'Flux'], skiprows=3)

time = curve_init['Time']
flux = curve_init['Flux']

line, = ax.plot(time, flux)
text = ax.text(-1.5, 1.0, '0')

def init():
    '''
    Cleans off the graph sometimes that's a lot more important than other times
    otherwise you might get a line stuck there -JL
    '''
    line.set_ydata([np.nan]*len(time))
    text.set_text('')
    
    plt.axhline(y=15, c='purple', linewidth=0.5, linestyle='--')
    
    for x in range(5):
        plt.axvline(x=8*x+3, c='r', linewidth=0.5)
    return line, 

def animate(i):
    filename = str(i)+'.txt'
    
    curve = pd.read_csv(filename, delimiter=' ',
              header=None, names=['Time', 'Time_Err', 'Flux'], skiprows=3)
    line.set_ydata(curve['Flux'])
    line.set_xdata(curve['Time'])
    
    xlim=ax.get_xlim()
    ylim=ax.get_ylim()
    
    text.set_position((xlim[0]+(xlim[1]-xlim[0])/100,ylim[0]+(ylim[1]-ylim[0])/100))
    text.set_text('dincl:'+str(i))
    
    return line,  text


ani = animation.FuncAnimation(fig, animate, init_func = init, interval = 500, blit=True, repeat=True, frames=45)