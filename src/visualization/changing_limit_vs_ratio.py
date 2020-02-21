# -*- coding: utf-8 -*-

import sys
sys.path.append("..") # Adds higher directory to python modules path.

import numpy as np

from curvefunc import load_curve_file_skip_footer, calc_alive_time
import matplotlib.pyplot as plt

curve = load_curve_file_skip_footer('../data/interim/curves/151_systems_0_inclination_curves/510-24.txt')

limit = 2.92

curve['time_diff'] = curve['Time'].diff()
curve_above = curve[curve['Flux'] > 2.92]
curve_above['index_diff'] = np.insert(np.diff(curve_above.index), 0, 0)





limits = np.linspace(2.7, 3, 100)
ratios = []
for limit in limits:
    alive, dead = calc_alive_time(curve, limit)
    ratio = alive/(alive+dead)
    ratios.append(ratio)
    

plt.plot(limits, ratios)
plt.xlabel('limit')
plt.ylabel('alive/dead ratio')


# curve.plot(x='Time', y='Flux')