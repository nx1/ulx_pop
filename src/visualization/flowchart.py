#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Sep 10 11:15:19 2019

@author: nk7g14
"""


import matplotlib.pyplot as plt
import matplotlib.patches as patches
import matplotlib.lines as mlines


N_dict = {'All': [36420, 2069, 34351],
 'ULX': [992, 913, 79],
 'Unbeamed_ULX': [765, 765, 0],
 'beamed_ULX': [227, 148, 79],
 'Beamed<45': [151, 88, 63],
 'Beamed>45': [76, 60, 16]}

fig,ax = plt.subplots(1)
xmax = 200
ymax = 100
ax.set_xlim([0,xmax])
ax.set_ylim([0,ymax])

MIDDLE = xmax/2
TOP = ymax


def DrawRect(xpos, ypos, xsize, ysize, ratio):
    rect1_xsize = ratio * xsize
    rect2_xsize = (1 - ratio) * xsize
    
    rect1_xpos_true = xpos - (xsize/2)
    rect1_ypos_true = ypos - (ysize/2)
    
    rect2_xpos_true = rect1_xpos_true + rect1_xsize
    rect2_ypos_true = rect1_ypos_true
    BH_COLOR = '#119cd8'
    NS_COLOR = '#3c3c3c'
    rect1 = patches.Rectangle((rect1_xpos_true, rect1_ypos_true),rect1_xsize,ysize,linewidth=1, 
                             edgecolor='black', facecolor=NS_COLOR)
    
    rect2 = patches.Rectangle((rect2_xpos_true, rect2_ypos_true),rect2_xsize,ysize,linewidth=1, 
                             edgecolor='black', facecolor=BH_COLOR)
    ax.add_patch(rect1)
    ax.add_patch(rect2)

xsize = 40
ysize = 8

r1 = N_dict['All'][1]/N_dict['All'][0]
r2 = N_dict['ULX'][1]/N_dict['ULX'][0]
r3 = N_dict['Unbeamed_ULX'][1]/N_dict['Unbeamed_ULX'][0]
r4 = N_dict['beamed_ULX'][1]/N_dict['beamed_ULX'][0]
r5 = N_dict['Beamed<45'][1]/N_dict['Beamed<45'][0]
r6 = N_dict['Beamed>45'][1]/N_dict['Beamed>45'][0]

#TOP RECTANGLE ALL SYSTEMS
DrawRect(MIDDLE, TOP-10, xsize, ysize, ratio=r1)
# ax.text(MIDDLE, TOP-10, N_dict['All'][0], fontsize=10, horizontalalignment='center',
        # verticalalignment='center', color='white')

# ax.arrow(MIDDLE, TOP-15, 0, -5, width=1)
# ax.text(MIDDLE+2, TOP-21, '$> 10^{39} \ \mathrm{erg \ s^{-1}}$', fontsize=10,)

#SECOND RECTANGLE ULXS
DrawRect(MIDDLE, TOP-30, xsize, ysize, ratio=r2)
# ax.text(MIDDLE, TOP-30, N_dict['ULX'][0], fontsize=10, horizontalalignment='center',
#         verticalalignment='center', color='white')

# ax.axvline(MIDDLE, ymin=58/ymax, ymax=65/ymax, c='black')


#RIGHT ROW3 RECTANGLE UNBEAMED ULXS
DrawRect(MIDDLE+35, TOP-50, xsize, ysize, ratio=r3)
# ax.text(MIDDLE+35, TOP-50, N_dict['Unbeamed_ULX'][0], fontsize=10, horizontalalignment='center',
#         verticalalignment='center', color='white')

#LEFT ROW3 RECTANGLE BEAMED ULXS
DrawRect(MIDDLE-35, TOP-50, xsize, ysize, ratio=r4)
# ax.text(MIDDLE-35, TOP-50, N_dict['beamed_ULX'][0], fontsize=10, horizontalalignment='center',
#         verticalalignment='center', color='white')

#LEFT ROW4 RECTANGLE BEAMED ULXS < 45
DrawRect(MIDDLE-70, TOP-70, xsize, ysize, ratio=r5)
# ax.text(MIDDLE-70, TOP-70, N_dict['Beamed<45'][0], fontsize=10, horizontalalignment='center',
#         verticalalignment='center', color='white')

#RIGHT ROW4 RECTANGLE BEAMED ULXS > 45
DrawRect(MIDDLE, TOP-70, xsize, ysize, ratio=r6)
# ax.text(MIDDLE, TOP-70, N_dict['Beamed>45'][0], fontsize=10, horizontalalignment='center',
#         verticalalignment='center', color='white')

plt.savefig('flowchart.eps', format='eps')