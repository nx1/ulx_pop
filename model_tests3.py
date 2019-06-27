#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jun 26 11:00:59 2019

@author: nk7g14
"""
import numpy as np
from numpy import sin, cos, tan, pi
import matplotlib.pyplot as plt
from matplotlib import animation

from matplotlib.lines import Line2D
import matplotlib.animation as animation


deg2rad = lambda deg : (2*np.pi/360) * deg

def CreateLineFromAngleDegrees(theta):
    theta_rad = deg2rad(theta)
    # print('theta: {}, theta_rad: {}'.format(theta,theta_rad))
    x = []
    y = []
    plt.xlim(-1,1)
    plt.ylim(-1,1)
    #QUADRANT 1
    if theta_rad <= np.pi/2:
        # print('q1')
        m = tan(theta_rad)   
        x = np.arange(-1,1,0.01)
        y = m*x
    #QUADRANT 2
    elif theta_rad > np.pi/2 and theta_rad <= np.pi:
        # print('q2')
        theta_q2 = np.pi - theta_rad #Angle between x axis and line in q2
        m = - tan(theta_q2)
        x = np.arange(-1,1,0.01)
        y = m*x
    #QUADRANT 3
    elif theta_rad > np.pi and theta_rad <= (3/4)*np.pi:
        # print('q3')
        theta_q3 = (3/4)*np.pi - theta_rad#Angle between y axis and line in q3
        m = 1/tan(theta_q3)
        x = np.arange(-1,1,0.01)
        y = m*x
    #QUADRANT 4
    elif theta_rad > (3/4)*np.pi and theta_rad <= 2*np.pi:
        # print('q4')
        theta_q4 = 2*np.pi - theta_rad #Angle between x axis and line in q4
        m = -tan(theta_q4)
        x = np.arange(-1,1,0.01)
        y = m*x
    else:
        print('theta: {}, theta_rad: {}'.format(theta,theta_rad))
        print('erm, something bad happened, we shouldnt be here')
    return x, y

        
def isVisible(eye_position, theta1, theta2):
    
    # a = round(theta1, 2)
    # b = round(theta2, 2)
    
    # print('Checking is visible for eye: {} t1: {} t2: {}'.format(eye_position, a, b))
    if theta2 < theta1: #The cone is overlapping the 0 axis
        #Check if eye is between 0 axis and theta 2
        if eye_position <= theta2:
            return 1
        #Check if eye is between theta 1 and 0 axis
        if eye_position >= theta1:
            return 1
    else: #The cone is not overlapping the 0 axis
        if eye_position > theta1 and eye_position < theta2:
            return 1
        else:
            return 0
    return 0 #This is sometimes activated, i'm not entirely sure on the logic but it works

        
def ObserverPosition(inclination_angle):
    inclination_angle_radians = deg2rad(inclination_angle)
    r = 0.5
    x = r * cos(inclination_angle_radians)
    y = r * sin(inclination_angle_radians)
    return x, y


def init():
    """
    The init function is used to reset the plot elements after every cycle.
    This prevents us plotting the same thing over the top of everything else
    The FuncAnimation requires that the init function return a list of artist objects
    an artist object is anything drawn using matplotlib's Artist module.
    """
    line1a.set_data([], [])
    line1b.set_data([], [])
    text1a.set_text('')
    text1b.set_text('')
    text1c.set_text('')
    
    line2.set_data([], [])
    artist_objects = (line1a, line1b, line2)
    return artist_objects

def animate(t):
    """
    The animate function is essentally provides the information for each
    frame of the animation.
    It is called by the matplotlib function 'FuncAnimation'
    and tells it what to do for each iteration of our variable, in this case t (time)
    """
    global xp1 #Setting global variables for the lightcurve plot
    global yp1 #Setting global variables for the lightcurve plot
    inclination_angle = 20
    opening_angle = 45
    swing_angle = 40
    period = 120

    theta1 = (swing_angle * sin(2*pi/period*t))%360
    theta2 = (theta1+opening_angle)%360
    
    x1, y1 = CreateLineFromAngleDegrees(theta1)
    x2, y2 = CreateLineFromAngleDegrees(theta2)
    
    ax1.collections.clear()
    ax1.fill_between(x1, y1, y2, facecolor='yellow', alpha=0.5)
    ax1.scatter(*ObserverPosition(inclination_angle), c='b')    #Plots observer position
    ax1.scatter(0,0, marker='*', c='b', s=5)
    
    vis = isVisible(inclination_angle, theta1, theta2)
    
    line1a.set_data(x1, y1)
    line1b.set_data(x2, y2)
    text1a.set_text('theta1: {}'.format(str(round(theta1,2))))
    text1b.set_text('theta2: {}'.format(str(round(theta2,2))))
    text1c.set_text('Opening angle: {}'.format(str(round(opening_angle,2))))
    
    xp1 = np.append(xp1, t)
    yp1 = np.append(yp1, vis)
    
    if len(yp1) >359:
        xp1 = np.zeros(0)
        yp1 = np.zeros(0)

    line2.set_data(xp1,yp1)
    ax2.set_xlim(0, 360) #Seting our xlimits
    ax2.set_ylim(0, 1.5) #Setting our ylimits

    artist_objects = (line1a, line1b, line2)
    return artist_objects





#Figure for plotting onto
fig = plt.figure(figsize=(8,4))


#ax1 will hold our spinning cone
ax1 = fig.add_subplot(1,2,1)
ax1.set_xlim(-1, 1) #Seting our xlimits
ax1.set_ylim(-1, 1) #Setting our ylimits
line1a, = ax1.plot([], [], lw=2, color='r') #Line 1a for one of the cone sides
line1b, = ax1.plot([], [], lw=2, color='b') #Line 1b for the other one
text1a = ax1.text(-0.9, 0.8, '', fontsize=7, color='r') #Text for displaying theta1
text1b = ax1.text(-0.9, 0.7, '', fontsize=7, color='b') #Text for displaying theta2
text1c = ax1.text(-0.9, 0.6, '', fontsize=7, color='b') #Text for displaying theta2


#ax2 will show a lightcurve created from our simulation.
ax2 = fig.add_subplot(1,2,2)
ax2.set_xlim(0, 360) #Seting our xlimits
ax2.set_ylim(0, 1.5) #Setting our ylimits
line2, = ax2.plot([], [], lw=2, color='b') #Line on second axis for our lightcurve


#Data placeholders
xp1 = np.zeros(0)
yp1 = np.zeros(0)


anim = animation.FuncAnimation(fig, animate, init_func=init,
                               frames=360, interval=20, blit=False)