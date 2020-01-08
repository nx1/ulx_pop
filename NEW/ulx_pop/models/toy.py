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
from tqdm import tqdm
from matplotlib.lines import Line2D
import matplotlib.animation as animation
plt.style.use('classic')

# Set up formatting for the movie files
# Writer = animation.writers['ffmpeg']
# writer = Writer(fps=20, metadata=dict(artist='Me'), bitrate=1800)

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
    text1d.set_text('')
    text1e.set_text('')
    
    line2.set_data([], [])
    line2b.set_data([], [])
    text2a.set_text('')
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
    inclination_angle = 350
    opening_angle = 30
    half_opening_angle = opening_angle / 2
    b = 1 - cos(deg2rad(half_opening_angle))
    swing_angle = 40
    period = 120

    theta1 = (swing_angle * sin(2*pi/period*t))%360
    theta2 = (theta1+opening_angle)%360
    
    vis = isVisible(inclination_angle, theta1, theta2)
    
    xp1 = np.append(xp1, t)
    yp1 = np.append(yp1, vis)
    alive_dead_ratio = np.sum(yp1) / len(yp1)
    
    if len(yp1) >359:
        xp1 = np.zeros(0)
        yp1 = np.zeros(0)
    
    x1, y1 = CreateLineFromAngleDegrees(theta1)
    x2, y2 = CreateLineFromAngleDegrees(theta2)
    
    ax1.collections.clear()
    ax1.fill_between(x1, y1, y2, facecolor='yellow', alpha=0.5)
    ax1.scatter(*ObserverPosition(inclination_angle), c='b')    #Plots observer position
    ax1.scatter(0,0, marker='*', c='b', s=5)
    
    line1a.set_data(x1, y1)
    line1b.set_data(x2, y2)
    text1a.set_text('theta1: {}'.format(str(round(theta1,2))))
    text1b.set_text('theta2: {}'.format(str(round(theta2,2))))
    text1c.set_text('Opening angle: {}'.format(str(round(opening_angle,2))))
    text1d.set_text('b: {}'.format(str(round(b,2))))
    text1e.set_text('Precession angle: {}'.format(swing_angle))
    
    line2.set_data(xp1,yp1)
    line2b.set_data(xp1, np.ones(len(xp1))*alive_dead_ratio)
    
    ax2.set_xlim(0, 360) #Seting our xlimits
    ax2.set_ylim(0, 1.5) #Setting our ylimits
    
    text2a.set_text('Alive/Dead: {}'.format(round(alive_dead_ratio,2)))
    artist_objects = (line1a, line1b, line2)
    return artist_objects


#Figure for plotting onto


def model(inclination_angle, opening_angle, swing_angle, period):
    '''
    This function contains the same code as the ones found in the 'animate'
    function with the plotting removed, it is used to investigate properties
    of multiple spinning cones.
    '''
    xp1 = np.zeros(0)
    yp1 = np.zeros(0)
    for t in range(360):
        half_opening_angle = opening_angle / 2
        b = 1 - cos(deg2rad(half_opening_angle))

        theta1 = (swing_angle * sin(2*pi/period*t))%360
        theta2 = (theta1+opening_angle)%360
        
        vis = isVisible(inclination_angle, theta1, theta2)
        
        xp1 = np.append(xp1, t)
        yp1 = np.append(yp1, vis)
        alive_dead_ratio = np.sum(yp1) / len(yp1)
        
        if len(yp1) >359:
            xp1 = np.zeros(0)
            yp1 = np.zeros(0)
        
        x1, y1 = CreateLineFromAngleDegrees(theta1)
        x2, y2 = CreateLineFromAngleDegrees(theta2)
    return alive_dead_ratio

 
def periodTest():
    '''
    Try randomising start times and try again
    '''
    inclination_angle = 10
    opening_angle = 20
    swing_angle = 20

    periods = np.arange(1,360,1)
    alive_deads = np.empty(len(periods))
    for i, period in tqdm(enumerate(periods)):
        alive_deads[i] = model(inclination_angle, opening_angle, swing_angle, period)
    
    return periods, alive_deads


def openingAngleTest():
    inclination_angle = 10
    swing_angle = 20
    period = 120

    opening_angles = np.arange(1,90,1)
    alive_deads = np.empty(len(opening_angles))
    for i, opening_angle in tqdm(enumerate(opening_angles)):
        alive_deads[i] = model(inclination_angle, opening_angle, swing_angle, period)
    
    return opening_angles, alive_deads


def swingAngleTest():
    inclination_angle = 10
    opening_angle = 20
    period = 120

    
    swing_angles = np.arange(0,45, 0.2)
    alive_deads = np.empty(len(swing_angles))
    for i, swing_angle in tqdm(enumerate(swing_angles)):
        alive_deads[i] = model(inclination_angle, opening_angle, swing_angle, period)
    
    return swing_angles, alive_deads


# opening_angles, alive_deads1 = openingAngleTest()
# periods, alive_deads2 = periodTest()
# swing_angles, alive_deads3 = swingAngleTest()

# plt.title('inclination_angle = 10, swing_angle = 20, period = 120')
# plt.xlabel('opening_angle')
# plt.ylabel('Alive/Dead Ratio')
# plt.plot(opening_angles, alive_deads1)


# plt.title('inclination_angle = 10, swing_angle = 20, opening_angle = 20')
# plt.xlabel('period')
# plt.ylabel('Alive/Dead Ratio')
# plt.plot(periods, alive_deads2)


# plt.title('inclination_angle = 10, opening_angle = 20, period = 120')
# plt.xlabel('swing_angle')
# plt.ylabel('Alive/Dead Ratio')
# plt.plot(swing_angles, alive_deads3)





def RandomInclination():
    r = np.random.randint(0,high=2)
    if r == 1:
        random = np.random.randint(0, high=90)
    else:
        random = np.random.randint(270, high=360)
    return random

def GetOpeningAngleFromDistribution():
    swing_angle = np.random.randint(1, high=20)
    return swing_angle

def GetSwingAngleFromDistribution():
    opening_angle = np.random.randint(5, high=20)
    return opening_angle

# =============================================================================
# ulx observation test
# =============================================================================
'''
NUMBER_OF_ULXS = 1000
period = 120
alive_deads = np.empty(NUMBER_OF_ULXS)
for N in tqdm(range(NUMBER_OF_ULXS)):
    inclination = RandomInclination()
    opening_angle = GetOpeningAngleFromDistribution()
    swing_angle = GetSwingAngleFromDistribution()
    
    alive_deads[N] = model(inclination, opening_angle, swing_angle, period)

plt.figure()
plt.hist(alive_deads)
plt.title(len(alive_deads))

alive_deads_nonzero = alive_deads[alive_deads != 0]
    
plt.figure()
plt.hist(alive_deads_nonzero)
plt.title(len(alive_deads_nonzero))

'''
#distribution of thetas
#distribution of swing angles
#how many of the 1000 can we see?
#of the ones we can see, how often do we see them?











fig = plt.figure(figsize=(10,5))


#ax1 will hold our spinning cone
ax1 = fig.add_subplot(1,2,1)
ax1.set_xlim(-1, 1) #Seting our xlimits
ax1.set_ylim(-1, 1) #Setting our ylimits
line1a, = ax1.plot([], [], lw=2, color='r') #Line 1a for one of the cone sides
line1b, = ax1.plot([], [], lw=2, color='b') #Line 1b for the other one
text1a = ax1.text(-0.9, 0.9, '', fontsize=10, color='r') #Text for displaying theta1
text1b = ax1.text(-0.9, 0.8, '', fontsize=10, color='b') #Text for displaying theta2
text1c = ax1.text(-0.9, 0.7, '', fontsize=10, color='b') #Text for displaying opening angle
text1d = ax1.text(-0.9, 0.6, '', fontsize=10, color='b') #Text for displaying b
text1e = ax1.text(-0.9, 0.5, '', fontsize=10, color='b') #Text for displaying swing angle


#ax2 will show a lightcurve created from our simulation.
ax2 = fig.add_subplot(1,2,2)
ax2.set_xlim(0, 360) #Seting our xlimits
ax2.set_ylim(0, 1.5) #Setting our ylimits
line2, = ax2.plot([], [], lw=2, color='b') #Line on second axis for our lightcurve
text2a = ax2.text(100, 1.2, '', fontsize=15, color='r') #Text for displaying alive/dead
line2b, = ax2.plot([], [], lw=2, color='r') #Line on second axis for alive/dead


#Data placeholders
xp1 = np.zeros(0)
yp1 = np.zeros(0)

anim = animation.FuncAnimation(fig, animate, init_func=init,
                               frames=360, interval=20, blit=False)

# anim.save('swing40open20.mp4', writer=writer)

'''
Next thing we want to do is:
    Simulate many precessing cones for many inclination angles and figure
    out the alive/dead time as a distribution of:
        -opening angle
        -Precession angle
        -period
        
    -We need a way of obtaining the alive/dead time for a system.
    This may be as simple as summing the array and dividing by the length of it.
    
   ''' 