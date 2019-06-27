"""
Created on Mon Jun 24 15:52:49 2019

@author: nk7g14
"""
import numpy as np
from numpy import sin, cos, tan
import matplotlib.pyplot as plt
from matplotlib import animation

from matplotlib.lines import Line2D
import matplotlib.animation as animation

deg2rad = lambda deg : (2*np.pi/360) * deg


def CaclulateNewTheta(theta, t, spin_period):
    theta = (2 * np.pi / (spin_period) * t)%(2*np.pi)
    return theta

def ObserverPosition(inclination_angle):
    inclination_angle_radians = deg2rad(inclination_angle)
    r = 0.5
    x = r * cos(inclination_angle_radians)
    y = r * sin(inclination_angle_radians)
    return x, y
    

def isVisible(eye_position, theta1, theta2):
    #The cone is overlapping the 0 axis
    a = round(theta1, 2)
    b = round(theta2, 2)
    
    print('Checking is visible for eye: {} t1: {} t2: {}'.format(eye_position, a, b))
    if theta2 < theta1:
        #Check if eye is between 0 axis and theta 2
        if eye_position <= theta2:
            return 1
        #Check if eye is between theta 1 and 0 axis
        if eye_position >= theta1:
            return 1
    else:    
        if eye_position > theta1 and eye_position < theta2:
            return 1
        else:
            return 0
    return 0


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
    elif theta_rad > (3/4)*np.pi and theta_rad < 2*np.pi:
        # print('q4')
        theta_q4 = 2*np.pi - theta_rad #Angle between x axis and line in q4
        m = -tan(theta_q4)
        x = np.arange(-1,1,0.01)
        y = m*x
    else:
        print('uwot????????????????????')
    return x, y


def init():
    line1.set_data([], [])
    line2.set_data([], [])
    return line1, line2

def animate(time):
    swing_angle = 30
    theta_min = -swing_angle
    theta_max = swing_angle
    
    opening_angle_degrees = 30
    
    theta = theta_max/2 * (np.sin(time * np.pi / 180) + theta_min + 1)
    x, y = CreateLineFromAngleDegrees(theta)
    x1, y1 = CreateLineFromAngleDegrees( (theta+opening_angle_degrees)%360 )
     
    line1.set_data(x, y)
    line2.set_data(x1, y1)
    
    return [line1, line2]




# eye_position = 15



# fig = plt.figure(figsize=(5,5))
# ax = plt.axes(xlim=(-1, 1), ylim=(-1, 1))
# PlotEye(eye_position)
# line1, = ax.plot([], [], lw=2, color='g')
# line2, = ax.plot([], [], lw=2, color='g')
# anim = animation.FuncAnimation(fig, animate, init_func=init,
#                                frames=360, interval=20, blit=True)





'''

time_range = 3
spin_period = 1    
opening_angle_degrees = 30
eye_position = 15    
theta1 = 0

x=[]    #Time
vis_arr = []    #Visible or not

time_range = np.arange(0,time_range,0.01)
for t in time_range:
    theta1 = CaclulateNewTheta(theta1, t, spin_period)
    theta2 = (theta1 + deg2rad(opening_angle_degrees))%(2*np.pi)
    vis_arr.append(isVisible(eye_position, theta1, theta2))

    x.append(t)
    
plt.plot(x,vis_arr)










test_opening_angles = np.arange(0,90,2)

for opening_angle_degrees in test_opening_angles:
    
    theta1 = 0
    theta2 = theta1 + opening_angle_degrees
    
    x=[]    #Time
    y0=[]   #Theta1
    y1=[]   #Theta2
    vis_arr = []    #Visible or not
    
    t_arr = np.arange(0,time_range,0.01)
    
    for t in t_arr:
        theta1 = CaclulateNewTheta(theta1, t)
        theta2 = (theta1 + deg2rad(opening_angle_degrees))%(2*np.pi)
        vis_arr.append(isVisible(eye_position, theta1, theta2))
        
        x.append(t)
        y0.append(theta1)
        y1.append(theta2)
        
    plt.plot(x,vis_arr)
    # plt.plot(x,y0)
    # plt.plot(x,y1)
    # print('opening_angle_degrees:',opening_angle_degrees)
    
    
    # plt.scatter(opening_angle_degrees, np.mean(vis_arr))

plt.show()
'''