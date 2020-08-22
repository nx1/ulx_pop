# -*- coding: utf-8 -*-
"""
Created on Tue Jul  7 14:31:12 2020

@author: norma
"""

#CGS UNITS
c    = 3E10          # Speed of Light in cm/s
M_sol = 1.989E33     # Mass of sun in g
R_sol = 6.957E10     #cm

yr   = 31557600     # 1 yr in Seconds
Myr  = yr * 1E6     # 1 Myr in Seconds

G_SI = 6.67408E-11  # N.(m^2)/(kg)^2 
G    = 6.674E-8     # (cm)^3 g^-1 s^-2


#GR constants
epsilon = 0.25 # Normally between 0.25 to 0.95 (check latex)
beta    = 1.4  # Velocity of the wind, distinct from the beta used in ulxlc
BH_ISCO = 1.25
BH_SPIN = 0.998
NS_ISCO = 6
NS_SPIN = 0.01