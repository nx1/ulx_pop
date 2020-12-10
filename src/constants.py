# -*- coding: utf-8 -*-
"""
Created on Tue Jul  7 14:31:12 2020

@author: norma
"""

#CGS UNITS
c    = 3E10          # Speed of Light in cm/s
M_sol = 1.989E33     # Mass of sun in g
R_sol = 6.957E10     # cm
L_sol = 3.828E33     # erg s^-1

sigma = 5.6704E-5    # erg cm^-2 s^-1 K^-4 Stefan-Boltzmann constant

yr   = 31557600      # 1 yr in Seconds
Myr  = yr * 1E6      # 1 Myr in Seconds

G    = 6.674E-8      # (cm)^3 g^-1 s^-2

eta = 1/12 # Accretion efficiency

#GR constants
epsilon = 0.25 # Normally between 0.25 to 0.95 (check latex)
beta    = 1.4  # Velocity of the wind, distinct from the beta used in ulxlc
BH_ISCO = 1.25
BH_SPIN = 0.998
NS_ISCO = 6
NS_SPIN = 0.01

#ULXLC Constants
period_default = 50.0
phase_default = 0.0
theta_default = 10.0
incl_default = 5.0
dincl_default = 20.0
beta_default = 0.3
dopulse_default = 0.0

params_default = [period_default,
                  phase_default,
                  theta_default,
                  incl_default,
                  dincl_default,
                  beta_default,
                  dopulse_default]

# StarTrack
stellar_dict = {0 : 'MS M < 0.7',
                1 : 'MS M > 0.7',
                2 : 'Hertzsprung gap star',
                3 : 'First giant branch star',
                4 : 'Core helium burning star',
                5 : 'Early asymptotic giant branch star',
                6 : 'Thermally pulsing asymptotic giant branch star',
                7 : 'Main sequence naked helium star (Wolf-Rayet star)',
                8 : 'Hertzsprung gap naked helium star',
                9 : 'Giant branch naked helium star',
                10 : 'Helium white dwarf',
                11 : 'Carbon/oxygen white dwarf',
                12 : 'Oxygen/neon white dwarf',
                13 : 'Neutron star',
                14 : 'Black hole',
                15 : 'Massless Remnant',
                16 : 'Hydrogen White Dwarf',
                17 : 'Hybrid white dwarf'}