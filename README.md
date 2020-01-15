ULX population study
====================

Population study on precessing ulxs

Project Organization
------------

    ├── LICENSE
    ├── Makefile           <- Makefile with commands like `make data` or `make train`
    ├── README.md          <- The top-level README for developers using this project.
    ├── data
    │   ├── external       <- Data from third party sources.
    │   ├── interim        <- Intermediate data that has been transformed.
    │   ├── processed      <- The final, canonical data sets for modeling.
    │   └── raw            <- The original, immutable data dump.
    │
    ├── docs               <- A default Sphinx project; see sphinx-doc.org for details
    │
    ├── models             <- Trained and serialized models, model predictions, or model summaries
    │
    ├── notebooks          <- Jupyter notebooks. Naming convention is a number (for ordering),
    │                         the creator's initials, and a short `-` delimited description, e.g.
    │                         `1.0-jqp-initial-data-exploration`.
    │
    ├── references         <- Data dictionaries, manuals, and all other explanatory materials.
    │
    ├── reports            <- Generated analysis as HTML, PDF, LaTeX, etc.
    │   └── figures        <- Generated graphics and figures to be used in reporting
    │
    ├── requirements.txt   <- The requirements file for reproducing the analysis environment, e.g.
    │                         generated with `pip freeze > requirements.txt`
    │
    ├── setup.py           <- makes project pip installable (pip install -e .) so src can be imported
    ├── src                <- Source code for use in this project.
    │   ├── __init__.py    <- Makes src a Python module
    │   │
    │   ├── data           <- Scripts to download or generate data
    │   │   └── make_dataset.py
    │   │
    │   ├── features       <- Scripts to turn raw data into features for modeling
    │   │   └── build_features.py
    │   │
    │   ├── models         <- Scripts to train models and then use trained models to make
    │   │   │                 predictions
    │   │   ├── predict_model.py
    │   │   └── train_model.py
    │   │
    │   └── visualization  <- Scripts to create exploratory and results oriented visualizations
    │       └── visualize.py
    │
    └── tox.ini            <- tox file with settings for running tox; see tox.testrun.org


--------

# ULX population estimate.

This github repo contains all the files required to perform a ulx population estimate using the xspec model ulxlc: http://www.sternwarte.uni-erlangen.de/~dauser/research/ulxlc/

The code currently relies on 3 files:
    data_mdot/test.py
    cuves/plotter.py
    ulxlc_code_v0.1/Batch.py
    
# Building the xspec model

The method by which the xspec command ulxlc is called requires first that the code be built using the xspec command:

    initpackage [model_name] [path to model .dat file] [path to model install location]
eg:

    initpackage ulxlc /home/nk7g14/Desktop/gitbox/ulx_pop/ulxlc_code_v0.1/lmodel_ulxlc.dat /home/nk7g14/Desktop/gitbox/ulx_pop/ulxlc_code_v0.1

# Program structure

The code currently consists of 3 files

## test.py

Currently this code takes in the raw files, combines them into one large data
frame, and calculates the required values for use in the simulation.

~~~~
FILE NAMES:
    eg: 'Z_0.002_tage_20.dat'
    'Z_0.002' = Metallicity = 0.002
    Metallicity range: 0.0002 - 0.02 
    'tage_20' = Age of system (time since ZAMS in Myr).
    Age range: 10 - 10000 
    
COLUMN NAMES:
mdot:
    mass transfer rate (in solar masses per Myr; to the disk, not to the accretor, 
    which is lower due to the mass loss in disk wind),
m:
    mass of a compact object (in solar masses; we assumed that these with a mass 
    lower than 2.5 Msun are NSs, whereas the heavier ones are BHs).
idum, iidd:
    identificators of the system in the simulations. They have no physical meaning,
    but will allow to find more info about some potentially odd, or interesting system.
Z:
    Metallicity | range: 0.0002 - 0.02 
tage:
    Age of system (time since ZAMS in Myr). | range: 10 - 10000
    
Mass Accretion rates:
    mdot_gs:
        Mass Accretion rate in grams per second
    mdot_ratio:
        Eddington Ratio
    MEdd:
        Mass Accretion rate for eddington.
        defined by the equation LEdd / n c^2 where n = 1/12 for BH and 0.2 for NS
        (See Grzegorz Wiktorowicz 2017 Appendix)
        
Luminosities:
    LEdd:
        Eddington Luminosity in erg/s from (1.2E38 * m)
    XLsph:
        A.R King method with log contribution 
    XLsph2:
        A.R King Method without log contribution
    LXtot:
        Grzegorz Wiktorowicz Method
    Lx:
        Beamed Luminosity from Lxtot/b
        In the 2018 paper by Greg they propose the following:
            Lx = Ledd * (1 + ln(mdot))  for mdot > 1
            Lx = Ledd * mdot            for mdot < 1
    ratio:
        Ratio of LXtot/Ledd
    ratio_beamed:
        Ratio of Lx/Ledd
        
Beaming:
    b:
        Beaming parameter set by the following:
            b = 1 for m_dot < 8.5
            b = 73/mdot**2 >= 8.5
            b = 3.2E-3 for >150
        
        Examples of b to theta (Full opening angle):
            b = 0.0032  | theta = 9.2 deg
            b = 0.01    | theta = 16.8 deg
            b = 0.022   | theta = 24 deg
            b = 0.077   | theta = 45 deg
            b = 0.1     | theta = 53 deg
            b = 0.2     | theta = 73.9 deg
            b = 0.3     | theta = 91 deg
            b = 1.0     | theta = 180 deg (No beaming)
            
    theta:
        cone opening angle in rad
        Calculated from b = 1 - cos(theta/2)
    theta_deg:
        cone opening angle in rad
    theta_half_deg:
        cone half opening angle in rad
ULXLC:
   norm - Normalisation
   period - Period of the light curve
   phase - Phase since the start of the lightcurve
   theta - Half opening angle
   incl - Inclination of system (90 is edge on)
   dincl - Presccesion angle
   beta - velocity of the outflow; given in units of the speed of light
   dopulse - 0 or 1 for full lightcurve or 1 pulse
   

1) Create light curve for given system for variety of dincl
   norm = 1 
   period = 10
   phase = 0
   theta = from df
   incl = random(0,90)
   dincl = [0, 5 , 10, 30, 45]
   beta = 0.5
   dopulse = 0 
   
2) Save light curve
3) Filename  = mdot,m,Z,tage,dincl
        
how to export the model from XSPEC:
    after plotting the model
    >ipl
    >wdata [filename]
    
    
1) use python to run xspec script
2) once you get to PLT> bit use python to run the command wdata
3) use python to repeat this for different parameters
4) Once you have all files perform the analysis in python


~~~~

## Batch.py

This code takes in the dataframe.csv file from test.py and then runs through each
of the required systems creating light curves for both inclination = 0 and inclination = random
using the xspec model [ulxlc](http://www.sternwarte.uni-erlangen.de/~dauser/research/ulxlc/ "ulxlc").


~~~~
Import and filter dataframe csv
Run through all systems
    for each system run through dincls
        MakeXCM()
        MakeInfo()
        Call Xspec
~~~~    

## plotter.py

This code takes all the output lightcurve files from batch.py in txt form and the associated
.info file which contains information about the system, then performs alive/dead time
analysis on each of the sytems.


~~~~
This script runs through the output curves generated from the xspec model 
ulxlc and performs alive/dead time analysis for each of them.

Lightcurves are stored in .txt files
Each lightcurve has an associated .info file stored in /info/

The filename format goes as:
    simulation_number - dincl - inclination.txt

where dincl is the precession angle

1) Load in all lightcurves into df_dict
2) Load in all info files into info_dict

For each simulation_number dincl combo there is a 0 inclination system and
a random inclination system.

3) Find the 0 inclination system and set the scaling factor, c, and the limiting
value, N_lim.

c = Lx / max(curves[key]['Flux']) #Scaling factor
N_lim = 1E39 / c                  #Limiting value
~~~~


~~~~


