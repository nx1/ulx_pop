project_name
==============================

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



# Run values
 	

~~~~
#TOTAL SUMS
#0
dincl = 5, ratio =  0.8317294459702234
dincl = 10, ratio =  0.8406542739290764
dincl = 15, ratio =  0.8360675776747886
dincl = 20, ratio =  0.8313453717813738
dincl = 25, ratio =  0.8256529008486989
dincl = 30, ratio =  0.8371486045261644
dincl = 35, ratio =  0.8382256430651508
dincl = 40, ratio =  0.83197472580942
dincl = 45, ratio =  0.8351058772117088

#1
dincl = 5, ratio =  0.8353364487126705
dincl = 10, ratio =  0.8396567212163679
dincl = 15, ratio =  0.8373005473723573
dincl = 20, ratio =  0.8338930974371417
dincl = 25, ratio =  0.8381568251703276
dincl = 30, ratio =  0.8344133578241867
dincl = 35, ratio =  0.834895887096725
dincl = 40, ratio =  0.8413651592607563
dincl = 45, ratio =  0.8381873248328369

#2
dincl = 5, ratio =  0.8370096090071556
dincl = 10, ratio =  0.8332802641896158
dincl = 15, ratio =  0.8384537557000523
dincl = 20, ratio =  0.8400171002175412
dincl = 25, ratio =  0.8390795036174923
dincl = 30, ratio =  0.8327522441916837
dincl = 35, ratio =  0.8341380592172919
dincl = 40, ratio =  0.8345432384202286
dincl = 45, ratio =  0.8358387671146716

#3
dincl = 5, ratio =  0.8331333228792229
dincl = 10, ratio =  0.8363463785142771
dincl = 15, ratio =  0.8345205901657279
dincl = 20, ratio =  0.8252856675769871
dincl = 25, ratio =  0.8467548500983522
dincl = 30, ratio =  0.842981615392678
dincl = 35, ratio =  0.827284771752153
dincl = 40, ratio =  0.8334447235240953
dincl = 45, ratio =  0.8340206975187066

#4
dincl = 5, ratio =  0.8390115730456478
dincl = 10, ratio =  0.8427490007504218
dincl = 15, ratio =  0.8368012337997748
dincl = 20, ratio =  0.8362223315972949
dincl = 25, ratio =  0.8331391317882217
dincl = 30, ratio =  0.8321833667824919
dincl = 35, ratio =  0.8283764492132099
dincl = 40, ratio =  0.8382321855228965
dincl = 45, ratio =  0.8355874206738593

#5
dincl = 5, ratio =  0.8218154822208887
dincl = 10, ratio =  0.8344779873673505
dincl = 15, ratio =  0.8390542000284258
dincl = 20, ratio =  0.8344839162751383
dincl = 25, ratio =  0.8352030352147007
dincl = 30, ratio =  0.8317972404742697
dincl = 35, ratio =  0.836655203549911
dincl = 40, ratio =  0.8373553123841916
dincl = 45, ratio =  0.8350817293459148

#6
dincl = 5, ratio =  0.83162743736234
dincl = 10, ratio =  0.8348237154576692
dincl = 15, ratio =  0.834634995411111
dincl = 20, ratio =  0.8366250096869557
dincl = 25, ratio =  0.8374347453029429
dincl = 30, ratio =  0.8330281975912883
dincl = 35, ratio =  0.8357450230346363
dincl = 40, ratio =  0.8395102169086934
dincl = 45, ratio =  0.83181355150354

#7
dincl = 5, ratio =  0.830963763579405
dincl = 10, ratio =  0.8417705360424537
dincl = 15, ratio =  0.8432165585401263
dincl = 20, ratio =  0.8330802819698838
dincl = 25, ratio =  0.835515052286546
dincl = 30, ratio =  0.8409213303619443
dincl = 35, ratio =  0.8319285325028283
dincl = 40, ratio =  0.8304982701846341
dincl = 45, ratio =  0.8373765649366446

#8
dincl = 5, ratio =  0.8407494146254205
dincl = 10, ratio =  0.8262362290130204
dincl = 15, ratio =  0.8274289396770118
dincl = 20, ratio =  0.8444285071840321
dincl = 25, ratio =  0.8458616296524615
dincl = 30, ratio =  0.8347044737185683
dincl = 35, ratio =  0.8344974583348533
dincl = 40, ratio =  0.8383964607851576
dincl = 45, ratio =  0.8380422364502151


#mean for alive/dead on the 240 beamed sources:

#0
0.30448171001025687
0.34137099890684913
0.3224126543891258
0.30289420336301176
0.2793653235079553
0.326880898708146
0.33133265800262335
0.3054955333456026
0.31843762580839596

#1
0.3193906546790381
0.337247781027654
0.32750892913907664
0.31342480274018547
0.33104821070402135
0.3155752123399718
0.31756966666646363
0.3443093249444591
0.3311742759757256

#2
0.3263063838962433
0.3108917586504119
0.3322755235602163
0.338737347565837
0.33486194828563476
0.3087092759922925
0.3144373114314734
0.31611205213694477
0.32146690407397593

#3
0.3102844012341212
0.3235650311923451
0.31601843935167545
0.2778474259848796
0.3665867137398556
0.3509906769564026
0.28611038990889887
0.31157152389959364
0.313952216410654

#4
0.33458116858867726
0.3500292031017432
0.3254450997057356
0.32305230393548534
0.31030841139131626
0.30635791603430024
0.2906226567479342
0.3313597001613054
0.32042800545195144

#5
0.2635039931796734
0.31584234778504855
0.3347573601174933
0.3158668539372383
0.3188392122207633
0.3047619272936485
0.3248415080062989
0.32773529118799216
0.31833781462978133

#6
0.30406007443100547
0.31727135722503247
0.3164913143659252
0.324716706706084
0.3280636139188306
0.309849883377325
0.3210794285431631
0.3366422298892659
0.3048293462146317
#7
0.301316889461541
0.3459848823088085
0.3519617752991889
0.31006516547551943
0.32012888278439033
0.34247483216270297
0.30530460101169016
0.29939285009648775
0.3278231350714643

#8
0.3417642471184051
0.2817764132538178
0.2867062839983152
0.35697116302733256
0.36289473589684124
0.3167784913700824
0.3159228277840605
0.33203870457865126
0.33057457732755624
~~~~


