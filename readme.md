# Code for ULX population estimate.

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

