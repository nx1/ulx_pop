# -*- coding: utf-8 -*-
"""
Created on Tue Jul 14 13:12:54 2020

@author: norma
"""

import numpy as np

rad2deg = 180.0 / np.pi
deg2deg = np.pi / 180.0


class Parameters:
    def __init__(self, period, phase, theta, incl, dincl, beta, dopulse):
        self.period = period            # time of period
        self.phase = phase              # going from 0 to 1
        self.theta = theta * deg2rad    # half of the opening angle
        self.incl = incl * deg2rad      # inclination to the system at phase zero
        self.dincl = dincl * deg2rad    # precession in angle
        self.beta = beta                # beta=v/c
        self.dopulse = dopulse          # pulse(=1) or lightcurve(=0)?
        self.check_per_model()

    def check_per_model(self):
        if self.beta<0:
            print(f" *** warning, value of beta={self.beta} < 0 not allowed. Resetting to beta=0.0");
            self.beta = 0.0
    
        if self.beta>0.999:
            print(f" *** warning, value of beta={self.beta} > 0.999 not allowed. Resettisng to beta=0.999")
            self.beta = 0.999
    
        if self.theta<0:
            print(f" *** warning, value of theta={self.theta} < 2 deg not allowed. Resetting to theta=0.0")
            self.theta = 2.0/rad2deg
    
        if self.theta>PI/4.0:
            print(f" *** warning, value of theta={self.theta} > 45.0 deg not allowed. Resetting to theta=45.0 deg")
            self.theta = PI/4.0
    
        if self.incl<0:
            print(" *** warning, value of incl={self.incl} < 0 deg not allowed. Resetting to incl=0.0 deg")
            self.incl = 0.0
    
        if self.incl>PI/2.0:
            print(f" *** warning, value of incl={self.incl} > 90 deg not allowed. Resetting to incl=90 deg")
            self.incl = PI/2.0
    
        if (self.dopulse>1) or (self.dopulse<0):
            print(f" *** warning, value of dopulse = {self.dopulse}, but can only be 0 (light curve) or 1 (pulse)")
            self.dopulse = 0

Parameters(10, 0, 10, 10, 10, 0.3, 1)

def binary_search(val, arr, n):
    """
    Binary search for to find interpolation interval
    *   - return value is the bin [ind,ind+1]
    *   - assume list is sorted ascending

    Parameters
    ----------
    val : float
    arr : float
    n : int
    """
    if (val < arr[0]) or (val > arr[n-1]):
        return -1;
    
    high = n-1
    low = 0
    
    while (high > low):
        mid = int((low+high) / 2)
        if (arr[mid] <= val):
            low = mid + 1
        else:
            high = mid;
    return low-1





class ulx_table:
    def __init__(self):
        self.ang = None 
        self.flux = None
        self.nang = 0
        self.ntheta = 0
        self.theta = NULL
        self.beta = NULL
        self.nbeta = 0




tab = ulx_table()

class ulxlc_error(Exception):
    pass

# calculate the flux array for the current value of theta
def interp_betTheta(tab, theta, beta, status):
    assert(tab!=NULL);

    ind_th = binary_search(theta, tab.theta, tab.ntheta)
    
    if ind_th == -1:
        print(f" *** error: values of theta={theta*rad2deg} not tabulated")
        # " failed to look up the value of theta in the loaded table"
        raise(ulxlc_error)
        return None
    
    ind_bet = binary_search(beta, tab.beta, tab.nbeta)
    
    if ind_bet == -1:
        print(f" *** error: values of beta={beta} not tabulated")
		# " failed to look up the value of theta in the loaded table"
        raise(ulxlc_error)
        return None

    # allocate memory for the flux array
    inp_flux = []

    # 2d-interpolate the flux now
    fac_th  = (theta - tab.theta[ind_th]) / (tab.theta[ind_th + 1] - tab.theta[ind_th])
    fac_bet = (beta  - tab.beta[ind_bet]) / (tab.beta[ind_bet + 1] - tab.beta[ind_bet])

    for ii in range(tab.nang):
        # \ means continue to next line
        inp_flux[ii] = (1-fac_bet)*(1-fac_th) * tab.flux[ind_bet][ind_th][ii] \
                        + (fac_bet)*(1-fac_th)   * tab.flux[ind_bet+1][ind_th][ii] \
                        + (1-fac_bet)*(fac_th)   * tab.flux[ind_bet][ind_th+1][ii] \
                        + (fac_bet)*(fac_th)     * tab.flux[ind_bet+1][ind_th+1][ii] 
    return inp_flux

def calc_orbit_ang(double phase, param):
    # need the sine distributed from [0,1]
    ang = 0.5*(np.sin(phase*2*np.pi+np.pi/2)+1)*(param.dincl*2) + param.incl - param.dincl
	return abs(ang) # return the absolute value

def get_flux_from_angle(ang, flux, ang_arr, nang, status):
    ang_min = 0.0
    ang_max = np.pi/2

    # allow angles between 90 to 180 degrees
	if (ang > np.pi/2) and (ang < np.pi):
		ang = np.pi - ang

	# check at the boundaries
	if (ang < ang_arr[0]) and (ang>=ang_min):
		return flux[0]
	else if (ang > ang_arr[nang-1]) and (ang <= ang_max):
		return flux[nang-1]
    
    # get the corresponding index
	ind = binary_search(ang, ang_arr, nang)

    # check if something went wrong (should not happen)
    if (ind == -1):
		print(f" *** error: angle {ang*180/np.pi} [deg] is not tabulated (should be in [{ang_arr[0]*180/PI},{ang_arr[nang-1]*180/PI}]")
		raise(ulxlc_eror)
		return -1.0

	# interpolate now the flux
    fac = (ang-ang_arr[ind]) / (ang_arr[ind+1]-ang_arr[ind])
    return (1-fac)*flux[ind]+fac*flux[ind+1]

# input:  t[nt+1], photar[nt]
def create_lc(t0, nt, photar, flux, tab, param, status):
	for ii in range(nt): # nt as photar[nt]
		# if we want only a single pulse
		phase=0.0

		if param.dopulse == 1:
			# take the phase in the middle again
			delta_t_phase_lo = (t0[ii] - t0[0]) / param.period - param.phase
			delta_t_phase_hi = (t0[ii+1] - t0[0]) / param.period - param.phase
			phase = 0.5*(delta_t_phase_lo+delta_t_phase_hi)
			if ( phase < 0) or (phase >= 1):
                continue
            else:
                # see which part in the period we are in (then time is between 0 and param.period)
                phase = (0.5*(t0[ii]+t0[ii+1])/param.period + param.phase)
                phase = phase - (double) ((int) phase) # just get everything behind the colon

		# now get the corresponding angle
		double ang = calc_orbit_ang(phase, param)
		photar[ii] = get_flux_from_angle(ang, flux, tab.ang, tab.nang, status)


# basic calculation of the model
# input:  t[nt+1], photar[nt]
# output: photar
def ulxlc_base(t, nt, photar, param_struct, status):
    # load the GLOBAL TABLE
	if (global_ulx_table==None):
		ulxlc_load_table(ULXTABLE_FILENAME, &global_ulx_table, status)
		CHECK_ULXLC_ERROR("reading the table failed", status);
		CHECK_STATUS_VOID(*status)
        
	assert(global_ulx_table!=None)

    # interpolate the flux for a given opening_angle/theta
	flux = interp_betTheta(global_ulx_table, param_struct.theta, param_struct.beta, status)
	CHECK_ULXLC_ERROR("reading of theta values from table failed",status)
	CHECK_STATUS_VOID(*status)

	create_lc(t, nt, photar, flux, global_ulx_table, param_struct, status)
	free(flux)
	CHECK_STATUS_VOID(*status)

    # free_ulxTable(&global_ulx_table)

def ulxlc_model(const t, nt, photar, parameter, n_parameter):
	status = EXIT_SUCCESS

	# fill in parameters
	inp_param param_struct;

	init_par_model(&param_struct, parameter, n_parameter, &status);
	CHECK_STATUS_RET(status, status)

	# call the function which calculates the light curve
	ulxlc_base(t, nt, photar, param_struct,&status);
	CHECK_STATUS_RET(status,status);

	return status;
}


void ulxlcmodxspec(const Real* energy, int Nflux, const Real* parameter,
		int spectrum, Real* flux, Real* fluxError, const char* init){

	int status = ulxlc_model(energy,Nflux,flux,parameter,NUM_PARAM);
	CHECK_ULXLC_ERROR("evaluation of the ulxlc model failed",&status);

	// XSPEC / ISIS want bin integrated values
	int ii;
	for (ii=0; ii < Nflux; ii++){
		flux[ii] = flux[ii] * (energy[ii+1] - energy[ii]);
	}


}
