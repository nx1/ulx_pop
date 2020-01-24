/*
   This file is part of the ULXLC model code.

   The ULXLC model is free software: you can redistribute it and/or modify it
   under the terms of the GNU General Public License as published by
   the Free Software Foundation, either version 3 of the License, or
   any later version.

   ULXLC is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
   GNU General Public License for more details.
   For a copy of the GNU General Public License see
   <http://www.gnu.org/licenses/>.

    Copyright 2016 Thomas Dauser, Remeis Observatory & ECAP
*/

#include "ulxlcmod.h"

/** global parameters, which can be used for several calls of the model */
ulx_table* global_ulx_table=NULL;

void check_par_model(inp_param* param){

	// convert from degree to radian
	double rad2deg = 180.0/PI;

	if (param->beta<0){
		printf(" *** warning, value of beta=%e < 0 not allowed. Resetting to beta=0.0 \n",param->beta);
		param->beta=0.0;
	}
	if (param->beta>0.999){
		printf(" *** warning, value of beta=%e > 0.999 not allowed. Resetting to beta=0.999 \n",param->beta);
		param->beta=0.999;
	}

	if (param->theta<0){
		printf(" *** warning, value of theta=%e < 2 deg not allowed. Resetting to theta=0.0 \n",param->theta*rad2deg);
		param->theta=2.0/rad2deg;
	}
	if (param->theta>PI/4.0){
		printf(" *** warning, value of theta=%e > 45.0 deg not allowed. Resetting to theta=45.0 deg \n",param->theta*rad2deg);
		param->theta=PI/4.0;
	}

	if (param->incl<0){
		printf(" *** warning, value of incl=%e < 0 deg not allowed. Resetting to incl=0.0 deg \n",param->incl*rad2deg);
		param->incl=0.0;
	}
	if (param->incl>PI/2.0){
		printf(" *** warning, value of incl=%e > 90 deg not allowed. Resetting to incl=90 deg \n",param->incl*rad2deg);
		param->incl=PI/2.0;
	}
	if (param->dopulse>1 || param->dopulse<0){
		printf(" *** warning, value of dopulse = %i, but can only be 0 (light curve) or 1 (pulse) \n",param->dopulse);
		param->dopulse=0;
	}

	return;
}

void init_par_model(inp_param* param, const double* inp_par, const int n_parameter, int* status){

	assert(n_parameter == NUM_PARAM);

	// convert from degree to radian
	double deg2rad = PI/180.0;

	param->period    = inp_par[0];
	param->phase     = inp_par[1];
	param->theta     = inp_par[2]*deg2rad;
	param->incl      = inp_par[3]*deg2rad;
	param->dincl     = inp_par[4]*deg2rad;
	param->beta      = inp_par[5];
	param->dopulse   = (int) inp_par[6];

	check_par_model(param);
}

// calculate the flux array for the current value of theta
static double* interp_betTheta(ulx_table* tab,double theta, double beta, int* status){

	assert(tab!=NULL);

	int ind_th = binary_search(theta, tab->theta, tab->ntheta);
	if (ind_th == -1){
		printf(" *** error: values of theta=%.2f not tabulated",theta*180/3.1415);
		ULXLC_ERROR(" failed to look up the value of theta in the loaded table ",status);
		return NULL;
	}

	int ind_bet = binary_search(beta, tab->beta, tab->nbeta);
	if (ind_bet == -1){
		printf(" *** error: values of beta=%e not tabulated",beta);
		ULXLC_ERROR(" failed to look up the value of beta in the loaded table ",status);
		return NULL;
	}

	//	allocate memory for the flux array
	double* inp_flux = (double*) malloc (tab->nang*sizeof(double));
	CHECK_MALLOC_RET_STATUS(inp_flux,status,NULL);

	// 2d-interpolate the flux now
	double fac_th  = (theta-tab->theta[ind_th]) / (tab->theta[ind_th+1]-tab->theta[ind_th]);
	double fac_bet = (beta-tab->beta[ind_bet]) / (tab->beta[ind_bet+1]-tab->beta[ind_bet]);
	int ii;
	for (ii=0;ii<tab->nang;ii++){
		inp_flux[ii] = (1-fac_bet)*(1-fac_th) * tab->flux[ind_bet][ind_th][ii]
				     + (fac_bet)*(1-fac_th)   * tab->flux[ind_bet+1][ind_th][ii]
				     + (1-fac_bet)*(fac_th)   * tab->flux[ind_bet][ind_th+1][ii]
				     + (fac_bet)*(fac_th)     * tab->flux[ind_bet+1][ind_th+1][ii];
	}

	return inp_flux;
}

static double calc_orbit_ang(double phase, inp_param param){
	// need the sine distributed from [0,1]
	double ang = 0.5*(sin(phase*2*PI+PI/2)+1)*(param.dincl*2) + param.incl - param.dincl;
	return fabs(ang); // return the absolute value
}

static double get_flux_from_angle(double ang, double* flux, double* ang_arr, int nang, int* status){

	double ang_min = 0.0;
	double ang_max = PI/2;

	// allow angles between 90 to 180 degrees
	if ( ang > PI/2 && ang < PI){
		ang = PI - ang;
	}

	// check at the boundaries
	if (ang < ang_arr[0] && ang>=ang_min ){
		return flux[0];
	} else if (ang > ang_arr[nang-1] && ang <= ang_max){
		return flux[nang-1];
	}

	// get the corresponding index
	int ind = binary_search(ang,ang_arr,nang);

	// check if something went wrong (should not happen)
	if (ind == -1){
		printf(" *** error: angle %e [deg] is not tabulated (shoud be in [%.2e,%.2e]\n",ang*180/PI,
				ang_arr[0]*180/PI,ang_arr[nang-1]*180/PI);
		ULXLC_ERROR("failed to get the flux for the desired phase",status);
		return -1.0;
	}

	// interpolate now the flux
	double fac = (ang-ang_arr[ind]) / (ang_arr[ind+1]-ang_arr[ind]);

	return (1-fac)*flux[ind]+fac*flux[ind+1];
}

// input:  t[nt+1], photar[nt]
static void create_lc(const double* t0,const int nt,double* photar,
		double* flux, ulx_table* tab,inp_param param,int* status){

	int ii;
	for (ii=0; ii<nt; ii++){ //nt as photar[nt]

		// if we want only a single pulse

		double phase=0.0;

		if (param.dopulse==1){
			// take the phase in the middle again
			double delta_t_phase_lo = (t0[ii] - t0[0])/ param.period - param.phase ;
			double delta_t_phase_hi = (t0[ii+1] - t0[0])/ param.period - param.phase ;
			phase=0.5*(delta_t_phase_lo+delta_t_phase_hi);
			if ( ( phase < 0 ) || (phase >= 1 ) ){
				continue;
			}

		} else {

			// see which part in the period we are in (then time is between 0 and param.period)
			phase = (0.5*(t0[ii]+t0[ii+1])/param.period + param.phase) ;
			phase = phase - (double) ((int) phase); // just get everything behind the colon

		}

		// now get the corresponding angle
		double ang = calc_orbit_ang(phase,param);
		photar[ii] = get_flux_from_angle(ang,flux,tab->ang,tab->nang,status);
		CHECK_STATUS_VOID(*status);
	}

	return;
}


// basic calculation of the model
// input:  t[nt+1], photar[nt]
// output: photar
static void ulxlc_base(const double* t, const int nt, double* photar, inp_param param_struct,int* status){

	// load the GLOBAL TABLE
	if (global_ulx_table==NULL){
		ulxlc_load_table(ULXTABLE_FILENAME, &global_ulx_table, status);
		CHECK_ULXLC_ERROR("reading the table failed",status);
		CHECK_STATUS_VOID(*status);
	}
	assert(global_ulx_table!=NULL);

	// interpolate the flux for a given opening_angle/theta
	double* flux = interp_betTheta(global_ulx_table, param_struct.theta,param_struct.beta,status);
	CHECK_ULXLC_ERROR("reading of theta values from table failed",status);
	CHECK_STATUS_VOID(*status);

	create_lc(t,nt,photar, flux,global_ulx_table,param_struct,status);
	free(flux);
	CHECK_STATUS_VOID(*status);

//	free_ulxTable(&global_ulx_table);

}

void free_ulxTable(void){
	ulx_table* tab = global_ulx_table;
	if (tab!=NULL){
		free(tab->theta);
		free(tab->ang);
		free(tab->beta);

		if(tab->flux!=NULL){
			int ii,jj;
			for (jj=0; jj<tab->nbeta; jj++){

				if (tab->flux[jj]!=NULL){
					for(ii=0;ii<tab->ntheta;ii++){
						free(tab->flux[jj][ii]);
					}
					free(tab->flux[jj]);
				}

			}
			free(tab->flux);
		}
		free(tab);
	}
	// set it to NULL to know the table is not there anymore
	tab = NULL;
}

/*void free_ulxTable(){
	ulx_table** tab = &global_ulx_table;
	if ((*tab)!=NULL){
		free((*tab)->theta);
		free((*tab)->ang);
		free((*tab)->beta);

		if((*tab)->flux!=NULL){
			int ii,jj;
			for (jj=0; jj<(*tab)->nbeta; jj++){

				if ((*tab)->flux[jj]!=NULL){
					for(ii=0;ii<(*tab)->ntheta;ii++){
						free((*tab)->flux[jj][ii]);
					}
					free((*tab)->flux[jj]);
				}

			}
			free((*tab)->flux);
		}
		free(*tab);
	}
	(*tab) = NULL;
} */



int ulxlc_model(const double* t, const int nt, double* photar, const double* parameter, const int n_parameter){

	int status = EXIT_SUCCESS;

	// fill in parameters
	inp_param param_struct;

	init_par_model(&param_struct,parameter,n_parameter,&status);
	CHECK_STATUS_RET(status,status);

	// call the function which calculates the light curve
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
