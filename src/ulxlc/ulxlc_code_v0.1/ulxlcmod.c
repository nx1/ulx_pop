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
	
	Modifcations 2020 by Norman Khan
*/



#include "ulxlcmod.h"
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include <sqlite3.h>
#include <stdbool.h> 

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
		// printf("t: %d \n", ii);
		// printf("flux: %f \n", flux);
		// printf("tab->ang %lf \n", tab->ang);
		// printf("tab->nang %lf \n", tab->nang);
		// printf("Photar[t]: %f \n", photar[ii]);
		// printf("status: %d \n", *status);
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


/*
##########################
Added functions start here
##########################
*/

// Auxilary (Things that any reasonable programming language should do)
double max(double* arr, int len){
	// Get the maximum value of an array
	// arr : array
	// len : length of array
	int i;
	double t;
	t = arr[0];
	for(i=1; i<len ;i++)
        {
		if(arr[i]>t)
			t=arr[i];
	}
	return(t);
}

double min(double* arr, int len){
	// Get the minimum value of an array
	// x : array
	// len : length of array
	int i;
	double t;
	t = arr[0];
	for(i=1; i<len ;i++)
        {
		if(arr[i]<t)
			t=arr[i];
	}
	return(t);
}


int classify_curve(double lc_ulx_lim, double lc_max_flux, double lc_min_flux){
	// 0 : dead
	// 1 : transient
	// 2 : alive
	int classification;
	
	if (lc_min_flux > lc_ulx_lim){
        classification = 2;
	}
    else if (lc_max_flux < lc_ulx_lim){
        classification = 0;
	}
    else if (lc_ulx_lim < lc_max_flux && lc_ulx_lim > lc_min_flux){
        classification = 1;
	}
    return classification;
}

int get_first_transient_cycle_index(int is_ulx[]){
	int diff;
	for(int i=1; i<8; i++){
		diff = is_ulx[i] - is_ulx[i-1];
		if (diff != 0){
			return i;
		}
	}
	return -1;
}



// some sort of sqlite thingy idk
static int callback(void *data, int argc, char **argv, char **azColName) {
   int i;
   
   fprintf(stderr, "%s: \n", (const char*)data);
   
   for(i = 0; i<argc; i++) {
	  printf("%d \n", i);
      printf("%s = %s\n", azColName[i], argv[i] ? argv[i] : "NULL");
   }
   printf("\n");
   return 0;
}


int main(int argc, char *argv[]){
	// argv[1] = run_id
	// argv[2] = size (system_N)
	// argv[3] = erass_system_period_cutoff
	
	if(argc!=4){
		printf("Incorrect number of parameters passed (%d / %d) \n", argc-1, 3); 
		exit(1);
	}

	double fmod(double x, double y);	// Modulo Operator
	srand(time(NULL));					// Initialise random seed

	char* run_id = argv[1];	// run identifier

	// ULXLC parameters
	double params[7];
	params[0] =	50.0;	// period
	params[1] =	0.0;	// phase
	params[2] =	NAN;	// theta
	params[3] =	NAN;	// incl
	params[4] =	NAN;	// dincl
	params[5] = 0.3;	// beta
	params[6] = 0;		// dopulse

	// Setup Curve
	int    lc_nt = 5000;         			// Length of time series
	float  lc_timestep = 0.01; 				// Timestep in seconds
	float  lc_duration = lc_timestep*lc_nt;	// Lightcurve duration
	double lc_t[lc_nt];		   				// Time array
	double lc_flux[lc_nt];     				// Photon detection Array (Flux)

	// lightcurve properties
	double lc_period = params[0];
	double lc_zero_incl_max_flux;
	double lc_max_flux;
	double lc_min_flux;
	double lc_flux_scaling_constant;
	double lc_P_wind_time_scaling_constant;
	double lc_P_sup_time_scaling_constant;
	double lc_ulx_lim=0;
	int    lc_classification; 		// 0 = Dead, 1 = Transient, 2 = Dead
	
	// System parameters
	int    system_N = atoi(argv[2]);
	int    system_row_id[system_N];
	double system_theta[system_N];
	int    system_inclination[system_N];
	int    system_dincl[system_N];
	double system_Lx[system_N];		// In units of 10^39 erg s^-1
	double system_P_wind_days[system_N];
	double system_P_sup_days[system_N];
	
	// Not Used
	int    system_mttype = 0; 				// Mass transfer type: 0 = None, 1 = Nuclear, 2 = Thermal, 3 = WD mass transfer
	
	// eRASS Curve Sampling
	int	   erass_system_period_cutoff = atoi(argv[3]);  // Treat sources that are over a wind period length as persistent sources
	int    erass_sample_interval = 30*6;	      // Sampling interval in days
	int    erass_sample_iterations = 10000;	      // Number of Monte Carlo cycle repeats
	
	int    erass_P_wind_start_time_days[erass_sample_iterations];		// Used to hold time of first observation
	int    erass_P_sup_start_time_days[erass_sample_iterations];		// Used to hold time of first observation
	
	int    erass_P_wind_sample_time_days[erass_sample_iterations][8]; 	// Sample times in days
	int    erass_P_sup_sample_time_days[erass_sample_iterations][8]; 	// Sample times in days
	
	double erass_P_wind_lc_sample_times[erass_sample_iterations][8]; 	// Sample times on the curve
	double erass_P_sup_lc_sample_times[erass_sample_iterations][8]; 	// Sample times on the curve
	
	int    erass_P_wind_lc_sample_indexs[erass_sample_iterations][8];	// Sample time index on the curve
	int    erass_P_sup_lc_sample_indexs[erass_sample_iterations][8];	// Sample time index on the curve
	
	double erass_P_wind_lc_sample_fluxes[erass_sample_iterations][8];	// Sample fluxes on the curve
	double erass_P_sup_lc_sample_fluxes[erass_sample_iterations][8];	// Sample fluxes on the curve
	
	double erass_P_wind_sample_fluxes[erass_sample_iterations][8];		// Sample fluxes in x10^39 erg s^-1
	double erass_P_sup_sample_fluxes[erass_sample_iterations][8];		// Sample fluxes in x10^39 erg s^-1

	int   erass_P_wind_is_ulx[erass_sample_iterations][8];				// Is the sampled time observed as a ULX?
	int   erass_P_sup_is_ulx[erass_sample_iterations][8];				// Is the sampled time observed as a ULX?
	
	double erass_1_ulx_prob;											// Probability of source being a ulx on first erass cycle
	
	double erass_P_wind_persistent_prob;								// probability of the source being persistent for the entire length of eRASS
	double erass_P_sup_persistent_prob;									// probability of the source being persistent for the entire length of eRASS

	double erass_P_wind_transient_prob[8];								// probability of the source being transient for the entire length of eRASS
	double erass_P_sup_transient_prob[8];								// probability of the source being transient for the entire length of eRASS


	// SQLite
	char filename_database[9];
	sprintf(filename_database, "ulxlc.db");
	sqlite3 *db;
	char *zErrMsg = 0;
	char *sql;
	int rc;
	const char* data = "Callback function called";
	sqlite3_stmt *res;

	/* Retrieve sim rows from sql table */
	sql = sqlite3_mprintf("SELECT system_row_id, theta_half_deg, Lx1, P_wind_days, P_sup_days, inclination, dincl from ERASS_MC_SAMPLED_SYSTEMS WHERE run_id=%Q", run_id);
	// sql = sqlite3_mprintf("SELECT * from ERASS_MC_SAMPLED_SYSTEMS WHERE run_id='%Q'", run_id);
	rc = sqlite3_open(filename_database, &db);
	rc = sqlite3_prepare_v2(db, sql, -1, &res, 0);
	int r = 0;
	while(sqlite3_step(res) == SQLITE_ROW){
		system_row_id[r] 	  = sqlite3_column_int(res, 0);
		system_theta[r] 	  = sqlite3_column_double(res, 1);
		system_Lx[r] 		  = sqlite3_column_double(res, 2);
		system_P_wind_days[r] = sqlite3_column_double(res, 3);
		system_P_sup_days[r]  = sqlite3_column_double(res, 4);
		system_inclination[r] = sqlite3_column_int(res, 5);
		system_dincl[r]	      = sqlite3_column_int(res, 6);
		r++;
	}
	sqlite3_finalize(res);
	sqlite3_close(db);

	// Populate time array
	for (int i=0; i<lc_nt; i++){
		lc_t[i] = lc_timestep*i;
	}

    // Loop over all systems
	for (int N=0;N<system_N;N++){
		params[2] =	system_theta[N];	// theta
		params[3] =	0;	                // incl
		params[4] =	system_dincl[N];    // dincl
		if (system_theta[N] > 45){
			continue;   // If system has opening angle > 45 deg (unbeamed) then don't perform any simulations
		}
		// Run the lightcurve model.
		ulxlc_model(lc_t, lc_nt, lc_flux, params, 7);
		
		// 0 Inclination is treated as maximum possible Lx (Looking straight down the windcone)
		if (params[3]==0){
			// Calculate Flux normalisation
			lc_zero_incl_max_flux = max(lc_flux, lc_nt);
			lc_flux_scaling_constant = system_Lx[N] / lc_zero_incl_max_flux; // Curve Normalisation constant
			lc_ulx_lim = 1 / lc_flux_scaling_constant;	//Units of 10^39 erg s^-1
		}
		
		params[3] =	system_inclination[N];	// incl
		
		// Run the lightcurve model.
		ulxlc_model(lc_t, lc_nt, lc_flux, params, 7);
		
		// Classify the lightcurve as alive/dead/transient
		lc_max_flux = max(lc_flux, lc_nt);
		lc_min_flux = min(lc_flux, lc_nt);
		lc_classification = classify_curve(lc_ulx_lim, lc_max_flux, lc_min_flux);


    	// Lightcurve is transient & P_wind or P_sup is below limit (and rolled during ourburst for lmxrb systems)
		if ((lc_classification == 1) & (system_P_wind_days[N] < erass_system_period_cutoff || system_P_sup_days[N] < erass_system_period_cutoff) ){ 
		
			int erass_1_N_ulx = 0;		// Number of times the system was observed as a ULX on the first cycle.
			
			int erass_P_wind_N_transient[8] = {0,0,0,0,0,0,0,0};	// Probability of lc being identified as transient at cycle i
			int erass_P_sup_N_transient[8] = {0,0,0,0,0,0,0,0};	    // Probability of lc being identified as transient at cycle i
			
			int erass_P_wind_N_persistent = 0;	// Times lc was being identified as persistent over eRASS
			int erass_P_sup_N_persistent = 0;	// Times lc was being identified as persistent over eRASS
			

			for (int j=0; j<erass_sample_iterations; j++){	// MC lightcurve sampling repeats
			
				// Calculate scaling constants
				lc_P_wind_time_scaling_constant = lc_period / system_P_wind_days[N]; // Time scaling constant
				lc_P_sup_time_scaling_constant  = lc_period / system_P_sup_days[N];  // Time scaling constant
				
				// Sampling in days
				erass_P_wind_start_time_days[j] = (rand() % ((int) system_P_wind_days[N] - 0 + 1)) + 0; // Random integer between 0 and P_wind_days
				erass_P_sup_start_time_days[j]  = (rand() % ((int) system_P_sup_days[N] - 0 + 1)) + 0;	// Random integer between 0 and system_P_sup_days
				
				// printf("eRASS Sampling: j = %d / %d \n", j, erass_sample_iterations);
				// printf("------------------------------- \n");
				// printf("erass_system_period_cutoff: %d days \n", erass_system_period_cutoff);
				// printf("erass_sample_interval: %d days \n", erass_sample_interval);
				// printf("erass_P_wind_start_time_days: %d days \n", erass_P_wind_start_time_days[j]);
				// printf("erass_P_sup_start_time_days: %d days \n", erass_P_sup_start_time_days[j]);
				// printf("system_P_wind_days %f \n", system_P_wind_days[N]);
				// printf("system_P_sup_days %f \n", system_P_sup_days[N]);
				// printf("lc_P_wind_time_scaling_constant %f \n", lc_P_wind_time_scaling_constant);
				// printf("lc_P_sup_time_scaling_constant %f \n", lc_P_sup_time_scaling_constant);
				// printf("---------------------------------------------------------------------------------------------------------------- \n");
				
			    // Sample the lightcurve at each erass cycle	
				for (int i=0; i<8; i++){
					erass_P_wind_sample_time_days[j][i] = erass_P_wind_start_time_days[j] + i * erass_sample_interval;
					erass_P_sup_sample_time_days[j][i]  = erass_P_sup_start_time_days[j] + i * erass_sample_interval;
					
					erass_P_wind_lc_sample_times[j][i] =  fmod(lc_P_wind_time_scaling_constant *erass_P_wind_sample_time_days[j][i], lc_period);	// Sampling points scaled for curve
					erass_P_sup_lc_sample_times[j][i]  = fmod(lc_P_sup_time_scaling_constant * erass_P_sup_sample_time_days[j][i], lc_period);		// Sampling points scaled for curve
					
					erass_P_wind_lc_sample_indexs[j][i] = binary_search(erass_P_wind_lc_sample_times[j][i], lc_t, lc_nt);	// Get corresponding indexs
					erass_P_sup_lc_sample_indexs[j][i]  = binary_search(erass_P_sup_lc_sample_times[j][i], lc_t, lc_nt); 	// Get corresponding indexs
					
					erass_P_wind_lc_sample_fluxes[j][i] = lc_flux[erass_P_wind_lc_sample_indexs[j][i]];			 			// Find corresponding fluxes on curve
					erass_P_sup_lc_sample_fluxes[j][i]  = lc_flux[erass_P_sup_lc_sample_indexs[j][i]];			 			// Find corresponding fluxes on curve
					
					// erass_P_wind_sample_fluxes[j][i] = lc_flux_scaling_constant * erass_P_wind_lc_sample_fluxes[j][i];		// Fluxes in 1x10^39 erg s^-1 units
					// erass_P_sup_sample_fluxes[j][i]  = lc_flux_scaling_constant * erass_P_sup_lc_sample_fluxes[j][i];		// Fluxes in 1x10^39 erg s^-1 units
					
					erass_P_wind_is_ulx[j][i] = erass_P_wind_lc_sample_fluxes[j][i] > lc_ulx_lim;							// True/False is sampled flux ulx?
					erass_P_sup_is_ulx[j][i] = erass_P_sup_lc_sample_fluxes[j][i] > lc_ulx_lim;							    // True/False is sampled flux ulx?
					
					// printf("i=%d \n", i);
					// printf("erass_P_wind_start_time_days: %d days 			| erass_P_sup_start_time_days: %d days \n", erass_P_wind_start_time_days[j], erass_P_sup_start_time_days[j]);
					// printf("erass_P_wind_sample_time_days[i]: %d days 		| erass_P_sup_sample_time_days[i]: %d days \n", erass_P_wind_sample_time_days[j][i], erass_P_sup_sample_time_days[j][i]);
					// printf("erass_P_wind_lc_sample_times[i]: %f units 	| erass_P_sup_lc_sample_times[i]: %f units \n", erass_P_wind_lc_sample_times[j][i], erass_P_sup_lc_sample_times[j][i]);
					// printf("erass_P_wind_lc_sample_indexs[i]: %d 			| erass_P_sup_lc_sample_indexs[i]: %d  \n", erass_P_wind_lc_sample_indexs[j][i], erass_P_sup_lc_sample_indexs[j][i]);
					// printf("erass_P_wind_lc_sample_fluxes[i]: %f 		| erass_P_sup_lc_sample_fluxes[i]: %f\n", erass_P_wind_lc_sample_fluxes[j][i], erass_P_sup_lc_sample_fluxes[j][i]);
					// printf("erass_P_wind_sample_fluxes[i]: %f x10^39 erg s^-1 | erass_P_sup_sample_fluxes[i]: %f x10^39 erg s^-1\n", erass_P_wind_sample_fluxes[j][i], erass_P_sup_sample_fluxes[j][i]);
					// printf("erass_P_wind_is_ulx[i] %d 				| erass_P_sup_is_ulx[i] %d   \n", erass_P_wind_is_ulx[j][i], erass_P_sup_is_ulx[j][i]);
					// printf("---------------------------------------------------------------------------------------------------------------- \n");
				}
			// printf("================================================================================================================ \n");
			
			

				erass_1_N_ulx = erass_1_N_ulx + erass_P_wind_is_ulx[j][0]; // This probability is equal for both P_wind and P_orb as uniform sampling in the first cycle is irrelevant of length.
				
				int ftci_P_wind;
				int ftci_P_sup;
				
				ftci_P_wind = get_first_transient_cycle_index(erass_P_wind_is_ulx[j]);		// Get the first cycle for which the source was observed as transent
				ftci_P_sup = get_first_transient_cycle_index(erass_P_sup_is_ulx[j]);		// Get the first cycle for which the source was observed as transent
				
				// Calculate final probabilities
				if (system_P_wind_days[N] < erass_system_period_cutoff){
					if (ftci_P_wind==-1){
						erass_P_wind_N_persistent++;
					}
					else{
						erass_P_wind_N_transient[ftci_P_wind]++;
					}
				}
				else{
					erass_P_wind_N_persistent = erass_sample_iterations;	// let every observation be a persistent one if system_P_wind_days > erass_system_period_cutoff
				}
				
				if (system_P_sup_days[N] < erass_system_period_cutoff){
					if (ftci_P_sup==-1){
						erass_P_sup_N_persistent++;
					}
					else{
						erass_P_sup_N_transient[ftci_P_sup]++;
					}
						// printf("%d \n", erass_P_wind_N_transient[1]);
					}
				else{
					erass_P_sup_N_persistent = erass_sample_iterations;	// let every observation be a persistent one if system_P_sup_days > erass_system_period_cutoff
				}
				
			}
		
			erass_1_ulx_prob = (double) erass_1_N_ulx / erass_sample_iterations;
			erass_P_wind_persistent_prob = (double) erass_P_wind_N_persistent / erass_sample_iterations;
			erass_P_sup_persistent_prob = (double) erass_P_sup_N_persistent / erass_sample_iterations;
			
			for (int i=0; i<8;i++){
				erass_P_wind_transient_prob[i] = (double) erass_P_wind_N_transient[i] / erass_sample_iterations;
				erass_P_sup_transient_prob[i] = (double) erass_P_sup_N_transient[i] / erass_sample_iterations;
			}
			
			
			
			// Insert TRANSIENT into sqlite table
			char* sql2 = sqlite3_mprintf("INSERT INTO TRANSIENT VALUES (" \
										"%d, %d, %d, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f, %Q);",
										system_row_id[N],
										system_dincl[N],
										system_inclination[N],
										erass_1_ulx_prob,
										erass_P_wind_transient_prob[1],
										erass_P_wind_transient_prob[2],
										erass_P_wind_transient_prob[3],
										erass_P_wind_transient_prob[4],
										erass_P_wind_transient_prob[5],
										erass_P_wind_transient_prob[6],
										erass_P_wind_transient_prob[7],
										erass_P_sup_transient_prob[1],
										erass_P_sup_transient_prob[2],
										erass_P_sup_transient_prob[3],
										erass_P_sup_transient_prob[4],
										erass_P_sup_transient_prob[5],
										erass_P_sup_transient_prob[6],
										erass_P_sup_transient_prob[7],
										erass_P_wind_persistent_prob,
										erass_P_sup_persistent_prob,
										run_id);

			rc = sqlite3_open(filename_database, &db);
			rc = sqlite3_exec(db, sql2, callback, 0, &zErrMsg);

			if( rc != SQLITE_OK ){
				fprintf(stderr, "SQL error: %s\n", zErrMsg);
				sqlite3_free(zErrMsg);
			} else {
				// fprintf(stdout, "TRANSIENT Records created successfully\n");
			}
			sqlite3_close(db);
		}
		
		// printf("\n");
		// printf("System Information: (ID %d) (iter = %d) \n", system_row_id, N);
		// printf("--------------------------------------- \n");
		// printf("system_theta: %f deg\n", system_theta);
		// printf("system_dincl: %d deg\n", system_dincl);
		// printf("system_inclination: %d deg\n", system_inclination);
		// printf("system_Lx: %f x 10^39 erg s^-1 \n", system_Lx);
		// printf("system_P_wind_days: %f days\n", system_P_wind_days);
		// printf("system_P_sup_days: %f days\n", system_P_sup_days);
		// printf("system_mttype: %d \n", system_mttype[0]);
		// 
		// printf("\n");
		// 
		// printf("Lightcurve Information: \n");
		// printf("----------------------- \n");
		// printf("lc_nt: %d \n", lc_nt);
		// printf("lc_timestep: %f seconds \n", lc_timestep);
		// printf("lc_duration: %f seconds \n", lc_duration);	
		// printf("lc_period: %f seconds \n", lc_period);
		// printf("lc_zero_incl_max_flux: %f units\n", lc_zero_incl_max_flux);
		// printf("lc_max_flux: %f units\n", lc_max_flux);
		// printf("lc_min_flux: %f units\n", lc_min_flux);
		// printf("lc_flux_scaling_constant: %f 10^39 erg s^-1\n", lc_flux_scaling_constant);
		// printf("lc_P_wind_time_scaling_constant: %f s\n", lc_P_wind_time_scaling_constant);
		// printf("lc_ulx_lim: %f units \n", lc_ulx_lim);
		// printf("lc_classification: %d \n", lc_classification);
		// 
		// printf("\n");
		// 
		// printf("@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@ \n");
		// printf("@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@ \n");
		// printf("@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@ \n");

		// (%d, %f, %d, %d, %f, %f, %f, %d, ?)
		// Insert CLASSIFICATIONS results into sqlite table
		sql = sqlite3_mprintf("INSERT INTO CLASSIFICATIONS VALUES (%d, %f, %d, %d, %f, %f, %f, %d, %Q);", system_row_id[N], system_theta[N], system_dincl[N], system_inclination[N], lc_min_flux, lc_max_flux, lc_flux_scaling_constant, lc_classification, run_id);
		
		/* Execute SQL statement */
		rc = sqlite3_open(filename_database, &db);
		rc = sqlite3_exec(db, sql, callback, 0, &zErrMsg);

		if( rc != SQLITE_OK ){
			fprintf(stderr, "SQL error: %s\n", zErrMsg);
			sqlite3_free(zErrMsg);
		} else {
			// fprintf(stdout, "CLASSIFICATIONS Records created successfully\n");
		}
		sqlite3_close(db);
	}

	return 0;
}
