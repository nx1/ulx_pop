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
    
    To compile as shared library:
    gcc -shared -Wl,-soname,ulxlc -o ulxlc_shared.so -fPIC ulxlcbase.c ulxlc_shared.c -L /home/x1/cfitsio/lib -lcfitsio -lm -lnsl
	
	Modifcations 2020 by Norman Khan
	params[0] =	50.0;	// period
	params[1] =	0.0;	// phase
	params[2] =	NAN;	// theta
	params[3] =	NAN;	// incl
	params[4] =	NAN;	// dincl
	params[5] = 0.3;	// beta
	params[6] = 0;		// dopulse
	
*/

#include "ulxlcmod.h"
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
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



// Auxilary (Things that any reasonable programming language should do)
double max(double* arr, int len){
	// Get the maximum value of an array
	// arr : array
	// len : length of array
	int i;
	double t;
	t = arr[0];
	for(i=1; i<len; i++){
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
	for(i=1; i<len; i++){
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



double calc_Lx_prec(double Lx, double lc_max_flux_zero_incl, double lc_flux){
    // Calculate randomly sampled luminosity from light curve
	// Lx : System Luminosity
	// lc_max_flux_zero_incl : ulxlc flux for 0 inclination lc
	// lc_flux : ulxlc flux to convert to erg s^-1
    double fsc =  Lx / lc_max_flux_zero_incl; // Flux scaling constant
    double Lx_prec = lc_flux*fsc;
    return Lx_prec;
}


void xlf_calc_L_prec(const double* t, const int nt, double* photar, int lc_idx[], double Lx_prec[], double lc_classification[], double Lx[], double theta[], double incl[], double dincl[], int N){
    // ULXLC parameters
	double parameter[7];
	parameter[0] =	50.0;	// period
	parameter[1] =	0.0;	// phase
	parameter[2] =	NAN;	// theta
	parameter[3] =	NAN;	// incl
	parameter[4] =	NAN;	// dincl
	parameter[5] = 0.3;	    // beta
	parameter[6] = 0;		// dopulse
	
    double L_ULX = 1E39;

	double lc_flux;
	double lc_max_flux_zero_incl;
    double lc_max_flux;
    double lc_min_flux;
    double fsc;
    double lc_ulx_lim;

    for(int n=0;n<N;n++){
        parameter[2] = theta[n];
        if (parameter[2] > 45){ // If opening angle > 45 then assume the same luminosity
            Lx_prec[n] = Lx[n];
            lc_classification[n] = 3;
        }
        else{
            parameter[4] = dincl[n]; //dincl
            parameter[3] = 0;       // inclination
            ulxlc_model(t, nt, photar, parameter, 7);
            lc_max_flux_zero_incl = max(photar, nt);
            parameter[3] = incl[n];
            ulxlc_model(t, nt, photar, parameter, 7);
            lc_max_flux = max(photar, nt);        // Maximum flux
            lc_min_flux = min(photar, nt);        // Minimum flux
            lc_flux = photar[lc_idx[n]];          // Flux at random point on lc
            fsc =  Lx[n] / lc_max_flux_zero_incl; // Flux scaling constant
			lc_ulx_lim = L_ULX / fsc;             // lc ULX lim

            Lx_prec[n] = lc_flux*fsc;             // New luminosity
            lc_classification[n] = classify_curve(lc_ulx_lim, lc_max_flux, lc_min_flux); // Curve classification
        }
        
    }
}


void lc_boost(const double* t, const int nt, double* photar, double c_arr[][6], int N_tot, int N_save_par){
    // ULXLC parameters
	double parameter[7];
	parameter[0] =	50.0;	// period
	parameter[1] =	0.0;	// phase
	parameter[2] =	NAN;	// theta
	parameter[3] =	NAN;	// incl
	parameter[4] =	NAN;	// dincl
	parameter[5] = 0.3;	    // beta
	parameter[6] = 0;		// dopulse

    int n=0;
    
    double lc_max;
    double lc_min;
    double lc_max_flux_zero_incl;
    double lc_boost;

    for(int dincl=0;dincl<46;dincl++){
        parameter[4] = dincl;
        for(int theta=0;theta<46;theta++){
            parameter[2] = theta;
            for(int incl=0;incl<91;incl++){
                parameter[3] = incl;
                
                ulxlc_model(t, nt, photar, parameter, 7);
                lc_max = max(photar, nt);
                lc_min = min(photar, nt);
                
                if(incl == 0){
                    lc_max_flux_zero_incl = lc_max;
                }
                lc_boost = lc_max_flux_zero_incl / lc_max;
                
                c_arr[n][0] = dincl;    // dincl
                c_arr[n][1] = theta;    // thetas
                c_arr[n][2] = incl;     // incl
                c_arr[n][3] = lc_min;   // lc_min
                c_arr[n][4] = lc_max;   // lc_max
                c_arr[n][5] = lc_boost; // lc_boost
                if(n%1000==0){
                    printf("%d \n", n);
                }
                n++;
            }
        }
    }
}



// Run ULXLC over a grid of inclinations (0-90) and dincl(0-45)
int grid_ulxlc_model(double theta[], double Lx[], int N, const double* t, const int nt, double* photar, const double* parameter, const int n_parameter){
    
    // Loop over systems
    for(int n=0; n<N; n++){
        theta[n];
        Lx[n];
        for(int incl=0;incl<91;incl++){
            // parameter[3] = incl;
            for(int dincl=0;dincl<46;dincl++){
                // parameter[4] = dincl;
                // ulxlc_model(t, nt, photar, parameter, n_parameter);
            }
        }
    }
}

void print_params(double parameter[], int DEBUG){
	if (DEBUG){
		printf( "period = %.2f \t phase = %.2f \t theta = %.2f \t incl = %.2f \t dincl = %.2f \t beta = %.2f \t dopulse = %.2f \n", parameter[0], parameter[1], parameter[2], parameter[3], parameter[4], parameter[5], parameter[6]);
	}
}


#define N_sys 500

struct system{
	double id;
	double theta;
	double inclination;
	double dincl;
	double Lx;
	double period;
	double phase;
};

struct MC_input{
	double s_id[N_sys];
	double s_theta[N_sys];
	double s_inclination[N_sys];
	double s_dincl[N_sys];
	double s_Lx[N_sys];
	double s_period[N_sys];
    double s_phase[N_sys];
};


struct MC_output{
	int N_alive;
	int N_dead;
	int N_transient;
	int N_alive_unsimmed;
	int N_alive_tot;
	
	int N_ulx[8];
	int N_not_ulx[8];
	int N_new[8];
	int N_dip[8];
	int N_delta_ulx[8];
	int N_transients[8];
};

int get_N_sys(){
	int n_sys;
	n_sys = N_sys;
	return n_sys;
}


void print_system(struct MC_input *inp, int i){
	printf("i=%d \t id=%.0f \t  theta=%.2f \t  incl=%.2f \t  dincl=%.2f \t  Lx=%.2e \t  period=%.2f \t  phase=%.2f \t \n", i, inp->s_id[i], inp->s_theta[i], inp->s_inclination[i], inp->s_dincl[i], inp->s_Lx[i], inp->s_period[i], inp->s_phase[i]);
}

void print_MC_input(struct MC_input *inp){
	printf("MC_input: \n");
	printf("---------: \n");
	for (int i=0; i<N_sys;i++){
		print_system(inp, i);
	}
}

void print_MC_output(struct MC_output *out){
	printf("MC_output: \n");
	printf("---------: \n");
	printf("N_alive=%d \n", out->N_alive);
	printf("N_dead=%d \n", out->N_dead);
	printf("N_transient=%d \n", out->N_transient);
	printf("N_alive_unsimmed=%d \n", out->N_alive_unsimmed);
	printf("N_alive_tot=%d \n", out->N_alive_tot);
	printf("\n");
	printf("Transient classification eRASS Evolution: \n");
	for (int i=0; i<8;i++){
		printf("i=%d \t N_ulx[i]=%d \t N_not_ulx[i]=%d \t N_new[i]=%d \t N_dip[i]=%d \t N_delta_ulx[i]=%d \t N_transients[i]=%d \n", i, out->N_ulx[i], out->N_not_ulx[i], out->N_new[i], out->N_dip[i], out->N_delta_ulx[i], out->N_transients[i]);
	}
}

int sim(struct MC_input *inp, struct MC_output *out){

	int DEBUG = 0; //0 or 1 for print statements
    
	double PERIOD_DEFAULT = 50.0;
	double BETA_DEFAULT = 0.3;
	double L_ULX = 1E39;
	
    // ULXLC
	// LC Parameters
	double lc_par[7];
	lc_par[0] = PERIOD_DEFAULT; // period
	lc_par[1] = NAN;			// phase
	lc_par[2] = NAN;		  	// theta
	lc_par[3] = 0;   	 	  	// incl
	lc_par[4] = NAN;		  	// dincl
	lc_par[5] = BETA_DEFAULT; 	// beta
	lc_par[6] = 0;			  	// dopulse

    // Setup Curve
	int    lc_nt = 5000;         			// Length of time series
	float  lc_timestep = 0.01; 				// Timestep in seconds
	float  lc_duration = lc_timestep*lc_nt;	// Lightcurve duration
	double lc_t[lc_nt];		   				// Time array
	double lc_flux[lc_nt];     				// Photon detection Array (Flux)

    // lightcurve properties
	double lc_period = lc_par[0];
	double lc_zero_incl_max_flux;
	double lc_max_flux;
	double lc_min_flux;
	double lc_max_L;
	double lc_min_L;
	double lc_flux_scaling_constant;
	double lc_P_time_scaling_constant; // c * time_in_days = time_on_curve
	double lc_ulx_lim;
	int    lc_classification; 		   // 0 = Dead, 1 = Transient, 2 = Dead, 3 = Alive with theta > 45
	
	// eRASS
	// Sample parameters
	int	   erass_s_period_cutoff = 99999;  		// Treat sources that are over a wind period length as persistent sources
	int    erass_N_cycles = 8;	        		// Number of cycles
	int    erass_sample_interval = 30*6;	    // Sampling interval in days
	
	//Curve Sampling
	int    erass_t_samp[erass_N_cycles];		// Sample time in days
	double erass_lc_samp[erass_N_cycles];		// Sample time in curve space
	int    erass_lc_samp_idx[erass_N_cycles];	// Sample indexs in curve space
	double erass_lc_samp_flux[erass_N_cycles];	// Sample flux in curve space
	double erass_lc_samp_L[erass_N_cycles];		// Sample luminosity in erg s^-1
	int    erass_lc_is_ulx[erass_N_cycles];		// Is sampled luminosity above 1e39?
	
	int erass_lc_obs_as_ulx;				    // Has the source been observed as a ulx?
	int erass_lc_first_transient_cycle;			// Cycle at which lightcurve was found to be transient (if not found then -1)
	int erass_lc_transient_classification; 		// 0: from moved up to ULX L        1: Moved down from ulx L
	
	// MC output
	// ULX Counting params
	out->N_alive     = 0;
	out->N_dead      = 0;
	out->N_transient = 0;
	out->N_alive_unsimmed = 0;
	out->N_alive_tot = 0;
	
	// Initialize counting arrays
	for (int i=0; i<erass_N_cycles; i++){
		out->N_ulx[i] = 0;
		out->N_not_ulx[i] = 0;
		out->N_new[i] = 0;
		out->N_dip[i] = 0;
		out->N_delta_ulx[i] = 0;
		out->N_transients[i] = 0;
	}
	
	// Initialize time array
	for (int i=0; i<lc_nt; i++){
		lc_t[i] = lc_timestep*i;
	}

    // DEBUGGING INFORMATION
    if (DEBUG){
        printf("DEBUG MODE ON \n");
		printf("------------- \n \n");
		
		printf("Setup Parameters: \n");
		printf("PERIOD_DEFAULT = %.2f \t BETA_DEFAULT = %.2f \n", PERIOD_DEFAULT, BETA_DEFAULT);
		printf("lc_nt = %d \t lc_timestep = %.2f \n", lc_nt, lc_timestep);
		printf("erass_s_period_cutoff = %d \t erass_sample_interval = %d \n", erass_s_period_cutoff, erass_sample_interval);
		printf("\n \n");
		
		print_MC_input(inp);
		print_MC_output(out);
	}
    
	
	
    // Loop over all systems
	for (int N=0;N<N_sys;N++){
		if (DEBUG){
			printf("====================== Loop Start N=%d ====================\n", N);
		}

		if (inp->s_theta[N] > 45){
			if (DEBUG){
				printf("inp->s_theta[N] = %.2f is > 45 , lc_classification=3 \n \n", inp->s_theta[N]);
			}
			lc_classification = 3;
		}
		
		else{
			lc_par[1] = inp->s_phase[N];	// phase
			lc_par[2] =	inp->s_theta[N];	// theta
			lc_par[3] =	0;	        // incl
			lc_par[4] =	inp->s_dincl[N]; // dincl
			
			print_params(lc_par, DEBUG);
			
			// Run the lightcurve model.
			ulxlc_model(lc_t, lc_nt, lc_flux, lc_par, 7);

			// 0 Inclination is treated as maximum possible Lx (Looking straight down the windcone)
			if (lc_par[3]==0.0){
				// Calculate Flux normalisation
				lc_zero_incl_max_flux = max(lc_flux, lc_nt);
				lc_flux_scaling_constant = inp->s_Lx[N] / lc_zero_incl_max_flux; // Curve Normalisation constant
				lc_ulx_lim = L_ULX / lc_flux_scaling_constant;
				
				if (DEBUG){
					printf("inp->s_Lx[N] = %.2e \t lc_zero_incl_max_flux = %.2f \t lc_flux_scaling_constant = %.2e \t lc_ulx_lim = %.2e \n \n", inp->s_Lx[N], lc_zero_incl_max_flux, lc_flux_scaling_constant, lc_ulx_lim);
				}
				
			}
			
			lc_par[3] =	inp->s_inclination[N];	// incl
			
			
			print_params(lc_par, DEBUG);

			
			// Run the lightcurve model.
			ulxlc_model(lc_t, lc_nt, lc_flux, lc_par, 7);
			
			// Classify the lightcurve as alive/dead/transient
			lc_max_flux = max(lc_flux, lc_nt);
			lc_min_flux = min(lc_flux, lc_nt);
			lc_classification = classify_curve(lc_ulx_lim, lc_max_flux, lc_min_flux);
			lc_max_L    =  lc_max_flux * lc_flux_scaling_constant;
			lc_min_L    =  lc_min_flux * lc_flux_scaling_constant;
			
			if (DEBUG){
				printf("lc_max_flux = %.2f \t lc_min_flux = %.2f \t lc_max_L = %.2e \t lc_min_L = %.2e  \t lc_classification = %d \n", lc_max_flux, lc_min_flux, lc_max_L, lc_min_L, lc_classification);
			}


			// Lightcurve is transient & P_wind or P_sup is below limit (and rolled during ourburst for lmxrb systems)
			if ((lc_classification == 1) & (inp->s_period[N] < erass_s_period_cutoff) ){ 

				// Calculate time scaling constant
				lc_P_time_scaling_constant = lc_period / inp->s_period[N]; // Time scaling constant
				
				erass_lc_first_transient_cycle = -1;
				erass_lc_transient_classification = -1;
				erass_lc_obs_as_ulx = 0;

				if (DEBUG){
					printf("lc_classification = %d (transient), performing erass sampling \n", lc_classification);
					printf("s_period[N] = %.2f \t lc_P_time_scaling_constant = %.2f \n", inp->s_period[N], lc_P_time_scaling_constant);
				}

				// Sample the lightcurve
				for (int c=0; c<erass_N_cycles; c++){
						erass_t_samp[c] = c * erass_sample_interval; 					 	   // Sample time in days
						erass_lc_samp[c] = erass_t_samp[c] * lc_P_time_scaling_constant; 	   // Sample time in curve space
						erass_lc_samp[c] = fmod(erass_lc_samp[c], lc_period); 		     	   // To account for overflowing times
						erass_lc_samp_idx[c] = binary_search(erass_lc_samp[c], lc_t, lc_nt);   // Find corresponding indexs on the curve
						erass_lc_samp_flux[c] = lc_flux[erass_lc_samp_idx[c]]; 				   // Find corresponding L on curve
						erass_lc_samp_L[c] = lc_flux_scaling_constant * erass_lc_samp_flux[c]; // Find corresponding L in ergs^-1
						erass_lc_is_ulx[c] = erass_lc_samp_flux[c] > lc_ulx_lim;			   // Is the source a ulx?

						// If we are past cycle 0 check transient status
						if ((c>0) & (erass_lc_first_transient_cycle==-1)){
							int prev = erass_lc_is_ulx[c] - erass_lc_is_ulx[c-1];
							if (prev!=0){ // If there was a change in ulx status, identify cycle as transient.
								erass_lc_first_transient_cycle = c;
								erass_lc_obs_as_ulx = 1;
							}
							if (prev == 1){ // Source moved upto ulx luminosity
								out->N_new[c]++;
							}
							else if (prev == -1){ // Source moved down from ulx luminosity
								out->N_dip[c]++;
							}
						}
						
						// Update other quantities
						if (erass_lc_is_ulx[c] == 0){
							out->N_not_ulx[c]++;
						}
						else if (erass_lc_is_ulx[c] == 1){
							out->N_ulx[c]++;
							if (c == 0){
								out->N_new[c]++;
							}
						}
						
						
						out->N_delta_ulx[c] = out->N_new[c] - out->N_dip[c];
						if (c>0)
							out->N_transients[c] = out->N_new[c] + out->N_dip[c];
						
						if (DEBUG){
							printf("Cycle %d / %d \t erass_t_samp[c] = %5d \t erass_lc_samp[c] = %.2f \t erass_lc_samp_idx[c] = %5d \t erass_lc_samp_flux[c] = %.2e \t erass_lc_samp_L[c] = %.2e \t erass_lc_is_ulx[c] = %d \n", c+1, erass_N_cycles, erass_t_samp[c], erass_lc_samp[c], erass_lc_samp_idx[c], erass_lc_samp_flux[c], erass_lc_samp_L[c], erass_lc_is_ulx[c]);
						}
				}
				
			}
		}
		
		// Add classifications to totals
		if      (lc_classification == 0){out->N_dead++;}
		else if (lc_classification == 1){out->N_transient++;}
		else if (lc_classification == 2){out->N_alive++;}
		else if (lc_classification == 3){out->N_alive_unsimmed++;}
		
		out->N_alive_tot = out->N_alive + out->N_alive_unsimmed;
		
	    if (DEBUG){	
		    printf("\n");
			print_system(inp, N);
			
		    printf("\n");
		    
		    printf("Lightcurve Information: \n");
		    printf("----------------------- \n");
		    printf("lc_nt: %d \n", lc_nt);
		    printf("lc_timestep: %.2f seconds \n", lc_timestep);
		    printf("lc_duration: %.2f seconds \n", lc_duration);	
		    printf("lc_period: %.2f seconds \n", lc_period);
		    printf("lc_zero_incl_max_flux: %.2f units\n", lc_zero_incl_max_flux);
		    printf("lc_max_flux: %.2f units\n", lc_max_flux);
		    printf("lc_min_flux: %.2f units\n", lc_min_flux);
		    printf("lc_flux_scaling_constant: %.2e  erg s^-1\n", lc_flux_scaling_constant);
		    printf("lc_P_time_scaling_constant: %.2e s\n", lc_P_time_scaling_constant);
		    printf("lc_ulx_lim: %.2f units \n", lc_ulx_lim);
		    printf("lc_classification: %d \n", lc_classification);
			printf("erass_lc_first_transient_cycle: %d \n", erass_lc_first_transient_cycle);
			
		    printf("\n");
			
			print_MC_output(out);

		    printf("====================== Loop End N=%d ======================\n\n\n", N);
        }
	}

	
	return 0;
}

