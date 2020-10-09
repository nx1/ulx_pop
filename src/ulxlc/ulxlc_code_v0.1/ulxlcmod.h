/*
   This file is part of the ULXLC model code.

   ULXLC is free software: you can redistribute it and/or modify it
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
#ifndef ULXLCMOD_H_
#define ULXLCMOD_H_

#define _GNU_SOURCE

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <assert.h>
#include <limits.h>
#include <float.h>
#include "fitsio.h"

#define NUM_PARAM 7

/** path to all RELXILL tables */
#define ULXLC_TABLE_PATH "."
#define ULXTABLE_FILENAME "table_ulx_v01.fits"

#define version_major 0
#define version_minor 4

/** input parameters */
typedef struct{
/** times are all given in seconds, angles in degrees */
	double period;   // time of period
	double phase;    // going from 0 to 1
	double theta;    // half of the opening angle
	double incl;     // inclination to the system at phase zero
	double dincl;    // precession in angle
	double beta;     // beta=v/c
	int dopulse;     // pulse(=1) or lightcurve(=0)?
} inp_param;

typedef struct{
	double* ang;
	int nang;

	int ntheta;
	double* theta;

	int nbeta;
	double* beta;

	double*** flux;
} ulx_table;

typedef double Real;


/***************/
/** functions **/
/***************/

#include "ulxlcbase.h"


// basic model function
int ulxlc_model(const double* t, const int nt, double* photar, const double* parameter, const int n_parameter);

// model function called by ISIS/Xspec
void ulxlcmodxspec(const Real* energy, int Nflux, const Real* parameter,int spectrum, Real* flux, Real* fluxError, const char* init);

void free_ulxTable(void);


#endif /* ULXLCMOD_H_ */
