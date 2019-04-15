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

#include "ulxlcbase.h"
#include "ulxlcmod.h"

ulx_table* new_ulxTable(int* status){
	ulx_table* tab = (ulx_table*) malloc (sizeof(ulx_table));
	CHECK_MALLOC_RET_STATUS(tab,status,tab);

	tab->ang=NULL;
	tab->flux=NULL;
	tab->nang=0;
	tab->ntheta=0;
	tab->theta=NULL;
	tab->beta=NULL;
	tab->nbeta=0;

	return tab;
}


/** read dimensions from the header of the FITS file   */
static int get_ulxtable_ndim(char* extname, fitsfile* fptr, int* status){

	int extver = 0;
	fits_movnam_hdu(fptr, BINARY_TBL, extname, extver ,status);
	if (*status!=EXIT_SUCCESS){
		printf(" *** error moving to extension %s\n",extname);
		return 0;
	}

	long ndim;
	if (fits_get_num_rows(fptr, &ndim, status)) return 0;

	// integer number of dimension is certainly enough
	return (int) ndim;
}

/** read one axis of the rel table from the FITS file   */
static void get_ulxtable_axis(int nrows, double** val, char* extname, char* colname, fitsfile* fptr, int* status){

	int extver = 0;
	fits_movnam_hdu(fptr, BINARY_TBL, extname, extver ,status);
	if (*status!=EXIT_SUCCESS){
		printf(" *** error moving to extension %s\n",extname);
		return;
	}

	// get the column id-number
	int colnum;
	if(fits_get_colnum(fptr, CASEINSEN, colname, &colnum, status)) return;

	// get the number of rows
	long n;
	if (fits_get_num_rows(fptr, &n, status)) return;

	if (nrows != n){
		ULXLC_ERROR("wrong dimension of at least one axis given in the ulx_table",status);
	}

	// allocate memory for the array
	*val=(double*)malloc(n*sizeof(double));
	CHECK_MALLOC_VOID_STATUS(*val,status);

    int anynul=0;
    double nullval=0.0;
    LONGLONG nelem = (LONGLONG) n;
    fits_read_col(fptr, TDOUBLE, colnum, 1, 1, nelem ,&nullval,*val, &anynul, status);

	return;
}

/** load the complete ulxlc table */
void ulxlc_load_table(char* filename, ulx_table** inp_tab, int* status){

	ulx_table* tab = (*inp_tab);
	fitsfile *fptr=NULL;

	char* fullfilename=NULL;

	do{ // Error handling loop
		if (tab != NULL){
			ULXLC_ERROR("relline table already loaded",status);
			break;
		}

		tab = new_ulxTable(status);
		CHECK_STATUS_BREAK(*status);

		// should be set by previous routine
		assert(tab!=NULL);


		// get the location of the table
		char* table_path;
		if (getenv("ULXLC_TABLE_PATH")!=NULL){
			table_path = strdup(getenv("ULXLC_TABLE_PATH"));
		} else {
			table_path = strdup(ULXLC_TABLE_PATH);
		}


		// get the full filename
		if (asprintf(&fullfilename, "%s/%s", table_path,filename) == -1){
			ULXLC_ERROR("failed to construct full path the ulx-table",status);
			break;
		}
		free(table_path);

		// open the file
		if (fits_open_table(&fptr, fullfilename, READONLY, status)) {
			CHECK_ULXLC_ERROR("opening of the ulx table failed",status);
			printf("    full path given: %s \n",fullfilename);
			break;
		}

		// first read the dimensions from the table
		tab->nang = get_ulxtable_ndim(EXTNAME_TAB_ANG, fptr, status);
		CHECK_ULXLC_ERROR("reading of dimensions of the table failed",status);

		tab->nbeta = get_ulxtable_ndim(EXTNAME_TAB_BET, fptr, status);
		CHECK_ULXLC_ERROR("reading of dimensions of the table failed",status);


		tab->ntheta = get_ulxtable_ndim(EXTNAME_TAB_DATA, fptr, status);
		CHECK_ULXLC_ERROR("reading of dimensions of the table failed",status);


		if (CHATTER>-1){
			printf("\n  Loading ULXLC table %s (dimensions [bet=%i, theta=%i, incl=%i])\n",fullfilename,
					tab->nbeta, tab->ntheta,tab->nang);
		}

		if (tab->beta != NULL){
			CHECK_ULXLC_ERROR("allocation problem: reading of angles from table failed",status);
		}

		tab->flux = (double***) malloc (tab->nbeta*sizeof(double**));
		CHECK_MALLOC_VOID_STATUS(tab->flux,status);
		int ii, jj;
		for (jj=0; jj<tab->nbeta; jj++){
			tab->flux[jj] = (double**) malloc (tab->ntheta*sizeof(double*));
			CHECK_MALLOC_VOID_STATUS(tab->flux[jj],status);
			for (ii=0; ii<tab->ntheta; ii++){
				tab->flux[jj][ii] = (double*) malloc (tab->nang*sizeof(double));
				CHECK_MALLOC_VOID_STATUS(tab->flux[jj][ii],status);
			}
		}

		// now read the values from the table
		// (1) first the axes
		get_ulxtable_axis(tab->nang, &(tab->ang), EXTNAME_TAB_ANG , "ang", fptr, status);
		CHECK_ULXLC_ERROR("reading of angles from table failed",status);

		get_ulxtable_axis(tab->nbeta, &(tab->beta), EXTNAME_TAB_BET , "beta", fptr, status);
		CHECK_ULXLC_ERROR("reading of beta values from table failed",status);

		get_ulxtable_axis(tab->ntheta, &(tab->theta), EXTNAME_TAB_DATA , "theta", fptr, status);
		CHECK_ULXLC_ERROR("reading of theta values from table failed",status);

		if (CHATTER==1){
			printf("  - Beta  values range [%.2f,%.2f] \n",tab->beta[0],tab->beta[tab->nbeta-1]);
			printf("  - Theta values range [%.1f,%.1f] \n",tab->theta[0]*180/3.1415,tab->theta[tab->ntheta-1]*180/3.1415);
			printf("  - Angle values range [%.1f,%.1f] \n",tab->ang[0]*180/3.1415,tab->ang[tab->nang-1]*180/3.1415);
		}

		// (2) and then the flux array
		int colnum_flux;
		char* colname_flux = NULL;

		for (jj=0; jj<tab->nbeta;jj++){

			// get the colname
			if (asprintf(&colname_flux, "fr_bet%.2fc", tab->beta[jj]) == -1){
				ULXLC_ERROR("failed to construct extname of the ulx-table",status);
				break;
			}

			// get the number and verify the column exists
			if(fits_get_colnum(fptr, CASEINSEN, colname_flux, &colnum_flux, status)){
				printf(" *** error: can not find column %s in the table %s \n",colname_flux,fullfilename);
				ULXLC_ERROR("failed to read the table",status);
				return;
			}
			free(colname_flux);

			int anynul=0;
			double nullval=0.0;
			assert(tab->flux[jj]!=NULL);
			LONGLONG nelem = (LONGLONG) tab->nang;

			for (ii=0; ii<tab->ntheta;ii++){
				if(fits_read_col(fptr, TDOUBLE, colnum_flux, ii+1, 1, nelem ,&nullval,tab->flux[jj][ii], &anynul, status)){
					CHECK_ULXLC_ERROR("reading data from table failed",status);
					return;
				}

			}
		}


	} while(0);

	if (*status==EXIT_SUCCESS){
		// assigne the value
		(*inp_tab) = tab;
	} else {
		printf(" *** error: problem loading the table \n");
		free_ulxTable();
		return;
	}
	free(fullfilename);

	if (fptr!=NULL) {fits_close_file(fptr,status);}

	return;
}


void ulxlc_error(const char* const func, const char* const msg, int* status){
	*status = EXIT_FAILURE;
	printf(" *** error in ulxlc (%s): %s!\n", func, msg);
}

void check_ulxlc_error(const char* const func, const char* const msg, int* status){
	if (*status!=EXIT_SUCCESS){
		*status = EXIT_FAILURE;
		printf(" *** error in ulxlc (%s): %s!\n", func, msg);
	}
}

void get_version_number(char** vstr, int* status){
	if (asprintf(vstr, "%i.%i", version_major, version_minor) == -1){
		ULXLC_ERROR("failed to get version number",status);
	}
}

/**  Binary search for to find interpolation interval
 *   - return value is the bin [ind,ind+1]
 *   - assume list is sorted ascending */
int binary_search(double val, double* arr, int n){

	if (val < arr[0] || val > arr[n-1]){
		return -1;
	}

	int high=n-1;
	int low=0;
	int mid;
	while (high > low) {
		mid=(low+high)/2;
		if (arr[mid] <= val) {
			low=mid+1;
		} else {
			high=mid;
		}
	}
	return low-1;
}
