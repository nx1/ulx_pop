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
#ifndef ULXLCBASE_H_
#define ULXLCBASE_H_

#include "ulxlcmod.h"

#define EXTNAME_TAB_ANG  "angles"
#define EXTNAME_TAB_DATA "data"
#define EXTNAME_TAB_BET "beta"

#define PI (3.141592653589793)

#define CHATTER 0

#define ULXLC_ERROR(msg,status) (ulxlc_error(__func__, msg,status))

#define CHECK_ULXLC_ERROR(msg,status) (check_ulxlc_error(__func__, msg,status))

#define CHECK_STATUS_RET(status, retval) \
  if (EXIT_SUCCESS!=status) return(retval);

#define CHECK_STATUS_VOID(status) \
  if (EXIT_SUCCESS!=status) return;

#define CHECK_STATUS_BREAK(status) \
  if (EXIT_SUCCESS!=status) break;

#define CHECK_MALLOC_VOID_STATUS(a,status) \
	if (NULL==a) { \
		ULXLC_ERROR("memory allocation failed\n",status); \
		return;\
	}

#define CHECK_MALLOC_RET_STATUS(a,status,retval) \
	if (NULL==a) { \
		ULXLC_ERROR("memory allocation failed\n",status); \
		return retval;\
	}

#define CHECK_STATUS_RET(status, retval) \
  if (EXIT_SUCCESS!=status) return(retval);



/***************/
/** functions **/
/***************/

/* load the table */
void ulxlc_load_table(char* fname,ulx_table** tab, int* status);

/* print ulxlc error message */
void ulxlc_error(const char* const func, const char* const msg, int* status);

/* check and print ulxlc error message */
void check_ulxlc_error(const char* const func, const char* const msg, int* status);

/* get the current version number */
void get_version_number(char** vstr, int* status);

/* binary search */
int binary_search(double val, double* arr, int n);

#endif /* ULXLCBASE_H_ */
