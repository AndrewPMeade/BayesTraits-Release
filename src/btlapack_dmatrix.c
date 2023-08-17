/*
*  BayesTriats 4.0
*
*  copyright 2022
*
*  Andrew Meade
*  School of Biological Sciences
*  University of Reading
*  Reading
*  Berkshire
*  RG6 6BX
*
* BayesTriats is free software: you can redistribute it and/or modify
* it under the terms of the GNU General Public License as published by
* the Free Software Foundation, either version 3 of the License, or
* (at your option) any later version.
* 
* This program is distributed in the hope that it will be useful,
* but WITHOUT ANY WARRANTY; without even the implied warranty of
* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
* GNU General Public License for more details.
* 
* You should have received a copy of the GNU General Public License
* along with this program.  If not, see <http://www.gnu.org/licenses/>
*
*/



#include <stdio.h>

#ifdef BTLAPACK

#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <time.h>

#include "btlapack_dmatrix.h"
#include "btlapack.h"
#include "btlin_alg.h"
#include "btlapack_dmatrix.h"

// Cholesky decomposition, blocked version
// function calls btlapack routines
int btlin_bcholeskyDMATRIX(DMATRIX* dm, double* det) {
	return btlin_bcholesky(dm->m, dm->nrows,det);
}

/*  *** Interafce with LAPACK ********* */

int btlapack_choleskyDMATRIX(DMATRIX* dm, double* det) {
	return btlapack_cholesky(dm->m, dm->nrows,det);
}

int btlapack_invcholeskyDMATRIX(DMATRIX* dm, double* det) {
	return btlapack_invcholesky(dm->m, dm->nrows,det);
}

int btlapack_choleskyUDMATRIX(DMATRIX* dm, double* det) {  // for testing
	return btlapack_choleskyU(dm->m, dm->nrows,det);
}

int btlapack_cholesky2DMATRIX(DMATRIX* dm, double* det) {
	return btlapack_cholesky2(dm->m, dm->nrows,det);
}


int btlapack_invluDMATRIX(DMATRIX* dm, double* det) {
	return btlapack_invlu(dm->m, dm->nrows,det);
}

int btlapack_ldlDMATRIX(DMATRIX* dm, double* det) {
	return btlapack_ldl(dm->m, dm->nrows, det);
}

int btlapack_invldlDMATRIX(DMATRIX* dm, double *det) {
	return btlapack_invldl(dm->m, dm->nrows, det);
}
#endif



