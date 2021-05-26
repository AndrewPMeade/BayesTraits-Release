
#include <stdio.h>

#ifdef BTLAPACK

#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <time.h>

#include "btlapack_dmatrix.h"
#include "btlapack.h"


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



