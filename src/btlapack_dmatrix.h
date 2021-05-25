#ifndef BTLAPACK_DMATRIX_H
#define BTLAPACK_DMATRIX_H

#ifdef BTLAPACK

#include "btlin_dmatrix.h"

// Blocked versions - hybrid
int btlin_bcholeskyDMATRIX(DMATRIX* m, double* det);


/* *** Interface with BTLAPACK *** */

int btlapack_choleskyDMATRIX(DMATRIX* m, double* det);
int btlapack_choleskyUDMATRIX(DMATRIX* m, double* det);  // for testing
int btlapack_cholesky2DMATRIX(DMATRIX* m, double* det); // unblocked
int btlapack_invcholeskyDMATRIX(DMATRIX* m, double* det);
int btlapack_invluDMATRIX(DMATRIX* dm, double* det);
int btlapack_ldlDMATRIX(DMATRIX* dm, double* det);
int btlapack_invldlDMATRIX(DMATRIX* dm, double *det);

#endif

#endif
