#ifndef BTOCL_DMATRIX_H
#define BTOCL_DMATRIX_H

#ifdef BTOCL

#include "btocl_runtime.h"
#include "btocl_runtime_kernels.h"
#include "btocl_lin.h"
#include "btlin_dmatrix.h"
#include "btlapack.h"

/* New interface */
int btocl_choleskyDMATRIX(cl_mem buffer, DMATRIX* dm,double* det,int cbs, int mbs);
int btocl_invcholeskyDMATRIX(cl_mem buffer, DMATRIX* dm,double* det,int cbs, int mbs);

/*  *** end new interface  ***/

//int btocl_choleskyDMATRIXbs(cl_mem buffer, BTLAPACK_DMATRIX* dm,double* det, int cbs, int mbs);

//int btocl_invdetcholeskyDMATRIX(cl_mem buffer, DMATRIX* dm,double* det,int cbs, int mbs);

//int btocl_choleskyDMATRIX(BTLAPACK_DMATRIX* m, double* det);

#endif // if BTOCL defined

#endif