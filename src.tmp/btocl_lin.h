
#ifndef BTOCL_LIN_H
#define BTOCL_LIN_H

#ifdef BTOCL

#include "btocl_runtime.h"
#include "btocl_kernels_bayestraits.h"
#include "btlin_basic.h"
#include "btlapack.h"

int btocl_writeLtoU(cl_mem buffer,int mat_idx, int mat_dim, int lda, int mbs);
int btocl_cholesky(cl_mem buffer, double* mat, int mat_dim, double* det,int cbs, int mbs);
int btocl_invcholesky(cl_mem buffer, double* mat, int mat_dim, double* det,int cbs, int mbs);
int btocl_invcholesky_pure(cl_mem buffer, double* mat, int mat_dim, double* det,int cbs, int mbs);
int btocl_ltri_inv(cl_mem buffer, double* mat, int mat_dim, double* det, int nb);

int btocl_ltri_LTbyLOld(cl_mem buffer,int n);
int btocl_ltri_LTbyL(cl_mem buffer,int n);
//int btocl_cholesky_unblocked(cl_mem buffer, int mat_idx, int mat_dim,double* det);

//int btocl_choleskyDMATRIX(BTLAPACK_DMATRIX* m, double* det);

// Kronecker Product with vector multiplication
int btocl_kronecker_vectmult(cl_mem vres, cl_mem v, cl_mem sigma, int sigma_dim, cl_mem mat, int mat_dim);
int btocl_kronecker_vectmult_one(cl_mem vres, cl_mem v, double sigma, cl_mem mat, int mat_dim);




#endif // if BTOCL defined

#endif