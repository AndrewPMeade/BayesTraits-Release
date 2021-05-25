#ifndef BTLIN_ALG_H
#define BTLIN_ALG_H

#include "btlin_basic.h"

#ifdef BTLAPACK
#include "btlapack.h"
#endif


int btlin_kronecker_vectmultOne(double* vres, double* v,double* sigma, int sigma_dim,double*  mat, int mat_dim);
int btlin_kronecker_vectmult(double* vres, double* v,double* sigma, int sigma_dim,double*  mat, int mat_dim);
int btlin_kronecker_vectmult2(double* vres, double* v,double* sigma, int sigma_dim,double*  mat, int mat_dim);
int btlin_kronecker_vectmult2Row(double* vres, double* v,double* sigma, int sigma_dim,double*  mat, int mat_dim);


int btlin_cholesky(double* m, int n,double* det);
int btlin_invcholesky(double* m, int n, double *det);

// computes determinant of lower triangular matrix
double btlin_ltri_det(double* a, int n, int lda);

#ifdef BTLAPACK
// blocked version
int btlin_bcholesky(double* m, int n, double* det);
int btlin_ltri_invhybrid(double* m, int n, double* det);
#endif

#endif