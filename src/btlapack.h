#ifndef BTLAPACK_H
#define BTLAPACK_H

#ifdef BTLAPACK

#include "btlin_basic.h"

// ****************** Interface with LAPACK and BLAS ****************************

// LU decomoposition of a general matrix
extern void dgetrf_(int* M, int *N, double* A, int* lda, int* IPIV, int* INFO);

// generate inverse of a matrix given its LU decomposition
extern void dgetri_(int* N, double* A, int* lda, int* IPIV, double* WORK, int* lwork, int* INFO);

// calculates the LDL factorisation of symmetric (indefinite) matrix A
extern void dsytrf_(char* uplo, int* N, double* A, int* lda, int* IPIV, double* WORK, int* LWORK, int* INFO);
// Computes the inverse of a real symmetric indefinite matrix A using the UDU/LDL
// factorization computed by dsytrf
extern void dsytri_(char* uplo, int* N, double* A, int* lda, int* IPIV, double* WORK, int* INFO);

// calculates the choleski factorisation of symmetric positive definite matrix A
extern void dpotrf_(char* uplo, int* N, double* A, int* lda, int* INFO);

// compute the Cholesky factorization of a real symmetric positive definite matrix A
// This is the unblocked version of the algorithm, calling Level 2 BLAS
extern void dpotf2_(char* uplo, int* N, double* A, int* lda, int* INFO);

// computes the inverse of the matrix given as input to dpotrf
extern void dpotri_(char* uplo, int* N, double* A, int* lda, int* INFO);

// *********** Lower Triangular Inverse  ****************
// Unblocked version
// SUBROUTINE DTRTI2( UPLO, DIAG, N, A, LDA, INFO )
extern void dtrti2_(char* uplo, char* diag, int* n, double* A, int* lda, int* INFO);

// Blocked version
// SUBROUTINE DTRTRI( UPLO, DIAG, N, A, LDA, INFO )
extern void dtrtri_(char* uplo, char* diag, int* n, double* A, int* lda, int* INFO);

// ************************************************************************

int btlapack_cholesky(double* m, int n, double* det);
int btlapack_invcholesky(double* m, int n, double* det);
int btlapack_choleskytU(double* m, int n, double* det);  // for testing
int btlapack_cholesky2(double* m, int n, double* det);  // unblocked

int btlapack_invlu(double* m, int n, double* det);
int btlapack_ldl(double* m, int n, double* det);
int btlapack_invldl(double* m, int n, double *det);
int btlapack_invldlW(double* m, double* w, int n, double *det);

// Inverse lower triangular
int btlapack_ltri_inv2(double* m, int n, int lda, double* det);  // unblocked
int btlapack_ltri_inv(double* m, int n, int lda, double* det);


// Auxiliary functions
// computes determinant of lower triangular matrix
double btlapack_ltri_det(double* a, int n, int lda);

// new naming/argument style - come back later
int btlapack_ducholesky(char uplo,double* m,int n,int lda, double* det);
int btlapack_ducholeskydet(char uplo,double* m,int n,int lda, double* det);

#endif  // if BTLAPACK defined

#endif
