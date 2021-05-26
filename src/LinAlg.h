/*	linalg.h
|
|	Prototypes for matrix-inversion and eigensystem functions
|
|	Copyright (c) 1998 by David L. Swofford, Smithsonian Institution.
|	All rights reserved.
|
|	NOTE: if ANSI function prototypes are not supported, define NO_PROTOTYPES
|		  before including this file.
*/

#define RC_COMPLEX_EVAL 2	/* code that complex eigenvalue obtained */

extern int  InvertMatrix (double **a, int n, double *col, int *indx, double **a_inv);
extern int  LUDecompose (double **a, int n, double *vv, int *indx, double *pd);
extern int  EigenRealGeneral (int n, double **a, double *v, double *vi, double **u, int *iwork, double *work);
extern int	MatrixDeterminant(double **a, int n, double *vv, int *indx, double *Det);
extern int	InvertMatrixAndDet (double **a, int n, double *col, int *indx, double **a_inv, double *Det);

#ifndef NO_ERROR
	#define NO_ERROR	0
#endif

#ifndef ERROR
	#define ERROR		1
#endif

#ifndef NO
	#define NO 0
#endif

#ifndef YES
	#define YES 1
#endif
