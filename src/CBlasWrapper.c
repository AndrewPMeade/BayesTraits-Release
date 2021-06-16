#include "TypeDef.h"
#ifdef BTLAPACK

#include <cblas.h>

/* Multiple Matrix by Vector, can be change for non squar matrixes*/
void MatrixVectProduct(double* M, double* V, double* U, int N)
{ 
	cblas_dgemv(CblasRowMajor, CblasNoTrans, N, N, 1.0, M, N, V, 1, 0, U, 1);
}

/* 
	Multiple Vector by Matrix, can be change for non squar matrixes
	Note: same code as above but using tranpos 
*/
void VectMatrixProduct(double* V, double* M, double* U, int N)
{
	cblas_dgemv(CblasRowMajor, CblasTrans, N, N, 1.0, M, N, V, 1, 0, U, 1);
}

/*
	sum of the product of two vectors
*/
double VectVectDotProduct(double *V1, double *V2, int N)
{
	return cblas_ddot(N, V1, 1, V2, 1);
}

#else
/* Multiple Matrix by Vector, can be change for non squar matrixes*/
void MatrixVectProduct(double* M, double* V, double* U, int N)
{ 
	int y, x;
	double Sum;

	for(y=0;y<N;y++)
	{
		Sum = 0;
		for(x=0;x<N;x++)
			Sum += V[x] * M[y * N + x];

		U[y] = Sum;
	}
}

/* 
Multiple Vector by Matrix, can be change for non squar matrixes
Note: same code as above but using tranpos 
*/
void VectMatrixProduct(double* V, double* M, double* U, int N)
{
	int		x,y;
	double	Sum;

	for(y=0;y<N;y++)
	{
		Sum = 0;
		for(x=0;x<N;x++)	
			Sum += M[x * N + y] * V[x];


		U[y] = Sum;
	}
}

double VectVectDotProduct(double *V1, double *V2, int N)
{
	double Ret;
	int x;

	Ret = 0;
	for(x=0;x<N;x++)
		Ret = Ret + (V1[x] * V2[x]);

	return Ret;
}

#endif