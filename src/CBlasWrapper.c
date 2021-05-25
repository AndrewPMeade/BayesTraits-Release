#include "TypeDef.h"
#ifdef BTLAPACK


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
	double Sum;
	int Index;


	Sum = cblas_dasum(N, V1, 1);
	printf("%f\n", Sum);

	Sum = 0;
	for(Index=0;Index<Sum;Index++)
	{
		
	}
	return cblas_ddot(N, V1, 1, V2, 1);
}

#endif