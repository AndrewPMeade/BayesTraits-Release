#ifndef CBLAS_WRAPPER_H
#define CBLAS_WRAPPER_H

#include "TypeDef.h"

void MatrixVectProduct(double* M, double* V, double* U, int N);
void VectMatrixProduct(double* V, double* M, double* U, int N);

double VectVectDotProduct(double *V1, double *V2, int N);

#endif