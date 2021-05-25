#ifndef PMATRIX_H
#define PMATRIX_H

#include "typedef.h"

int	CreateMSAMatrix(INVINFO *InvInfo, RATES* Rates, TREES* Trees);
int	CreateMSAMatrixCoVar(INVINFO *InvInfo, RATES* Rates, TREES* Trees);

int	CreateDEPAMatrix(INVINFO* InvInfo, double *R, RATES* Rates, TREES* Trees);
int	CreateDEPAMatrixCoVar(INVINFO *InvInfo, RATES* Rates, TREES* Trees);

int	CreateInDEPAMatrix(INVINFO* InvInfo, double *R, RATES* Rates, TREES* Trees);
int	CreateInDEPAMatrixCoVar(INVINFO *InvInfo, RATES* Rates, TREES* Trees);

int	CreateDepCVAMatrix(INVINFO *InvInfo, double *R, RATES* Rates, TREES* Trees);


double	CreatFullPMatrix(double t, INVINFO	*InvInfo, MATRIX *Mat, TREES* Trees, MATRIX *A, double *Et);
double	Create4SPMat(double t, INVINFO *InvInfo, MATRIX *Mat, TREES* Trees, MATRIX *A, double *Et);
double	Create2SPMat(double t, INVINFO *InvInfo, MATRIX *Mat, TREES* Trees, MATRIX *A, double *Et);

#endif
