#ifndef PMATRIX_H
#define PMATRIX_H

#include "typedef.h"

int	CreateMSAMatrix(INVINFO *InvInfo, int NOS, double *Rates, double *Pis);
int	CreateMSAMatrixCoVar(INVINFO *InvInfo, RATES* Rates, TREES* Trees, double *RateP, double *Pi);


int	CreateDEPAMatrix(INVINFO* InvInfo, RATES* Rates, TREES* Trees, double *RateP);
int	CreateDEPAMatrixCoVar(INVINFO *InvInfo, RATES* Rates, TREES* Trees, double *RateP);

int	CreateInDEPAMatrix(INVINFO* InvInfo, RATES* Rates, TREES* Trees, double *RateP);
int	CreateInDEPAMatrixCoVar(INVINFO *InvInfo, RATES* Rates, TREES* Trees, double *RateP);

int	CreateDepCVAMatrix(INVINFO *InvInfo, RATES* Rates, TREES* Trees, double *R);


double	CreatFullPMatrix(double t, INVINFO	*InvInfo, MATRIX *Mat, int NOS, int ThrNo);
double	Create4SPMat(double t, INVINFO *InvInfo, MATRIX *Mat, int ThrNo);
double	Create2SPMat(double t, INVINFO *InvInfo, MATRIX *Mat, int ThrNo);

#endif
