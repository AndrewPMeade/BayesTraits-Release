#ifndef PMATRIX_H
#define PMATRIX_H

#include "typedef.h"

void	CreateMSAMatrix(INVINFO *InvInfo, int NOS, double *Rates, double *Pis);
void	CreateMSAMatrixCoVar(INVINFO *InvInfo, RATES* Rates, TREES* Trees, double *RateP, double *Pi);


void	CreateDEPAMatrix(INVINFO* InvInfo, RATES* Rates, TREES* Trees, double *RateP);
void	CreateDEPAMatrixCoVar(INVINFO *InvInfo, RATES* Rates, TREES* Trees, double *RateP);

void	CreateInDEPAMatrix(INVINFO* InvInfo, RATES* Rates, TREES* Trees, double *RateP);
void	CreateInDEPAMatrixCoVar(INVINFO *InvInfo, RATES* Rates, TREES* Trees, double *RateP);

void	CreateDepCVAMatrix(INVINFO *InvInfo, RATES* Rates, TREES* Trees, double *R);


double	CreatFullPMatrix(double t, INVINFO	*InvInfo, MATRIX *Mat, int NOS, int ThrNo);
double	Create4SPMat(double t, INVINFO *InvInfo, MATRIX *Mat, int ThrNo);
double	Create2SPMat(double t, INVINFO *InvInfo, MATRIX *Mat, int ThrNo);

int		InvMat(INVINFO	*InvInfo, int NoStates);

#endif
