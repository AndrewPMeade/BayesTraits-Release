/*
*  BayesTriats 3.0
*
*  copyright 2017
*
*  Andrew Meade
*  School of Biological Sciences
*  University of Reading
*  Reading
*  Berkshire
*  RG6 6BX
*
* BayesTriats is free software: you can redistribute it and/or modify
* it under the terms of the GNU General Public License as published by
* the Free Software Foundation, either version 3 of the License, or
* (at your option) any later version.
* 
* This program is distributed in the hope that it will be useful,
* but WITHOUT ANY WARRANTY; without even the implied warranty of
* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
* GNU General Public License for more details.
* 
* You should have received a copy of the GNU General Public License
* along with this program.  If not, see <http://www.gnu.org/licenses/>
*
*/



#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

#include "TypeDef.h"
#include "GenLib.h"
#include "LinAlg.h"

int	InvMat(INVINFO	*InvInfo, int NoStates)
{
	int				Ret;
	int				*iwork;
	double			*work;
	double			*vi;
	

	iwork = (int*)InvInfo->TempVect2;
	work = InvInfo->TempVect3;
	vi = InvInfo->TempVect4;

	CopyMatrix(InvInfo->Q, InvInfo->A);

//	PrintMatrix(InvInfo->A, "A = ", stdout);exit(0);

	Ret = EigenRealGeneral(NoStates, InvInfo->A->me, InvInfo->val, vi, InvInfo->vec->me, iwork, work);


	if(Ret != NO_ERROR)
		return ERROR;

	CopyMatrix(InvInfo->Q, InvInfo->vec);
	Ret = InvertMatrix(InvInfo->Q->me, NoStates, work, iwork, InvInfo->inv_vec->me);

	if(Ret != NO_ERROR)
	{
		PrintMatrix(InvInfo->inv_vec, "inv vec", stdout);
		printf("%s::%d Inver Err\n", __FILE__, __LINE__);
		exit(0);
		return ERROR;
	}

	return NO_ERROR;
}


void CreateMSAMatrix(INVINFO *InvInfo, int NOS, double *Rates, double *Pis)
{
	int		Outter,Inner, RPos;
	double	Tot;
	MATRIX *A;
	
	A = InvInfo->A;

	RPos = 0;
	for(Outter=0;Outter<NOS;Outter++)
	{
		Tot = 0;
		for(Inner=0;Inner<NOS;Inner++)
		{
			if(Inner != Outter)
			{
				A->me[Outter][Inner] = Rates[RPos] * Pis[Inner];
				Tot += A->me[Outter][Inner];
				RPos++;
			}
		}

		A->me[Outter][Outter] = -Tot;
	}
}

void	SetUpCoVarMatrix(MATRIX *A, RATES* Rates, TREES* Trees)
{
	int		Outter;
	int		Inner;
	int		NOS;
	int		RPos=0;

	for(Outter=0;Outter<Trees->NoStates;Outter++)
		for(Inner=0;Inner<Trees->NoStates;Inner++)
			A->me[Outter][Inner] = 0;

	NOS = Trees->NoStates / 2;

	for(Outter=0;Outter<NOS;Outter++)
		for(Inner=0;Inner<NOS;Inner++)
			A->me[Outter][Inner] = 0;

	Outter = 0;
	Inner = NOS;
	for(Outter=0;Outter<NOS;Outter++, Inner++)
		A->me[Outter][Inner] = Rates->OffToOn;

	Outter = NOS;
	Inner = 0;
	for(Outter=NOS;Outter<Trees->NoStates;Outter++, Inner++)
		A->me[Outter][Inner] = Rates->OnToOff;
}

void CreateMSAMatrixCoVar(INVINFO *InvInfo, RATES* Rates, TREES* Trees, double *RateP, double *Pi)
{
	int		Outter,Inner, NOS, RPos;
	double	Tot;
	MATRIX *A;

	A = InvInfo->A;
	RPos = 0;
	NOS = Trees->NoStates / 2;

	SetUpCoVarMatrix(A, Rates, Trees);

	for(Outter=NOS;Outter<Trees->NoStates;Outter++)
	{
		Tot = 0;
		for(Inner=NOS;Inner<Trees->NoStates;Inner++)
		{
			if(Inner != Outter)
			{
				A->me[Outter][Inner] = RateP[RPos] * Pi[Inner];
				Tot += A->me[Outter][Inner];
				RPos++;
			}
			else
				A->me[Outter][Inner] = 0;
		}
	}

	/*  Set man diagonal to -row */
	for(Outter=0;Outter<Trees->NoStates;Outter++)
	{
		Tot = 0;
		for(Inner=0;Inner<Trees->NoStates;Inner++)
			Tot += A->me[Outter][Inner];
		A->me[Outter][Outter] = -Tot;
	}
}

void	CreateDEPAMatrixCoVar(INVINFO *InvInfo, RATES* Rates, TREES* Trees, double *RateP)
{
	int		Inner,Outter;
	double	Tot;
	MATRIX *A;

	A = InvInfo->A;

	SetUpCoVarMatrix(A, Rates, Trees);

/*	A->me[4][4] = -(Rates->FullRates[0] + Rates->FullRates[1]); */
	A->me[4][5] = RateP[0];
	A->me[4][6] = RateP[1];
	A->me[4][7] = 0;

	A->me[5][4] = RateP[2];
/*	A->me[5][5] = -(Rates->FullRates[2] + Rates->FullRates[3]); */
	A->me[5][6] = 0;
	A->me[5][7] = RateP[3];

	A->me[6][4] = RateP[4];
	A->me[6][5] = 0;
/*	A->me[6][6] = -(Rates->FullRates[4] + Rates->FullRates[5]); */
	A->me[6][7] = RateP[5];

	A->me[7][4] = 0;
	A->me[7][5] = RateP[6];
	A->me[7][6] = RateP[7];
/*	A->me[7][7] = -(Rates->FullRates[6] + Rates->FullRates[7]); */

	for(Outter=0;Outter<Trees->NoStates;Outter++)
	{
		Tot = 0;
		for(Inner=0;Inner<Trees->NoStates;Inner++)
			Tot += A->me[Outter][Inner];
		A->me[Outter][Outter] = -Tot;
	}

}

void	CreateInDEPAMatrixCoVar(INVINFO *InvInfo, RATES* Rates, TREES* Trees, double *RateP)
{
	int		Inner,Outter;
	double	Alpha1, Beta1, Alpha2, Beta2;
	double	Tot;
	MATRIX *A;

	A = InvInfo->A;

	SetUpCoVarMatrix(A, Rates, Trees);

	Alpha1	= RateP[0];
	Beta1	= RateP[1];
	Alpha2	= RateP[2];
	Beta2	= RateP[3];

	A->me[4][5] = Alpha2;
	A->me[4][6] = Alpha1;
	A->me[4][7] = 0;

	A->me[5][4] = Beta2;
	A->me[5][6] = 0;
	A->me[5][7] = Alpha1;

	A->me[6][4] = Beta1;
	A->me[6][5] = 0;
	A->me[6][7] = Alpha1;

	A->me[7][4] = 0;
	A->me[7][5] = Beta1;
	A->me[7][6] = Beta2;

	for(Outter=0;Outter<Trees->NoStates;Outter++)
	{
		Tot = 0;
		for(Inner=0;Inner<Trees->NoStates;Inner++)
			Tot += A->me[Outter][Inner];
		A->me[Outter][Outter] = -Tot;
	}

}

void	CreateInDEPAMatrix(INVINFO* InvInfo, RATES* Rates, TREES* Trees, double *RateP)
{
	double	Alpha1, Beta1, Alpha2, Beta2;
	MATRIX *A;
	
	A = InvInfo->A;

	Alpha1	= RateP[0];
	Beta1	= RateP[1];
	Alpha2	= RateP[2];
	Beta2	= RateP[3];

	A->me[0][0] = -(Alpha2 + Alpha1);
	A->me[0][1] = Alpha2;
	A->me[0][2] = Alpha1;
	A->me[0][3] = 0;

	A->me[1][0] = Beta2;
	A->me[1][1] = -(Beta2 + Alpha1);
	A->me[1][2] = 0;
	A->me[1][3] = Alpha1;

	A->me[2][0] = Beta1;
	A->me[2][1] = 0;
	A->me[2][2] = -(Beta1 + Alpha2);
	A->me[2][3] = Alpha2;

	A->me[3][0] = 0;
	A->me[3][1] = Beta1;
	A->me[3][2] = Beta2;
	A->me[3][3] = -(Beta1 + Beta2);
}

/* Include Pis into the calcl*/
/*
int	CreateDEPAMatrix(INVINFO* InvInfo, double *R, RATES* Rates, TREES* Trees)
{
	double **Mat;
	double *Pi;

	Pi = Rates->Pis;
	Mat = InvInfo->A->me;
	
	Mat[0][1] = R[0] * Pi[1];
	Mat[0][2] = R[1] * Pi[2];
	Mat[0][3] = 0;
	
	Mat[1][0] = R[2] * Pi[0];
	Mat[1][2] = 0;
	Mat[1][3] = R[3] * Pi[3];

	Mat[2][0] = R[4] * Pi[0];
	Mat[2][1] = 0;
	Mat[2][3] = R[5] * Pi[3];

	Mat[3][0] = 0;
	Mat[3][1] = R[6] * Pi[1];
	Mat[3][2] = R[7] * Pi[2];

	Mat[0][0] = -(Mat[0][1] + Mat[0][2]);
	Mat[1][1] = -(Mat[1][0] + Mat[1][3]);
	Mat[2][2] = -(Mat[2][0] + Mat[2][3]);
	Mat[3][3] = -(Mat[3][1] + Mat[3][2]);

	return PreCalc(InvInfo, Trees, Rates);
}
*/


void	CreateDEPAMatrix(INVINFO* InvInfo, RATES* Rates, TREES* Trees, double *RateP)
{
	MATRIX *A;
	
	A = InvInfo->A;

	A->me[0][0] = -(RateP[0] + RateP[1]);
	A->me[0][1] = RateP[0];
	A->me[0][2] = RateP[1];
	A->me[0][3] = 0;

	A->me[1][0] = RateP[2];
	A->me[1][1] = -(RateP[2] + RateP[3]);
	A->me[1][2] = 0;
	A->me[1][3] = RateP[3];

	A->me[2][0] = RateP[4];
	A->me[2][1] = 0;
	A->me[2][2] = -(RateP[4] + RateP[5]);
	A->me[2][3] = RateP[5];

	A->me[3][0] = 0;
	A->me[3][1] = RateP[6];
	A->me[3][2] = RateP[7];
	A->me[3][3] = -(RateP[6] + RateP[7]);

}


void	SetADiag(MATRIX *A)
{
	int i,j;
	double Total;

	for(i=0;i<A->NoOfRows;i++)
	{
		Total = 0;
		A->me[i][i] = 0;
		for(j=0;j<A->NoOfCols;j++)
		{
			Total += A->me[i][j];
		}
		A->me[i][i] = -Total;
	}
}

void	CreateDepCVAMatrix(INVINFO *InvInfo, RATES* Rates, TREES* Trees, double *R)
{
	double Alpha1, Beta1, Alpha2, Beta2;
	double q12,q13,q21,q24,q31,q34,q42,q43;
	double 	qDI00, qDI01, qDI10, qDI11;
	double 	qID00, qID01, qID10, qID11;
	double	qDI, qID;
	MATRIX *A;

	A = InvInfo->A;

	FillMatrix(A, 0);

	Alpha1	= R[0];
	Beta1	= R[1];
	Alpha2	= R[2];
	Beta2	= R[3];

	q12		= R[4];
	q13		= R[5];
	q21		= R[6];
	q24		= R[7];
	q31		= R[8];
	q34		= R[9];
	q42		= R[10];
	q43		= R[11];
	
	qDI		= R[12];
	qID		= R[13];

	qDI00	= qDI;
	qDI01	= qDI;
	qDI10	= qDI;
	qDI11	= qDI;

	qID00	= qID;
	qID01	= qID;
	qID10	= qID;
	qID11	= qID;
	
	A->me[0][1] = Alpha2;
	A->me[0][2] = Alpha1;

	A->me[1][0] = Beta2;
	A->me[1][3] = Alpha1;

	A->me[2][0] = Beta1;
	A->me[2][3] = Alpha2;

	A->me[3][1] = Beta1;
	A->me[3][2] = Beta2;

	A->me[4][5] = q12;
	A->me[4][6] = q13;

	A->me[5][4] = q21;
	A->me[5][7] = q24;

	A->me[6][4] = q31;
	A->me[6][7] = q34;

	A->me[7][5] = q42;
	A->me[7][6] = q43;

	A->me[4][0] = qDI00;
	A->me[5][1] = qDI01;
	A->me[6][2] = qDI10;
	A->me[7][3] = qDI11;

	A->me[0][4] = qID00;
	A->me[1][5] = qID01;
	A->me[2][6] = qID10;
	A->me[3][7] = qID11;

	SetADiag(A);
}

double	Create2SPMat(double t, INVINFO *InvInfo, MATRIX *Mat, int ThrNo)
{
	double  t1, t2;
	double	*Val;
	double	**Vec, **InvVec, **M, **Am, *Et;
		
	Am		= InvInfo->As[ThrNo]->me;
	Et		= InvInfo->Ets[ThrNo];

	Val		= InvInfo->val;
	Vec		= InvInfo->vec->me;
	InvVec	= InvInfo->inv_vec->me;
	M		= Mat->me;
	

	Et[0] = exp(t*InvInfo->val[0]);
	Et[1] = exp(t*InvInfo->val[1]);

	Am[0][0] = Vec[0][0] * Et[0];
	Am[0][1] = Vec[0][1] * Et[1];

	Am[1][0] = Vec[1][0] * Et[0];
	Am[1][1] = Vec[1][1] * Et[1];

	t1 = 0.0;

		t2=1.0;

		M[0][0] =	Am[0][0] * InvVec[0][0] +
					Am[0][1] * InvVec[1][0];
		t2-=M[0][0];
		if(M[0][0] < 0)
			return 1000;

	
		M[0][1] =	Am[0][0] * InvVec[0][1] +
					Am[0][1] * InvVec[1][1];
		t2-=M[0][1];
		if(M[0][1] < 0)
			return 1000;
		
		t1 += t2 * t2;

		t2=1.0;

		M[1][0] =	Am[1][0] * InvVec[0][0] +
					Am[1][1] * InvVec[1][0];
		t2-=M[1][0];
		if(M[1][0] < 0)
			return 1000;

	
		M[1][1] =	Am[1][0] * InvVec[0][1] +
					Am[1][1] * InvVec[1][1];
		t2-=M[1][1];
		if(M[1][1] < 0)
			return 1000;
		
		t1 += t2 * t2;
	
//	PrintMatrix(Mat, "P=", stdout);exit(0);
	return t1;
}

double	Create4SPMat(double t, INVINFO *InvInfo, MATRIX *Mat, int ThrNo)
{
	int		i;
	double  t1, t2;
	double	*Val;
	double	**Vec, **InvVec, **M, **Am, *Et;
	
	Val		= InvInfo->val;
	Vec		= InvInfo->vec->me;
	InvVec	= InvInfo->inv_vec->me;
	M		= Mat->me;
	Am		= InvInfo->As[ThrNo]->me;
	Et		= InvInfo->Ets[ThrNo];


	Et[0] = exp(t*InvInfo->val[0]);
	Et[1] = exp(t*InvInfo->val[1]);
	Et[2] = exp(t*InvInfo->val[2]);
	Et[3] = exp(t*InvInfo->val[3]);

	Am[0][0] = Vec[0][0] * Et[0];
	Am[0][1] = Vec[0][1] * Et[1];
	Am[0][2] = Vec[0][2] * Et[2];
	Am[0][3] = Vec[0][3] * Et[3];

	Am[1][0] = Vec[1][0] * Et[0];
	Am[1][1] = Vec[1][1] * Et[1];
	Am[1][2] = Vec[1][2] * Et[2];
	Am[1][3] = Vec[1][3] * Et[3];

	Am[2][0] = Vec[2][0] * Et[0];
	Am[2][1] = Vec[2][1] * Et[1];
	Am[2][2] = Vec[2][2] * Et[2];
	Am[2][3] = Vec[2][3] * Et[3];

	Am[3][0] = Vec[3][0] * Et[0];
	Am[3][1] = Vec[3][1] * Et[1];
	Am[3][2] = Vec[3][2] * Et[2];
	Am[3][3] = Vec[3][3] * Et[3];


	t1 = 0.0;
	for(i=0;i<4;i++)
	{
		t2=1.0;

		M[i][0] =	Am[i][0] * InvVec[0][0] +
					Am[i][1] * InvVec[1][0] +
					Am[i][2] * InvVec[2][0] +
					Am[i][3] * InvVec[3][0];
		t2-=M[i][0];
		if(M[i][0] < 0)
			return 1000;

	
		M[i][1] =	Am[i][0] * InvVec[0][1] +
					Am[i][1] * InvVec[1][1] +
					Am[i][2] * InvVec[2][1] +
					Am[i][3] * InvVec[3][1];
		t2-=M[i][1];
		if(M[i][1] < 0)
			return 1000;

		M[i][2]	=	Am[i][0] * InvVec[0][2] +
					Am[i][1] * InvVec[1][2] +
					Am[i][2] * InvVec[2][2] +
					Am[i][3] * InvVec[3][2];
		t2-=M[i][2];
		if(M[i][2] < 0)
			return 1000;

		M[i][3]	=	Am[i][0] * InvVec[0][3] +
					Am[i][1] * InvVec[1][3] +
					Am[i][2] * InvVec[2][3] +
					Am[i][3] * InvVec[3][3];

		t2-=M[i][3];

		if(M[i][3] < 0)
			return 1000;
		
		t1 += t2 * t2;
	}

	return t1;
}
/*
double	CreateDiscretePMat(double t, INVINFO *InvInfo, MATRIX *Mat, TREES* Trees, MATRIX *A, double *Et)
{
	int		NOS, i;
	double  t1, t2;
	double	*Val;
	double	**Vec, **InvVec, **M, **Am;
	
	NOS		= Trees->NoStates;
	Val		= InvInfo->val;
	Vec		= InvInfo->vec->me;
	InvVec	= InvInfo->inv_vec->me;
	M		= Mat->me;
	Am		= A->me;


	Et[0] = exp(t*InvInfo->val[0]);
	Et[1] = exp(t*InvInfo->val[1]);
	Et[2] = exp(t*InvInfo->val[2]);
	Et[3] = exp(t*InvInfo->val[3]);

	Am[0][0] = Vec[0][0] * Et[0];
	Am[0][1] = Vec[0][1] * Et[1];
	Am[0][2] = Vec[0][2] * Et[2];
	Am[0][3] = Vec[0][3] * Et[3];

	Am[1][0] = Vec[1][0] * Et[0];
	Am[1][1] = Vec[1][1] * Et[1];
	Am[1][2] = Vec[1][2] * Et[2];
	Am[1][3] = Vec[1][3] * Et[3];

	Am[2][0] = Vec[2][0] * Et[0];
	Am[2][1] = Vec[2][1] * Et[1];
	Am[2][2] = Vec[2][2] * Et[2];
	Am[2][3] = Vec[2][3] * Et[3];

	Am[3][0] = Vec[3][0] * Et[0];
	Am[3][1] = Vec[3][1] * Et[1];
	Am[3][2] = Vec[3][2] * Et[2];
	Am[3][3] = Vec[3][3] * Et[3];


	t1 = 0.0;
	for(i=0;i<NOS;i++)
	{
		t2=1.0;

		M[i][0]=0.0;
		M[i][0] += Am[i][0] * InvVec[0][0];
		M[i][0] += Am[i][1] * InvVec[1][0];
		M[i][0] += Am[i][2] * InvVec[2][0];
		M[i][0] += Am[i][3] * InvVec[3][0];
			

		t2-=M[i][0];

		if(M[i][0] < 0)
			return 1000;
	
		M[i][1]=0.0;
		M[i][1] += Am[i][0] * InvVec[0][1];
		M[i][1] += Am[i][1] * InvVec[1][1];
		M[i][1] += Am[i][2] * InvVec[2][1];
		M[i][1] += Am[i][3] * InvVec[3][1];
			

		t2-=M[i][1];

		if(M[i][1] < 0)
			return 1000;

		M[i][2]=0.0;
		M[i][2] += Am[i][0] * InvVec[0][2];
		M[i][2] += Am[i][1] * InvVec[1][2];
		M[i][2] += Am[i][2] * InvVec[2][2];
		M[i][2] += Am[i][3] * InvVec[3][2];
			

		t2-=M[i][2];

		if(M[i][2] < 0)
			return 1000;

		M[i][3]=0.0;
		M[i][3] += Am[i][0] * InvVec[0][3];
		M[i][3] += Am[i][1] * InvVec[1][3];
		M[i][3] += Am[i][2] * InvVec[2][3];
		M[i][3] += Am[i][3] * InvVec[3][3];
			

		t2-=M[i][3];

		if(M[i][3] < 0)
			return 1000;
		
		t1 += t2 * t2;
	}

	return t1;
}
*/

void	SetIDMatrix(MATRIX *Mat)
{
	int i, S;

	S = Mat->NoOfCols * Mat->NoOfRows;
	for(i=0;i<S;i++)
		Mat->me[0][i] = 0;

	for(i=0;i<Mat->NoOfCols;i++)
		Mat->me[i][i] = 1.0;
}

double	CreatFullPMatrix(double t, INVINFO	*InvInfo, MATRIX *Mat, int NOS, int ThrNo)
{
	int		i, j, k;
	double	t1, t2;
	double	*Val;
	double	**Vec, **InvVec, **M, **Am, *Et;

	Am		= InvInfo->As[ThrNo]->me;
	Et		= InvInfo->Ets[ThrNo];


	Val		= InvInfo->val;
	Vec		= InvInfo->vec->me;
	InvVec	= InvInfo->inv_vec->me;
	M		= Mat->me;
	
	if(t < MIN_BL)
	{
		SetIDMatrix(Mat);
		return 0;
	}	
	
	for(i=0;i<NOS;i++)
		Et[i] = exp(t*Val[i]);

	for(i=0;i<NOS;i++)
		for(j=0;j<NOS;j++)
			Am[i][j] = Vec[i][j] * Et[j];

	t1 = 0.0;
	for(i=0;i<NOS;i++)
	{
		t2=1.0;
		for(j=0;j<NOS;j++)
		{
			M[i][j]=0.0;

			for(k=0;k<NOS;k++)
				M[i][j] += Am[i][k] * InvVec[k][j];

			t2-=M[i][j];

			if(M[i][j] < 0)
			{
				SetIDMatrix(Mat);
				return 0;
			}	
			//	return 1000;
		}
		t1 += t2 * t2;
	}

	if(t1 > 0.001)
		return t1;


	return t1;
} 
/*
double	CreatFullPMatrix(double t, INVINFO	*InvInfo, MATRIX *Mat, TREES* Trees, MATRIX *A, double *Et)
{
	int		i, j, k;
	double	t1, t2;
	int		NOS;
	
	NOS		= Trees->NoStates;
	
	for(i=0;i<NOS;i++)
		Et[i] = exp(t*InvInfo->val[i]);
	


	for(i=0;i<NOS;i++)
		for(j=0;j<NOS;j++)
			A->me[i][j] = InvInfo->vec->me[i][j] * Et[j];

	t1 = 0.0;
	for(i=0;i<NOS;i++)
	{
		t2=1.0;
		for(j=0;j<NOS;j++)
		{
			Mat->me[i][j]=0.0;
			for(k=0;k<NOS;k++)
				Mat->me[i][j] += A->me[i][k] * InvInfo->inv_vec->me[k][j];

			t2-=Mat->me[i][j];

			if(Mat->me[i][j] < 0)
				return 1000;
		}
		t1 += t2 * t2;
	}

	return t1;
}
*/
