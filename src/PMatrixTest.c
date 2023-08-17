/*
*  BayesTriats 4.0
*
*  copyright 2022
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

#include "TypeDef.h"
#include "GenLib.h"
#include "Matrix.h"
#include "LinAlg.h"
#include "Rates.h"

extern INVINFO*	AllocInvInfo(int NOS);
extern void FreeInvInfo(INVINFO* InvInfo);

double*	GetRates(void)
{
	double*	Ret;
	int		Index;
	
	Ret = (double*)SMalloc(sizeof(double) * 8);

	for(Index=0;Index<8;Index++)
	{
		do
		{
			printf("%s:\t", DEPPRAMS[Index]);
			scanf("%lf", &Ret[Index]);
		} while(Ret[Index] < 0);
	}

	for(Index=0;Index<8;Index++)
		printf("%s\t%f\n", DEPPRAMS[Index], Ret[Index]);

	return Ret;
}

void	LocPreCalc(INVINFO	*InvInfo)
{
	int		Ret;
	int		NoStates;
	
	static int		*iwork=NULL;
	static double	*work=NULL;
	static double	*vi=NULL;

	NoStates = 4;

	if(iwork == NULL)
	{
		iwork = (int*)malloc(sizeof(int) * NoStates);
		work = (double*)malloc(sizeof(double) * NoStates);
		vi = (double*)malloc(sizeof(double) * NoStates);
	}

	CopyMatrix(InvInfo->Q, InvInfo->A);

	Ret = EigenRealGeneral(NoStates, InvInfo->A->me, InvInfo->val, vi, InvInfo->vec->me, iwork, work);

	CopyMatrix(InvInfo->Q, InvInfo->vec);
	Ret = InvertMatrix(InvInfo->Q->me, NoStates, work, iwork, InvInfo->inv_vec->me);

	if(Ret == ERROR)
	{
		printf("%s::%d Inver Err\n", __FILE__, __LINE__);
		exit(0);
	}
}

void CreateLocDEPAMatrix(MATRIX *A, double* Rates)
{
	A->me[0][0] = -(Rates[0] + Rates[1]);
	A->me[0][1] = Rates[0];
	A->me[0][2] = Rates[1];
	A->me[0][3] = 0;

	A->me[1][0] = Rates[2];
	A->me[1][1] = -(Rates[2] + Rates[3]);
	A->me[1][2] = 0;
	A->me[1][3] = Rates[3];

	A->me[2][0] = Rates[4];
	A->me[2][1] = 0;
	A->me[2][2] = -(Rates[4] + Rates[5]);
	A->me[2][3] = Rates[5];

	A->me[3][0] = 0;
	A->me[3][1] = Rates[6];
	A->me[3][2] = Rates[7];
	A->me[3][3] = -(Rates[6] + Rates[7]);
}

double	LocCreatFullPMatrix(double t, MATRIX *Mat, INVINFO	*InvInfo)
{
	int		i, j, k;
	double	t1, t2;
	static  MATRIX	*A=NULL;
	static	double	*Et=NULL;
	int		NOS;

	NOS		= 4;
	
	if(A==NULL)
	{
		A = AllocMatrix(NOS, NOS);
		Et = (double*)malloc(sizeof(double)* NOS);
	}

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

void	PMatrixTest(void)
{
	MATRIX *P;
	MATRIX *A;
	double	t;
	double	*Rates;
	INVINFO	*InvInfo;
	
	P = AllocMatrix(4,4);
	A = AllocMatrix(4,4);

	Rates = GetRates();
	printf("T:\t");
	scanf("%lf", &t);
	
	CreateLocDEPAMatrix(A, Rates);
	InvInfo = AllocInvInfo(4);
	CopyMatrix(InvInfo->A, A);
	LocPreCalc(InvInfo);

	LocCreatFullPMatrix(t, P, InvInfo);

	printf("\tP(0,0)\tP(0,1)\tP(1,0)\tP(1,1)\n");
	
	printf("P(0,0)\t%f\t%f\t%f\t%f\n", P->me[0][0], P->me[0][1], P->me[0][2], P->me[0][3]);
	printf("P(0,1)\t%f\t%f\t%f\t%f\n", P->me[1][0], P->me[1][1], P->me[1][2], P->me[1][3]);
	printf("P(1,0)\t%f\t%f\t%f\t%f\n", P->me[2][0], P->me[2][1], P->me[2][2], P->me[2][3]);
	printf("P(1,1)\t%f\t%f\t%f\t%f\n", P->me[3][0], P->me[3][1], P->me[3][2], P->me[3][3]);
	
	FreeInvInfo(InvInfo);
	FreeMatrix(A);
	FreeMatrix(P);
	free(Rates);

	exit(0);
}
