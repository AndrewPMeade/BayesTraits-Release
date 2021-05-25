#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

#include "typedef.h"
#include "genlib.h"
#include "linalg.h"

int	DB = FALSE;

double	CreatFullPMatrix(double t, INVINFO	*InvInfo, MATRIX *Mat, TREES* Trees, MATRIX *A, double *Et);

int	PreCalc(INVINFO	*InvInfo, TREES *Trees, RATES *Rates)
{
	int				Ret;
	int				NoOfStates;
	int				*iwork;
	double			*work;
	double			*vi;
	
	
	NoOfStates = Trees->NoOfStates;

	iwork = (int*)InvInfo->TempVect2;
	work = InvInfo->TempVect3;
	vi = InvInfo->TempVect4;

	CopyMatrix(InvInfo->Q, InvInfo->A);

	if(DB == TRUE)
	{
		PrintMatrix(InvInfo->A, "A Matrix", stdout);		
	}

//	PrintMatrix(InvInfo->A, "A = ", stdout);exit(0);

	Ret = EigenRealGeneral(NoOfStates, InvInfo->A->me, InvInfo->val, vi, InvInfo->vec->me, iwork, work);
	
	if(DB == TRUE)
		PrintMatrix(InvInfo->vec, "egi Vec Matrix", stdout);

	if(Ret != NO_ERROR)
	{
	/*	for(Ret=0;Ret<Rates->NoOfRates;Ret++)
			printf("%d\t%f\n", Ret, Rates->Rates[Ret]);
		printf("\n\n\n");
		PrintMatrix(InvInfo->Q, "QMat", stdout);
		PrintMathematicaMatrix(InvInfo->Q, "Q", stdout);
		exit(0); 
		*/
		// TODO Phoneim remove 
		return ERROR;
	}

	CopyMatrix(InvInfo->Q, InvInfo->vec);
	Ret = InvertMatrix(InvInfo->Q->me, NoOfStates, work, iwork, InvInfo->inv_vec->me);
	if(Ret != NO_ERROR)
		PrintMatrix(InvInfo->inv_vec, "inv vec", stdout);

	if(Ret != NO_ERROR)
	{
		printf("%s::%d Inver Err\n", __FILE__, __LINE__);
		exit(0);
		return ERROR;
	}

	if(DB == TRUE)
	{
		CreatFullPMatrix(1, InvInfo, Trees->PList[0], Trees, Trees->InvInfo->A, Trees->InvInfo->TempVect1);
		PrintMatrix(Trees->PList[0], "My P", stdout);
	}

	return NO_ERROR;
}


int	CreateMSAMatrix(INVINFO *InvInfo, RATES* Rates, TREES* Trees)
{
	int		Outter,Inner, RPos;
	double	Tot;
	MATRIX *A;
	
	A = InvInfo->A;

	RPos = 0;
	for(Outter=0;Outter<Trees->NoOfStates;Outter++)
	{
		Tot = 0;
		for(Inner=0;Inner<Trees->NoOfStates;Inner++)
		{
			if(Inner != Outter)
			{
				A->me[Outter][Inner] = Rates->FullRates[RPos] * Rates->Pis[Inner];
				Tot += A->me[Outter][Inner];
				RPos++;
			}
		}

		A->me[Outter][Outter] = -Tot;
	}

	return PreCalc(InvInfo, Trees, Rates);
}

void	SetUpCoVarMatrix(MATRIX *A, RATES* Rates, TREES* Trees)
{
	int		Outter;
	int		Inner;
	int		NOS;
	int		RPos=0;

	for(Outter=0;Outter<Trees->NoOfStates;Outter++)
		for(Inner=0;Inner<Trees->NoOfStates;Inner++)
			A->me[Outter][Inner] = 0;

	NOS = Trees->NoOfStates / 2;

	for(Outter=0;Outter<NOS;Outter++)
		for(Inner=0;Inner<NOS;Inner++)
			A->me[Outter][Inner] = 0;

	Outter = 0;
	Inner = NOS;
	for(Outter=0;Outter<NOS;Outter++, Inner++)
		A->me[Outter][Inner] = Rates->OffToOn;

	Outter = NOS;
	Inner = 0;
	for(Outter=NOS;Outter<Trees->NoOfStates;Outter++, Inner++)
		A->me[Outter][Inner] = Rates->OnToOff;
}

int	CreateMSAMatrixCoVar(INVINFO *InvInfo, RATES* Rates, TREES* Trees)
{
	int		Outter,Inner, NOS, RPos;
	double	Tot;
	MATRIX *A;

	A = InvInfo->A;
	RPos = 0;
	NOS = Trees->NoOfStates / 2;

	SetUpCoVarMatrix(A, Rates, Trees);

	for(Outter=NOS;Outter<Trees->NoOfStates;Outter++)
	{
		Tot = 0;
		for(Inner=NOS;Inner<Trees->NoOfStates;Inner++)
		{
			if(Inner != Outter)
			{
				A->me[Outter][Inner] = Rates->FullRates[RPos] * Rates->Pis[Inner];
				Tot += A->me[Outter][Inner];
				RPos++;
			}
			else
				A->me[Outter][Inner] = 0;
		}
	}

	/*  Set man diagonal to -row */
	for(Outter=0;Outter<Trees->NoOfStates;Outter++)
	{
		Tot = 0;
		for(Inner=0;Inner<Trees->NoOfStates;Inner++)
			Tot += A->me[Outter][Inner];
		A->me[Outter][Outter] = -Tot;
	}

	return PreCalc(InvInfo, Trees, Rates);
}

int	CreateDEPAMatrixCoVar(INVINFO *InvInfo, RATES* Rates, TREES* Trees)
{
	int		Inner,Outter;
	double	Tot;
	MATRIX *A;

	A = InvInfo->A;

	SetUpCoVarMatrix(A, Rates, Trees);

/*	A->me[4][4] = -(Rates->FullRates[0] + Rates->FullRates[1]); */
	A->me[4][5] = Rates->FullRates[0];
	A->me[4][6] = Rates->FullRates[1];
	A->me[4][7] = 0;

	A->me[5][4] = Rates->FullRates[2];
/*	A->me[5][5] = -(Rates->FullRates[2] + Rates->FullRates[3]); */
	A->me[5][6] = 0;
	A->me[5][7] = Rates->FullRates[3];

	A->me[6][4] = Rates->FullRates[4];
	A->me[6][5] = 0;
/*	A->me[6][6] = -(Rates->FullRates[4] + Rates->FullRates[5]); */
	A->me[6][7] = Rates->FullRates[5];

	A->me[7][4] = 0;
	A->me[7][5] = Rates->FullRates[6];
	A->me[7][6] = Rates->FullRates[7];
/*	A->me[7][7] = -(Rates->FullRates[6] + Rates->FullRates[7]); */

	for(Outter=0;Outter<Trees->NoOfStates;Outter++)
	{
		Tot = 0;
		for(Inner=0;Inner<Trees->NoOfStates;Inner++)
			Tot += A->me[Outter][Inner];
		A->me[Outter][Outter] = -Tot;
	}

	return PreCalc(InvInfo, Trees, Rates);

}

int	CreateInDEPAMatrixCoVar(INVINFO *InvInfo, RATES* Rates, TREES* Trees)
{
	int		Inner,Outter;
	double	Alpha1, Beta1, Alpha2, Beta2;
	double	Tot;
	MATRIX *A;

	A = InvInfo->A;

	SetUpCoVarMatrix(A, Rates, Trees);

	Alpha1	= Rates->FullRates[0];
	Beta1	= Rates->FullRates[1];
	Alpha2	= Rates->FullRates[2];
	Beta2	= Rates->FullRates[3];

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

	for(Outter=0;Outter<Trees->NoOfStates;Outter++)
	{
		Tot = 0;
		for(Inner=0;Inner<Trees->NoOfStates;Inner++)
			Tot += A->me[Outter][Inner];
		A->me[Outter][Outter] = -Tot;
	}

	return PreCalc(InvInfo, Trees, Rates);

}

int	CreateInDEPAMatrix(INVINFO* InvInfo, double *R, RATES* Rates, TREES* Trees)
{
	double	Alpha1, Beta1, Alpha2, Beta2;
	MATRIX *A;
	
	A = InvInfo->A;

	Alpha1	= R[0];
	Beta1	= R[1];
	Alpha2	= R[2];
	Beta2	= R[3];

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

	return PreCalc(InvInfo, Trees, Rates);
}


int	CreateDEPAMatrix(INVINFO* InvInfo, double *R, RATES* Rates, TREES* Trees)
{
	MATRIX *A;
	
	A = InvInfo->A;

	A->me[0][0] = -(R[0] + R[1]);
	A->me[0][1] = R[0];
	A->me[0][2] = R[1];
	A->me[0][3] = 0;

	A->me[1][0] = R[2];
	A->me[1][1] = -(R[2] + R[3]);
	A->me[1][2] = 0;
	A->me[1][3] = R[3];

	A->me[2][0] = R[4];
	A->me[2][1] = 0;
	A->me[2][2] = -(R[4] + R[5]);
	A->me[2][3] = R[5];

	A->me[3][0] = 0;
	A->me[3][1] = R[6];
	A->me[3][2] = R[7];
	A->me[3][3] = -(R[6] + R[7]);


	return PreCalc(InvInfo, Trees, Rates);
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

int	CreateDepCVAMatrix(INVINFO *InvInfo, double *R, RATES* Rates, TREES* Trees)
{
	double Alpha1, Beta1, Alpha2, Beta2;
	double q12,q13,q21,q24,q31,q34,q42,q43;
	double 	qDI00, qDI01, qDI10, qDI11;
	double 	qID00, qID01, qID10, qID11;
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
	
	qDI00	= R[12];
	qDI01	= R[13];
	qDI10	= R[14];
	qDI11	= R[15];

	qID00	= R[16];
	qID01	= R[17];
	qID10	= R[18];
	qID11	= R[19];
	
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

	return PreCalc(InvInfo, Trees, Rates);
}

double	Create2SPMat(double t, INVINFO *InvInfo, MATRIX *Mat, TREES* Trees, MATRIX *A, double *Et)
{
	int		NOS;
	double  t1, t2;
	double	*Val;
	double	**Vec, **InvVec, **M, **Am;
	
	NOS		= Trees->NoOfStates;
	Val		= InvInfo->val;
	Vec		= InvInfo->vec->me;
	InvVec	= InvInfo->inv_vec->me;
	M		= Mat->me;
	Am		= A->me;


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

	return t1;
}

double	Create4SPMat(double t, INVINFO *InvInfo, MATRIX *Mat, TREES* Trees, MATRIX *A, double *Et)
{
	int		NOS, i;
	double  t1, t2;
	double	*Val;
	double	**Vec, **InvVec, **M, **Am;
	
	NOS		= Trees->NoOfStates;
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
	
	NOS		= Trees->NoOfStates;
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

double	CreatFullPMatrix(double t, INVINFO	*InvInfo, MATRIX *Mat, TREES* Trees, MATRIX *A, double *Et)
{
	int		i, j, k;
	double	t1, t2;
	int		NOS;
	double	*Val;
	double	**Vec, **InvVec, **M, **Am;

	Val		= InvInfo->val;
	Vec		= InvInfo->vec->me;
	InvVec	= InvInfo->inv_vec->me;
	M		= Mat->me;
	Am		= A->me;
	
	NOS		= Trees->NoOfStates;
	
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
				return 1000;
		}
		t1 += t2 * t2;
	}

	return t1;
} 
/*
double	CreatFullPMatrix(double t, INVINFO	*InvInfo, MATRIX *Mat, TREES* Trees, MATRIX *A, double *Et)
{
	int		i, j, k;
	double	t1, t2;
	int		NOS;
	
	NOS		= Trees->NoOfStates;
	
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
