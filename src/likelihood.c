#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "typedef.h"
#include "data.h"
#include "likelihood.h"
#include "matrix.h"
#include "genlib.h"
#include "linalg.h"
#include "rates.h"
#include "continuous.h"
#include "gamma.h"
#include "trees.h"
#include "praxis.h"
#include "RandLib.h"
#include "contrasts.h"
#include "threaded.h"
#include "BigLh.h"


int	DB = FALSE;

#ifndef isnan
	extern int isnan(double x);
#endif

double	CreatFullPMatrix(double t, INVINFO	*InvInfo, MATRIX *Mat, TREES* Trees, MATRIX *A, double *Et);

int		IsNum(double n)
{
	if(isnan(n) == 1)
		return FALSE;

	if(n == n + 1)
		return FALSE;

	if(n != n)
		return FALSE;

	if(n == ERRLH)
		return FALSE;

	if(n == -ERRLH)
		return FALSE;
	
	return TRUE;
}

double  AddLog(double a, double b)
{
  /************************************/
  /* addlog(log(p),log(q)) = log(p+q) */
  /************************************/
  double          x;

  x = .5 * (a - b);
  return (log(cosh(x) * 2.0) + a - x);
}


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

//	PrintMatrix(A, "A = ", stdout);
//	exit(0);

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
	double 	qDI,qID;
	int		i;
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

	A->me[4][0] = qDI;
	A->me[5][1] = qDI;
	A->me[6][2] = qDI;
	A->me[7][3] = qDI;

	A->me[0][4] = qID;
	A->me[1][5] = qID;
	A->me[2][6] = qID;
	A->me[3][7] = qID;

	SetADiag(A);

	return PreCalc(InvInfo, Trees, Rates);
}

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

INVINFO*	AllocInvInfo(int NOS)
{
	INVINFO	*Ret;
	int		Index;

	Ret = (INVINFO*)malloc(sizeof(INVINFO));
	if(Ret==NULL)
		MallocErr();

	Ret->vec		= AllocMatrix(NOS, NOS);
	Ret->inv_vec	= AllocMatrix(NOS, NOS);
	Ret->Q			= AllocMatrix(NOS, NOS);
	Ret->A			= AllocMatrix(NOS, NOS);
	Ret->TempA		= AllocMatrix(NOS, NOS);

	Ret->val = (double*)malloc(sizeof(double)*NOS);
	if(Ret->val == NULL)
		MallocErr();

	Ret->TempVect1	= (double*)malloc(sizeof(double)*NOS);
	if(Ret->TempVect1 == NULL)
		MallocErr();
	Ret->TempVect2	= (double*)malloc(sizeof(double)*NOS);
	if(Ret->TempVect2 == NULL)
		MallocErr();
	Ret->TempVect3	= (double*)malloc(sizeof(double)*NOS);
	if(Ret->TempVect3 == NULL)
		MallocErr();
	Ret->TempVect4	= (double*)malloc(sizeof(double)*NOS);
	if(Ret->TempVect4 == NULL)
		MallocErr();
	
	Ret->NoThreads = GetMaxThreads();
	Ret->Ets = (double**)malloc(sizeof(double*) * Ret->NoThreads);
	Ret->As = (MATRIX**)malloc(sizeof(MATRIX*) * Ret->NoThreads);
	if((Ret->Ets == NULL) || (Ret->As == NULL))
		MallocErr();
	for(Index=0;Index<Ret->NoThreads;Index++)
	{
		Ret->As[Index] = AllocMatrix(NOS, NOS);
		Ret->Ets[Index] = (double*)malloc(sizeof(double) * NOS );
		if(Ret->Ets[Index] == NULL)
			MallocErr();
	}



	return Ret;
}

void	FreeInvInfo(INVINFO* InvInfo)
{
	int Index;

	FreeMatrix(InvInfo->vec);
	FreeMatrix(InvInfo->inv_vec);
	FreeMatrix(InvInfo->Q);
	FreeMatrix(InvInfo->A);
	
	FreeMatrix(InvInfo->TempA);
	free(InvInfo->TempVect1);

	free(InvInfo->TempVect2);
	free(InvInfo->TempVect3);
	free(InvInfo->TempVect4);
	
	free(InvInfo->val);

	for(Index=0;Index<InvInfo->NoThreads;Index++)
	{
		FreeMatrix(InvInfo->As[Index]);
		free(InvInfo->Ets[Index]);
	}
	free(InvInfo->As);
	free(InvInfo->Ets);

	free(InvInfo);
}

void	AllocLHInfo(TREES *Trees, OPTIONS *Opt)
{
	int	NOS, Index;
	
	NOS = Trees->NoOfStates;

	Trees->InvInfo = AllocInvInfo(NOS);

	Trees->PList = (MATRIX**)malloc(sizeof(MATRIX*) * Trees->MaxNodes);
	if(Trees->PList == NULL)
		MallocErr();

	for(Index=0;Index<Trees->MaxNodes;Index++)
		Trees->PList[Index] = AllocMatrix(NOS, NOS);

}
/*
void	SumLikeDep(NODE N, TREES *Trees, int SiteNo, double Kappa, int *Err)
{
	int		Inner;
	int		Outter;
	double	Lr;
	double	Ll;
	double	Len;
	double	Val;

	if(N->Left->Tip == FALSE)
		SumLikeDep(N->Left, Trees, SiteNo, Kappa, Err);

	if(N->Right->Tip == FALSE)
		SumLikeDep(N->Right, Trees, SiteNo, Kappa, Err);

	Len = N->Left->Length;
	if(Kappa != -1)
		Len = pow(Len, Kappa);	

	Val = CreatFullPMatrix(Len, Trees->PLeft, Trees);
	if(Val > 0.001)
		(*Err) = TRUE;

	Len = N->Right->Length;
	if(Kappa != -1)
		Len = pow(Len , Kappa);

	Val = CreatFullPMatrix(Len, Trees->PRight, Trees);
	if(Val > 0.001)
		(*Err) = TRUE;


	for(Outter=0;Outter<Trees->NoOfStates;Outter++)
	{
		Ll = 0;
		Lr = 0;

		for(Inner=0;Inner<Trees->NoOfStates;Inner++)
		{
			Ll += N->Left->Partial[SiteNo][Inner] * Trees->PLeft->me[Outter][Inner];
			Lr += N->Right->Partial[SiteNo][Inner] * Trees->PRight->me[Outter][Inner];
		}

		N->Partial[SiteNo][Outter] = Ll * Lr;
	}
}
*/

void	SetFossilStates(NODE N, int SiteNo, int s00, int s01, int s10, int s11)
{
#ifdef BIG_LH
	FossilDepLhBig(N, SiteNo, s00, s01, s10, s11);
	return;
#endif

	if(s00 == 0)
		N->Partial[SiteNo][0] = 0;

	if(s01 == 0)
		N->Partial[SiteNo][1] = 0;

	if(s10 == 0)
		N->Partial[SiteNo][2] = 0;

	if(s11 == 0)
		N->Partial[SiteNo][3] = 0;
}
/*
X 	Likilhood values unchanged
-	Likilhood set to zero


Symbol	0,0	0,1	1,0	1,1
0		X	-	-	-
1		-	X	-	-
2		-	-	X	-
3		-	-	-	X
				
10		X	X	-	-
11		X	-	X	-
12		X	-	-	X
13		-	X	X	-
14		-	X	-	X
15		-	-	X	X
				
20		X	X	X	-
21		X	X	-	X
22		X	-	X	X
23		-	X	X	X
*/

void	FossilLh(NODE N, TREES *Trees, int SiteNo)
{
	int	Index;

	/* Are we using the expanded discite fossil states? */
	if(N->FossilState < Trees->NoOfStates)
	{
		#ifdef BIG_LH
			FossilLhBig(N, Trees, SiteNo);
			return;
		#endif

		for(Index=0;Index<Trees->NoOfStates;Index++)
		{
			if(Index != N->FossilState)
				N->Partial[SiteNo][Index] = 0;
		}
	}
	else
	{
		switch(N->FossilState)
		{
			case 10:
				SetFossilStates(N, SiteNo, 1, 1, 0, 0);
			break;

			case 11:
				SetFossilStates(N, SiteNo, 1, 0, 1, 0);
			break;

			case 12:
				SetFossilStates(N, SiteNo, 1, 0, 0, 1);
			break;

			case 13:
				SetFossilStates(N, SiteNo, 0, 1, 1, 0);
			break;

			case 14:
				SetFossilStates(N, SiteNo, 0, 1, 0, 1);
			break;

			case 15:
				SetFossilStates(N, SiteNo, 0, 0, 1, 1);
			break;

			case 20:
				SetFossilStates(N, SiteNo, 1, 1, 1, 0);
			break;

			case 21:
				SetFossilStates(N, SiteNo, 1, 1, 0, 1);
			break;

			case 22:
				SetFossilStates(N, SiteNo, 1, 0, 1, 1);
			break;

			case 23:
				SetFossilStates(N, SiteNo, 0, 1, 1, 1);
			break;
		}
	}

}

double	CreatFullAP(double T, double Mue, int K, MATRIX *Mat)
{
	double	Hit;
	double	Miss;
	double	temp;
	double	dK;
	int		i,j;

	dK = (double)K;
	temp = exp(-(K * Mue * T));

	Hit = (dK-1)/dK;
	Hit = Hit * temp;
	Hit = (1/dK) + Hit;

	Miss= (1/dK)*temp;
	Miss= (1/dK)-Miss;

	for(i=0;i<K;i++)
	{
		for(j=0;j<K;j++)
			Mat->me[i][j] = Miss;
		Mat->me[i][i] = Hit;
	}

	return 0;
}
/*
void	SumLikeMultiStateAP(NODE N, TREES *Trees, int SiteNo, double Kappa, int* Err, double BLMult, double Rate)
{
	int		Inner;
	int		Outter;
	double	Lr;
	double	Ll;
	double	Len;
	double	Val;
	int		NOS;

	if(N->Left->Tip == FALSE)
		SumLikeMultiStateAP(N->Left, Trees, SiteNo, Kappa, Err, BLMult, Rate);

	if(N->Right->Tip == FALSE)
		SumLikeMultiStateAP(N->Right, Trees, SiteNo, Kappa, Err, BLMult, Rate);

	if(Trees->NOSPerSite == TRUE)
		NOS = Trees->NOSList[SiteNo];
	else
		NOS = Trees->NoOfStates;

	Len = N->Left->Length;
	if(Kappa != -1)
		Len = pow(Len, Kappa);
	Len = Len * BLMult;
	Val = CreatFullAP(Len, Rate, NOS, Trees->PLeft);

	Len = N->Right->Length;
	if(Kappa != -1)
		Len = pow(Len, Kappa);
	Len = Len * BLMult;
	Val = CreatFullAP(Len, Rate, NOS, Trees->PRight);
	
	for(Outter=0;Outter<NOS;Outter++)
	{
		Ll = 0;
		Lr = 0;

		for(Inner=0;Inner<NOS;Inner++)
		{
			Ll += N->Left->Partial[SiteNo][Inner] * Trees->PLeft->me[Outter][Inner];
			Lr += N->Right->Partial[SiteNo][Inner] * Trees->PRight->me[Outter][Inner];
		}

		N->Partial[SiteNo][Outter] = Ll * Lr;
	}

	if(N->FossilState != -1)
		FossilLh(N, Trees, SiteNo);
}

void	SumLikeMultiState(NODE N, TREES *Trees, int SiteNo, double Kappa, int* Err, double BLMult)
{
	int		Inner, Outter, NIndex;
	double	Lr, Ll, Lh;
	double	Len;
	double	Val;

	if(N->Left->Tip == FALSE)
		SumLikeMultiState(N->Left, Trees, SiteNo, Kappa, Err, BLMult);

	if(N->Right->Tip == FALSE)
		SumLikeMultiState(N->Right, Trees, SiteNo, Kappa, Err, BLMult);

	Len = N->Left->Length;
	if(Kappa != -1)
		Len = pow(Len, Kappa);
	
	Len = Len * BLMult;

	Val = CreatFullPMatrix(Len, Trees->PLeft, Trees);
	if(Val > 0.001)
	{
		(*Err) = TRUE;
		return;
	}

	Len = N->Right->Length;
	if(Kappa != -1)
		Len = pow(Len, Kappa);

	Len = Len * BLMult;
	Val = CreatFullPMatrix(Len, Trees->PRight, Trees);
	if(Val > 0.001)
	{
		(*Err) = TRUE;
		return;
	}

	for(Outter=0;Outter<Trees->NoOfStates;Outter++)
	{
		Ll = 0;
		Lr = 0;

		for(Inner=0;Inner<Trees->NoOfStates;Inner++)
		{
			Ll += N->Left->Partial[SiteNo][Inner] * Trees->PLeft->me[Outter][Inner];
			Lr += N->Right->Partial[SiteNo][Inner] * Trees->PRight->me[Outter][Inner];
		}

		N->Partial[SiteNo][Outter] = Ll * Lr;
	}

	if(N->FossilState != -1)
		FossilLh(N, Trees, SiteNo);
}
*/

void	SumLikeMultiState(NODE N, TREES *Trees, int SiteNo,  int Rec)
{
	int		Inner, Outter, NIndex;
	double	Lh;
	double	**Mat, **Partial;

	if(Rec == TRUE)
	{
		for(NIndex=0;NIndex<N->NoNodes;NIndex++)
		{
			if(N->NodeList[NIndex]->Tip == FALSE)
				SumLikeMultiState(N->NodeList[NIndex], Trees, SiteNo, Rec);
		}
	}

	for(Outter=0;Outter<Trees->NoOfStates;Outter++)
	{
		N->Partial[SiteNo][Outter] = 1;

		for(NIndex=0;NIndex<N->NoNodes;NIndex++)
		{
			Mat = Trees->PList[N->NodeList[NIndex]->ID]->me;
			Partial = N->NodeList[NIndex]->Partial;

			Lh = 0;
			for(Inner=0;Inner<Trees->NoOfStates;Inner++)
				Lh += Partial[SiteNo][Inner] * Mat[Outter][Inner];
			
			N->Partial[SiteNo][Outter] *= Lh;
		}
	}

	if(N->FossilState != -1)
		FossilLh(N, Trees, SiteNo);
}


void	SumLikeRModel(NODE N, TREES *Trees, int SiteNo, RATES *Rates)
{
	return;

/* Need to fix for polytomes */
/*
	int		Inner;
	int		Outter;
	double	Lr;
	double	Ll;
	double	LP[2];
	double	RP[2];
	double	PiT;


	if(N->Left->Tip == FALSE)
		SumLikeRModel(N->Left, Trees, SiteNo, Rates);

	if(N->Right->Tip == FALSE)
		SumLikeRModel(N->Right, Trees, SiteNo, Rates);

	for(Outter=0;Outter<Trees->NoOfStates;Outter++)
	{
		Ll = 0;
		Lr = 0;
		PiT = 0;

		for(Inner=0;Inner<Trees->NoOfStates;Inner++)
		{
			if(Inner != Outter)
			{
				Ll += N->Left->Partial[SiteNo][Inner];
				Lr += N->Right->Partial[SiteNo][Inner];
				PiT += Rates->Pis[Inner];
			}
		}
	
		Ll = Ll / (Trees->NoOfStates-1);
		Lr = Lr / (Trees->NoOfStates-1);

		LP[0] = exp(-(Rates->FullRates[0]*N->Left->Length*PiT));
		LP[1] = 1 - LP[0];

		RP[0] = exp(-(Rates->FullRates[0]*N->Right->Length*PiT));
		RP[1] = 1 - RP[0];

		Ll = LP[1] * Ll;
		Lr = RP[1] * Lr;

		Ll += LP[0] * N->Left->Partial[SiteNo][Outter]; 
		Lr += RP[0] * N->Right->Partial[SiteNo][Outter];

		N->Partial[SiteNo][Outter] = Ll * Lr;
	}

	*/
}

void	CreatIndepPMatrix(double t, MATRIX *Mat, double Alpha, double Beta)
{
	double	Temp;
	double	Body;

	Body = -((Alpha+Beta)*t);
	Body = 1 - exp(Body);

	Temp = Alpha / (Alpha + Beta);
	Mat->me[0][1] = Temp * Body;
	Mat->me[0][0] = 1 - Mat->me[0][1];

	Temp = Beta / (Alpha + Beta);
	Mat->me[1][0] = Temp * Body;
	Mat->me[1][1] = 1 - Mat->me[1][0];

}
/*
void	SumLikeInDep(NODE N, TREES *Trees, double *Rates, int SiteNo)
{
	int		Inner;
	int		Outter;
	double	Lr;
	double	Ll;


	if(N->Left->Tip == FALSE)
		SumLikeInDep(N->Left, Trees, Rates, SiteNo);

	if(N->Right->Tip == FALSE)
		SumLikeInDep(N->Right, Trees, Rates, SiteNo);

	
		CreatIndepPMatrix(N->Left->Length, Trees->PLeft, Rates[0], Rates[1]);
	CreatIndepPMatrix(N->Right->Length, Trees->PRight, Rates[0], Rates[1]);

	for(Outter=0;Outter<Trees->NoOfStates;Outter++)
	{
		Ll = 0;
		Lr = 0;

		for(Inner=0;Inner<Trees->NoOfStates;Inner++)
		{
			Ll += N->Left->Partial[SiteNo][Inner] * Trees->PLeft->me[Outter][Inner];
			Lr += N->Right->Partial[SiteNo][Inner] * Trees->PRight->me[Outter][Inner];
		}

		N->Partial[SiteNo][Outter] = Ll * Lr;
	}
}
*/
/*
void	SumLikeInDep(NODE N, TREES *Trees, double *Rates, int SiteNo)
{
	int		Inner, Outter, Index;
	double	Lh;

	for(Index=0;Index<N->NoNodes;Index++)
		if(N->NodeList[Index]->Tip == FALSE)
			SumLikeInDep(N->NodeList[Index], Trees, Rates, SiteNo);

	for(Index=0;Index<N->NoNodes;Index++)
	{
		CreatIndepPMatrix(N->NodeList[Index]->Length, Trees->PList[Index], Rates[0], Rates[1]);
	}


	for(Outter=0;Outter<Trees->NoOfStates;Outter++)
	{
		N->Partial[SiteNo][Outter] = 1;

		for(Index=0;Index<N->NoNodes;Index++)
		{
			Lh = 0;
			for(Inner=0;Inner<Trees->NoOfStates;Inner++)
				Lh += N->NodeList[Index]->Partial[SiteNo][Inner] * Trees->PList[Index]->me[Outter][Inner];
			N->Partial[SiteNo][Outter] *= Lh;
		}
	}
}
*/

/* Only for first site */

void	PrintTipData(TREES* Trees, int TreeNo)
{
	int		NIndex;
	TREE	*Tree=NULL;
	NODE	N;

	Tree = Trees->Tree[TreeNo];

	for(NIndex=0;NIndex<Tree->NoNodes;NIndex++)
	{
		N = Tree->NodeList[NIndex];
		if(N->Tip==TRUE)
		{
			printf("%d\t", N->TipID);

			printf("%f\t", N->Partial[0][0]);
			printf("%f\t", N->Partial[0][1]);
			printf("%f\t", N->Partial[0][2]);
			printf("%f\t", N->Partial[0][3]);

			printf("\n");
		}
	}
}

int		SetUpAMatrix(RATES* Rates, TREES *Trees, OPTIONS *Opt)
{
	int	Err;
	HETERO *Hetero;

 	if(Opt->UseRModel == TRUE)
		return NO_ERROR;

	if(Opt->Model == MULTISTATE)
	{
		if(Trees->UseCovarion == FALSE)
			Err = CreateMSAMatrix(Trees->InvInfo, Rates, Trees);
		else
			Err = CreateMSAMatrixCoVar(Trees->InvInfo, Rates, Trees);
	}

	if(Opt->Model == DESCDEP)
	{
		if(Trees->UseCovarion == FALSE)
			Err = CreateDEPAMatrix(Trees->InvInfo, Rates->FullRates, Rates, Trees);
		else
			Err = CreateDEPAMatrixCoVar(Trees->InvInfo, Rates, Trees);
	}

	if(Opt->Model == DESCINDEP)
	{
		if(Trees->UseCovarion == FALSE)
			Err = CreateInDEPAMatrix(Trees->InvInfo, Rates->FullRates, Rates, Trees);
		else
			Err = CreateInDEPAMatrixCoVar(Trees->InvInfo, Rates, Trees);
	}

	if(Opt->Model == DESCCV)
	{
		Err = CreateDepCVAMatrix(Trees->InvInfo, Rates->FullRates, Rates, Trees);
	}

	if(Opt->Model == DESCHET)
	{
		Hetero = Rates->Hetero;
		if(Trees->UseCovarion == FALSE)
		{
			Err = CreateInDEPAMatrix(Hetero->ModelInv[0], Rates->FullRates, Rates, Trees);

			
			Err += CreateDEPAMatrix(Hetero->ModelInv[1], &Rates->FullRates[4], Rates, Trees);
		}
		else
		{
			printf("CV not supported\n");
			exit(0);
		}
	}

	if(Err > 1)
		Err = 1;
	return Err;
}

void	SetGammaBlank(RATES* Rates, OPTIONS* Opt)
{
	NODE	N;
	TREE	*Tree;
	TREES	*Trees;
	int		NIndex;
	int		SiteIndex;
	int		SIndex;
	int		NOS;

	Trees	= Opt->Trees;
	Tree	= Trees->Tree[Rates->TreeNo];
	
	for(NIndex=0;NIndex<Tree->NoNodes;NIndex++)
	{
		N = Tree->NodeList[NIndex];

		if(N->Tip == FALSE)
		{
			for(SiteIndex=0;SiteIndex<Trees->NoOfSites;SiteIndex++)
			{
				if(Trees->NOSPerSite == FALSE)
					NOS = Trees->NoOfStates;
				else
					NOS = Trees->NOSList[SiteIndex];
							
				for(SIndex=0;SIndex<NOS;SIndex++)
					N->GammaPartial[SiteIndex][SIndex] = 0;
			}
		}
	}
}

void	SetUpGamma(RATES* Rates, OPTIONS* Opt)
{
	static double	*RateW;
	
	SetGammaBlank(Rates, Opt);

	if(Rates->LastGamma == Rates->Gamma)
		return;

	if(RateW == NULL)
	{
		RateW = (double*)malloc(sizeof(double) * Opt->GammaCats);
		if(RateW == NULL)
			MallocErr();
	}

	DiscreteGamma(RateW, Rates->GammaMults, Rates->Gamma, Rates->Gamma, Rates->GammaCats, 0);
	Rates->LastGamma = Rates->Gamma;
}

void	ProcessGamma(RATES *Rates, TREES* Trees)
{
	int		NIndex;
	TREE	*Tree;
	NODE	N;
	double	Weight;
	int		SIndex;
	int		SiteIndex;
	int		NOS;



	Tree	= Trees->Tree[Rates->TreeNo];
	Weight	= (double)1 / Rates->GammaCats;
	NOS = Trees->NoOfStates;

	for(NIndex=0;NIndex<Tree->NoNodes;NIndex++)
	{
		N = Tree->NodeList[NIndex];

		if(N->Tip == FALSE)
		{
			for(SiteIndex=0;SiteIndex<Trees->NoOfSites;SiteIndex++)
			{
				if(Trees->NOSPerSite == TRUE)
					NOS = Trees->NOSList[SiteIndex];

				for(SIndex=0;SIndex<NOS;SIndex++)
					N->GammaPartial[SiteIndex][SIndex] += N->Partial[SiteIndex][SIndex] * Weight;
			}
		}
	}
}

void	FinishUpGamma(RATES* Rates, OPTIONS* Opt, TREES* Trees)
{
	int		NIndex;
	TREE	*Tree;
	NODE	N;
	int		SIndex;
	int		SiteIndex;
	int		NOS;

	Tree	= Trees->Tree[Rates->TreeNo];

	for(NIndex=0;NIndex<Tree->NoNodes;NIndex++)
	{
		N = Tree->NodeList[NIndex];
		
		if(N->Tip == FALSE)
		{
			for(SiteIndex=0;SiteIndex<Trees->NoOfSites;SiteIndex++)
			{
				if(Trees->NOSPerSite == FALSE)
					NOS = Trees->NoOfStates;
				else
					NOS = Trees->NOSList[SiteIndex];

				for(SIndex=0;SIndex<NOS;SIndex++)
					N->Partial[SiteIndex][SIndex] = N->GammaPartial[SiteIndex][SIndex];
			}
		}
	}
}

void	SetDiscEstDataTaxa(TAXA *Taxa, char S1, char S2)
{
	if((S1 == UNKNOWNSTATE) || (S2 == UNKNOWNSTATE))
	{
		free(Taxa->DesDataChar[0]);
		Taxa->DesDataChar[0] = SetDescUnknownStates(S1, S2);
		return;
	}

	if((S1 == '0') && (S2 == '0'))
		Taxa->DesDataChar[0][0] = '0';

	if((S1 == '0') && (S2 == '1'))
		Taxa->DesDataChar[0][0] = '1';

	if((S1 == '1') && (S2 == '0'))
		Taxa->DesDataChar[0][0] = '2';

	if((S1 == '1') && (S2 == '1'))
		Taxa->DesDataChar[0][0] = '3';
}

void SetDiscEstData(RATES* Rates, TREES *Trees, OPTIONS *Opt)
{
	TREE	*Tree;
	int		NIndex;
	int		SIndex;
	int		MDPos;
	NODE	N;
	TAXA	*Taxa;

	MDPos = 0;	
	Tree = Trees->Tree[Rates->TreeNo];

	for(NIndex=0;NIndex<Tree->NoNodes;NIndex++)
	{
		N = Tree->NodeList[NIndex];
		if(N->Tip == TRUE)
		{
			if(N->Taxa->EstData == TRUE)
			{
				Taxa = N->Taxa;
				if(Opt->Model == MULTISTATE)
				{
					for(SIndex=0;SIndex<Trees->NoOfSites;SIndex++)
					{
						if(N->Taxa->EstDataP[SIndex] == TRUE)
							Taxa->DesDataChar[SIndex][0] = Trees->SymbolList[Rates->EstDescData[MDPos++]];
					}
				}
				else
				{
					if((N->Taxa->EstDataP[0] == TRUE) && (N->Taxa->EstDataP[1] == TRUE))
						SetDiscEstDataTaxa(Taxa, '0'+Rates->EstDescData[MDPos++], '0'+Rates->EstDescData[MDPos++]);
					else
					{
						if(N->Taxa->EstDataP[0] == TRUE)
							SetDiscEstDataTaxa(Taxa, '0'+Rates->EstDescData[MDPos++], Taxa->RealData[1]);

						if(N->Taxa->EstDataP[1] == TRUE)
							SetDiscEstDataTaxa(Taxa, Taxa->RealData[0], '0'+Rates->EstDescData[MDPos++]);
					}
				}

				SetNodeTipData(Opt, N, Tree, Trees);
			}
		}
	}
}

void	LongProc(int x)
{
	int i,j;
	double	d;

	
	for(i=0;i<x;i++)
		for(j=0;j<x;j++)
		{
			d = 999;
			d = pow(d, 99);
		}

}

int		SetStdPMatrix(INVINFO *InvInfo, TREES *Trees, NODE N, MATRIX *P, double Gamma, double Kappa)
{
	double Len;
	double ErrVal;
	int		ThrNo;

	Len = N->Length;
	if(Kappa != -1)
		Len = pow(Len, Kappa);
		
	Len = Len * Gamma;

	ThrNo = GetThreadNo();
	
	ErrVal = CreatFullPMatrix(Len, InvInfo, P, Trees, InvInfo->As[ThrNo], InvInfo->Ets[ThrNo]);
	if(ErrVal > 0.001)
		return TRUE;

	return FALSE;
}

int		SetAnalyticalPMatrix(TREES *Trees, NODE N, MATRIX *P, double Rate, double Gamma, double Kappa)
{
	double Len, ErrVal;
	
	Len = N->Length;
	if(Kappa != -1)
		Len = pow(Len, Kappa);
		
	Len = Len * Gamma;

	ErrVal = CreatFullAP(Len, Rate, Trees->NOSPerSite, P);			

	return FALSE;
}


int		SetRPMatrix(TREES *Trees, NODE N, MATRIX *P, double Gamma, double Kappa)
{
	return FALSE;
	/* Needs to be fiex for R model as its not fixed form polytomie code. */
}




int		SetAllPMatrix(RATES* Rates, TREES *Trees, OPTIONS *Opt, double Gamma, double Kappa)
{
	int NIndex, Err, NoErr, PMatNo;
	TREE *Tree;
	NODE N;
	
	NoErr = 0;
	Err = FALSE;
	Tree = Trees->Tree[Rates->TreeNo];
	/* No need to compute P for root */

	
#ifdef THREADED
	#pragma omp parallel for private(N, Err) num_threads(Opt->Cores)
#endif	
	for(NIndex=1;NIndex<Tree->NoNodes;NIndex++)
	{
		N = Tree->NodeList[NIndex];
		Err = FALSE;
		if(NoErr == 0)
		{

			if(Opt->UseRModel == TRUE)
				Err = SetRPMatrix(Trees, N, Trees->PList[N->ID], Gamma, Kappa);
			else
			{
				if(Opt->AnalyticalP == FALSE)
				{
					if(Opt->Model != DESCHET)
						Err = SetStdPMatrix(Trees->InvInfo, Trees, N, Trees->PList[N->ID], Gamma, Kappa);
					else
					{
						PMatNo = Rates->Hetero->MList[NIndex];
						Err = SetStdPMatrix(Rates->Hetero->ModelInv[PMatNo], Trees, N, Trees->PList[N->ID], Gamma, Kappa);
					}
				}
				else
					Err = SetAnalyticalPMatrix(Trees, N, Trees->PList[N->ID], Rates->FullRates[0], Gamma, Kappa);
			}
		}

		if(Err == TRUE)
			NoErr++;
	}

	if(NoErr > 0)
		return TRUE;
	
	return FALSE;
}

void	RunNodeGroup(int GroupNo, int Parallel, RATES* Rates, TREE *Tree, TREES *Trees, OPTIONS *Opt, int SiteNo)
{
	int NIndex;

	if(Parallel == FALSE)
	{
		for(NIndex=0;NIndex<Tree->NoPNodes[GroupNo];NIndex++)
		{
			#ifndef BIG_LH
				SumLikeMultiState(Tree->PNodes[GroupNo][NIndex], Trees, SiteNo, FALSE);
			#else
				LhBigLh(Tree->PNodes[GroupNo][NIndex], Trees, Opt->Precision, SiteNo);
			#endif
		}

		return;
	}

#ifdef THREADED
	#pragma omp parallel for num_threads(Opt->Cores)
#endif
	for(NIndex=0;NIndex<Tree->NoPNodes[GroupNo];NIndex++)
	{
		#ifndef BIG_LH
			SumLikeMultiState(Tree->PNodes[GroupNo][NIndex], Trees, SiteNo, FALSE);
		#else
			LhBigLh(Tree->PNodes[GroupNo][NIndex], Trees, Opt->Precision, SiteNo);
		#endif
	}
}

void	SumLhLiner(RATES* Rates, TREES *Trees, OPTIONS *Opt, int SiteNo)
{
	int	GIndex, NIndex;
	TREE	*Tree;
	
	Tree = Trees->Tree[Rates->TreeNo];
	for(GIndex=0;GIndex<Tree->NoPGroups;GIndex++)
	{
#ifndef THREADED
		RunNodeGroup(GIndex, FALSE, Rates, Tree, Trees, Opt, SiteNo);
#else
		if(Tree->NoPNodes[GIndex] >= Trees->InvInfo->NoThreads * MIN_NODES_PER_PROC)
			RunNodeGroup(GIndex, TRUE, Rates, Tree, Trees, Opt, SiteNo);
		else
			RunNodeGroup(GIndex, FALSE, Rates, Tree, Trees, Opt, SiteNo);
#endif
	}
}

double	CombineLh(RATES* Rates, TREES *Trees, OPTIONS *Opt)
{
	int SiteNo, Index, NOS;
	double Sum, SiteLh, Ret;
	TREE	*Tree;

	Tree = Trees->Tree[Rates->TreeNo];

	Ret = 0;
	for(SiteNo=0;SiteNo<Trees->NoOfSites;SiteNo++)
	{
		if(Trees->NOSPerSite == TRUE)
		{
			NOS = Trees->NOSList[SiteNo];
			for(Index=0;Index<NOS;Index++)
				Rates->Pis[Index]  = (double)1/NOS;
		}
		else
			NOS = Trees->NoOfStates;

#ifndef BIG_LH
		Sum = 0;
		for(Index=0;Index<NOS;Index++)
			Sum += Tree->Root->Partial[SiteNo][Index] * Rates->Pis[Index];

		SiteLh = 0;
		for(Index=0;Index<NOS;Index++)
		{
			if(Tree->Root->Partial[SiteNo][Index] / Sum < 0)
				return ERRLH;

			SiteLh += Tree->Root->Partial[SiteNo][Index] * Rates->Pis[Index];
		}

		if(IsNum(log(SiteLh)) == FALSE)
			return ERRLH;

		Ret += log(SiteLh);
#else
		Ret += CombineBigLh(Rates, Trees, Opt, SiteNo, NOS);
#endif
	}

	return Ret;
}


double	Likelihood(RATES* Rates, TREES *Trees, OPTIONS *Opt)
{
	double	Ret;
	int		SiteNo;
	TREE	*Tree;
	int		Err=NO_ERROR;
	int		GammaCat;
	double	RateMult;
//	int		Index, NOS;
//	double	SiteLh, Sum;


	MapRates(Rates, Opt);

	if(Opt->Model == CONTRASTM)
		return CalcContrastLh(Opt, Trees, Rates);
	
	if(Opt->DataType == CONTINUOUS)
		return LHRandWalk(Opt, Trees, Rates);

	Tree = Trees->Tree[Rates->TreeNo];

	if(Rates->UseEstData == TRUE)
		SetDiscEstData(Rates, Trees, Opt);

	if(Opt->AnalyticalP == FALSE)
	{
		Err = SetUpAMatrix(Rates, Trees, Opt);

		if(Err == ERROR)
			return ERRLH;
	}

	if(Opt->UseGamma == TRUE)
		SetUpGamma(Rates, Opt);

	for(GammaCat=0;GammaCat<Rates->GammaCats;GammaCat++)
	{

		if(Opt->UseGamma == FALSE)
			RateMult = 1;
		else
			RateMult = Rates->GammaMults[GammaCat];

		Err = SetAllPMatrix(Rates, Trees, Opt, RateMult, Rates->Kappa);
		if(Err == TRUE)
			return ERRLH;

		for(SiteNo=0;SiteNo<Trees->NoOfSites;SiteNo++)
		{
			SumLhLiner(Rates, Trees, Opt, SiteNo);
//			Err = SumLh(Rates, Trees, Opt, SiteNo);

			if(Err == TRUE)
				return ERRLH;
		}

		if(Opt->UseGamma == TRUE)
			ProcessGamma(Rates, Trees);
	}

	if(Opt->UseGamma == TRUE)
		FinishUpGamma(Rates, Opt, Trees);

	Ret = CombineLh(Rates, Trees, Opt);

	if((Ret > 0) || (Ret == Ret+1) || (Ret != Ret) || (Ret == ERRLH))
		return ERRLH;

	return Ret;
}

double	LhPraxisCon(PRAXSTATE* PState, double *List)
{
	int	Index=0;

	if(PState->Opt->EstKappa == TRUE)
	{
		PState->Rates->Kappa = List[Index];
		if(PState->Rates->Kappa > MAX_KAPPA)
			PState->Rates->Kappa = MAX_KAPPA;
		if(PState->Rates->Kappa < MIN_KAPPA)
			PState->Rates->Kappa = MIN_KAPPA;

		Index++;
	}

	if(PState->Opt->EstDelta == TRUE)
	{
		PState->Rates->Delta = List[Index];
		if(PState->Rates->Delta > MAX_DELTA)
			PState->Rates->Delta = MAX_DELTA;
		if(PState->Rates->Delta < MIN_DELTA)
			PState->Rates->Delta = MIN_DELTA;

		Index++;
	}

	if(PState->Opt->EstLambda == TRUE)
	{
		PState->Rates->Lambda = List[Index];
		if(PState->Rates->Lambda > MAX_LAMBDA)
			PState->Rates->Lambda = MAX_LAMBDA;
		if(PState->Rates->Lambda < MIN_LAMBDA)
			PState->Rates->Lambda = MIN_LAMBDA;
		
		Index++;
	}

	return LHRandWalk(PState->Opt, PState->Trees, PState->Rates);
}

double	LhPraxis(void* P, double *List)
{
//	static	Called;
	int		Index;
	double		Ret;
	double		*TRates;
	PRAXSTATE	*PState;

	PState = (PRAXSTATE*)P;

	if(PState->Opt->DataType == CONTINUOUS)
		return -LhPraxisCon(PState, List);

	for(Index=0;Index<PState->Rates->NoOfRates;Index++)
		if((List[Index] < MINRATE) || (IsNum(List[Index]) == FALSE))
			List[Index] = MINRATE;

	TRates = PState->Rates->Rates;
	PState->Rates->Rates = List;

	Ret = Likelihood(PState->Rates, PState->Trees, PState->Opt);

	PState->Rates->Rates = TRates;

	return -Ret;
}
