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
#include <string.h>
#include <math.h>

#include "TypeDef.h"
#include "Data.h"
#include "Likelihood.h"
#include "Matrix.h"
#include "GenLib.h"
#include "Rates.h"
#include "Continuous.h"
#include "Gamma.h"
#include "Trees.h"
#include "Praxis.h"
#include "RandLib.h"
#include "Contrasts.h"
#include "Threaded.h"
#include "BigLh.h"
#include "PMatrix.h"
#include "LinAlg.h"
#include "ModelFile.h"
#include "QuadDouble.h"
#include "FatTail.h"
#include "TransformTree.h"
#include "VarRates.h"
#include "LocalTransform.h"
#include "DistData.h"
#include "TimeSlices.h"
#include "Landscape.h"

#ifdef BTOCL
	#include "btocl_discrete.h"
#endif

#ifdef CLIK_P
	#include <cilk/cilk.h>
	#include <cilk/cilk_api.h>
#endif

int		IsNum(double n)
{
	if(isnan(n) != 0)
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

	Ret->val = (double*)SMalloc(sizeof(double)*NOS);
	Ret->TempVect1	= (double*)SMalloc(sizeof(double)*NOS);
	Ret->TempVect2	= (double*)SMalloc(sizeof(double)*NOS);
	Ret->TempVect3	= (double*)SMalloc(sizeof(double)*NOS);
	Ret->TempVect4	= (double*)SMalloc(sizeof(double)*NOS);
	
	Ret->NoThreads = GetMaxThreads();
	Ret->Ets = (double**)SMalloc(sizeof(double*) * Ret->NoThreads);
	Ret->As = (MATRIX**)SMalloc(sizeof(MATRIX*) * Ret->NoThreads);

	for(Index=0;Index<Ret->NoThreads;Index++)
	{
		Ret->As[Index] = AllocMatrix(NOS, NOS);
		Ret->Ets[Index] = (double*)SMalloc(sizeof(double) * NOS );
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
	int	Index, NOS, NoPatterns;
	
	NOS = Trees->NoStates;
	NoPatterns = Opt->NoPatterns + 1;

	Trees->InvInfo = (INVINFO**)SMalloc(sizeof(INVINFO*) * NoPatterns);

	for(Index=0;Index<NoPatterns;Index++)
		Trees->InvInfo[Index] = AllocInvInfo(NOS);

	Trees->PList = AllocMultiMatrixLinMem(Trees->MaxNodes, NOS, NOS);
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


void	CheckBigLh(NODE N, int SiteNo, TREES *Trees)
{
	int Index, NOS, UnderFlow;


	NOS = Trees->NoStates;

	UnderFlow = FALSE;
	for(Index=0;Index<NOS;Index++)
	{
		if(N->Partial[SiteNo][Index] < LH_UNDER_FLOW)
			UnderFlow = TRUE;
	}

	if(UnderFlow == FALSE)
		return;

	N->NoUnderFlow++;
	for(Index=0;Index<NOS;Index++)
		N->Partial[SiteNo][Index] = N->Partial[SiteNo][Index] / LH_UNDER_FLOW;
}

void	SumLikeMultiState(NODE N, OPTIONS *Opt, TREES *Trees, int SiteNo)
{
	int		Inner, Outter, NIndex;
	double	Lh;
	double	**Mat, **Partial;
	double *CPart;

	CPart = N->Partial[SiteNo];

	for(Outter=0;Outter<Trees->NoStates;Outter++)
	{
		N->Partial[SiteNo][Outter] = 1;

		for(NIndex=0;NIndex<N->NoNodes;NIndex++)
		{
			Mat = Trees->PList[N->NodeList[NIndex]->ID]->me;
			Partial = N->NodeList[NIndex]->Partial;

			Lh = 0;
			for(Inner=0;Inner<Trees->NoStates;Inner++)
				Lh += Partial[SiteNo][Inner] * Mat[Outter][Inner];

			N->Partial[SiteNo][Outter] *= Lh;
		}
	}
	
	CheckBigLh(N, SiteNo, Trees);

	if(N->FossilMask != NULL)
		FossilLh(N, Opt, Trees, SiteNo);
}

void	RecSumLikeMultiState(NODE N, OPTIONS *Opt, TREES *Trees, int SiteNo)
{
	int Index;
	
	if(N->Tip == TRUE)
		return;


	for(Index=0;Index<N->NoNodes-1;Index++)
	{
//		cilk_spawn RecSumLikeMultiState(N->NodeList[Index], Opt, Trees, SiteNo);
	}

	RecSumLikeMultiState(N->NodeList[Index], Opt, Trees, SiteNo);
//	cilk_sync;

	SumLikeMultiState(N, Opt, Trees, SiteNo);
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

	for(Outter=0;Outter<Trees->NoStates;Outter++)
	{
		Ll = 0;
		Lr = 0;
		PiT = 0;

		for(Inner=0;Inner<Trees->NoStates;Inner++)
		{
			if(Inner != Outter)
			{
				Ll += N->Left->Partial[SiteNo][Inner];
				Lr += N->Right->Partial[SiteNo][Inner];
				PiT += Rates->Pis[Inner];
			}
		}
	
		Ll = Ll / (Trees->NoStates-1);
		Lr = Lr / (Trees->NoStates-1);

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


/*
int		SetUpAMatrixOld(RATES* Rates, TREES *Trees, OPTIONS *Opt)
{
	int		Err;
	HETERO *Hetero;

 	if(Opt->UseRModel == TRUE)
		return NO_ERROR;

	if(Opt->Model == M_MULTISTATE)
	{
		if(Trees->UseCovarion == FALSE)
			Err = CreateMSAMatrix(Trees->InvInfo, Rates, Trees);
		else
			Err = CreateMSAMatrixCoVar(Trees->InvInfo, Rates, Trees);
	}

	if(Opt->Model == M_DESCDEP)
	{
		if(Trees->UseCovarion == FALSE)
			Err = CreateDEPAMatrix(Trees->InvInfo, Rates->FullRates, Rates, Trees);
		else
			Err = CreateDEPAMatrixCoVar(Trees->InvInfo, Rates, Trees);
	}

	if(Opt->Model == M_DESCINDEP)
	{
		if(Trees->UseCovarion == FALSE)
			Err = CreateInDEPAMatrix(Trees->InvInfo, Rates->FullRates, Rates, Trees);
		else
			Err = CreateInDEPAMatrixCoVar(Trees->InvInfo, Rates, Trees);
	}

	if(Opt->Model == M_DESCCV)
	{
		Err = CreateDepCVAMatrix(Trees->InvInfo, Rates->FullRates, Rates, Trees);
	}

	if(Opt->Model == M_DESCHET)
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
*/

void		SetUpAMatrix(MODEL Model, RATES *Rates, TREES *Trees, int NOS, INVINFO *InvInfo, double *RateP, double *Pi)
{

	if(Model == M_MULTISTATE)
	{
		if(Trees->UseCovarion == FALSE)
			CreateMSAMatrix(InvInfo, NOS, RateP, Pi);
		else
			CreateMSAMatrixCoVar(InvInfo, Rates, Trees, RateP, Pi);
	}
	
	if(Model == M_DISC_DEP)
	{
		if(Trees->UseCovarion == FALSE)
			CreateDEPAMatrix(InvInfo, Rates, Trees, RateP);
		else
			CreateDEPAMatrixCoVar(InvInfo, Rates, Trees, RateP);
	}

	if(Model == M_DISC_INDEP)
	{
		if(Trees->UseCovarion == FALSE)
			CreateInDEPAMatrix(InvInfo, Rates, Trees, RateP);
		else
			CreateInDEPAMatrixCoVar(InvInfo, Rates, Trees, RateP);
	}

	if(Model == M_DISC_CV)
		CreateDepCVAMatrix(InvInfo, Rates, Trees, RateP);

	if(Model == M_DISC_HET)
	{
		if(Trees->UseCovarion == FALSE)
		{
			CreateInDEPAMatrix(Rates->Hetero->ModelInv[0], Rates, Trees, RateP);
			CreateDEPAMatrix(Rates->Hetero->ModelInv[1], Rates, Trees, &RateP[4]);
		}
		else
			exit(0);
	}
}



int		SetInvMat(MODEL Model, RATES *Rates, int NOS, INVINFO *InvInfo)
{
	if(Model != M_DISC_HET)
		return InvMat(InvInfo, NOS);

	if(InvMat(Rates->Hetero->ModelInv[0], NOS) == ERROR)
		return ERROR;

	return InvMat(Rates->Hetero->ModelInv[1], NOS);
}

double CaclNormConst(MATRIX *A, double *Pi)
{
	double Ret;
	int Index;

	Ret = 0;
	for(Index=0;Index<A->NoOfCols;Index++)
		Ret += Pi[Index] * A->me[Index][Index];

	return -1.0 / Ret;
}

void NormAMatrix(double NormC, MATRIX *A)
{
	int Index, No;

	No = A->NoOfCols * A->NoOfRows;

	for(Index=0;Index<No;Index++)
		A->me[0][Index] = A->me[0][Index] * NormC;
}

void	NormaliseAMatrix(RATES *Rates, MATRIX *A)
{
	int Index;

	Index = 0;
	Rates->NormConst = CaclNormConst(A, Rates->Pis);
	NormAMatrix(Rates->NormConst, A);


	// Will need to normalise the rates, not just the matrix. 

//	for(Index=0;Index<Rates->NoOfFullRates;Index++)
//		Rates->FullRates[Index] = Rates->FullRates[Index] * Rates->NormConst;

//	PrintMatrix(A, "A=", stdout);
//	exit(0);
}

int		SetUpAllAMatrix(RATES *Rates, TREES *Trees, OPTIONS *Opt)
{
	int PIndex, Pos, Err;

	Pos = 0;

	for(PIndex=0;PIndex<Rates->NoPatterns;PIndex++)
	{
		SetUpAMatrix(Opt->Model, Rates, Trees, Trees->NoStates, Trees->InvInfo[PIndex], &Rates->FullRates[Pos], Rates->Pis);

//		PrintMatrix( Trees->InvInfo[PIndex]->A, "A=", stdout);fflush(stdout);

		if(Opt->NormQMat == TRUE)
			NormaliseAMatrix(Rates, Trees->InvInfo[PIndex]->A);


		Err = SetInvMat(Opt->Model, Rates, Trees->NoStates, Trees->InvInfo[PIndex]); 

		if(Err == ERROR)
			return Err;
		
		Pos += Opt->DefNoRates;
	}

	return 0;
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
			for(SiteIndex=0;SiteIndex<Trees->NoSites;SiteIndex++)
			{
				if(Trees->NOSPerSite == FALSE)
					NOS = Trees->NoStates;
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
	double	*RateW;
	
	SetGammaBlank(Rates, Opt);

	RateW = (double*)SMalloc(sizeof(double) * Opt->GammaCats);

	DiscreteGamma(RateW, Rates->GammaMults, Rates->Gamma, Rates->Gamma, Rates->GammaCats, 0);

	free(RateW);
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
	NOS = Trees->NoStates;

	for(NIndex=0;NIndex<Tree->NoNodes;NIndex++)
	{
		N = Tree->NodeList[NIndex];

		if(N->Tip == FALSE)
		{
			for(SiteIndex=0;SiteIndex<Trees->NoSites;SiteIndex++)
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
			for(SiteIndex=0;SiteIndex<Trees->NoSites;SiteIndex++)
			{
				if(Trees->NOSPerSite == FALSE)
					NOS = Trees->NoStates;
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
				if(Opt->Model == M_MULTISTATE)
				{
					for(SIndex=0;SIndex<Trees->NoSites;SIndex++)
					{
						if(N->Taxa->EstDataP[SIndex] == TRUE)
							Taxa->DesDataChar[SIndex][0] = Trees->SymbolList[Rates->EstDescData[MDPos++]];
					}
				}
				else
				{
					if((N->Taxa->EstDataP[0] == TRUE) && (N->Taxa->EstDataP[1] == TRUE))
					{
						SetDiscEstDataTaxa(Taxa, '0'+Rates->EstDescData[MDPos], '0'+Rates->EstDescData[MDPos+1]);
						MDPos += 2;
					}
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

int		SetStdPMatrix(RATES *Rates, INVINFO *InvInfo, TREES *Trees, NODE N, MATRIX *P, double Gamma)
{
	double Len;
	double ErrVal;
	int		ThrNo;

	Len = N->Length * Gamma * Rates->GlobablRate;

	ThrNo = GetThreadNo();
	
	switch(Trees->NoStates)
	{
		case 2:
			ErrVal = Create2SPMat(Len, InvInfo, P, ThrNo);
		break;

		case 4:
			ErrVal = Create4SPMat(Len, InvInfo, P, ThrNo);
		break;

		default:
			ErrVal = CreatFullPMatrix(Len, InvInfo, P, Trees->NoStates, ThrNo);
		break;
	}
		
	if(ErrVal > 0.001)
		return TRUE;

	return FALSE;
}

int		SetAnalyticalPMatrix(TREES *Trees, NODE N, MATRIX *P, double Rate, double Gamma)
{
	double Len, ErrVal;
	
	Len = N->Length * Gamma;
	
	ErrVal = CreatFullAP(Len, Rate, Trees->NoStates, P);			

	return FALSE;
}


int		SetRPMatrix(TREES *Trees, NODE N, MATRIX *P, double Gamma)
{
	return FALSE;
	/* Needs to be fiex for R model as its not fixed form polytomie code. */
}

int		SetAllPMatrix(RATES* Rates, TREES *Trees, OPTIONS *Opt, double Gamma)
{
	int NIndex, Err, NoErr, PMatNo;
	TREE *Tree;
	NODE N;
	INVINFO *InvInfo;
	
	NoErr = 0;
	Err = FALSE;
	Tree = Trees->Tree[Rates->TreeNo];
	
#ifdef OPENMP_THR
	#pragma omp parallel for private(N, Err)
#endif	
	for(NIndex=1;NIndex<Tree->NoNodes;NIndex++)
	{
		N = Tree->NodeList[NIndex];
		Err = FALSE;
		if(NoErr == 0)
		{
			if(Opt->UseRModel == TRUE)
				Err = SetRPMatrix(Trees, N, Trees->PList[N->ID], Gamma);
			else
			{
				if(Opt->AnalyticalP == TRUE)
					Err = SetAnalyticalPMatrix(Trees, N, Trees->PList[N->ID], Rates->FullRates[0], Gamma);
				else
				{
					InvInfo = Trees->InvInfo[N->PatternNo];

					if(Opt->Model == M_DISC_HET)
					{
						PMatNo = Rates->Hetero->MList[NIndex];
						InvInfo = Rates->Hetero->ModelInv[PMatNo];
					}

					Err = SetStdPMatrix(Rates, InvInfo, Trees, N, Trees->PList[N->ID], Gamma);
				}
			}
		}

		if(Err == TRUE)
			NoErr++;
	}

	if(NoErr > 0)
		return TRUE;
	
	return FALSE;
}

void	RunNodeGroup(int GroupNo, RATES* Rates, TREE *Tree, TREES *Trees, OPTIONS *Opt, int SiteNo)
{
	int NIndex;

#ifdef OPENMP_THR
	#pragma omp parallel for
#endif
	for(NIndex=0;NIndex<Tree->NoFNodes[GroupNo];NIndex++)
	{
		#ifdef BIG_LH
			LhBigLh(Tree->FNodes[GroupNo][NIndex], Opt, Trees, Opt->Precision, SiteNo);
		#else
			#ifdef QUAD_DOUBLE
				NodeLhQuadDouble(Tree->FNodes[GroupNo][NIndex], Opt, Trees, SiteNo);
			#else
				SumLikeMultiState(Tree->FNodes[GroupNo][NIndex], Opt, Trees, SiteNo);
			#endif
		#endif
	}
}

void	SumLhLiner(RATES* Rates, TREES *Trees, OPTIONS *Opt, int SiteNo)
{
	int	GIndex, NIndex;
	TREE	*Tree;

	NIndex = 0;

	Tree = Trees->Tree[Rates->TreeNo];

#ifdef CLIK_P
	RecSumLikeMultiState(Tree->Root, Opt, Trees, SiteNo);
	return;
#endif

	
	for(GIndex=0;GIndex<Tree->NoFGroups;GIndex++)
		RunNodeGroup(GIndex, Rates, Tree, Trees, Opt, SiteNo);
}

double AddBigLh(RATES *Rates, TREES *Trees, OPTIONS *Opt)
{
	int NoUnderFlow, Index;
	TREE	*Tree;

	Tree = Trees->Tree[Rates->TreeNo];

	NoUnderFlow = 0;
	for(Index=0;Index<Tree->NoNodes;Index++)
		NoUnderFlow += Tree->NodeList[Index]->NoUnderFlow;

	return NoUnderFlow * log(LH_UNDER_FLOW);
}

double	CombineLh(RATES* Rates, TREES *Trees, OPTIONS *Opt)
{
	int SiteNo, Index, NOS;
	double Sum, SiteLh, Ret;
	TREE	*Tree;

	Tree = Trees->Tree[Rates->TreeNo];

	Ret = 0;
	for(SiteNo=0;SiteNo<Trees->NoSites;SiteNo++)
	{
		if(Trees->NOSPerSite == TRUE)
		{
			NOS = Trees->NOSList[SiteNo];
			for(Index=0;Index<NOS;Index++)
				Rates->Pis[Index]  = (double)1.0/NOS;
		}
		else
			NOS = Trees->NoStates;

#ifdef BIG_LH
		Ret += CombineBigLh(Rates, Trees, Opt, SiteNo, NOS);
#else
	#ifdef QUAD_DOUBLE
		Ret += CombineQuadDoubleLh(Rates, Trees, Opt, SiteNo, NOS);
	#else
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
	#endif
#endif
	}

	Ret += AddBigLh(Rates, Trees, Opt);

	return Ret;
}

void	LhTransformTree(RATES* Rates, TREES *Trees, OPTIONS *Opt)
{
	if(Opt->ModelType == MT_CONTINUOUS)
		return;

	if(Opt->UseLandscape == TRUE)
		MapLandscape(Opt, Trees, Rates);
	
	if(NeedToReSetBL(Opt, Rates) == TRUE)
	{
		SetUserBranchLength(Trees->Tree[Rates->TreeNo]);

		if(Rates->TimeSlices != NULL)
			ApplyTimeSlices(Rates, Trees);

//		PrintTreeBL(Trees->Tree[Rates->TreeNo]); exit(0);			

		TransformTree(Opt, Trees, Rates, NORMALISE_TREE_CON_SCALING);

		ApplyLocalTransforms(Rates, Trees, Opt, NORMALISE_TREE_CON_SCALING);

		if(Rates->VarRates != NULL)
			VarRatesTree(Opt, Trees, Rates, NORMALISE_TREE_CON_SCALING);
		
		
	}
//	ScaleTrees(Trees, 0.0000000001);
//	ScaleTrees(Trees, 0.000000001);
//	SaveTrees("DTest.trees", Trees); exit(0);
}

int		ValidLh(double LH, MODEL_TYPE MT)
{
	if(LH == LH+1 || LH != LH || LH == ERRLH)
		return FALSE;

	if(LH > 0 && MT == MT_DISCRETE)
		return FALSE;

	if(LH == ERRLH)
		return FALSE;

	return TRUE;
}

void	ZeroNoUnderFlow(TREE *Tree)
{
	int Index;

	for(Index=0;Index<Tree->NoNodes;Index++)
		Tree->NodeList[Index]->NoUnderFlow = 0;
}

double	Likelihood(RATES* Rates, TREES *Trees, OPTIONS *Opt)
{
	double	Ret;
	int		SiteNo;
	TREE	*Tree;
	int		Err;
	int		GammaCat;
	double	RateMult;
		

 	if(Rates->ModelFile == NULL)
		MapRates(Rates, Opt);
	else
		MapModelFile(Opt, Rates);
	
	LhTransformTree(Rates, Trees, Opt);
	
	if(Opt->NoLh == TRUE)
		return -1.0;

	if(Rates->AutoAccept == TRUE || Rates->CalcLh == FALSE)
		return Rates->Lh;
	
	if(Opt->UseDistData == TRUE && Opt->ModelType != MT_FATTAIL)
		SetTreeDistData(Rates, Opt, Trees);
	
	if(Opt->ModelType == MT_FATTAIL)
		return CalcTreeStableLh(Opt, Trees, Rates);

	if(Opt->ModelType == MT_CONTRAST)
		return CalcContrastLh(Opt, Trees, Rates);
	
	if(Opt->ModelType == MT_CONTINUOUS)
		return LHRandWalk(Opt, Trees, Rates);
	
	Tree = Trees->Tree[Rates->TreeNo];

	if(Rates->UseEstData == TRUE)
		SetDiscEstData(Rates, Trees, Opt);

	if(Opt->AnalyticalP == FALSE)
	{
		Err = SetUpAllAMatrix(Rates, Trees, Opt);

		if(Err == ERROR)
			return ERRLH;
	}

	if(Opt->UseGamma == TRUE)
		SetUpGamma(Rates, Opt);

	ZeroNoUnderFlow(Tree);

	for(GammaCat=0;GammaCat<Rates->GammaCats;GammaCat++)
	{
		if(Opt->UseGamma == FALSE)
			RateMult = 1;
		else
			RateMult = Rates->GammaMults[GammaCat];

#ifndef BTOCL
		Err = SetAllPMatrix(Rates, Trees, Opt, RateMult);
#else
		// use GPU for setting PMatrix
		Err = btocl_SetAllPMatrix(Rates, Trees, Opt, RateMult);
		// Use CPU
		//Err = SetAllPMatrix(Rates, Trees, Opt, RateMult);
#endif
		//}

			
		if(Err == TRUE) 
		{
			return ERRLH;
		}
		
//#ifdef BTOCL_DSC
//		btdebug_enter("partiallh");
//		for(SiteNo=0;SiteNo<100;SiteNo++){
		//printf("NoSites = %d\n",Trees->NoOfSites);
//		btocl_computePartialLh(Rates, Trees, Opt); // all sites
		
		//for(SiteNo=0;SiteNo<Trees->NoOfSites;SiteNo++)
		//{
		//	SumLhLiner(Rates, Trees, Opt, SiteNo);
			//Err = SumLh(Rates, Trees, Opt, SiteNo);
		
		//	if(Err == TRUE)
		//		return ERRLH;
		//}
		
//		}
//		btdebug_exit("partiallh");
//#else  

		for(SiteNo=0;SiteNo<Trees->NoSites;SiteNo++)
		{
			SumLhLiner(Rates, Trees, Opt, SiteNo);
//			Err = SumLh(Rates, Trees, Opt, SiteNo);

			if(Err == TRUE)
				return ERRLH;
		}
//#endif	

		if(Opt->UseGamma == TRUE)
			ProcessGamma(Rates, Trees);
	}
	
	// exit(0);

	if(Opt->UseGamma == TRUE)
		FinishUpGamma(Rates, Opt, Trees);

	Ret = CombineLh(Rates, Trees, Opt);
		
	if(ValidLh(Ret, Opt->ModelType) == FALSE)
		return ERRLH;

	return Ret;
}
