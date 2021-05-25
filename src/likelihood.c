#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "typedef.h"
#include "data.h"
#include "likelihood.h"
#include "matrix.h"
#include "genlib.h"
#include "Rates.h"
#include "continuous.h"
#include "gamma.h"
#include "trees.h"
#include "praxis.h"
#include "RandLib.h"
#include "contrasts.h"
#include "Threaded.h"
#include "BigLh.h"
#include "pmatrix.h"
#include "linalg.h"
#include "modelfile.h"
#include "QuadDouble.h"
#include "FatTail.h"
#include "TransformTree.h"
#include "VarRates.h"
#include "LocalTransform.h"
#include "DistData.h"

#ifdef BTOCL
	#include "btocl_discrete.h"
#endif


double	CreatFullPMatrix(double t, INVINFO	*InvInfo, MATRIX *Mat, TREES* Trees, MATRIX *A, double *Et);

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
	int	NOS;
	
	NOS = Trees->NoOfStates;

	Trees->InvInfo = AllocInvInfo(NOS);

	Trees->PList = AllocMultiMatrixLinMem(Trees->MaxNodes, NOS, NOS);
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

void	SumLikeMultiState(NODE N, OPTIONS *Opt, TREES *Trees, int SiteNo,  int Rec)
{
	int		Inner, Outter, NIndex;
	double	Lh;
	double	**Mat, **Partial;

	if(Rec == TRUE)
	{
		for(NIndex=0;NIndex<N->NoNodes;NIndex++)
		{
			if(N->NodeList[NIndex]->Tip == FALSE)
				SumLikeMultiState(N->NodeList[NIndex], Opt, Trees, SiteNo, Rec);
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

	if(N->FossilMask != NULL)
		FossilLh(N, Opt, Trees, SiteNo);
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

	if(RateW == NULL)
	{
		RateW = (double*)malloc(sizeof(double) * Opt->GammaCats);
		if(RateW == NULL)
			MallocErr();
	}

	DiscreteGamma(RateW, Rates->GammaMults, Rates->Gamma, Rates->Gamma, Rates->GammaCats, 0);
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
				if(Opt->Model == M_MULTISTATE)
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

int		SetStdPMatrix(INVINFO *InvInfo, TREES *Trees, NODE N, MATRIX *P, double Gamma)
{
	double Len;
	double ErrVal;
	int		ThrNo;

	Len = N->Length * Gamma;

	ThrNo = GetThreadNo();
	
	switch(Trees->NoOfStates)
	{
		case 2:
			ErrVal = Create2SPMat(Len, InvInfo, P, Trees, InvInfo->As[ThrNo], InvInfo->Ets[ThrNo]);
		break;

		case 4:
			ErrVal = Create4SPMat(Len, InvInfo, P, Trees, InvInfo->As[ThrNo], InvInfo->Ets[ThrNo]);
		break;

		default:
			ErrVal = CreatFullPMatrix(Len, InvInfo, P, Trees, InvInfo->As[ThrNo], InvInfo->Ets[ThrNo]);
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
	
	

	ErrVal = CreatFullAP(Len, Rate, Trees->NoOfStates, P);			

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
	
	NoErr = 0;
	Err = FALSE;
	Tree = Trees->Tree[Rates->TreeNo];
	/* No need to compute P for root */

	
#ifdef OPENMP_THR
	#pragma omp parallel for private(N, Err) num_threads(Opt->Cores)
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
				if(Opt->AnalyticalP == FALSE)
				{
					if(Opt->Model != M_DESCHET)
						Err = SetStdPMatrix(Trees->InvInfo, Trees, N, Trees->PList[N->ID], Gamma);
					else
					{
						PMatNo = Rates->Hetero->MList[NIndex];
						Err = SetStdPMatrix(Rates->Hetero->ModelInv[PMatNo], Trees, N, Trees->PList[N->ID], Gamma);
					}
				}
				else
					Err = SetAnalyticalPMatrix(Trees, N, Trees->PList[N->ID], Rates->FullRates[0], Gamma);
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
		for(NIndex=0;NIndex<Tree->NoFNodes[GroupNo];NIndex++)
		{
			#ifdef BIG_LH
				LhBigLh(Tree->FNodes[GroupNo][NIndex], Opt, Trees, Opt->Precision, SiteNo);
			#else
				#ifdef QUAD_DOUBLE
					NodeLhQuadDouble(Tree->FNodes[GroupNo][NIndex], Opt, Trees, SiteNo);
				#else
					SumLikeMultiState(Tree->FNodes[GroupNo][NIndex], Opt, Trees, SiteNo, FALSE);	
				#endif
			#endif
		}

		return;
	}

#ifdef OPENMP_THR
	#pragma omp parallel for num_threads(Opt->Cores)
#endif
	for(NIndex=0;NIndex<Tree->NoFNodes[GroupNo];NIndex++)
	{
		#ifdef BIG_LH
			LhBigLh(Tree->FNodes[GroupNo][NIndex], Opt, Trees, Opt->Precision, SiteNo);
		#else
			#ifdef QUAD_DOUBLE
				NodeLhQuadDouble(Tree->FNodes[GroupNo][NIndex], Opt, Trees, SiteNo);
			#else
				SumLikeMultiState(Tree->FNodes[GroupNo][NIndex], Opt, Trees, SiteNo, FALSE);
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
	
	for(GIndex=0;GIndex<Tree->NoFGroups;GIndex++)
	{
#ifndef OPENMP_THR
		RunNodeGroup(GIndex, FALSE, Rates, Tree, Trees, Opt, SiteNo);
#else
		if(Opt->Cores == 1)
			RunNodeGroup(GIndex, FALSE, Rates, Tree, Trees, Opt, SiteNo);
		else
			RunNodeGroup(GIndex, TRUE, Rates, Tree, Trees, Opt, SiteNo);
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
				Rates->Pis[Index]  = (double)1.0/NOS;
		}
		else
			NOS = Trees->NoOfStates;

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

//		printf("site:\t%d\t%f\t%f\n", SiteNo, log(SiteLh), Rates->Rates[0]);fflush(stdout);
		if(IsNum(log(SiteLh)) == FALSE)
			return ERRLH;

		Ret += log(SiteLh);
	#endif
#endif
	}

	return Ret;
}

void	LhTransformTree(RATES* Rates, TREES *Trees, OPTIONS *Opt)
{
	if(Opt->ModelType == MT_CONTINUOUS)
		return;
	
	if(NeedToReSetBL(Opt, Rates) == TRUE)
	{
		ReSetBranchLength(Trees->Tree[Rates->TreeNo]);
				
		TransformTree(Opt, Trees, Rates, NORMALISE_TREE_CON_SCALING);

		ApplyLocalTransforms(Rates, Trees, Opt, NORMALISE_TREE_CON_SCALING);

		if(Rates->VarRates != NULL)
			VarRatesTree(Opt, Trees, Rates, NORMALISE_TREE_CON_SCALING);
		
	//	SaveTrees("DTest.trees", Trees); exit(0);
	}
}

int		ValidLh(double LH)
{
	if(LH > 0 || LH == LH+1 || LH != LH || LH == ERRLH)
		return FALSE;

	return TRUE;
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

		for(SiteNo=0;SiteNo<Trees->NoOfSites;SiteNo++)
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

	
	if(ValidLh(Ret) == FALSE)
		return ERRLH;

	return Ret;
}
