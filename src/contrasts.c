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
#include "Contrasts.h"
#include "Continuous.h"
#include "Praxis.h"
#include "RandLib.h"
#include "VarRates.h"
#include "Threaded.h"
#include "Trees.h"
#include "Part.h"
#include "Matrix.h"
#include "LinAlg.h"
#include "TransformTree.h"
#include "Likelihood.h"
#include "ContrastsStd.h"
#include "Priors.h"
#include "RJDummy.h"


#ifdef CLIK_P
	#include <cilk/cilk.h>
#endif

void RetSetConTraitData(TREE* Tree, int NoSites)
{
	int Index, SIndex;
	NODE Node;
	CONDATA* Con;
	
	for(Index=0;Index<Tree->NoNodes;Index++)
	{
		Node = Tree->NodeList[Index];
		if(Node->Tip == TRUE)
		{
			Con = Node->ConData;
			for(SIndex=0;SIndex<NoSites;SIndex++)
				Con->Contrast[0]->Data[SIndex] = Node->Taxa->ConData[SIndex];

//			memcpy(Con->Contrast[0]->Data, Node->Taxa->ConData, sizeof(double) * NoSites);
		}
	}
}

//void	RatesToConVals(OPTIONS *Opt, RATES *Rates, CONTRASTR* ConR);

int			GetNoContrastSites(OPTIONS *Opt, TREES *Trees)
{
	if(Opt->RJDummy == TRUE)
		return GetMaxDummy(Opt, Trees) + Trees->NoSites;

	return Trees->NoSites;
}



CONTRAST*	AllocContrastMem(int NoSites)
{
	CONTRAST* Ret;

	Ret = (CONTRAST*)SMalloc(sizeof(CONTRAST));

	Ret->Data = (double*)SMalloc(sizeof(double) * NoSites);
	Ret->Cont = (double*)SMalloc(sizeof(double) * NoSites);
//	Ret->Var  = (double*)SMalloc(sizeof(double) * NoSites);
//	Ret->Err  = (double*)SMalloc(sizeof(double) * NoSites);

	return Ret;
}

void	AllocContrast(NODE N, TREES *Trees, int NoSites)
{
	CONDATA		*ConData;
	CONTRAST*	Ret;
	int			Index, SIndex, NoC;
			
	ConData = (CONDATA*)SMalloc(sizeof(CONDATA));

	if(N->Tip == TRUE)
		NoC = 1;
	else
		NoC = N->NoNodes - 1;

	ConData->NoContrast = NoC;

	ConData->Contrast = (CONTRAST**)SMalloc(sizeof(CONTRAST*) * NoC);

	for(Index=0;Index<NoC;Index++)
	{
		Ret = AllocContrastMem(NoSites);
		Ret->Var = 0;
		Ret->Err = 0;

		for(SIndex=0;SIndex<NoSites;SIndex++)
		{
			Ret->Cont[SIndex]		= 0;
			Ret->Data[SIndex]		= -1;

			if(N->Tip == TRUE)
			{
				if(SIndex < Trees->NoSites)
					Ret->Data[SIndex] = N->Taxa->ConData[SIndex];
			}
			
		}

		ConData->Contrast[Index] = Ret;
	}

	ConData->GVar = (double*)SMalloc(sizeof(double) * NoSites);
	ConData->SumLogVar = (double*)SMalloc(sizeof(double) * NoSites);

	for(SIndex=0;SIndex<NoSites;SIndex++)
	{
		ConData->GVar[SIndex] = 0;
		ConData->SumLogVar[SIndex] = 0;
	}

	N->ConData = ConData;
}

int		GetTotalContrasts(TREE *Tree)
{
	int Ret, NIndex;
	NODE N;

	Ret = 0;
	for(NIndex=0;NIndex<Tree->NoNodes;NIndex++)
	{
		N = Tree->NodeList[NIndex];
		if(N->Tip == FALSE)
			Ret += N->ConData->NoContrast;
	}

	return Ret;
}

void	InitContrastTree(OPTIONS *Opt, TREES* Trees, int TNo, int NoSites)
{
	TREE *Tree;
	int NIndex;
	NODE N;

	Tree = Trees->Tree[TNo];
	for(NIndex=0;NIndex<Tree->NoNodes;NIndex++)
	{
		N = Tree->NodeList[NIndex];
		AllocContrast(N, Trees, NoSites);
	}

	Tree->NoContrast = GetTotalContrasts(Tree);

	SetNodesVPos(Trees, Tree);
}

void	InitContrastAll(OPTIONS *Opt, TREES* Trees)
{
	int TIndex;
	int	NoSites;

	CheckZeroTaxaBL(Trees);
	SetTreesDistToRoot(Trees);

	NoSites = GetNoContrastSites(Opt, Trees);

	for(TIndex=0;TIndex<Trees->NoTrees;TIndex++)
		InitContrastTree(Opt, Trees, TIndex, NoSites);

//	TransformContrastTreeFixed(Opt, Trees);
}

void	FreeContrastS(CONTRAST* C)
{
	free(C->Cont);
	free(C->Data);
	free(C);
}

void	FreeContrast(OPTIONS *Opt, TREES* Trees, int TreeNo)
{
	TREE *Tree;
	NODE N;
	CONDATA *ConData;
	int NIndex, CIndex;

	Tree = Trees->Tree[TreeNo];
	for(NIndex=0;NIndex<Tree->NoNodes;NIndex++)
	{
		N = Tree->NodeList[NIndex];
		ConData = N->ConData;
		
		for(CIndex=0;CIndex<ConData->NoContrast;CIndex++)
			FreeContrastS(ConData->Contrast[CIndex]);

		free(ConData->GVar);
		free(ConData->SumLogVar);
		free(ConData->Contrast);
		free(ConData);
	}
}

void	FreeAllContrast(OPTIONS *Opt, TREES* Trees)
{
	int TIndex;

	for(TIndex=0;TIndex<Trees->NoTrees;TIndex++)
		FreeContrast(Opt, Trees, TIndex);
}

CONTRAST*	CopyC(CONTRAST *Con, int NoSites)
{
	int Index;
	CONTRAST* Ret;

	Ret = AllocContrastMem(NoSites);
	
	Ret->Err = Con->Err;
	Ret->Var = Con->Var;

	for(Index=0;Index<NoSites;Index++)
	{
		Ret->Cont[Index] = Con->Cont[Index];
		Ret->Data[Index] = Con->Data[Index];
	}

	return Ret;
}

void	AddPolyContrast(CONTRAST *C0, CONTRAST *Dest, NODE Add, int NoSites)
{
	CONTRAST	*C1;
	double		t;
	double		l0, l1;
	int			SIndex;

	C1 = Add->ConData->Contrast[0];
	for(SIndex=0;SIndex<NoSites;SIndex++)
	{
		l0 = C0->Err;
		l1 = Add->Length + C1->Err;

		t = (l0 * C1->Data[SIndex]) +  (l1 * C0->Data[SIndex]);
		t = t / (l0 + l1);

		Dest->Data[SIndex] = t;
		Dest->Cont[SIndex] = C0->Data[SIndex] - C1->Data[SIndex];

		Dest->Err = (l0 * l1) / (l0 + l1);
		Dest->Var = l0 + l1;		
	}
}



void	PrintConNodeDx(NODE N, TREES *Trees)
{
	CONTRAST *C, *AncC;
	NODE Anc;
	double dx;
	int Index;

	C = N->ConData->Contrast[0];


	for(Index=0;Index<N->NoNodes;Index++)
	{
		Anc = N->NodeList[Index];

		AncC = Anc->ConData->Contrast[0];


		printf("%d\t", N->NoNodes);

		dx = C->Data[0] - AncC->Data[0];

		printf("%f\t", dx);
		printf("%f\t", Anc->Length);
		printf("%d\t", N->Part->NoTaxa);



//		PrintPart(stdout, Trees, N->Part);

		printf("\n");
	}
	
}


void	PrintPrintContrastDx(RATES *Rates, TREES *Trees)
{
	TREE *Tree;
	NODE N;
	int NIndex;

	Tree = Trees->Tree[0];

	printf("NoNodes\tdx\tLen\tNoTaxa\n");

	for(NIndex=0;NIndex<Tree->NoNodes;NIndex++)
	{
		N = Tree->NodeList[NIndex];
		if(N->Tip == FALSE)
			PrintConNodeDx(N, Trees);

	}

	exit(0);
}

void	PrintNodeContrast(NODE N, int NoSites)
{
	int Index;
	int CIndex;
	CONDATA *ConData;

	ConData = N->ConData;
	
//	printf("Con\t");
	// Data	Cont	Var	Err	u	y	v	V

	if(N->Tip == TRUE)
		return;

	for(CIndex=0;CIndex<ConData->NoContrast;CIndex++)
	{
		printf("%d\t%d\t", N->ID, ConData->NoContrast);

		for(Index=0;Index<NoSites;Index++)
		{ 
			printf("%12.12f\t", ConData->Contrast[CIndex]->Data[Index]); 
			printf("%12.12f\t", ConData->Contrast[CIndex]->Cont[Index]);
			printf("%12.12f\t", ConData->Contrast[CIndex]->Var);
			printf("%12.12f\t", ConData->Contrast[CIndex]->Err);


			printf("%12.12f\t", N->Length);
			printf("%12.12f\t", N->UserLength);

			printf("%12.12f\t", GetNodeHeight(N));
						
			printf("%d\t", NoTaxa(N));
			RecPRintNodeTaxa(N, ',');

			printf("\t");
		}

		printf("\n");
	}
}

void	PrintContrast(RATES *Rates, TREES *Trees)
{
	TREE *Tree;
	int Index;
	NODE N;
	CONTRASTR *CR; 

//	PrintPrintContrastDx(Rates, Trees);

	CR = Rates->Contrast;

	Tree = Trees->Tree[Rates->TreeNo];

	for(Index=0;Index<CR->NoSites;Index++)
		printf("NodeID\tNo Node Con\tData\tContrast\tVar\tError\tLen\tUser Length\tHeight\tNo Taxa\tTaxa"); 
	printf("\n");

	for(Index=0;Index<Tree->NoNodes;Index++)
	{
		N = Tree->NodeList[Index];
		PrintNodeContrast(N, CR->NoSites);
	}

//	exit(0);
}


void	CalcNodeContrast(NODE N, int NoSites)
{
	CONTRAST	*C, *C0, *C1;
	double		t;
	double		l0, l1;
	NODE		N0, N1;
	int			Index, SIndex;

	if(N->Tip == TRUE)
		return;

	N0 = N->NodeList[0];
	N1 = N->NodeList[1];

	C = N->ConData->Contrast[0];
	C0 = N0->ConData->Contrast[0];
	C1 = N1->ConData->Contrast[0];

	for(SIndex=0;SIndex<NoSites;SIndex++)
	{
		l0 = N0->Length + C0->Err;
		l1 = N1->Length + C1->Err;

		t = (l0 * C1->Data[SIndex]) +  (l1 * C0->Data[SIndex]);
		t = t / (l0 + l1);

		C->Data[SIndex] = t;
		C->Cont[SIndex] = C0->Data[SIndex] - C1->Data[SIndex];

		C->Err = (l0 * l1) / (l0 + l1);
		C->Var = l0 + l1;
	}

	for(Index=1;Index<N->ConData->NoContrast;Index++)
		AddPolyContrast(N->ConData->Contrast[Index-1], N->ConData->Contrast[Index], N->NodeList[Index+1], NoSites);

	C = N->ConData->Contrast[0];
	N->ConData->Contrast[0] = N->ConData->Contrast[N->ConData->NoContrast-1];
	N->ConData->Contrast[N->ConData->NoContrast-1] = C;

}




void	RecCalcContrast(NODE N, int NoSites)
{
/*	CONTRAST	*C, *C0, *C1;
	double		t;
	double		l0, l1;
	NODE		N0, N1; 
	int			SIndex; */
	int			Index;

	if(N->Tip == TRUE)
		return;

	for(Index=0;Index<N->NoNodes;Index++)
//#ifdef OPENMP_THR
//		if(N->NodeList[Index]->Visited == FALSE)
//#endif
			RecCalcContrast(N->NodeList[Index], NoSites);


	CalcNodeContrast(N, NoSites);
}



void	PrintContrasts(TREES* Trees, RATES* Rates)
{
	TREE	*Tree;
	int		Index, CIndex, SIndex;
	NODE	N;
	CONTRAST	*C;


	Tree = Trees->Tree[Rates->TreeNo];

	for(Index=0;Index<Tree->NoNodes;Index++)
	{
		N = Tree->NodeList[Index];

		for(CIndex=0;CIndex<N->ConData->NoContrast;CIndex++)
		{
			printf("Node Contrasts\t");
			RecPrintNode(N);
			printf("\t");
			C = N->ConData->Contrast[CIndex];

			for(SIndex=0;SIndex<Trees->NoSites;SIndex++)
				printf("%f\t%f\t%f\t%f\t", C->Cont[SIndex], C->Data[SIndex], C->Err, C->Var);

			printf("\n");
		}
	}
}

void	CalcContrastP(OPTIONS *Opt, TREES* Trees, RATES* Rates)
{
	TREE	*Tree;
	int		PIndex, NIndex;
		
	Tree = Trees->Tree[Rates->TreeNo];
/*
	#pragma omp parallel for num_threads(Opt->Cores)
	for(NIndex=0;NIndex<Tree->NoPNodes;NIndex++)
		RecCalcContrast(Tree->PNodes[NIndex], Trees->NoSites);

	RecCalcContrast(Tree->Root, Trees->NoSites);

	return;
*/
	for(PIndex=0;PIndex<Tree->NoParallelGroups;PIndex++)
	{
		if(Tree->ParallelGroupSize[PIndex] < 10)
		{
			for(NIndex=0;NIndex<Tree->ParallelGroupSize[PIndex];NIndex++)
				CalcNodeContrast(Tree->ParallelNodes[PIndex][NIndex], Trees->NoSites);
		}
		else
		{
#ifdef OPENMP_THR
			#pragma omp parallel for num_threads(Opt->Cores)
#endif
			for(NIndex=0;NIndex<Tree->ParallelGroupSize[PIndex];NIndex++)
				CalcNodeContrast(Tree->ParallelNodes[PIndex][NIndex], Trees->NoSites);
		}
	}
}

void	CalcContrast(OPTIONS *Opt, TREES* Trees, RATES* Rates)
{
	TREE	*Tree;
	
	Tree = Trees->Tree[Rates->TreeNo];

#ifdef OPENMP_THR
	CalcContrastP(Opt, Trees, Rates);
#else
	RecCalcContrast(Tree->Root, Rates->Contrast->NoSites);
#endif	
		
}

void RecCaclNodeGVar(NODE N, int SiteNo, double *GlobalVar, double *SumLogVar, int *NoCont)
{
	int			NIndex, CIndex;
	CONTRAST	*Con;

	if(N->Tip == TRUE)
		return;

	for(NIndex=0;NIndex<N->NoNodes;NIndex++)
		RecCaclNodeGVar(N->NodeList[NIndex], SiteNo, GlobalVar, SumLogVar, NoCont);

	for(CIndex=0;CIndex<N->ConData->NoContrast;CIndex++)
	{
		Con =  N->ConData->Contrast[CIndex];
		*GlobalVar += (Con->Cont[SiteNo] * Con->Cont[SiteNo]) / Con->Var;
		*SumLogVar += log(Con->Var);
		(*NoCont)++;
	}
}
/*
void	RecIntNode(NODE N, int SiteNo, double *Alpha, double *Sigma, double *Lh)
{
	double GlobalVar, SumLogVar, T1;
	int NoCont;

	SumLogVar = GlobalVar = 0;
	NoCont = 0;

	RecCaclNodeGVar(N, SiteNo, &GlobalVar, &SumLogVar, &NoCont);

	
	SumLogVar += log(N->ConData->Contrast[0]->Err[SiteNo]);

	T1 = GlobalVar;
	GlobalVar = GlobalVar / (NoCont+1);
	*Lh = (NoCont+1) * log(6.28318530717958647692528676655900576839 * GlobalVar);
	*Lh += SumLogVar + (T1 / GlobalVar);
	*Lh *= -0.5;

	*Alpha = N->ConData->Contrast[0]->Data[SiteNo];
	*Sigma = GlobalVar;
}
*/
void	RecIntNode(NODE N, int SiteNo, double *Alpha)
{
	*Alpha = N->ConData->Contrast[0]->Data[SiteNo];
}

double	CalcContrastCoVar(TREE *Tree, int S1, int S2)
{
	int			Index, CIndex;
	CONTRAST	*Con;
	NODE		N;
	double		Ret;

	Ret = 0;
	
	for(Index=0;Index<Tree->NoNodes;Index++)
	{
		N = Tree->NodeList[Index];
		if(N->Tip == FALSE)
		{
			for(CIndex=0;CIndex<N->ConData->NoContrast;CIndex++)
			{
				Con = N->ConData->Contrast[CIndex];

				Ret += (Con->Cont[S1] * Con->Cont[S2]) / Con->Var;
			}
		}
	}

	Ret = Ret / Tree->NoContrast;

	return Ret;
}



int	CalcMLContrastSigma(OPTIONS *Opt, TREES* Trees, RATES* Rates)
{
	TREE		*Tree;
	CONTRASTR	*ConRates;
	MATRIX		*Sig;
	int			x,y;
	
	ConRates = Rates->Contrast;
	Tree = Trees->Tree[Rates->TreeNo];

	Sig = Rates->Contrast->SigmaMat;
	for(x=0;x<Trees->NoSites;x++)
	{
		for(y=x;y<Trees->NoSites;y++)
		{
			Sig->me[x][y] = CalcContrastCoVar(Tree, x, y);
			Sig->me[y][x] = Sig->me[x][y];

			if(IsNum(Sig->me[x][y]) == FALSE)
				return FALSE;
		}
	}

//	PrintMatrix(Sig, "sig=", stdout);

	if(Opt->TestCorrel == FALSE)
	{
		for(x=0;x<Trees->NoSites;x++)
		{
			for(y=x+1;y<Trees->NoSites;y++)
			{
				Sig->me[x][y] = 0;
				Sig->me[y][x] = 0;
			}
		}
	}

	return TRUE;

//	PrintMatrix(Sig, "Sig = ", stdout);exit(0);

	// Used to fix sigma
//	Sig->me[0][0] = 0.003;
}

void	GetCalcSiteLhParams(TREE *Tree, int SiteNo, double *GlobalVar, double *SumLogVar)
{
	int			Index, CIndex;
	CONTRAST	*Con;
	NODE		N;

	*GlobalVar = 0;
	*SumLogVar = 0;
	
	for(Index=0;Index<Tree->NoNodes;Index++)
	{
		N = Tree->NodeList[Index];
		if(N->Tip == FALSE)
		{
			for(CIndex=0;CIndex<N->ConData->NoContrast;CIndex++)
			{
				Con = N->ConData->Contrast[CIndex];

				*GlobalVar += (Con->Cont[SiteNo] * Con->Cont[SiteNo]) / Con->Var;
				*SumLogVar += log(Con->Var);
			}
		}
	}
}

double	CalcSiteLh(OPTIONS *Opt, TREES* Trees, RATES* Rates, int SiteNo)
{
	TREE		*Tree;
	CONTRASTR	*ConRates;
	double		GlobalVar;
	double		SumLogVar;
	double		T1, Ret;
	int			NoCon;


	ConRates = Rates->Contrast;
	Tree = Trees->Tree[Rates->TreeNo];

	NoCon = 0;
	
	NoCon = Tree->NoContrast;

	GetCalcSiteLhParams(Tree, SiteNo, &GlobalVar, &SumLogVar);

	SumLogVar += log(Tree->Root->ConData->Contrast[0]->Err);

	T1 = GlobalVar;
	
	GlobalVar = GlobalVar / (NoCon+1);
	
	Ret = (NoCon+1) * log(6.28318530717958647692528676655900576839 * GlobalVar);
	Ret += SumLogVar + (T1 / GlobalVar);
	Ret *= -0.5;
	
	Rates->Contrast->Alpha[SiteNo] = Tree->Root->ConData->Contrast[0]->Data[SiteNo];
	Rates->Contrast->SigmaMat->me[SiteNo][SiteNo] = T1 / NoCon;

	return Ret;
}

double	CaclAllSiteLhContMult(double *u, double **InvC, double *Temp, int NoSites)
{
	int x,y;
	double Ret;

	Ret = 0;
	for(x=0;x<NoSites;x++)
	{
		Temp[x] = 0;
		for(y=0;y<NoSites;y++)
		{
			Temp[x] += u[y] * InvC[x][y];
		}
		Ret += Temp[x] * u[x];
	}

	return Ret;
}


double	CaclAllSiteLhContMultS2(double *u, double **InvC)
{
	double Ret;

	Ret = ((u[0] * InvC[0][0]) + (u[1] * InvC[0][1])) * u[0];
	Ret += ((u[0] * InvC[1][0]) + (u[1] * InvC[1][1])) * u[1];
	
	return Ret;
}
double	CaclAllSiteLhCont(TREE *Tree, RATES *Rates, double *SumLV)
{
	int			Index, CIndex;
	CONTRAST	*Con;
	NODE		N;
	double		Ret, Temp, NoSites;
	MATRIX		*InvSig;

	InvSig = Rates->Contrast->SigmaInvInfo->Inv;
	
	NoSites = Rates->Contrast->SigmaInvInfo->Inv->NoOfCols;

	*SumLV = 0;
	Ret = 0;
	for(Index=0;Index<Tree->NoNodes;Index++)
	{
		N = Tree->NodeList[Index];
		if(N->Tip == FALSE)
		{
			for(CIndex=0;CIndex<N->ConData->NoContrast;CIndex++)
			{
				Con = N->ConData->Contrast[CIndex];

				*SumLV += log(Con->Var);

				if(NoSites == 1)
					Temp = Con->Cont[0] * InvSig->me[0][0] * Con->Cont[0];

				if(NoSites == 2)
					Temp = CaclAllSiteLhContMultS2(Con->Cont, InvSig->me);

				if(NoSites > 2)
					Temp = CaclAllSiteLhContMult(Con->Cont, InvSig->me, Rates->Contrast->SigmaInvInfo->TempD, Rates->Contrast->SigmaInvInfo->Inv->NoOfCols);

				Ret += Temp / Con->Var;

		//		Ret += log(Con->Var[0]);

/*				printf("(*====================================*)\n");
				printf("u={{%f}, {%f}};\n", Con->Cont[0], Con->Cont[1]);
				PrintMathematicaMatrix(Rates->Contrast->SigmaInvInfo->Inv, "CI=", stdout);
				printf("Transpose[u].CI.u\n");
				printf("(* Res\t%f *)\n", Temp); */
			//	exit(0);
			}
		}
	}
	return Ret;
}

void	GetConStdMLAlpha(TREES *Trees, RATES *Rates)
{
	TREE *Tree;
	int Index;

	Tree = Trees->Tree[Rates->TreeNo];

	for(Index=0;Index<Trees->NoSites;Index++)
		Rates->Contrast->Alpha[Index] = Tree->Root->ConData->Contrast[0]->Data[Index];
}

double	CalcAllSiteLh(OPTIONS *Opt, TREES* Trees, RATES* Rates, double AlphaErr)
{
	TREE		*Tree;
	CONTRASTR	*ConRates;
	double		T1, Ret, SumLV, GlobalVar;
	int			NoCon;
	int			N, K;

	ConRates = Rates->Contrast;
	Tree = Trees->Tree[Rates->TreeNo];

	NoCon = Tree->NoContrast; 
	
	GlobalVar = CaclAllSiteLhCont(Tree, Rates, &SumLV);

	SumLV += log(Tree->Root->ConData->Contrast[0]->Err);
	
	N = NoCon + 1;
	K = Trees->NoSites;
	T1 = log(Rates->Contrast->SigmaInvInfo->Det);
	
	Ret = N * K * 1.83787706640935 + K * SumLV + N * log(Rates->Contrast->SigmaInvInfo->Det) + GlobalVar;

//	printf("\tPLh:\t%d\t%d\t%f\t%f\t%f\t%f\n", N, K, SumLV, GlobalVar, log(Rates->Contrast->SigmaInvInfo->Det), AlphaErr);

//	AlphaErr = (N * GlobalVar) + AlphaErr;
//	AlphaErr = AlphaErr / Rates->Contrast->Sigma->me[0][0];

	Ret += AlphaErr;

	Ret = -0.5 * Ret;

	if(Opt->Analsis == ANALML)
		GetConStdMLAlpha(Trees, Rates);
	
	return Ret;

// For Rob
/*	Ret = SumLV + GlobalVar;
	Ret += (NoCon + 1) * log(Rates->Contrast->SigmaInvInfo->Det);
	Ret += Trees->NoOfSites * (NoCon + 1) * 0.798179868;

	Ret = Ret * - 0.5;
	return Ret;

	*/
	
	T1 = GlobalVar;

	GlobalVar = GlobalVar / (NoCon + 1);

	Ret = (NoCon+1) * log(6.283185307 * GlobalVar);
	Ret += SumLV + (T1 / GlobalVar);

	Ret *= -0.5;
	
	printf("LH:\t%f\n", Ret);exit(0);

	return Ret;

}
/*
void	CalcContLh(OPTIONS *Opt, TREES* Trees, RATES* Rates)
{
	int Index;

	Rates->Lh = 0;
	for(Index=0;Index<Trees->NoOfSites;Index++)
		Rates->Lh += CalcSiteLh(Opt, Trees, Rates, Index);
} */

double	CalcAlphaErr(NODE Node, double* EstAlpha, MATRIX *Sig2, int NoSites)
{
	int Index;
	double Ret, SErr;
	CONTRAST *Con;

	Con = Node->ConData->Contrast[0];

	Ret = 0;
	for(Index=0;Index<NoSites;Index++)
	{
		SErr = (EstAlpha[Index] - Con->Data[Index]) * (EstAlpha[Index] - Con->Data[Index]);

		Ret += SErr / (Con->Err * Sig2->me[Index][Index]);
	}

	return Ret;
}

double NewCalcSig(OPTIONS *Opt, TREES* Trees, RATES* Rates, double Alpha)
{
	TREE		*Tree;
	CONTRASTR	*ConRates;
	double		Var;
	CONTRAST	*Con;
	NODE		N;
	int			CIndex, Index;


	ConRates = Rates->Contrast;
	Tree = Trees->Tree[Rates->TreeNo];
//	Alpha = ConRates->Alpha[0];
	
//	Alpha = ConRates->Alpha[0] + 0.1;
//	Alpha = ConRates->Alpha[0] + 0.1;

	Var = 0;
	
	for(Index=0;Index<Tree->NoNodes;Index++)
	{
		N = Tree->NodeList[Index];
		if(N->Tip == FALSE)
		{
			for(CIndex=0;CIndex<N->ConData->NoContrast;CIndex++)
			{
				Con = N->ConData->Contrast[CIndex];

				Var += ((Con->Cont[0] - Alpha) * (Con->Cont[0] - Alpha)) / Con->Var;
			}
		}
	}

	Var = Var / Tree->NoContrast;
	Var = Var - Alpha;

	return Var;
}

void	TestCalcSig(OPTIONS *Opt, TREES* Trees, RATES* Rates)
{
	TREE		*Tree;
	CONTRASTR	*ConRates;
	double		Alpha, Var;
	CONTRAST	*Con;
	NODE		N;
	int			CIndex, Index;
	
	
	for(Alpha = -10; Alpha < 10; Alpha +=0.1)
	{
		Var = NewCalcSig(Opt, Trees, Rates, Alpha);
		printf("%f\t%f\n", Alpha, Var);
	}

	exit(0);


	ConRates = Rates->Contrast;
	Tree = Trees->Tree[Rates->TreeNo];
	Alpha = ConRates->Alpha[0];
	
//	Alpha = ConRates->Alpha[0] + 0.1;
//	Alpha = ConRates->Alpha[0] + 0.1;

	Var = 0;
	
	for(Index=0;Index<Tree->NoNodes;Index++)
	{
		N = Tree->NodeList[Index];
		if(N->Tip == FALSE)
		{
			for(CIndex=0;CIndex<N->ConData->NoContrast;CIndex++)
			{
				Con = N->ConData->Contrast[CIndex];

				Var += ((Con->Cont[0] - Alpha) * (Con->Cont[0] - Alpha)) / Con->Var;

			}
		}
	}

	Var = Var / Tree->NoContrast;
	Var = Var - Alpha;
	printf("Var:\t%f\n", Var);
	exit(0);
}

double CalcContLh(OPTIONS *Opt, TREES* Trees, RATES* Rates)
{
	TREE		*Tree;
	CONTRASTR	*ConRates;
	double		Ret, AlphaErr;
	int			NoSites;

	ConRates = Rates->Contrast;
	Tree = Trees->Tree[Rates->TreeNo];

	NoSites = ConRates->NoSites;
	
	AlphaErr = 0;
	
	if(CalcMLContrastSigma(Opt, Trees, Rates) == FALSE)
		return ERRLH;

	if(Opt->Analsis == ANALMCMC)
		AlphaErr = CalcAlphaErr(Tree->Root, Rates->Contrast->Alpha, ConRates->SigmaMat, NoSites);
	
	memcpy(ConRates->SigmaInvVec, ConRates->SigmaMat->me[0], sizeof(double) * NoSites * NoSites);

//	printf("%f\t%f\n", Rates->OU, Rates->Kappa);
//	PrintMatrix(ConRates->SigmaMat, "Sig=", stdout);fflush(stdout);

	if(Matrix_Invert(ConRates->SigmaMat, ConRates->SigmaInvInfo) == ERROR)
	{
		printf("Sig invert error\n");
		return ERRLH;
	}

	memcpy(ConRates->SigmaMat->me[0], ConRates->SigmaInvVec, sizeof(double) * NoSites * NoSites);

	Ret = CalcAllSiteLh(Opt, Trees, Rates, AlphaErr);

//	TestCalcSig(Opt, Trees, Rates);

	return Ret;  
	
// Old Methods give same Lh; but no covar.
/*	int Index;
	Rates->Lh = 0;
	for(Index=0;Index<Trees->NoOfSites;Index++)
		Rates->Lh += CalcSiteLh(Opt, Trees, Rates, Index);*/
}

double CalcContrastMCMCSiteLh(OPTIONS *Opt, TREES* Trees, RATES* Rates, int SiteNo)
{
	TREE		*Tree;
	CONTRASTR	*ConRates;
	double		GlobalVar;
	double		SumLogVar;
	double		T1;
	int			NoCon;
	double		Ret;

	ConRates = Rates->Contrast;
	Tree = Trees->Tree[Rates->TreeNo];
	
	NoCon = Tree->NoContrast;
	GetCalcSiteLhParams(Tree, SiteNo, &GlobalVar, &SumLogVar);

//	GlobalVar = T->Root->ConData->GVar[SiteNo];
//	SumLogVar = T->Root->ConData->SumLogVar[SiteNo];
	
	SumLogVar += log(Tree->Root->ConData->Contrast[0]->Err);

	T1 = GlobalVar;
	GlobalVar = GlobalVar / NoCon;

	Ret = (NoCon+1) * log(6.283185307 * ConRates->SigmaMat->me[SiteNo][SiteNo]);

//	Alpha Err no longer avalable, can be calculated.
//	T2 = ((NoCon+1) * GlobalVar) + ConRates->AlphaErr[SiteNo];
//	Ret += SumLogVar + (T2 / ConRates->Sigma->me[SiteNo][SiteNo]);
//	Ret *= -0.5;

	Rates->Contrast->Alpha[SiteNo] = Tree->Root->ConData->Contrast[0]->Data[SiteNo];
	Rates->Contrast->SigmaMat->me[SiteNo][SiteNo] = GlobalVar;

	return Ret;
}

void	CalcContrastMCMC(OPTIONS *Opt, TREES* Trees, RATES* Rates)
{
	int SIndex;

	Rates->Lh = 0;

	for(SIndex=0;SIndex<Trees->NoSites;SIndex++)
		Rates->Lh += CalcContrastMCMCSiteLh(Opt, Trees, Rates, SIndex);
}

double	CalcMCMCAlpha(NODE N, double Alpha, int SiteNo)
{
	double Ret;
	CONTRAST	*Con;

	Con = N->ConData->Contrast[0];

	Ret = (Alpha - Con->Data[SiteNo]) * (Alpha - Con->Data[SiteNo]);
	Ret = Ret / Con->Err;	

	return Ret;
}
/*
double*		GetConList(TREE *Tree, int NoCon, int SiteNo)
{
	double *Ret;
	int		NIndex, CIndex, Pos;
	NODE	N;

	Ret = (double*)malloc(sizeof(double) * NoCon);
	if(Ret == NULL)
		MallocErr();
	Pos = 0;

	for(NIndex=0;NIndex<Tree->NoNodes;NIndex++)
	{
		N = Tree->NodeList[NIndex];
		if(N->Tip == FALSE)
		{
			for(CIndex=0;CIndex<N->ConData->NoContrast;CIndex++)
		//		Ret[Pos++] = N->ConData->Contrast[CIndex]->Cont[SiteNo];
				Ret[Pos++] = N->ConData->Contrast[CIndex]->Cont[SiteNo] / sqrt(N->ConData->Contrast[CIndex]->Var[SiteNo]);
		}
	}

	return Ret;
}
*/
double CaclRegSigma(TREE *Tree, double *Beta, int NoSites)
{
	double	Ret, Con, Var, Uy, Ux;
	int		NIndex, CIndex, SIndex;
	NODE	N;

	Ret = 0;
	
	for(NIndex=0;NIndex<Tree->NoNodes;NIndex++)
	{
		N = Tree->NodeList[NIndex];
		if(N->Tip == FALSE)
		{
			for(CIndex=0;CIndex<N->ConData->NoContrast;CIndex++)
			{
				Var = N->ConData->Contrast[CIndex]->Var;

				Uy = N->ConData->Contrast[CIndex]->Cont[0] / sqrt(Var);
				
				Ux = 0;
//				for(SIndex=0;SIndex<NoSites;SIndex++)
//					Ux += (N->ConData->Contrast[CIndex]->Cont[SIndex+1] / sqrt(Var)) * Beta[SIndex];

				for(SIndex=1;SIndex<NoSites;SIndex++)
					Ux += (N->ConData->Contrast[CIndex]->Cont[SIndex] / sqrt(Var)) * Beta[SIndex-1];

				Con = Uy - Ux;

				Con = Con * Con;
								
				Ret += Con;
			}
		}
	}

	Ret = Ret / (Tree->NoContrast + 1.0);
	return Ret;
}

double	CaclRegAlphaLh(OPTIONS *Opt, TREES* Trees, RATES* Rates, double Alpha)
{
	double AlphaErr;
	CONTRAST	*Con;
	TREE		*Tree;


	Tree = Trees->Tree[Rates->TreeNo];

	Con = Tree->Root->ConData->Contrast[0];

	AlphaErr = Con->Data[0] - Alpha;
	AlphaErr = AlphaErr  * AlphaErr;
	AlphaErr = AlphaErr / (Con->Err * Con->Var);

	return AlphaErr;
}

double	CalcRegLh(OPTIONS *Opt, TREES* Trees, RATES* Rates, double Alpha, double *Beta)
{
	int			Index, CIndex, SIndex, NoSites;
	TREE		*T;
	NODE		N;
	CONTRAST	*Con;
	CONTRASTR	*ConRates;
	double		GlobalVar;
	double		SumLogVar;
	double		T1, Ret;

	ConRates = Rates->Contrast;
	T = Trees->Tree[Rates->TreeNo];

	GlobalVar = 0;
	SumLogVar = 0;
//	NoSites = Trees->NoOfSites;
	NoSites = Rates->Contrast->NoSites;

	for(Index=0;Index<T->NoNodes;Index++)
	{
		N = T->NodeList[Index];
		if(N->Tip == FALSE)
		{
			for(CIndex=0;CIndex<N->ConData->NoContrast;CIndex++)
			{
				Con = N->ConData->Contrast[CIndex];

				T1 = 0;
				for(SIndex=1;SIndex<NoSites;SIndex++)
					T1 += Con->Cont[SIndex] * Beta[SIndex-1];

				T1 = Con->Cont[0] - T1;
				T1 = T1 * T1;

				GlobalVar += T1 / Con->Var;

				SumLogVar += log(Con->Var);
			}
		}
	}
		
//	GlobalVar = T->Root->ConData->GVar[SiteNo];
//	SumLogVar = T->Root->ConData->SumLogVar[SiteNo];
	
//	exit(0);
	SumLogVar += log(T->Root->ConData->Contrast[0]->Err);

	T1 = GlobalVar;

	
	GlobalVar = GlobalVar / (T->NoContrast+1);
	ConRates->GlobalVar = GlobalVar;

	Ret = (T->NoContrast+1) * log(6.283185307 * GlobalVar);
	Ret += SumLogVar + (T1 / GlobalVar);
	
	Ret *= -0.5;
	
	return Ret;
}

//	./Seq/Mammal-1.trees ./Seq/MamData.txt < in.txt > sout.txt
//	./Seq/ContrastReg/Freckleton.trees ./Seq/ContrastReg/Freckleton.txt < in.txt > sout.txt 
//	./Seq/ContrastReg/RegTree.trees ./Seq/ContrastReg/Data.txt < in.txt > sout.txt 
//	./Seq/ContrastReg/RegTree.trees ./Seq/ContrastReg/GoodC.txt < in.txt > sout.txt 

void	SetRegUyMatrix(TREE *Tree, MATRIX *Uy)
{
	int		NIndex, CIndex, Pos;
	NODE	N;

	Pos = 0;

	for(NIndex=0;NIndex<Tree->NoNodes;NIndex++)
	{
		N = Tree->NodeList[NIndex];
		if(N->Tip == FALSE)
		{
			for(CIndex=0;CIndex<N->ConData->NoContrast;CIndex++)
				Uy->me[0][Pos++] = N->ConData->Contrast[CIndex]->Cont[0] / sqrt(N->ConData->Contrast[CIndex]->Var);
		}
	}

	return ;
}

void	SetRegUxMatrix(TREE *Tree, MATRIX *Ux, int NoSites)
{
	int		NIndex, CIndex, Pos, SIndex;
	NODE	N;

		Pos = 0;

	for(SIndex=1;SIndex<NoSites;SIndex++)
	{
		Pos = 0;
		for(NIndex=0;NIndex<Tree->NoNodes;NIndex++)
		{
			N = Tree->NodeList[NIndex];
			if(N->Tip == FALSE)
			{
				for(CIndex=0;CIndex<N->ConData->NoContrast;CIndex++)
				{
					Ux->me[Pos][SIndex-1] = N->ConData->Contrast[CIndex]->Cont[SIndex] / sqrt(N->ConData->Contrast[CIndex]->Var);
					Pos++;
				}
			}
		}
	}
}

/*		Caulate the ML Beta Values using.

	(Ux'Ux)^-1(Ux'Uy)
*/

void	CaclRegBeta(OPTIONS *Opt, TREES* Trees, RATES* Rates)
{
	REG_BETA_SPACE *BSpace;
	int		NoCont, Index, Err;
	TREE	*Tree;

	Tree = Trees->Tree[Rates->TreeNo];
	NoCont = Tree->NoContrast;
	
	BSpace = Rates->Contrast->RegSapce;

	SetRegUxMatrix(Tree, BSpace->Ux, Trees->NoSites);

	Transpose(BSpace->Ux, BSpace->TUx);

	MatrixMult(BSpace->TUx, BSpace->Ux, BSpace->Prod1);
	
	SetRegUyMatrix(Tree, BSpace->Uy);
	
//	Mathematica Code
//	PrintMathematicaMatrix(BSpace->Ux, "Ux=", stdout);
//	PrintMathematicaMatrix(BSpace->Uy, "Uy=", stdout);
//	printf("Inverse[Transpose[Ux].Ux].(Transpose[Ux].Uy)\n");
	
	if(Trees->NoSites == 2)
	{
		if(Opt->TestCorrel == FALSE)
			BSpace->InvUx->me[0][0] = 1;
		else
			BSpace->InvUx->me[0][0] = 1.0 / BSpace->Prod1->me[0][0];
	}
	else
	{
		if(Opt->TestCorrel == FALSE)
			SetIdentityMatrix(BSpace->InvUx);
		else
		{
			

			Err = InvertMatrix(BSpace->Prod1->me, Trees->NoSites-1, BSpace->TempDVect, BSpace->TempIVect, BSpace->InvUx->me);

			if(Err != NO_ERROR)
			{
				printf("Matrix singular: %s::%d\n", __FILE__, __LINE__);
				exit(0);
			}
		}
	}

	MatrixMult(BSpace->TUx, BSpace->Uy, BSpace->Prod2);

	MatrixMult(BSpace->InvUx, BSpace->Prod2, BSpace->Prod3);

	for(Index=0;Index<Trees->NoSites-1;Index++)
		Rates->Contrast->RegBeta[Index] = BSpace->Prod3->me[0][Index];
}


// Method use pre Rob Code
/*
double	CaclStdContrastLh(OPTIONS *Opt, TREES* Trees, RATES* Rates)
{
	CONTRASTR	*Con;

	Con = Rates->Contrast;

	// if using ML paramter estiamte for MCMC.
#ifndef CONTRAST_ML_PARAM
	if(Opt->Analsis == ANALML)
#endif
	{
		CalcContLh(Opt, Trees, Rates);
				
		return Rates->Lh;
	}
		
	CalcContrastMCMC(Opt, Trees, Rates);
	
	return Rates->Lh;
} 
*/

double	CalcRegAlpha(TREE *Tree, CONTRASTR *CR, int NoSites)
{
	double Ret;
	int Index;
	CONTRAST	*C;

	C = Tree->Root->ConData->Contrast[0];
	Ret = 0;
	for(Index=1;Index<NoSites;Index++)
		Ret += C->Data[Index] * CR->RegBeta[Index-1];
	
	Ret = C->Data[0] - Ret;

// When normlised z scores. 
//	return 0;
	return Ret;
}

double CaclRegContrastLh(OPTIONS *Opt, TREES* Trees, RATES* Rates)
{
	CONTRASTR *CR;
	int		NoCon;
	TREE	*Tree;
	double	V0, Ret;
	int		NoSites;

	Tree = Trees->Tree[Rates->TreeNo];
	V0 = Tree->Root->ConData->Contrast[0]->Err;
	NoCon = Tree->NoContrast;

	CR = Rates->Contrast;

	if(Opt->Analsis == ANALML)
		CaclRegBeta(Opt, Trees, Rates);

	NoSites = Rates->Contrast->NoSites;

	CR->RegAlpha = CalcRegAlpha(Tree, CR, NoSites);
//	CR->RegSigma = CaclRegSigma(Tree, CR->RegBeta, NoSites);
		
		Ret = CalcRegLh(Opt, Trees, Rates, CR->RegAlpha, CR->RegBeta);
/*
	printf("LH\t%f\t%d\tAlpha:\t%f\tBeta\t%f\tSigma2\t%f\tV0\t%f\n", Ret, Rates->TreeNo, CR->RegAlpha, CR->RegBeta[0], CR->RegSigma, V0);

	printf("Lh:\t%f\n", Ret);
	printf("Alpha:\t%f\n", CR->RegAlpha);

	for(NoCon=0;NoCon<Trees->NoOfSites-1;NoCon++)
		printf("Beta %d\t%f\n", NoCon+1, CR->RegBeta[NoCon]);

	exit(0); */
	return Ret;
}

void	DummCodePreLh(OPTIONS *Opt, TREES* Trees, RATES* Rates)
{
	CONTRASTR *CRates;
	int NoRSize;

	NoRSize = Trees->NoSites - 1;

	CRates = Rates->Contrast;

	BuildDummyCodeBeta(Rates->RJDummy);
	MapDummyCodes(Trees, Rates);

	memcpy(&CRates->RegBeta[NoRSize], Rates->RJDummy->DummyBeta, sizeof(double) * Rates->RJDummy->NoDummyCode);
}


double	CalcContrastLh(OPTIONS *Opt, TREES* Trees, RATES* Rates)
{
	double Lh;

//	TransformTreeDelta(Trees->Tree[Rates->TreeNo]->Root, 60, TRUE);
//	TransContNodeKappa(Trees->Tree[Rates->TreeNo]->Root, 0.1, TRUE);
//	TransformTreeOU(Trees->Tree[Rates->TreeNo]->Root, 0.000000000000001, TRUE);
//	TransformTreeOU(Trees->Tree[Rates->TreeNo]->Root, 0.25, TRUE);
//	SaveTrees("DTest.trees", Trees); exit(0);
//	TransformTreeLambda(Trees->Tree[Rates->TreeNo]->Root, 0.5, TRUE);
	
	if(Opt->RJDummy == TRUE)
		DummCodePreLh(Opt, Trees, Rates);

	CalcContrast(Opt, Trees, Rates);	
	
	if(Opt->Model == M_CONTRAST)
		Lh = CaclStdContrastLh(Opt, Trees, Rates);

	if(Opt->Model == M_CONTRAST_CORREL)
		Lh = CalcContLh(Opt, Trees, Rates);
	
	if(Opt->Model == M_CONTRAST_REG)
		Lh = CaclRegContrastLh(Opt, Trees, Rates);
		
	if(ValidLh(Lh, Opt->ModelType) == FALSE)
		Lh = ERRLH;

//	PrintContrast(Rates, Trees); exit(0);

	return Lh;
}

CONTRASTR*	AllocContrastRates(OPTIONS *Opt, RATES *Rates)
{
	CONTRASTR*	Ret;

	Ret = (CONTRASTR*)SMalloc(sizeof(CONTRASTR));

	Ret->Alpha		= NULL;
	Ret->Sigma		= NULL;
	Ret->SigmaMat	= NULL;
	Ret->SigmaInvInfo=NULL;
	Ret->SigmaInvVec= NULL;
	Ret->RegBeta	= NULL;

	Ret->RegAlpha	=	-1;
	Ret->RegSigma	=	-1;

	Ret->RegSapce	=	NULL;

	return Ret;
}

void		StdConMCMCLHTest(OPTIONS *Opt, RATES *Rates)
{
	TREES *Trees;
	CONTRASTR *ConR;

	ConR = Rates->Contrast;

	Trees = Opt->Trees;

	CalcContrast(Opt, Trees, Rates);
	
	GetConStdMLAlpha(Trees, Rates);
	CalcMLContrastSigma(Opt, Trees, Rates);

	return;
}

void		InitStdContrastRates(OPTIONS *Opt, RATES *Rates, CONTRASTR* ConRates)
{
	int NoSites;

	NoSites = Opt->Trees->NoSites;
	ConRates->SigmaMat = AllocMatrix(NoSites, NoSites);
	ConRates->SigmaInvInfo = CreatMatInvertInfo(NoSites);


	ConRates->Alpha = (double*)SMalloc(sizeof(double) * NoSites);
	ConRates->SigmaInvVec = (double*)SMalloc(sizeof(double) * NoSites * NoSites);
	

	if(Opt->Analsis == ANALMCMC)
		StdConMCMCLHTest(Opt, Rates);
}

void		InitFullContrastRates(OPTIONS *Opt, RATES *Rates, CONTRASTR* ConRates)
{
	ConRates->Alpha = (double*)SMalloc(sizeof(double) * Opt->Trees->NoSites);
	ConRates->Sigma = (double*)SMalloc(sizeof(double) * Opt->Trees->NoSites);


	if(Opt->Analsis == ANALMCMC)
	{
		CalcContrast(Opt, Opt->Trees, Rates);
		CaclStdContrastLhML(Opt, Opt->Trees, Rates);
	}
}

int			GetMaxNoContrasts(TREES *Trees)
{
	int Ret, Index;

	Ret = Trees->Tree[0]->NoContrast;

	for(Index=1;Index<Trees->NoTrees;Index++)
	{
		if(Trees->Tree[Index]->NoContrast > Ret)
			Ret = Trees->Tree[Index]->NoContrast;
	}

	return Ret;
}

REG_BETA_SPACE*	InitRegBetaSpace(int NoSites, int NoCont)
{
	REG_BETA_SPACE* Ret;
	
	Ret = (REG_BETA_SPACE*)SMalloc(sizeof(REG_BETA_SPACE));
	
	Ret->Uy		= AllocMatrix(NoCont, 1);
	Ret->Ux		= AllocMatrix(NoCont, NoSites);
	Ret->TUx	= AllocMatrix(NoSites, NoCont);
	Ret->InvUx	= AllocMatrix(NoSites, NoSites);

	Ret->Prod1	= AllocMatrix(NoSites, NoSites);
	Ret->Prod2	= AllocMatrix(NoSites, 1);
	Ret->Prod3	= AllocMatrix(NoSites, 1);
	
	Ret->TempDVect	= (double*)malloc(sizeof(double) * NoSites);
	Ret->TempIVect	= (int*)malloc(sizeof(int) * NoSites);
	if((Ret->TempDVect == NULL) || (Ret->TempIVect == NULL))
		MallocErr();

	return Ret;
}


void		FreeRegBetaSpace(REG_BETA_SPACE* RSpace)
{
	FreeMatrix(RSpace->Uy);
	FreeMatrix(RSpace->Ux);
	FreeMatrix(RSpace->TUx);
	FreeMatrix(RSpace->InvUx);

	FreeMatrix(RSpace->Prod1);
	FreeMatrix(RSpace->Prod2);
	FreeMatrix(RSpace->Prod3);

	free(RSpace->TempDVect);
	free(RSpace->TempIVect);

	free(RSpace);
}

void		InitRegContrastRates(OPTIONS *Opt, RATES *Rates, CONTRASTR* ConRates)
{
	int NoSites;
	int	NoContrasts;

	NoSites = GetNoContrastSites(Opt, Opt->Trees);

	if(Opt->RJDummy == TRUE)
		Rates->RJDummy = CreatRJDummyCode(Opt, Opt->Trees);

	ConRates->RegBeta = (double*)SMalloc(sizeof(double) * NoSites);

	NoContrasts = GetMaxNoContrasts(Opt->Trees);
	
	ConRates->RegSapce = InitRegBetaSpace(Opt->Trees->NoSites-1, NoContrasts);

	if(Opt->Analsis == ANALMCMC)
	{
		CalcContrast(Opt, Opt->Trees, Rates);

//		PrintContrast(Rates, Opt->Trees);exit(0);
		
		CaclRegBeta(Opt, Opt->Trees, Rates);
		ConRates->RegAlpha = CalcRegAlpha(Opt->Trees->Tree[0], ConRates, Opt->Trees->NoSites);
	}
}

void FreeContrastRates(RATES *Rates)
{
	CONTRASTR*	CR;

	CR = Rates->Contrast;

	if(CR->Alpha != NULL)
		free(CR->Alpha);

	if(CR->Sigma != NULL)
		free(CR->Sigma);
	
	if(CR->SigmaMat != NULL)
		FreeMatrix(CR->SigmaMat);

	if(CR->RegBeta != NULL)
		free(CR->RegBeta);

	if(Rates->Contrast->RegSapce != NULL)
		FreeRegBetaSpace(Rates->Contrast->RegSapce);

	if(CR->SigmaInvInfo != NULL)
		FreeMatInvertInfo(CR->SigmaInvInfo);

	if(CR->SigmaInvVec != NULL)
		free(CR->SigmaInvVec);
	
	free(CR);
}

void	MapConValsToRates(OPTIONS *Opt, RATES *Rates, CONTRASTR* ConR)
{
	int NoSites;

	NoSites = Opt->Trees->NoSites;
	
	if(Opt->Model == M_CONTRAST_CORREL)
	{
		memcpy(Rates->Rates, ConR->Alpha, sizeof(double) * NoSites);
		return;
	}

	if(Opt->Model == M_CONTRAST)
	{
		memcpy(Rates->Rates, ConR->Alpha, sizeof(double) * NoSites);
		memcpy(&Rates->Rates[NoSites], ConR->Sigma, sizeof(double) * NoSites);
		return;
	}

	if(Opt->Model == M_CONTRAST_REG)
	{
	//	Rates->Rates[0] = ConR->RegAlpha;
		
		memcpy((void*)&Rates->Rates[0], ConR->RegBeta, sizeof(double) * (NoSites-1));
		return;
	}
}

void	MapRatesToConVals(OPTIONS *Opt, RATES *Rates, CONTRASTR* ConR)
{
	int NoSites;
	int Index;

//	NoSites = Opt->Trees->NoOfSites;
	NoSites = ConR->NoSites;
	
	if(Opt->Model == M_CONTRAST_CORREL)
		memcpy(ConR->Alpha, Rates->Rates, sizeof(double) * NoSites);
	
	if(Opt->Model == M_CONTRAST_REG)
	{
		// Reg Alpha will be set to ML values from Beta's
		ConR->RegAlpha = 0;

		if(Opt->TestCorrel == FALSE)
			for(Index=0;Index<Opt->Trees->NoSites - 1;Index++)
				Rates->Rates[Index] = 0.0;

		memcpy(ConR->RegBeta, (void*)Rates->Rates, sizeof(double) * (Opt->Trees->NoSites - 1));
	}

	if(Opt->Model == M_CONTRAST)
	{
		memcpy(ConR->Alpha, Rates->Rates, sizeof(double) * NoSites);
		memcpy(ConR->Sigma, &Rates->Rates[NoSites], sizeof(double) * NoSites);
	}
}



CONTRASTR*	CreatContrastRates(OPTIONS *Opt, RATES *Rates)
{
	CONTRASTR*	Ret;

	Ret = AllocContrastRates(Opt, Rates);

	Ret->NoSites = Opt->Trees->NoSites;

	Rates->Contrast = Ret;

	if(Opt->Model == M_CONTRAST_CORREL)
		InitStdContrastRates(Opt, Rates, Ret);
	
	if(Opt->Model == M_CONTRAST_REG)
		InitRegContrastRates(Opt, Rates, Ret);

	if(Opt->Model == M_CONTRAST)
		InitFullContrastRates(Opt, Rates, Ret);

	if(Opt->Analsis == ANALMCMC)
		MapConValsToRates(Opt, Rates, Ret);
	
	
	return Ret;
}

double	ChangeContrastRate(double Rate, double Dev, RANDSTATES *RS)
{
	double Ret;
	
//	do
//	{
		Ret = (RandDouble(RS) * Dev) - (Dev / 2.0); 
		Ret += Rate;
//	} while(Ret <= 0);

	if(Ret < 0)
		Ret = -Ret;

	return Ret;
}


void	MutateRegContrastRates(OPTIONS *Opt, TREES* Trees, RATES* Rates, SCHEDULE*	Shed)
{
}

void	MutateContrastRates(OPTIONS *Opt, TREES* Trees, RATES* Rates, SCHEDULE*	Shed)
{
	int Pos;
	double Dev;

	Shed->PNo = RandUSInt(Rates->RS) % Shed->NoParm;

	Pos = Shed->PNo;

	Shed->CurrentAT = Shed->RateDevATList[Pos];
	Dev = Shed->CurrentAT->CDev;

/* Do not use, it casues the SS to go mad, i don't know why. 
	if(Opt->Model == M_CONTRAST && Pos % 2 == 1)
	{
		Scale = exp(Dev * (RandDouble(Rates->RS) - 0.5));
		Rates->LnHastings = (Rates->Rates[Pos] * Scale) / Rates->Rates[Pos];
		Rates->Rates[Pos] = Rates->Rates[Pos] * Scale;
	}
	else*/
		Rates->Rates[Pos] += (RandDouble(Rates->RS) * Dev) - (Dev / 2.0);
}

void	CopyContrastRatesStd(RATES *R1, RATES* R2, CONTRASTR *C1, CONTRASTR	*C2, int NoSites)
{
	CopyMatrix(C1->SigmaMat, C2->SigmaMat);
	memcpy(C1->Alpha, C2->Alpha, sizeof(double) * NoSites);
}

void	CopyContrastRatesFull(RATES *R1, RATES* R2, CONTRASTR *C1, CONTRASTR *C2, int NoSites)
{
	memcpy(C1->Alpha, C2->Alpha, sizeof(double) * NoSites);
	memcpy(C1->Sigma, C2->Sigma, sizeof(double) * NoSites);
}

void	CopyContrastRatesReg(RATES *R1, RATES* R2, CONTRASTR *C1, CONTRASTR	*C2, int NoSites)
{
	memcpy(C1->RegBeta, C2->RegBeta, sizeof(double) * (NoSites - 1));
	C1->GlobalVar = C2->GlobalVar;
}

void	CopyContrastRates(OPTIONS *Opt, RATES* R1, RATES* R2, int NoSites)
{
	CONTRASTR	*C1, *C2;

	C1 = R1->Contrast;
	C2 = R2->Contrast;

	C1->NoSites = C2->NoSites;

	if(Opt->Analsis == ANALMCMC)
		memcpy(R1->Rates, R2->Rates, sizeof(double) * R1->NoOfRates);

	if(Opt->Model == M_CONTRAST_CORREL)
		CopyContrastRatesStd(R1, R2, C1, C2, NoSites);

	if(Opt->Model == M_CONTRAST_REG)
		CopyContrastRatesReg(R1, R2, C1, C2, NoSites);
	
	if(Opt->Model == M_CONTRAST)
		CopyContrastRatesFull(R1, R2, C1, C2, NoSites);
} 

double		DataToZScore(double X, double Mean, double SD)
{
	return (X - Mean) / SD;
}

void		ConvertDataToZScore(OPTIONS *Opt, TREES *Trees, RATES *Rates)
{
	double Temp;
	int Index, NIndex, SIndex;
	TREE	*Tree;
	NODE	N;

	for(Index=0;Index<Trees->NoTrees;Index++)
	{
		Tree = Trees->Tree[Index];
		for(NIndex=0;NIndex<Tree->NoNodes;NIndex++)
		{
			N = Tree->NodeList[NIndex];

			if(N->Tip == TRUE)
			{
				for(SIndex=0;SIndex<Trees->NoSites;SIndex++)
				{
					Temp = DataToZScore(N->Taxa->ConData[SIndex], Trees->PMean[SIndex], Trees->PSD[SIndex]);
					N->ConData->Contrast[0]->Data[SIndex] = Temp;
				}
			}
		}
	}
}

void		NormaliseReg(OPTIONS *Opt, TREES *Trees, RATES *Rates)
{
	int Index;	
	CONTRASTR *CR;

	return;
	
	CalcContrast(Opt, Opt->Trees, Rates);	
	
	Trees->PMean = (double*)SMalloc(sizeof(double) * Trees->NoSites);
	Trees->PSD	 = (double*)SMalloc(sizeof(double) * Trees->NoSites);

	CR = Rates->Contrast;

	CR->Alpha = (double*)SMalloc(sizeof(double) * Trees->NoSites);
	CR->Sigma = (double*)SMalloc(sizeof(double) * Trees->NoSites);


	CaclStdContrastLhML(Opt, Trees, Rates);

	for(Index=0;Index<Trees->NoSites;Index++)
	{
		Trees->PMean[Index] = CR->Alpha[Index];
		Trees->PSD[Index] = sqrt(CR->Sigma[Index]);
//		printf("%d\t%f\t%f\n", Index, Trees->PMean[Index], Trees->PSD[Index]);
	}
//	exit(0);
	ConvertDataToZScore(Opt, Trees, Rates);
}
