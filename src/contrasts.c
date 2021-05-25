#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

#include "typedef.h"
#include "genlib.h"
#include "contrasts.h"
#include "praxis.h"
#include "rand.h"
#include "phyloplasty.h"

CONTRAST*	AllocContrastMem(int NoSites)
{
	CONTRAST* Ret;

	Ret = (CONTRAST*)malloc(sizeof(CONTRAST));
	if(Ret == NULL)
		MallocErr();

	Ret->Data = (double*)malloc(sizeof(double) * NoSites);
	Ret->Cont = (double*)malloc(sizeof(double) * NoSites);
	Ret->Var  = (double*)malloc(sizeof(double) * NoSites);
	Ret->Err  = (double*)malloc(sizeof(double) * NoSites);

	if( (Ret->Data == NULL) || 
		(Ret->Cont == NULL) || 
		(Ret->Var == NULL) || 
		(Ret->Err == NULL))
		MallocErr();

	return Ret;
}


void	AllocContrast(NODE N, TREES *Trees)
{
	CONTRAST*	Ret;
	int			Index, SIndex;
	
	if(N->Tip == TRUE)
		N->NoContrast = 1;
	else
		N->NoContrast = N->NoNodes - 1;

	N->Contrast = (CONTRAST**)malloc(sizeof(CONTRAST*) * N->NoContrast);
	if(N->Contrast == NULL)
		MallocErr();

	for(Index=0;Index<N->NoContrast;Index++)
	{
		Ret = AllocContrastMem(Trees->NoOfSites);

		for(SIndex=0;SIndex<Trees->NoOfSites;SIndex++)
		{
			if(N->Tip == TRUE)
				Ret->Data[SIndex] = N->Taxa->ConData[SIndex];
			else
				Ret->Data[SIndex] = -1;

			Ret->Cont[SIndex]		= 0;
			Ret->Var[SIndex]		= 0;
			Ret->Err[SIndex]		= 0;
		}

		N->Contrast[Index] = Ret;
	}
}


void	InitContrastTree(OPTIONS *Opt, TREES* Trees, int TNo)
{
	TREE *Tree;
	int NIndex;
	NODE N;

	Tree = &Trees->Tree[TNo];
	for(NIndex=0;NIndex<Tree->NoNodes;NIndex++)
	{
		N = Tree->NodeList[NIndex];
		AllocContrast(N, Trees);
	}

}

void	InitContrastAll(OPTIONS *Opt, TREES* Trees)
{
	int TIndex;

	for(TIndex=0;TIndex<Trees->NoOfTrees;TIndex++)
		InitContrastTree(Opt, Trees, TIndex);
}

void	FreeContrastS(CONTRAST* C)
{
	free(C->Cont);
	free(C->Data);
	free(C->Err);
	free(C->Var);
	free(C);
}

void	FreeContrast(OPTIONS *Opt, TREES* Trees, int TreeNo)
{
	TREE *Tree;
	NODE N;
	int NIndex, CIndex;

	Tree = &Trees->Tree[TreeNo];
	for(NIndex=0;NIndex<Tree->NoNodes;NIndex++)
	{
		N = Tree->NodeList[NIndex];

		for(CIndex=0;CIndex<N->NoContrast;CIndex++)
			FreeContrastS(N->Contrast[CIndex]);
		free(N->Contrast);
	}
}

void	FreeAllContrast(OPTIONS *Opt, TREES* Trees)
{
	int TIndex;

	for(TIndex=0;TIndex<Trees->NoOfTrees;TIndex++)
		FreeContrast(Opt, Trees, TIndex);
}

CONTRAST*	CopyC(CONTRAST *Con, int NoSites)
{
	int Index;
	CONTRAST* Ret;

	Ret = AllocContrastMem(NoSites);

	for(Index=0;Index<NoSites;Index++)
	{
		Ret->Cont[Index] = Con->Cont[Index];
		Ret->Data[Index] = Con->Data[Index];
		Ret->Err[Index] = Con->Err[Index];
		Ret->Var[Index] = Con->Var[Index];
	}

	return Ret;
}

void	AddPolyContrast(CONTRAST *C0, CONTRAST *Dest, NODE Add, int NoSites)
{
	CONTRAST	*C1;
	double		t;
	double		l0, l1;
	int			SIndex;

	C1 = Add->Contrast[0];
	for(SIndex=0;SIndex<NoSites;SIndex++)
	{
		l0 = C0->Err[SIndex];
		l1 = Add->Length + C1->Err[SIndex];

		t = (l0 * C1->Data[SIndex]) +  (l1 * C0->Data[SIndex]);
		t = t / (l0 + l1);

		Dest->Data[SIndex] = t;
		Dest->Cont[SIndex] = C0->Data[SIndex] - C1->Data[SIndex];

		Dest->Err[SIndex] = (l0 * l1) / (l0 + l1);
		Dest->Var[SIndex] = l0 + l1;		
	}
}

void	RecCalcContrast(NODE N, int NoSites)
{
	CONTRAST	*C, *C0, *C1;
	double		t;
	double		l0, l1;
	NODE		N0, N1;
	int			Index, SIndex;

	if(N->Tip == TRUE)
		return;

	for(Index=0;Index<N->NoNodes;Index++)
		RecCalcContrast(N->NodeList[Index], NoSites);
	
	if(N->NoNodes == 1)
	{
		printf("nodes with only one descente caouse problems. \n");
		exit(0);
	}

	N0 = N->NodeList[0];
	N1 = N->NodeList[1];

	C = N->Contrast[0];
	C0 = N0->Contrast[0];
	C1 = N1->Contrast[0];

	for(SIndex=0;SIndex<NoSites;SIndex++)
	{
		l0 = N0->Length + C0->Err[SIndex];
		l1 = N1->Length + C1->Err[SIndex];

		t = (l0 * C1->Data[SIndex]) +  (l1 * C0->Data[SIndex]);
		t = t / (l0 + l1);
		
		C->Data[SIndex] = t;
		C->Cont[SIndex] = C0->Data[SIndex] - C1->Data[SIndex];

		C->Err[SIndex] = (l0 * l1) / (l0 + l1);
		C->Var[SIndex] = l0 + l1;

		if((C->Data[SIndex] != C->Data[SIndex]) ||
			(C->Cont[SIndex] != C->Cont[SIndex]) ||
			(C->Err[SIndex] != C->Err[SIndex]) ||
			(C->Var[SIndex] != C->Var[SIndex]))
			printf("err\n");
	}

	for(Index=1;Index<N->NoContrast;Index++)
		AddPolyContrast(N->Contrast[Index-1], N->Contrast[Index], N->NodeList[Index+1], NoSites);
		
	C = N->Contrast[0];
	N->Contrast[0] = N->Contrast[N->NoContrast-1];
	N->Contrast[N->NoContrast-1] = C;
}

void	RecPrintNode(NODE N)
{
	int Index;

	if(N->Tip == TRUE)
	{
		printf("%s ", N->Taxa->Name);
		return;
	}

	for(Index=0;Index<N->NoNodes;Index++)
		RecPrintNode(N->NodeList[Index]);
}


void	PrintContrasts(TREES* Trees, RATES* Rates)
{
	TREE	*Tree;
	int		Index, CIndex;
	NODE	N;
	CONTRAST	*C;


	Tree = &Trees->Tree[Rates->TreeNo];

	for(Index=0;Index<Tree->NoNodes;Index++)
	{
		N = Tree->NodeList[Index];

		for(CIndex=0;CIndex<N->NoContrast;CIndex++)
		{
			printf("Node Contrasts\t");
			RecPrintNode(N);
			printf("\t");
			C = N->Contrast[CIndex];

			printf("%f\t%f\t%f\t%f\t", C->Cont, C->Data, C->Err, C->Var);

			printf("\n");
		}
	}
}

void	CalcContrast(TREES* Trees, RATES* Rates)
{
	TREE	*Tree;

	Tree = &Trees->Tree[Rates->TreeNo];
	RecCalcContrast(Tree->Root, Trees->NoOfSites);
}

double	CalcSiteLh(OPTIONS *Opt, TREES* Trees, RATES* Rates, int SiteNo)
{
	int			Index;
	TREE		*T;
	NODE		N;
	CONTRAST	*Con;
	CONTRASTR	*ConRates;
	double		GlobalVar;
	double		SumLogVar;
	double		T1, Ret;
	int			NoCon, CIndex;

	ConRates = Rates->Contrast;
	T = &Trees->Tree[Rates->TreeNo];

	NoCon = 0;
	GlobalVar = 0;
	SumLogVar = 0;

	for(Index=0;Index<T->NoNodes;Index++)
	{
		N = T->NodeList[Index];
		if(N->Tip == FALSE)
		{
			for(CIndex=0;CIndex<N->NoContrast;CIndex++)
			{
				Con = N->Contrast[CIndex];

				GlobalVar += (Con->Cont[SiteNo] * Con->Cont[SiteNo]) / Con->Var[SiteNo];

//				printf("NL\t%f\t%f\t%f\n", N->Length, Con->Cont[SiteNo], Con->Var[SiteNo]);

				SumLogVar += log(Con->Var[SiteNo]);
				NoCon++;
			}
		}
	}

//	exit(0);
	SumLogVar += log(T->Root->Contrast[0]->Err[SiteNo]);

	T1 = GlobalVar;
	GlobalVar = GlobalVar / (NoCon+1);
	Ret = (NoCon+1) * log(6.28318530717958647692528676655900576839 * GlobalVar);
	Ret += SumLogVar + (T1 / GlobalVar);
	Ret *= -0.5;

	Rates->Contrast->Alpha[SiteNo] = T->Root->Contrast[0]->Data[SiteNo];
	Rates->Contrast->Sigma[SiteNo] = GlobalVar;

	return Ret;
}

void	CalcContLh(OPTIONS *Opt, TREES* Trees, RATES* Rates)
{
	int			Index;

	Rates->Lh = 0;
	for(Index=0;Index<Trees->NoOfSites;Index++)
		Rates->Lh += CalcSiteLh(Opt, Trees, Rates, Index);
}

double CalcContrastMCMCSiteLh(OPTIONS *Opt, TREES* Trees, RATES* Rates, int SiteNo)
{
	int			Index;
	TREE		*T;
	NODE		N;
	CONTRAST	*Con;
	CONTRASTR	*ConRates;
	double		GlobalVar;
	double		SumLogVar;
	double		T1, T2;
	int			NoCon, CIndex;
	double		Ret;

	ConRates = Rates->Contrast;
	T = &Trees->Tree[Rates->TreeNo];

	NoCon = 0;
	GlobalVar = 0;
	SumLogVar = 0;

	for(Index=0;Index<T->NoNodes;Index++)
	{
		N = T->NodeList[Index];
		if(N->Tip == FALSE)
		{
			for(CIndex=0;CIndex<N->NoContrast;CIndex++)
			{
				Con = N->Contrast[CIndex];
				GlobalVar += (Con->Cont[SiteNo] * Con->Cont[SiteNo]) / Con->Var[SiteNo];

				SumLogVar += log(Con->Var[SiteNo]);
				NoCon++;
			}
		}
	}

	SumLogVar += log(T->Root->Contrast[0]->Err[SiteNo]);

	T1 = GlobalVar;
	GlobalVar = GlobalVar / NoCon;

	Ret = (NoCon+1) * log(6.28318530717958647692528676655900576839 * ConRates->EstSigma[SiteNo]);

	T2 = ((NoCon+1) * GlobalVar) + ConRates->AlphaErr[SiteNo];
	Ret += SumLogVar + (T2 / ConRates->EstSigma[SiteNo]);
	Ret *= -0.5;

	Rates->Contrast->Alpha[SiteNo] = T->Root->Contrast[0]->Data[SiteNo];
	Rates->Contrast->Sigma[SiteNo] = GlobalVar;

	return Ret;
}

void	CalcContrastMCMC(OPTIONS *Opt, TREES* Trees, RATES* Rates)
{
	int SIndex;

	Rates->Lh = 0;

	for(SIndex=0;SIndex<Trees->NoOfSites;SIndex++)
		Rates->Lh += CalcContrastMCMCSiteLh(Opt, Trees, Rates, SIndex);
}

double	CalcMCMCAlpha(NODE N, double Alpha, int SiteNo)
{
	double Ret;
	CONTRAST	*Con;

	Con = N->Contrast[0];

	Ret = (Alpha - Con->Data[SiteNo]) * (Alpha - Con->Data[SiteNo]);
	Ret = Ret / Con->Err[SiteNo];	

	return Ret;
}

double	CalcContrastLh(OPTIONS *Opt, TREES* Trees, RATES* Rates)
{
	CONTRASTR	*Con;
	int			SIndex;

	Con = Rates->Contrast;

	if(Opt->UsePhyloPlasty == TRUE)
		Plasty(Opt, Trees, Rates);

	CalcContrast(Trees, Rates);	

	if(Opt->Analsis == ANALML)
	{
		CalcContLh(Opt, Trees, Rates);
		//	Use ML values for MCMC
		if(Opt->Analsis == ANALMCMC)
		{
			for(SIndex=0;SIndex<Trees->NoOfSites;SIndex++)
			{
				Con->EstAlpha[SIndex] = Con->Alpha[SIndex];
				Con->EstSigma[SIndex] = Con->Sigma[SIndex];
			}
		}

//		exit(0);
		return Rates->Lh;
	}

	for(SIndex=0;SIndex<Trees->NoOfSites;SIndex++)
		Con->AlphaErr[SIndex] = CalcMCMCAlpha(Trees->Tree[0].Root, Con->EstAlpha[SIndex], SIndex);

	CalcContrastMCMC(Opt, Trees, Rates);
	
//	PrintContrasts(Trees, Rates);
//	printf("Lh\t%f\n", Rates->Lh);
//	exit(0);

	return Rates->Lh;
}


CONTRASTR*	AllocContrastRates(OPTIONS *Opt, RATES *Rates)
{
	CONTRASTR*	Ret;
	int			Index;

	Ret = (CONTRASTR*)malloc(sizeof(CONTRASTR));
	if(Ret == NULL)
		MallocErr();

	Ret->Alpha = (double*)malloc(sizeof(double) * Opt->Trees->NoOfSites);
	Ret->Sigma = (double*)malloc(sizeof(double) * Opt->Trees->NoOfSites);

	if((Ret->Alpha == NULL) || (Ret->Sigma == NULL))
		MallocErr();

	Ret->EstAlpha = NULL;
	Ret->EstSigma = NULL;

	if(Opt->Analsis == ANALMCMC)
	{
		Ret->EstAlpha = (double*)malloc(sizeof(double) * Opt->Trees->NoOfSites);
		Ret->EstSigma = (double*)malloc(sizeof(double) * Opt->Trees->NoOfSites);
		Ret->AlphaErr = (double*)malloc(sizeof(double) * Opt->Trees->NoOfSites);

		Rates->Contrast = Ret;
		CalcContrast(Opt->Trees, Rates);
		CalcContLh(Opt, Opt->Trees, Rates);
	
		for(Index=0;Index<Opt->Trees->NoOfSites;Index++)
		{
			Ret->EstAlpha[Index] = Ret->Alpha[Index];	
			Ret->EstSigma[Index] = Ret->Sigma[Index];
		}
	}
	
	return Ret;
}

double	MLContrastPraxis(void* P, double *List)
{
	int			Index;
	PRAXSTATE	*PState;

	TREES		*Trees;
	OPTIONS		*Opt;
	RATES		*Rates;

	PState = (PRAXSTATE*)P;

	Trees	= PState->Trees;
	Opt		= PState->Opt;
	Rates	= PState->Rates;
	
	for(Index=0;Index<Trees->NoOfSites;Index++)
		Rates->Contrast->Sigma[Index] = List[Index];

	CalcContLh(PState->Opt, PState->Trees, PState->Rates);

	return -PState->Rates->Lh;
}

double	ChangeContrastRate(double Rate, double Dev, RANDSTATES *RS)
{
	double Ret;
	
	do
	{
		Ret = (GenRandState(RS) * Dev) - (Dev / 2.0); 
		Ret += Rate;
	} while(Ret <= 0);

	return Ret;
}

void	MutateContrastRates(OPTIONS *Opt, TREES* Trees, RATES* Rates)
{
	int			Index;
	CONTRASTR*	Con;
	double		Dev;

	Con = Rates->Contrast;

	if(GenRandState(Rates->RandStates) < 0.5)
	{
		Dev = Opt->RateDevList[0];

		for(Index=0;Index<Trees->NoOfSites;Index++)
			Con->EstAlpha[Index] += (GenRandState(Rates->RandStates) * Dev) - (Dev / 2.0);
	}	
	else
	{
		Dev = Opt->RateDevList[1];
		for(Index=0;Index<Trees->NoOfSites;Index++)
			Con->EstSigma[Index] = ChangeContrastRate(Con->EstSigma[Index], Dev, Rates->RandStates);
	} 
}

void	CopyContrastRates(RATES* R1, RATES* R2, int NoSites)
{
	int			Index;
	CONTRASTR	*C1, *C2;

	C1 = R1->Contrast;
	C2 = R2->Contrast;

	for(Index=0;Index<NoSites;Index++)
	{
		C1->Alpha[Index] = C2->Alpha[Index];
		C1->Sigma[Index] = C2->Sigma[Index];
	}

	if(C2->EstAlpha == NULL)
		return;

	for(Index=0;Index<NoSites;Index++)
	{
		C1->EstAlpha[Index] = C2->EstAlpha[Index];
		C1->EstSigma[Index] = C2->EstSigma[Index];
		C1->AlphaErr[Index] = C2->AlphaErr[Index];
	}
} 
