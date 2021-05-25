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


CONTRAST*	AllocContrast(NODE N)
{
	CONTRAST*	Ret;

	Ret = (CONTRAST*)malloc(sizeof(CONTRAST));
	if(Ret == NULL)
		MallocErr();

	if(N->Tip == TRUE)
		Ret->Data = N->Taxa->ConData[0];
	else
		Ret->Data = -1;

	Ret->Cont = 0;
	Ret->Var = 0;
	Ret->Err = 0;

	return Ret;
}

void	InitContrastTree(OPTIONS *Opt, TREES* Trees, int TNo)
{
	TREE *Tree;
	int NIndex;
	NODE N;

	Tree = &Trees->Tree[TNo];
	for(NIndex=0;NIndex<Trees->NoOfNodes;NIndex++)
	{
		N = &Tree->NodeList[NIndex];
		N->Contrast = AllocContrast(N);
	}

}

void	InitContrastAll(OPTIONS *Opt, TREES* Trees)
{
	int TIndex;

	for(TIndex=0;TIndex<Trees->NoOfTrees;TIndex++)
		InitContrastTree(Opt, Trees, TIndex);

}

void	FreeContrast(OPTIONS *Opt, TREES* Trees, int TreeNo)
{
	TREE *Tree;
	NODE N;
	int NIndex;

	Tree = &Trees->Tree[TreeNo];
	for(NIndex=0;NIndex<Trees->NoOfNodes;NIndex++)
	{
		N = &Tree->NodeList[NIndex];
		free(N->Contrast);
		N->Contrast = NULL;
	}
}

void	FreeAllContrast(OPTIONS *Opt, TREES* Trees)
{
	int TIndex;

	for(TIndex=0;TIndex<Trees->NoOfTrees;TIndex++)
		FreeContrast(Opt, Trees, TIndex);
}


void	RecCalcContrast(NODE N)
{
	CONTRAST	*C, *C1, *C2;
	double		t;
	double		l1, l2;

	if(N->Tip == TRUE)
		return;

	RecCalcContrast(N->Left);
	RecCalcContrast(N->Right);

	C = N->Contrast;
	C1 = N->Left->Contrast;
	C2 = N->Right->Contrast;

	l1 = N->Left->Length + C1->Err;
	l2 = N->Right->Length + C2->Err;

	t = (l1 * C2->Data) +  (l2 * C1->Data);
	t = t / (l1 + l2);

	C->Data = t;
	C->Cont = C1->Data - C2->Data;

	C->Err = (l1 * l2) / (l1 + l2);
	C->Var = l1 + l2;
}

void	RecPrintNode(NODE N)
{
	if(N->Tip == TRUE)
	{
		printf("%s ", N->Taxa->Name);
		return;
	}

	RecPrintNode(N->Left);
	RecPrintNode(N->Right);
}


void	PrintContrasts(TREES* Trees, RATES* Rates)
{
	TREE	*Tree;
	int		Index;
	NODE	N;
	CONTRAST	*C;


	Tree = &Trees->Tree[Rates->TreeNo];

	for(Index=0;Index<Trees->NoOfNodes;Index++)
	{
		N = &Tree->NodeList[Index];
		printf("Node Contrasts\t");
		RecPrintNode(N);
		printf("\t");
		C = N->Contrast;

		printf("%f\t%f\t%f\t%f\t", C->Cont, C->Data, C->Err, C->Var);

		printf("\n");
	}


}

void	CalcContrast(TREES* Trees, RATES* Rates)
{
	TREE	*Tree;

	Tree = &Trees->Tree[Rates->TreeNo];
	RecCalcContrast(Tree->Root);
}

void	CalcContLh(OPTIONS *Opt, TREES* Trees, RATES* Rates)
{
	int			Index;
	TREE		*T;
	NODE		N;
	CONTRAST	*Con;
	CONTRASTR	*ConRates;
	double		GlobalVar;
	double		SumLogVar;
	double		T1;
	int			NoCon;

	ConRates = Rates->Contrast;
	T = &Trees->Tree[Rates->TreeNo];

	NoCon = 0;
	GlobalVar = 0;
	SumLogVar = 0;

	for(Index=0;Index<Trees->NoOfNodes;Index++)
	{
		N = &T->NodeList[Index];
		if(N->Tip == FALSE)
		{
			Con = N->Contrast;
			GlobalVar += (Con->Cont * Con->Cont) / Con->Var;

			SumLogVar += log(Con->Var);
			NoCon++;
		}
	}

	SumLogVar += log(T->Root->Contrast->Err);

	

	T1 = GlobalVar;
	GlobalVar = GlobalVar / (NoCon+1);
	Rates->Lh = (NoCon+1) * log(6.28318530717958647692528676655900576839 * GlobalVar);
	Rates->Lh += SumLogVar + (T1 / GlobalVar);
	Rates->Lh *= -0.5;

	Rates->Contrast->Alpha[0] = T->Root->Contrast->Data;
	Rates->Contrast->Sigma[0] = GlobalVar;
}

void CalcContrastMCMC(OPTIONS *Opt, TREES* Trees, RATES* Rates)
{
	int			Index;
	TREE		*T;
	NODE		N;
	CONTRAST	*Con;
	CONTRASTR	*ConRates;
	double		GlobalVar;
	double		SumLogVar;
	double		T1, T2;
	int			NoCon;

	ConRates = Rates->Contrast;
	T = &Trees->Tree[Rates->TreeNo];

	NoCon = 0;
	GlobalVar = 0;
	SumLogVar = 0;

	for(Index=0;Index<Trees->NoOfNodes;Index++)
	{
		N = &T->NodeList[Index];
		if(N->Tip == FALSE)
		{
			Con = N->Contrast;
			GlobalVar += (Con->Cont * Con->Cont) / Con->Var;

			SumLogVar += log(Con->Var);
			NoCon++;
		}
	}

	SumLogVar += log(T->Root->Contrast->Err);

	T1 = GlobalVar;
	GlobalVar = GlobalVar / NoCon;
	Rates->Lh = (NoCon+1) * log(6.28318530717958647692528676655900576839 * ConRates->EstSigma[0]);

	T2 = ((NoCon+1) * GlobalVar) + ConRates->AlphaErr;
	Rates->Lh += SumLogVar + (T2 / ConRates->EstSigma[0]);
	Rates->Lh *= -0.5;

	Rates->Contrast->Alpha[0] = T->Root->Contrast->Data;
	Rates->Contrast->Sigma[0] = GlobalVar;
}

double	CalcMCMCAlpha(NODE N, double Alpha)
{
	double Ret;
	CONTRAST	*Con;

	Con = N->Contrast;

	Ret = (Alpha - Con->Data) * (Alpha - Con->Data);
	Ret = Ret / Con->Err;	

	return Ret;
}

double	CalcContrastLh(OPTIONS *Opt, TREES* Trees, RATES* Rates)
{
	CONTRASTR	*Con;

	Con = Rates->Contrast;

	if(Opt->UsePhyloPlasty == TRUE)
		Plasty(Opt, Trees, Rates);

	CalcContrast(Trees, Rates);	
	
//	Use ML values for MCMC
	if(Opt->Analsis == ANALML)
	{
		CalcContLh(Opt, Trees, Rates);
		if(Opt->Analsis != ANALMCMC)
		{
			Con->EstAlpha[0] = Con->Alpha[0];
			Con->EstSigma[0] = Con->Sigma[0];
		}
		return Rates->Lh;
	}

	Con->AlphaErr = CalcMCMCAlpha(Trees->Tree[0].Root, Con->EstAlpha[0]);
	CalcContrastMCMC(Opt, Trees, Rates);

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
	}

}