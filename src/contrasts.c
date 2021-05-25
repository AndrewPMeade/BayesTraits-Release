#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

#include "typedef.h"
#include "genlib.h"
#include "contrasts.h"
#include "praxis.h"

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

	
/* Common Varince */
//	C->Var = l1 + l2;
//	C->Var = C->Cont / sqrt(C->Var);

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
	GlobalVar = GlobalVar / NoCon;
	Rates->Lh = (NoCon+1) * log(6.28318530717958647692528676655900576839 * GlobalVar);
	Rates->Lh += SumLogVar + (T1 / GlobalVar);
	Rates->Lh *= -0.5;

	Rates->Contrast->Alpha[0] = T->Root->Contrast->Data;
	Rates->Contrast->Sigma[0] = GlobalVar;
}


double	CalcContrastLh(OPTIONS *Opt, TREES* Trees, RATES* Rates)
{
	int i;

	PrintTime(stdout);
	
	for(i=0;i<1000000;i++)
	{
		CalcContrast(Trees, Rates);
		CalcContLh(Opt, Trees, Rates);

		if(i%10000 == 0)
		{
			printf("i\t%d\n", i);
			fflush(stdout);
		}
	}

	PrintTime(stdout);
	return Rates->Lh;

	/*

	Sigma = ConRates->Sigma[0];

	Lh = 0;
	P = Trees->NoOfSites;
	Const = CalcContrastConst(Trees->NoOfSites, Sigma);

	for(Index=0;Index<Trees->NoOfNodes;Index++)
	{
		N = &T->NodeList[Index];
		if(N->Tip == FALSE)
		{
			C = N->Contrast;
			CR = N->Right->Contrast;
			CL = N->Left->Contrast;

			SumBL = N->Left->Length + CL->Err + N->Right->Length + CR->Err;
	
			Temp = C->Data * C->Data;

			Temp = Temp / (Sigma * Sigma * SumBL);
			Temp = exp(-(Temp * 0.5));
			

			T2 = pow(SumBL, Trees->NoOfSites/2.0) * Const;
			T2 = 1.0 / T2;

			Lh += log(Temp * T2);
		}
	}

	Rates->Lh = Lh; 
	return Lh; */
}


CONTRASTR*	AllocContrastRates(OPTIONS *Opt, RATES *Rates)
{
	CONTRASTR*	Ret;

	Ret = (CONTRASTR*)malloc(sizeof(CONTRASTR));
	if(Ret == NULL)
		MallocErr();

	Ret->Alpha = (double*)malloc(sizeof(double) * Opt->Trees->NoOfSites);
	Ret->Sigma = (double*)malloc(sizeof(double) * Opt->Trees->NoOfSites);

	if((Ret->Alpha == NULL) || (Ret->Sigma == NULL))
		MallocErr();
	
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

void	MLContrast(OPTIONS *Opt, TREES *Trees, RATES *Rates)
{
	PRAXSTATE	*PState;
	double		*TempV;
	int			Index;
	double		Sig;
	TREE		*Tree;
	
/*	double		Sig;

	Rates->TreeNo = 45;
	InitContinusTree(Opt, Trees, Rates->TreeNo);
	for(Sig = 0.0001;Sig<10;Sig+=0.01)
	{
		Rates->Contrast->Sigma[0] = Sig;
		CalcContrastLh(Opt, Trees, Rates);
		printf("%f\t%f\n", Sig, Rates->Lh);
	}

	exit(0);
*/

//	printf("=================== Tree %d ===================\n", Rates->TreeNo+1);

	Tree = &Trees->Tree[Rates->TreeNo];

	CalcContrast(Trees, Rates);
	CalcContLh(Opt, Trees, Rates);

	return;

	
//	PrintContrasts(Trees,Rates);
//	exit(0);

	TempV = (double*)malloc(sizeof(double) * Trees->NoOfSites);
	if(TempV == NULL) MallocErr();
	for(Index=0;Index<Trees->NoOfSites;Index++)
		TempV[Index] = Rates->Contrast->Sigma[Index];

	TempV[0] = 3; 
		
	PState = IntiPraxis(MLContrastPraxis, TempV, Trees->NoOfSites, 0, 1, 4, 1000);
//	PState = IntiPraxis(LhPraxis, TempV, Rates->NoOfRates, 0, 1, 4, 5000);

 	PState->Trees	= Trees;
	PState->Opt		= Opt;
	PState->Rates	= Rates;

	praxis(PState);

	for(Index=0;Index<Trees->NoOfSites;Index++)
		Rates->Contrast->Sigma[Index] = TempV[Index];

	CalcContrast(Trees, Rates);
	CalcContLh(Opt, Trees, Rates);
	FreePracxStates(PState);


	for(Sig=0.0;Sig<10;Sig+=0.01)
	{
		Rates->Contrast->Sigma[0] = Sig;
		CalcContLh(Opt, Trees, Rates);
	}
}