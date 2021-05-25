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
	int TIndex;
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
	C->Var = l1 + l2;

	C->Var = C->Cont / sqrt(C->Var);

	C->Err = (l1 * l2) / (l1 + l2);
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
	double	Alpha;
	int		i;

	Tree = &Trees->Tree[Rates->TreeNo];
//	Tree = &Trees->Tree[0];

//	PrintTime(stdout);
//	printf("\n\n");

	RecCalcContrast(Tree->Root);
/*	Alpha = Tree->Root->Contrast->Data;
	printf("%d\tRoot Val\t%f\t%f\n", i, Alpha, Alpha + Tree->Root->Contrast->Var);
	printf("\n\n");
	PrintTime(stdout);
	exit(0);  */
}

double	CalcContrastConst(double P, double Sig)
{
	return pow(2 * 3.141592653589793238462643383279502884197, P / 2.0) * pow(Sig, P);
}


void	CalcContLh(OPTIONS *Opt, TREES* Trees, RATES* Rates)
{
	double		Lh;
	int			Index;
	TREE		*T;
	NODE		N;
	CONTRAST	*C, *CL, *CR;
	double		Temp, T2;
	CONTRASTR	*ConRates;
	double		Sigma;
	double		SumBL, P, Const;

	ConRates = Rates->Contrast;
	T = &Trees->Tree[Rates->TreeNo];

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
			Temp = exp(-0.5 * Temp);
			

			T2 = pow(SumBL, Trees->NoOfSites/2.0) * Const;
			T2 = 1.0 / T2;

			Lh += log(Temp * T2);
		}
	}

	Rates->Lh = Lh; 
}

void	RecPrintConLh(OPTIONS *Opt, TREES* Trees, RATES* Rates)
{
	double		Lh;
	int			Index;
	TREE		*T;
	NODE		N;
	CONTRAST	*C, *CL, *CR;
	double		Temp, T2;
	CONTRASTR	*ConRates;
	double		Sigma;
	double		SumBL, P, Const;

	ConRates = Rates->Contrast;
	T = &Trees->Tree[Rates->TreeNo];

	Sigma = ConRates->Sigma[0];


	Lh = 0;
	P = (double)Trees->NoOfSites;
	Const = CalcContrastConst(Trees->NoOfSites, Sigma);

	printf("Sigma\t%f\n", Sigma);
	printf("P\t%f\n", P);
	printf("Constatn (2Pi)^(p/2)*Sigma^P\t%f\n", Const);

	for(Index=0;Index<Trees->NoOfNodes;Index++)
	{
		N = &T->NodeList[Index];
		if(N->Tip == FALSE)
		{
			printf("Node Lh\t");
			RecPrintNode(N);
			printf("\t");
			
			C = N->Contrast;
			CR = N->Right->Contrast;
			CL = N->Left->Contrast;

			SumBL = N->Left->Length + CL->Err + N->Right->Length + CR->Err;

			printf("SumBL\t%f\t", SumBL);
	
			Temp = C->Data * C->Data;

			printf("Data\t%f\tData^2\t%f\t", C->Data, Temp);

			Temp = Temp / (Sigma * Sigma * SumBL);
			printf("Data/(Sig^2 * SumBL)\t%f\t", Temp);

			Temp = exp(-0.5 * Temp);

			printf("T1 = Exp[All]\t%f\t", Temp);
		
			T2 = pow(SumBL, Trees->NoOfSites/2.0) * Const;
			printf("T2 = Const * (SumBL^p/2)\t%f\t", T2);

			T2 = 1.0 / T2;

			printf("T2 = 1 / T2\t%f\t", T2);

			printf("T2 * T1\t%f\t", Temp * T2);
			printf("log(T2 * T1)\t%f\t", log(Temp * T2));

			Lh += log(Temp * T2);

			printf("\n");
		}
	}

	Rates->Lh = Lh; 
}

double	CalcContrastLh(OPTIONS *Opt, TREES* Trees, RATES* Rates)
{
	CalcContrast(Trees, Rates);
	CalcContLh(Opt, Trees, Rates);
	
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

	RecPrintConLh(Opt, Trees, Rates);
	
	printf("Contrast Lh Res\t%d\t%f\t%f\t%f", Rates->TreeNo+1, Rates->Lh, Tree->Root->Contrast->Data,  Rates->Contrast->Sigma[0]);

	for(Sig=0.0;Sig<10;Sig+=0.01)
	{
		Rates->Contrast->Sigma[0] = Sig;
		CalcContLh(Opt, Trees, Rates);
		printf("%f\t%f\n", Sig, Rates->Lh);
	}

	exit(0);

//	printf("\n\n\n\n\n\n");
//	exit(0);
}