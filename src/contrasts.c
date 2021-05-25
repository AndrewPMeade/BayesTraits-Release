#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

#include "typedef.h"
#include "genlib.h"
#include "contrasts.h"

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

void	InitContrast(OPTIONS *Opt, TREES* Trees)
{
	TREE *Tree;
	NODE N;
	int TIndex;
	int NIndex;

	for(TIndex=0;TIndex<Trees->NoOfTrees;TIndex++)
	{
		Tree = &Trees->Tree[TIndex];
		for(NIndex=0;NIndex<Trees->NoOfNodes;NIndex++)
		{
			N = &Tree->NodeList[NIndex];
			N->Contrast = AllocContrast(N);
		}
	}
}


void	FreeContrast(OPTIONS *Opt, TREES* Trees)
{
	TREE *Tree;
	NODE N;
	int TIndex;
	int NIndex;

	for(TIndex=0;TIndex<Trees->NoOfTrees;TIndex++)
	{
		Tree = &Trees->Tree[TIndex];
		for(NIndex=0;NIndex<Trees->NoOfNodes;NIndex++)
		{
			N = &Tree->NodeList[NIndex];
			free(N->Contrast);
		}
	}
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

void	CalcContrast(TREES* Trees, RATES* Rates)
{
	TREE *Tree;
	double	Alpha;
	int		i;

//	Tree = &Trees->Tree[Rates->TreeNo];
	Tree = &Trees->Tree[0];

	PrintTime(stdout);
	printf("\n\n");

	RecCalcContrast(Tree->Root);


	Alpha = Tree->Root->Contrast->Data;

/*
	printf("%d\tRoot Val\t%f\t%f\n", i, Alpha, Alpha + Tree->Root->Contrast->Var);
	printf("\n\n");
	PrintTime(stdout);
	exit(0); */
}



double	CalcContrastLh(TREES* Trees, RATES* Rates)
{
	double	Lh;



}