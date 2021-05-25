#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "part.h"
#include "typedef.h"
#include "genlib.h"
#include "continuous.h"
#include "trees.h"

void	FreePart(PART *Part)
{
	if(Part != NULL)
	{
		free(Part->Taxa);
		free(Part);
	}
}

PART*	CreatPart(int NoTaxa)
{
	PART *Part;

	Part = (PART*)malloc(sizeof(PART));
	if(Part == NULL)
		MallocErr();

	Part->NoTaxa = NoTaxa;
	Part->Taxa = (int*)malloc(sizeof(int) * NoTaxa);
	if(Part->Taxa == NULL)
		MallocErr();
	
	return Part;
}

int		TaxaInPart(int TaxaIndex, PART *P)
{
	int Size, Index;
	int *TList;

	TList = P->Taxa;
	Size = P->NoTaxa;

	for(Index=0;Index<Size;Index++)
	{
		if(TList[Index] == TaxaIndex)
			return TRUE;

		if(TList[Index] > TaxaIndex)
			return FALSE;
	}

	return FALSE;
}

int		PartEqual(PART *A, PART *B)
{
	int Index;

	if(A->NoTaxa != B->NoTaxa)
		return FALSE;

	for(Index=0;Index<A->NoTaxa;Index++)
	{
		if(A->Taxa[Index] != B->Taxa[Index])
			return FALSE;
	}

	return TRUE;
}
/* is B a subset of A */
int		PartSubSet(PART *A, PART *B)
{
	int Index, SIndex;

	SIndex = 0;
	if(B->NoTaxa > A->NoTaxa)
		return FALSE;

	for(Index=0;Index<A->NoTaxa;Index++)
	{
		if(A->Taxa[Index] == B->Taxa[SIndex])
			SIndex++;

		if(SIndex == B->NoTaxa)
			return TRUE;
	}

	return FALSE;
}



int PartCompID(const void *av, const void *bv)
{
	int *a, *b;

	a = (int*)av;
	b = (int*)bv;

	if(*a > *b)
		return 1;

	if(*a < *b)
		return -1;

	return 0;
}
/*
void	PrintPart(PART *Part)
{
	int Index;

	printf("%d\t[", Part->NoTaxa);

	for(Index=0;Index<Part->NoTaxa-1;Index++)
		printf("%d,", Part->Taxa[Index]);
	printf("%d]\n", Part->Taxa[Index]);
}
*/
void	SetIntPart(NODE N)
{
	int NIndex, PIndex, Pos;
	PART *Part, *DPart;

	Part = N->Part;

	Pos = 0;
	for(NIndex=0;NIndex<N->NoNodes;NIndex++)
	{
		DPart = N->NodeList[NIndex]->Part;
		for(PIndex=0;PIndex<DPart->NoTaxa;PIndex++)
			Part->Taxa[Pos++] = DPart->Taxa[PIndex];
	}

	qsort(Part->Taxa, Part->NoTaxa, sizeof(int), PartCompID);
}

void	SetTreePart(NODE N)
{
	int Index, NoTaxa;

	if(N->Tip == TRUE)
	{
		N->Part = CreatPart(1);
		N->Part->Taxa[0] = N->Taxa->No;
		return;
	}

	NoTaxa = 0;
	for(Index=0;Index<N->NoNodes;Index++)
	{
		SetTreePart(N->NodeList[Index]);
		NoTaxa += N->NodeList[Index]->Part->NoTaxa;
	}

	N->Part = CreatPart(NoTaxa);
	SetIntPart(N);
}

void	FreeParts(TREES *Trees)
{
	int TIndex, NIndex;
	TREE *T;

	for(TIndex=0;TIndex<Trees->NoTrees;TIndex++)
	{
		T = Trees->Tree[TIndex];

	
		for(NIndex=0;NIndex<T->NoNodes;NIndex++)
		{
			FreePart(T->NodeList[NIndex]->Part);
			T->NodeList[NIndex]->Part = NULL;
		}
	}
}

void	SetParts(TREES *Trees)
{
	int TIndex;

	for(TIndex=0;TIndex<Trees->NoTrees;TIndex++)
		SetTreePart(Trees->Tree[TIndex]->Root);
}



void	GetPartDiff(PART *Ans, PART *Cur, PART *Diff)
{
	int Index;
	//TaxaInPart

	Diff->NoTaxa = 0;

	for(Index=0;Index<Ans->NoTaxa;Index++)
	{
		if(TaxaInPart(Ans->Taxa[Index], Cur) == FALSE)
			Diff->Taxa[Diff->NoTaxa++] = Ans->Taxa[Index];
	}
}

void	PrintPart(FILE *Str, TREES *Trees, PART *Part)
{
	int Index, ID;
	TAXA *T;

	fprintf(Str, "Part:\t%d\t", Part->NoTaxa);
	for(Index=0;Index<Part->NoTaxa;Index++)
	{
		ID = Part->Taxa[Index];
		T = Trees->Taxa[ID];
		fprintf(Str, "%s\t", T->Name);
	}
}

void	PrintPartTaxaOnly(FILE *Str, TREES *Trees, PART *Part)
{
	int Index, ID;
	TAXA *T;

	for(Index=0;Index<Part->NoTaxa;Index++)
	{
		ID = Part->Taxa[Index];
		T = Trees->Taxa[ID];
		fprintf(Str, "%s\t", T->Name);
	}
}


PART*	CreatePart(TREES *Trees, int NoTaxa, char **TaxaList)
{
	PART	*Ret;
	TAXA	*Taxa;
	int		Index;

	Ret = CreatPart(NoTaxa);

	for(Index=0;Index<NoTaxa;Index++)
	{
		Taxa = GetTaxaFromName(TaxaList[Index], Trees->Taxa, Trees->NoTaxa);

		if(Taxa == NULL)
		{
			printf("Invalid taxa name:\t%s\n", TaxaList[Index]); 
			exit(0);
		}
		
		Ret->Taxa[Index] = Taxa->No;
	}

	qsort(Ret->Taxa, Ret->NoTaxa, sizeof(int), PartCompID);

	for(Index=0;Index<NoTaxa-1;Index++)
	{
		if(Ret->Taxa[Index] == Ret->Taxa[Index+1])
		{
			printf("Taxa %s included multiple times.\n", Trees->Taxa[Ret->Taxa[Index]]->Name);
			exit(0);
		}
	}

	return Ret;
}

NODE	PartGetMRCA(TREE *Tree, PART *Part)
{
	NODE Ret, N;
	int Index;

	Ret = Tree->Root;

	for(Index=0;Index<Tree->NoNodes;Index++)
	{
		N = Tree->NodeList[Index];
		if(PartSubSet(N->Part, Part) == TRUE)
		{
			if(N->Part->NoTaxa < Ret->Part->NoTaxa)
				Ret = N;
		}
	}

	return Ret;
}