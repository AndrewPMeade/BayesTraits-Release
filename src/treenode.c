#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

#include "typedef.h"
#include "trees.h"
#include "genlib.h"
/*
void	SetNodePartition(NODE N, int* List, int *No)
{
	int Index;

	if(N->Tip == TRUE)
	{
		if(List != NULL)
			List[*No] = N->TipID;
		(*No)++;
		return;
	}

	for(Index=0;Index<N->NoNodes;Index++)
		SetNodePartition(N->NodeList[Index], List, No);
}

int PartComp(const void *av, const void *bv)
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

void	SetPart(NODE N)
{
	int	No;

	if(N->Tip == TRUE)
	{
		N->Part = NULL;
		N->PSize= 0;
		return;
	}

	N->PSize = 0;
	
	SetNodePartition(N, NULL, &N->PSize);

	N->Part = (int*)malloc(sizeof(int) * (N->PSize));
	if(N->Part == NULL)
		MallocErr();

	No = 0;
	SetNodePartition(N, N->Part, &No);

	qsort(N->Part, N->PSize, sizeof(int), PartComp);
}

int		IsPartEqual(int *Part1, int Len1, int* Part2, int Len2)
{
	int	Index;

	if(Len1 != Len2)
		return FALSE;

	for(Index=0;Index<Len1;Index++)
		if(Part1[Index] != Part2[Index])
			return FALSE;

	return TRUE;
}

void	FreePartitions(TREES *Trees)
{
	int	TIndex;
	int	NIndex;
	TREE	*Tree;
	NODE	N;

	for(TIndex=0;TIndex<Trees->NoOfTrees;TIndex++)
	{
		Tree = Trees->Tree[TIndex];

		for(NIndex=0;NIndex<Tree->NoNodes;NIndex++)
		{
			N = Tree->NodeList[NIndex];
			if(N->Part != NULL)
				free(N->Part);
			N->Part = NULL;
			N->PSize = 0;
		}
	}
}

void	SetPartitions(TREES *Trees)
{
	int	TIndex;
	int	NIndex;
	TREE	*Tree;
	NODE	N;

	for(TIndex=0;TIndex<Trees->NoOfTrees;TIndex++)
	{
		Tree = Trees->Tree[TIndex];

		for(NIndex=0;NIndex<Tree->NoNodes;NIndex++)
		{
			N = Tree->NodeList[NIndex];
			if(N->Tip == FALSE)
			{
				if(N->Part != NULL)
					free(N->Part);
				SetPart(N);
			}
		}
	}	
}

int		IsPartitionSumSet(int *Part, int Len, int* SubSet, int SubSetLen)
{
	int	Index;
	int	SIndex;

	if(SubSetLen > Len)
		return FALSE;

	SIndex=0;

	for(Index=0;Index<Len;Index++)
	{
		if(Part[Index] == SubSet[SIndex])
			SIndex++;

		if(SIndex == SubSetLen)
			return TRUE;
	}


	return FALSE;
}

NODE	FindNode(RECNODE RNode, TREE *Tree, int *Depth, int NoOfNodes)
{
	NODE	Ret;
	int		NIndex;
	NODE	TempNode;
	
	Ret = NULL;

	for(NIndex=0;NIndex<NoOfNodes;NIndex++)
	{
		TempNode = Tree->NodeList[NIndex];
		
		if(IsPartEqual(TempNode->Part, TempNode->PSize, RNode->TaxaID, RNode->NoOfTaxa) == TRUE)
		{
			*Depth = TempNode->PSize;
			return TempNode;
		}

		if(RNode->NodeType != NODEREC)
		{
			if(IsPartitionSumSet(TempNode->Part, TempNode->PSize, RNode->TaxaID, RNode->NoOfTaxa) == TRUE)
			{
				if(Ret == NULL)
					Ret = TempNode;
				else
				{
					if(Ret->PSize > TempNode->PSize)
						Ret = TempNode;
				}
			}
		}
	}
	if(Ret != NULL)
		*Depth = Ret->PSize;
	else
		*Depth = 0;
	return Ret;
}

void	SetTaxaIDList(RECNODE RNode)
{
	int	Index;

	if(RNode->TaxaID != NULL)
		free(RNode->TaxaID);

	RNode->TaxaID = (int*)malloc(sizeof(int) * RNode->NoOfTaxa);
	if(RNode->TaxaID == NULL)
		MallocErr();

	for(Index=0;Index<RNode->NoOfTaxa;Index++)
		RNode->TaxaID[Index] = RNode->Taxa[Index]->No;
	
	qsort(RNode->TaxaID, RNode->NoOfTaxa, sizeof(int), PartComp);
}

void	SetRecNodes(RECNODE RNode, OPTIONS *Opt)
{
	NODE	N;
	TREES	*Trees;
	int		TIndex;
	int		Depth;

	SetTaxaIDList(RNode);

	Trees = Opt->Trees;
	RNode->Hits = 0;
	for(TIndex=0;TIndex<Trees->NoOfTrees;TIndex++)
	{
		Depth = 0;
		N = FindNode(RNode, Trees->Tree[TIndex], &Depth, Trees->Tree[TIndex]->NoNodes);
		RNode->TreeNodes[TIndex] = N;

		if(N != NULL)
		{
			if((RNode->NodeType == MRCA) ||(RNode->NodeType == FOSSIL))
				RNode->Hits += Depth;

			if(RNode->NodeType == NODEREC)
				RNode->Hits++;
		}
	}
}

 */