/*
*  BayesTriats 4.0
*
*  copyright 2022
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

#include "Part.h"
#include "TypeDef.h"
#include "GenLib.h"
#include "Continuous.h"
#include "Trees.h"

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

	Part = (PART*)SMalloc(sizeof(PART));

	Part->NoTaxa = NoTaxa;
	Part->Taxa = (int*)SMalloc(sizeof(int) * NoTaxa);
	
	Part->PartID = -1;
	Part->Freq = 0;
	Part->Prob = 0.0;

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

size_t	GetMaxNoParts(TREES *Trees)
{
	size_t Ret;
	int TIndex;

	Ret = 0;

	for(TIndex=0;TIndex<Trees->NoTrees;TIndex++)
		Ret += Trees->Tree[TIndex]->NoNodes;
	
	return Ret;
}

void	FreeTreeParts(TREES *Trees)
{
	size_t Index;

	if(Trees->PartList == NULL)
		return;

	for(Index=0;Index<Trees->NoParts;Index++)
		FreePart(Trees->PartList[Index]);

	free(Trees->PartList);
	Trees->PartList = NULL;
	Trees->NoParts = 0;
}

PART*	PartInList(size_t NoParts, PART **PList, PART *Part)
{
	size_t Index;

	for(Index=0;Index<NoParts;Index++)
	{
		if(PartEqual(Part, PList[Index]) == TRUE)
			return PList[Index];
	}

	return NULL;
}

void	ProcessTreeParts(size_t *NoParts, PART **PList, TREE *Tree)
{
	int NIndex;
	NODE N;
	PART *ListP;

	for(NIndex=0;NIndex<Tree->NoNodes;NIndex++)
	{
		N = Tree->NodeList[NIndex];

		ListP = PartInList(*NoParts, PList, N->Part);
		if(ListP != NULL)
			FreePart(N->Part);
		else
		{
			ListP = N->Part;
			PList[*NoParts] = ListP;
			(*NoParts)++;
		}
		
		ListP->Freq++;
		N->Part = ListP;
	}
}

void	SetPartInfo(TREES *Trees)
{
	int TIndex;
	size_t MaxNodes, PIndex;
	

	MaxNodes = GetMaxNoParts(Trees);
	Trees->PartList = (PART**)SMalloc(sizeof(PART*) * MaxNodes);

	Trees->NoParts = 0;
	for(TIndex=0;TIndex<Trees->NoTrees;TIndex++)
		ProcessTreeParts(&Trees->NoParts, Trees->PartList, Trees->Tree[TIndex]);

	for(PIndex=0;PIndex<Trees->NoParts;PIndex++)
	{
		Trees->PartList[PIndex]->PartID = PIndex;
		Trees->PartList[PIndex]->Prob = Trees->PartList[PIndex]->Freq / (double)Trees->NoTrees;
	}
}


void	SetParts(TREES *Trees)
{
	int TIndex;

	FreeTreeParts(Trees);

	for(TIndex=0;TIndex<Trees->NoTrees;TIndex++)
		SetTreePart(Trees->Tree[TIndex]->Root);

	SetPartInfo(Trees);

//	PrintParts(stdout, Trees);
//	exit(0);

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

	fprintf(Str, "%zu\t%d\t%f\t", Part->PartID, Part->Freq, Part->Prob);
	fprintf(Str, "\t%d\t", Part->NoTaxa);
	
	for(Index=0;Index<Part->NoTaxa;Index++)
	{
		ID = Part->Taxa[Index];
		T = Trees->Taxa[ID];
		fprintf(Str, "%s\t", T->Name);
	}
}

void	PrintParts(FILE *Str, TREES *Trees)
{
	int PIndex;

	for(PIndex=0;PIndex<Trees->NoParts;PIndex++)
	{
		PrintPart(Str, Trees, Trees->PartList[PIndex]);
		fprintf(Str, "\n");
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
/*
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
*/

NODE	PartGetMRCA(TREE *Tree, PART *Part)
{
	NODE Node;

	Node = GetTreeTaxaNode(Tree, Part->Taxa[0]);

	if(Part->NoTaxa == 1)
		return Node;

	while(Node != NULL)
	{
		if(PartSubSet(Node->Part, Part) == TRUE)
			return Node;

		Node = Node->Ans;
	}

	return NULL;
}


PART*	GetPartFromID(size_t ID, TREES *Trees)
{
	int Index;


	for(Index=0;Index<Trees->NoParts;Index++)
		if(Trees->PartList[Index]->PartID == ID)
			return Trees->PartList[Index];

	return NULL;
}