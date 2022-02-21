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

#include "TypeDef.h"
#include "GenLib.h"
#include "PTrees.h"
#include "Trees.h"
#include "Part.h"

typedef struct
{
	int		Score;

	NODE	*TList;
	int		NoTList;

	NODE	*PList;
	int		NoPList;
} PTINFO;


void	FindParallelPoints(TREES *Trees, TREE *Tree, int Cores);
void	SetFlatternedNodes(TREES *Trees, TREE *T);
void	SetPNodes(TREE *Tree);

void	SetPTrees(OPTIONS *Opt, TREES *Trees)
{
	int TIndex;
	TREE *Tree;

	 for(TIndex=0;TIndex<Trees->NoTrees;TIndex++)
	 {
		 Tree = Trees->Tree[TIndex];
	 	 SetFlatternedNodes(Trees, Tree);
		 FindParallelPoints(Trees, Tree, Opt->Cores);
//		FindParallelPoints(Trees, Tree, 4);
		 SetPNodes(Tree);
	 }		
}


 void	InitPNodes(TREE *T)
 {
	 int Index;
	 NODE N;

	 for(Index=0;Index<T->NoNodes;Index++)
	 {
		 N = T->NodeList[Index];
		 N->Visited = N->Tip;
	 }
 }

 /* 
	NODE		**PNodes;
	int			*NoPNodes;
	int			NoPGroups;
*/

 int	ValidPNode(NODE N)
 {
	 int Index;

	if(N->Visited == TRUE)
		return FALSE;

	 for(Index=0;Index<N->NoNodes;Index++)
	 {
		if(N->NodeList[Index]->Visited == FALSE)
			return FALSE;
	 }
	 return TRUE;
 }

 NODE*	GetNextPList(TREE *T, int *Size)
 {
	 NODE	*Ret, *List, N;
	 int	No, Index;

	 List = (NODE*)malloc(sizeof(NODE) * T->NoNodes);
	 if(List == NULL)
		 MallocErr();

	 No = 0;
	 for(Index=0;Index<T->NoNodes;Index++)
	 {
		 N = T->NodeList[Index]; 

		 if(ValidPNode(N) == TRUE)
			 List[No++] = N;
	 }

	 if(No != 0)
	 {
		 Ret = (NODE*)malloc(sizeof(NODE) * No);
		 if(Ret == NULL)
			 MallocErr();

		 for(Index=0;Index<No;Index++)
		 {
			 Ret[Index] = List[Index];
			 Ret[Index]->Visited = TRUE;
		 }
	 }
	 else
		 Ret = NULL;

	 free(List);

	 *Size = No;

	 return Ret;
 }

 void	RecPrintNodeBelow(NODE N)
 {
	 int Index;
	 if(N->Tip == TRUE)
	 {
		 printf("%s,", N->Taxa->Name);
		 return;
	 }

	 for(Index=0;Index<N->NoNodes;Index++)
		 RecPrintNodeBelow(N->NodeList[Index]);
 }

 void	PrintPNodes(TREE *T)
 {
	 int Index,j;

	 printf("No Groups\t%d\n", T->NoParallelGroups);
	 printf("Size\t");
	 for(Index=0;Index<T->NoParallelGroups;Index++)
		 printf("%d\t", T->ParallelGroupSize[Index]);
	 printf("\n");

	 for(Index=0;Index<T->NoParallelGroups;Index++)
	 {
		 printf("%d\t%d\n", Index, T->ParallelGroupSize[Index]);
		 for(j=0;j<T->ParallelGroupSize[Index];j++)
		 {
			printf("\t%d\t", T->ParallelNodes[Index][j]->Part->NoTaxa);
			RecPrintNodeBelow(T->ParallelNodes[Index][j]);
			printf("\n");
		 }
	 }
 }

 void	SetFlatternedNodes(TREES *Trees, TREE *T)
 {
	 int	*NoNodes;
	 int	NoGroups, Size;
	 NODE	**PNodes, *PN;


	 NoNodes = (int*)malloc(sizeof(int) * T->NoNodes);
	 PNodes = (NODE**)malloc(sizeof(NODE*) * T->NoNodes);

	 if((NoNodes == NULL) || (PNodes == NULL))
		 MallocErr();

	  NoGroups = 0;
	
	  InitPNodes(T);

	 do
	 {
		 PN = GetNextPList(T, &Size);
		 if(PN != NULL)
		 {
			 NoNodes[NoGroups] = Size;
			 PNodes[NoGroups] = PN;
			 NoGroups++;
		 }
	 }while(PN != NULL);

	 T->NoParallelGroups = NoGroups;
	 T->ParallelGroupSize = (int*)SMalloc(sizeof(int) * NoGroups);
	 T->ParallelNodes = (NODE**)SMalloc(sizeof(NODE*) * NoGroups);

	 memcpy(T->ParallelGroupSize, NoNodes, sizeof(int) * NoGroups);
	 memcpy(T->ParallelNodes, PNodes, sizeof(NODE*) * NoGroups);
	 
	 free(PNodes);
	 free(NoNodes);
 }

 void	SetTreePNodes(TREES *Trees)
 {
//	 int TIndex;
//	 TREE *Tree;

	 printf("Info\n");
 }

 NODE	FindNodeOfSize(TREE *Tree, int SizeNode)
 {
	 int Index;

	 for(Index=0;Index<Tree->NoNodes;Index++)
	 {
		 if(Tree->NodeList[Index]->Part->NoTaxa == SizeNode)
			 return Tree->NodeList[Index];
	 }

	 return NULL;
 }


 /* 
 typedef struct
{
	int		Score;

	NODE	*NList;
	int		NoNList;

	NODE	*CList;
	int		NoCList;
} PTINFO;
*/
 PTINFO*	AllocPTInfo(TREE *Tree)
 {
	 PTINFO* Ret;
	 int	 Index;

	 Ret = (PTINFO*)malloc(sizeof(PTINFO));
	 if(Ret == NULL)
		 MallocErr();

	 Ret->Score = -1;
	 Ret->NoTList = 0;
	 Ret->NoPList = 0;

	 Ret->PList = (NODE*)malloc(sizeof(NODE) * Tree->NoNodes);
	 Ret->TList = (NODE*)malloc(sizeof(NODE) * Tree->NoNodes);

	 if((Ret->PList == NULL) || (Ret->TList == NULL))
		 MallocErr();

	 for(Index=0;Index<Tree->NoNodes;Index++)
	 {
		 Ret->PList[Index] = NULL;
		 Ret->TList[Index] = NULL;
	 }

	 return Ret;
 }

 void	FreePTInfo(PTINFO *PTInfo)
 {
	 free(PTInfo->PList);
	 free(PTInfo->TList);
	 
	 free(PTInfo);
 }


 void	SetBirdPoints(TREES *Trees, TREE *Tree)
 {
	 Tree->NoPNodes = 4;

	 Tree->PNodes = (NODE*)malloc(sizeof(NODE) * Tree->NoPNodes);
	 if(Tree->PNodes == NULL)
		 MallocErr();

	 Tree->PNodes[0] = FindNodeOfSize(Tree, 3810);
	 Tree->PNodes[1] = FindNodeOfSize(Tree, 1146);
	 Tree->PNodes[2] = FindNodeOfSize(Tree, 980);
	 Tree->PNodes[3] = FindNodeOfSize(Tree, 816);
}

void	RecSetTPList(NODE N, PTINFO* Info, int MinSize)
{
	int Index;
	
	if(N->Tip == TRUE)
		return;

	if((N->Visited == FALSE) && (N->Part->NoTaxa >= MinSize))
		Info->TList[Info->NoTList++] = N;

	for(Index=0;Index<N->NoNodes;Index++)
		RecSetTPList(N->NodeList[Index], Info, MinSize);
}


NODE	GetNextPNode(PTINFO *PTInfo, int TSize, int Last)
{
	int Index, Score, NScore;
	NODE N, Ret;


	Ret = NULL;
	for(Index=0;Index<PTInfo->NoTList;Index++)
	{
		N = PTInfo->TList[Index];
		if(Ret == NULL)
		{
			Ret = N;
			if(Last == TRUE)
				Score = Ret->Part->NoTaxa;
			else
				Score = abs(Ret->Part->NoTaxa - TSize);
		}
		else
		{
			if(Last == TRUE)
			{
				if(N->Part->NoTaxa > Score)
				{
					Ret = N;
					Score = N->Part->NoTaxa;
				}
			}
			else
			{
				NScore = abs(N->Part->NoTaxa - TSize);
				if(NScore < Score)
				{
					Ret = N;
					Score = NScore;
				}
			}
		}
	}

	return Ret;
}

void	SetSelectedPNode(NODE N)
{
	SetVisitedNode(N, TRUE);
	while(N != NULL)
	{
		N->Visited = TRUE;
		N = N->Ans;
	}
}

void	PrintPPartInfo(TREES *Trees, TREE *Tree)
{
	int PTaxa,Index;
	double Pct;
	NODE N;

	PTaxa = 0;
	for(Index=0;Index<Tree->NoPNodes;Index++)
	{
		N = Tree->PNodes[Index];
		PTaxa += N->Part->NoTaxa;
		
		PrintPart(stdout, Trees, N->Part);
		printf("\n");
	}

	Pct = (double)PTaxa / Trees->NoTaxa;
	Pct = Pct * 100;

	printf("No P Taxa\t%f\n", Pct);
	printf("No S Taxa\t%f\n", 100.0 - Pct);
	fflush(stdout);
}

void	FindParallelPoints(TREES *Trees, TREE *Tree, int Cores)
{
	PTINFO *PTInfo;
	int		Index, OptNSize;
	NODE	NN;
	int		Last;

//	Cores = 4;
//	Cores = 1;
	OptNSize = Trees->NoTaxa / Cores;

	PTInfo = AllocPTInfo(Tree);

	SetVisitedTree(Tree, FALSE);

	RecSetTPList(Tree->Root, PTInfo, 2);

	Last = FALSE;
	for(Index=0;Index<Cores;Index++)
	{
		if(Index == Cores-1)
			Last = TRUE;

		NN = GetNextPNode(PTInfo, OptNSize, Last);
		if(NN == NULL)
			Index = Cores;
		else
		{

			PTInfo->PList[PTInfo->NoPList++] = NN;

			SetSelectedPNode(NN);
			PTInfo->NoTList = 0;
			RecSetTPList(Tree->Root, PTInfo, 2);
		}
	}
	
	Tree->PNodes = (NODE*)malloc(sizeof(NODE) * PTInfo->NoPList);
	if(Tree->PNodes == NULL)
		MallocErr();

	memcpy(Tree->PNodes, PTInfo->PList, sizeof(NODE) * PTInfo->NoPList);
	Tree->NoPNodes = PTInfo->NoPList;

//	PrintPPartInfo(Trees, Tree);
	FreePTInfo(PTInfo);
//	exit(0);
}

void	SetPNodes(TREE *Tree)
{
	int Index;

	SetVisitedTree(Tree, FALSE);

	if(Tree->NoPNodes > 1)
	{
		for(Index=0;Index<Tree->NoPNodes;Index++)
			Tree->PNodes[Index]->Visited = TRUE;
	}
}