#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "typedef.h"
#include "genlib.h"
#include "continuous.h"
#include "ckappa.h"

void	InitCKappaTree(TREES* Trees, TREE* Tree);
void	KappaVarCoVar(TREES* Trees, TREE* Tree);

void	InitCKappaTree(TREES* Trees, TREE* Tree)
{
	CONVAR		*ConVar;
	TAXADIST	*TaxaDist;
	int			Index;
	int			TIndex;

	ConVar = Tree->ConVars;

	ConVar->TaxaDist = (TAXADIST*)malloc(sizeof(TAXADIST) * Trees->NoOfTaxa);
	if(ConVar->TaxaDist == NULL)
		MallocErr();

	TaxaDist = ConVar->TaxaDist;

	for(TIndex=0;TIndex<Trees->NoOfTaxa;TIndex++)
	{
		TaxaDist[TIndex].NoOfBLVect = (int*)malloc(sizeof(int) * Trees->NoOfTaxa);
		if(TaxaDist == NULL)
			MallocErr();
	} 

	for(TIndex=0;TIndex<Trees->NoOfTaxa;TIndex++)
	{
		TaxaDist[TIndex].TToTPath = (double**)malloc(sizeof(double*) * Trees->NoOfTaxa);
		if(TaxaDist->TToTPath == NULL)
			MallocErr();
	}


	for(TIndex=0;TIndex<Trees->NoOfTaxa;TIndex++)
	{
		for(Index=0;Index<Trees->NoOfTaxa;Index++)
		{
			TaxaDist[TIndex].NoOfBLVect[Index] = 0;
			TaxaDist[TIndex].TToTPath[Index] = NULL;
		}
	}
}

/*
typedef struct
{
	int		*NoOfBLVect;
	double	**TToTPath;
}  TAXADIST;
*/

void	FreeCKappaTree(CONVAR *ConVar, int NoOfTaxa)
{
	TAXADIST	*TaxaDist;
	int			TIndex;
	int			Index;

	TaxaDist = ConVar->TaxaDist;

	for(TIndex=0;TIndex<NoOfTaxa;TIndex++)
	{
		for(Index=0;Index<NoOfTaxa;Index++)
		{
			if(TaxaDist[TIndex].TToTPath[Index] != NULL)
				free(TaxaDist[TIndex].TToTPath[Index]);
		}
		free(TaxaDist[TIndex].TToTPath);
		free(TaxaDist[TIndex].NoOfBLVect);	
	}

	free(ConVar->TaxaDist);	
}

void	NodesBelowNode(TREE *Tree, NODE N, int *Count)
{
	if(Tree->Root == N)
		return;

	(*Count)++;

	NodesBelowNode(Tree, N->Ans, Count);
}

void	BlankVisited(TREES* Trees, TREE* Tree)
{
	int	NIndex;
	for(NIndex=0;NIndex<Tree->NoNodes;NIndex++)
		Tree->NodeList[NIndex]->Visited = FALSE;
}

void	SetUpTrace(NODE N, TREE *Tree)
{
	N->Visited = TRUE;
	if(N == Tree->Root)
		return;

	SetUpTrace(N->Ans, Tree);
}

void	KappaTaxaVarCoVar(TREES* Trees, TREE* Tree, TAXA* Taxa, int TNo)
{
	int			NodesBlow;
	int			TIndex;
	NODE		MyNode;
	NODE		TNode;
	NODE		CNode;
	TAXA		*CTaxa;
	TAXADIST	*TaxaDist;
	int			BLIndex;

	MyNode = TaxaToNode(Trees, Tree, Taxa);
	TaxaDist = &Tree->ConVars->TaxaDist[TNo];
	
	for(TIndex=0;TIndex<Trees->NoOfTaxa;TIndex++)
	{
		CTaxa = &Trees->Taxa[TIndex];
		
		BlankVisited(Trees, Tree);
		SetUpTrace(MyNode, Tree);

		TNode = TaxaToNode(Trees, Tree, CTaxa);
		while(TNode->Visited == FALSE)
			TNode = TNode->Ans;
	
		CNode = TNode;
		NodesBlow = 0;
		NodesBelowNode(Tree, TNode, &NodesBlow);

		TaxaDist->NoOfBLVect[TIndex] = NodesBlow;

		if(NodesBlow != 0)
		{
			TaxaDist->TToTPath[TIndex] = (double*)malloc(sizeof(double) * NodesBlow);
			if(TaxaDist->TToTPath[TIndex] == NULL)
				MallocErr();

			BLIndex = 0;
			while(CNode != Tree->Root)
			{
				TaxaDist->TToTPath[TIndex][BLIndex] = CNode->Length;
				BLIndex++;
				CNode = CNode->Ans;
			}
		}
		else
			TaxaDist->TToTPath[TIndex] = NULL;

	}
}



void	KappaVarCoVar(TREES* Trees, TREE* Tree)
{
	int	TIndex;

	for(TIndex=0;TIndex<Trees->NoOfTaxa;TIndex++)
	{
		KappaTaxaVarCoVar(Trees, Tree, &Trees->Taxa[TIndex], TIndex);
	}
}

void	InitCKappa(OPTIONS* Opt, TREES* Trees)
{
	int		TIndex;

	for(TIndex=0;TIndex<Trees->NoOfTrees;TIndex++)
	{
		InitCKappaTree(Trees, &Trees->Tree[TIndex]);
		KappaVarCoVar(Trees, &Trees->Tree[TIndex]);
	}
}

double	 CaclKappaVarCoVar(TREE *Tree, int T1, int T2, double Kappa)
{
	int			Index;
	int			Len;
	TAXADIST	*TaxaDist;
	double		Ret=0;

	TaxaDist = &Tree->ConVars->TaxaDist[T1];
	Len = TaxaDist->NoOfBLVect[T2];

	for(Index=0;Index<Len;Index++)
		Ret += pow(TaxaDist->TToTPath[T2][Index], Kappa);

	return Ret;
}

void	MakeKappaV(TREES* Trees, TREE* Tree, double Kappa)
{
	int		x,y;
	CONVAR	*ConVar;
	double	Temp;

	ConVar = Tree->ConVars;

	for(x=0;x<Trees->NoOfTaxa;x++)
	{
		for(y=x;y<Trees->NoOfTaxa;y++)
		{
			Temp = CaclKappaVarCoVar(Tree, x, y, Kappa);

			ConVar->V->me[x][y] = Temp;
			ConVar->V->me[y][x] = Temp;

		}
	}
}
