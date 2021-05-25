#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "typedef.h"
#include "genlib.h"
#include "phyloplasty.h"

int			FindNoPPLists(int NoTaxa)
{
	int Ret;

	Ret = NoTaxa * NoTaxa;
	Ret = (Ret - NoTaxa) / 2;
	Ret = Ret + NoTaxa;

	return Ret;
}

NODE	GetTNodeFromID(char* TName, TREES *Trees, TREE *Tree)
{
	int Index;
	NODE N;
	for(Index=0;Index<Trees->NoOfNodes;Index++)
	{
		N = &Tree->NodeList[Index];
		if(N->Tip == TRUE)
		{
			if(strcmp(TName, N->Taxa->Name) == 0)
				return N;
		}
	}

	return NULL;
}

int	NoNodesBelow(TREE *Tree, NODE N)
{
	int Ret;

	Ret = 0;

	while(N != Tree->Root)
	{
		N = N->Ans;
		Ret++;
	}

	return Ret;
}

PPCOVARPAIR*	MakePPVar(TREES *Trees, TREE *Tree, int TaxaID)
{
	PPCOVARPAIR* Ret;
	NODE Taxa;
	int		Index;

	Ret = (PPCOVARPAIR*)malloc(sizeof(PPCOVARPAIR));
	if(Ret == NULL)
		MallocErr();

	Taxa = GetTNodeFromID(Trees->Taxa[TaxaID].Name, Trees, Tree);

	Ret->x = TaxaID;
	Ret->y = TaxaID;

	Ret->Sum = 0;
	Ret->No = NoNodesBelow(Tree, Taxa);

	if(Ret->No == 0)
	{
		Ret->NList = NULL;
		return Ret;
	}

	Ret->NList = (NODE*)malloc(sizeof(NODE) * Ret->No);
	if(Ret->NList == NULL)
		MallocErr();

	Index=0;
	while(Taxa != Tree->Root)
	{
		Ret->NList[Index] = Taxa;
		Index++;
		Taxa = Taxa->Ans;
	}
	
	for(Index=0;Index<Ret->No;Index++)
		Ret->Sum += Ret->NList[Index]->Length;

	return Ret;
}

NODE	FindCommonNode(TREES *Trees, TREE *Tree, NODE NX, NODE NY)
{
	int Index;

	for(Index=0;Index<Trees->NoOfNodes;Index++)
		Tree->NodeList[Index].Visited = FALSE;
	
	while(NX!=Tree->Root)
	{
		NX = NX->Ans;
		NX->Visited = TRUE;
	}

	do
	{
		NY = NY->Ans;
	} while(NY->Visited != TRUE);

	return NY;
}

PPCOVARPAIR*	MakePPCoVar(TREES *Trees, TREE *Tree, int XID, int YID)
{
	PPCOVARPAIR* Ret;
	NODE		XTaxa;
	NODE		YTaxa;
	NODE		CNode;
	int			Index;

	Ret = (PPCOVARPAIR*)malloc(sizeof(PPCOVARPAIR));
	if(Ret == NULL)
		MallocErr();

	XTaxa = GetTNodeFromID(Trees->Taxa[XID].Name, Trees, Tree);
	YTaxa = GetTNodeFromID(Trees->Taxa[YID].Name, Trees, Tree);
	CNode = FindCommonNode(Trees, Tree, XTaxa, YTaxa);

	Ret->x = XID;
	Ret->y = YID;

	Ret->Sum = 0;

	Ret->No = NoNodesBelow(Tree, CNode);

	if(Ret->No == 0)
	{
		Ret->NList = NULL;
		return Ret;
	}

	Ret->NList = (NODE*)malloc(sizeof(NODE) * Ret->No);
	if(Ret->NList == NULL)
		MallocErr();
	
	Index=0;
	while(CNode != Tree->Root)
	{
		Ret->NList[Index] = CNode;
		Index++;
		CNode = CNode->Ans;
	}

	return Ret;
}


void		BuildPhyloPlastyConVar(TREES *Trees, TREE *Tree, PPCOVARV *PPCV)
{
	int x,y;
	int	Pos;

	Pos = 0;
	for(x=0;x<Trees->NoOfTaxa;x++)
	{
		for(y=x;y<Trees->NoOfTaxa;y++)
		{
			if(x == y)
				PPCV->List[Pos] = MakePPVar(Trees, Tree, x);
			else
				PPCV->List[Pos] = MakePPCoVar(Trees, Tree, x, y);
			Pos++;
		}			
	}
}

PPCOVARV*	InitPhyloPlastyConVar(TREES	*Trees, TREE *Tree)
{
	PPCOVARV	*Ret;
	
	Ret = (PPCOVARV*)malloc(sizeof(PPCOVARV));
	if(Ret == NULL)
		MallocErr();

	Ret->No =FindNoPPLists(Trees->NoOfTaxa);

	Ret->List = (PPCOVARPAIR**)malloc(sizeof(PPCOVARPAIR*) * Ret->No);
	if(Ret->List == NULL)
		MallocErr();

	BuildPhyloPlastyConVar(Trees, Tree, Ret);

	return Ret;
}

void	CalcPSize(PPCOVARPAIR* Pair)
{
	int Index;

	Pair->Sum = 0;
	for(Index=0;Index<Pair->No;Index++)
		Pair->Sum += Pair->NList[Index]->Length;
}

void	MapPhyloPlastyToV(TREES	*Trees, TREE *Tree)
{
	int			Index;
	CONVAR		*ConVar;
	PPCOVARV	*PPCoVarV;
	int			x,y;
	double		Dist;

	ConVar = Tree->ConVars;

	PPCoVarV = ConVar->PPCoVarV;

	for(Index=0;Index<PPCoVarV->No;Index++)
		CalcPSize(PPCoVarV->List[Index]);

	for(Index=0;Index<PPCoVarV->No;Index++)
	{
		x = PPCoVarV->List[Index]->x;
		y = PPCoVarV->List[Index]->y;
		Dist = PPCoVarV->List[Index]->Sum;

		ConVar->V->me[x][y] = Dist;
		ConVar->V->me[y][x] = Dist;
	}
}