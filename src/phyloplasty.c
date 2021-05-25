#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "typedef.h"
#include "genlib.h"
#include "phyloplasty.h"
#include "matrix.h"
#include "rand.h"
#include "likelihood.h"

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

	CopyMatrix(ConVar->TrueV, ConVar->V);
}


void		PrintNode(NODE N)
{
	if(N->Tip == TRUE)
	{
		printf("%s\t", N->Taxa->Name);
		return;
	}

	PrintNode(N->Left);
	PrintNode(N->Right);
}

PHYLOPLASTY*	CreatPhyloPlasty(OPTIONS *Opt, RATES *Rates)
{
	TREES		*Trees;
	TREE		*Tree;
	PHYLOPLASTY	*Ret;
	int			Index;
	
	Trees = Opt->Trees;
	Tree  = &Trees->Tree[0];

	Ret = (PHYLOPLASTY*)malloc(sizeof(PHYLOPLASTY));
	if(Ret == NULL)
		MallocErr();

	Ret->NoBL	 = Trees->NoOfNodes;
	
	Ret->RealBL= (double*)malloc(sizeof(double) * Ret->NoBL);
	if(Ret->RealBL == NULL)
		MallocErr();
	
	for(Index=0;Index<Ret->NoBL;Index++)
		Ret->RealBL[Index] = Tree->NodeList[Index].Length;
	
	Ret->NodeList	= NULL;
	Ret->NoNodes	= 0;

	return Ret;
}

void	MapPPNodeToTree(NODE N, double Scale)
{
	N->Length =  N->Length * Scale;
	if(N->Tip == TRUE)
		return;

	MapPPNodeToTree(N->Left, Scale);
	MapPPNodeToTree(N->Right, Scale);
}

void	MapPPToTree(TREES* Trees, RATES* Rates)
{
	TREE		*Tree;
	PHYLOPLASTY	*PP;
	int			Index;

	PP = Rates->PhyloPlasty;
	Tree = &Trees->Tree[Rates->TreeNo];

	for(Index=0;Index<PP->NoBL;Index++)
		Tree->NodeList[Index].Length = PP->RealBL[Index];

	for(Index=0;Index<PP->NoNodes;Index++)
		MapPPNodeToTree(PP->NodeList[Index]->Node, PP->NodeList[Index]->Scale);
}

void	SetPhyloPlastyV(TREES* Trees, RATES* Rates)
{
	TREE	*Tree;
	int		Index;

	Tree = &Trees->Tree[Rates->TreeNo];

	MapPPToTree(Trees, Rates);
//	CalcPVarCoVar(Trees, Tree); 
	MapPhyloPlastyToV(Trees, Tree);
}

void			FindScaleNodeNTaxa(NODE N, PPSCALENODE* SNode)
{
	if(N->Tip == TRUE)
	{
		SNode->NoTaxa++;
		return;
	}

	FindScaleNodeNTaxa(N->Left, SNode);
	FindScaleNodeNTaxa(N->Right, SNode);
}

PPSCALENODE*	InitPPScaleNode(NODE N, double Scale)
{
	PPSCALENODE*	Ret;

	Ret = (PPSCALENODE*)malloc(sizeof(PPSCALENODE));
	if(Ret == NULL)
		MallocErr();

	Ret->Node	= N;
	Ret->Scale	= Scale;

	Ret->NoTaxa = 0;
	FindScaleNodeNTaxa(N, Ret);
	
	return Ret;
}

void	AddPPNode(RATES *Rates, NODE N, double Scale)
{
	PPSCALENODE*	SNode;
	PPSCALENODE**	NList;
	PHYLOPLASTY*	PP;

	PP = Rates->PhyloPlasty;
	SNode = InitPPScaleNode(N, Scale);

	NList = (PPSCALENODE**)malloc(sizeof(PPSCALENODE*) * (PP->NoNodes + 1));
	if(NList == NULL)
		MallocErr();

	if(PP->NoNodes > 1)
	{
		memcpy(NList, PP->NodeList, sizeof(PPSCALENODE*) * PP->NoNodes);
		free(PP->NodeList);
	}

	NList[PP->NoNodes] = SNode;
	PP->NodeList = NList;	
	PP->NoNodes++;
}

void	BlankPP(PHYLOPLASTY *PP)
{
	int				Index;

	for(Index=0;Index<PP->NoNodes;Index++)
		free(PP->NodeList[Index]);

	free(PP->NodeList);
	PP->NodeList	= NULL;
	PP->NoNodes		= 0;
}

void	CopyPhyloPlasty(PHYLOPLASTY *A, PHYLOPLASTY *B)
{
	int Index;

	BlankPP(A);

	A->NodeList = (PPSCALENODE**)malloc(sizeof(PPSCALENODE*) * B->NoNodes);
	if(A->NodeList == NULL)
		MallocErr();

	for(Index=0;Index<B->NoNodes;Index++)
	{
		A->NodeList[Index] = (PPSCALENODE*)malloc(sizeof(PPSCALENODE));
		if(A->NodeList[Index] == NULL)
			MallocErr();

		A->NodeList[Index]->Node	= B->NodeList[Index]->Node;
		A->NodeList[Index]->Scale	= B->NodeList[Index]->Scale;
		A->NodeList[Index]->NoTaxa	= B->NodeList[Index]->NoTaxa;
	}

	A->NoNodes = B->NoNodes;
}

void	PhyloPlasyMove(RATES *Rates)
{
	
}