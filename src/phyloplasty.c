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

int		NoFixedTo1(PHYLOPLASTY *PP)
{
	int	Ret;
	int Index;

	Ret = 0;
	for(Index=1;Index<PP->NoRates;Index++)
		if(PP->Rates[Index] == 1)
			Ret++;

	return Ret;
}

int		PhyloPlasyMoveChoice(RATES *Rates)
{
	PHYLOPLASTY *PP;
	int			No1;

	PP = Rates->PhyloPlasty;
	No1 = NoFixedTo1(PP);
	
	do
	{
/*		if((GenRandState(Rates->RandStates) < 0.75) &&
			(PP->NoDiffRates > 0))
				return 0;
		
*/		if((GenRandState(Rates->RandStates) < 0.5) &&
			(PP->NoDiffRates > 0))
				return 1;
		
		if(No1 > 0)
			return 2;			

	} while(1==1);

	return -1;
}

int		GetPPRateMutPos(RATES *Rates)
{
	int Index;
	int No;
	PHYLOPLASTY *PP;

	PP = Rates->PhyloPlasty;

	No = rand() % PP->NoDiffRates;

	for(Index=0;Index<PP->NoRates;Index++)
	{
		if(PP->Rates[Index] != 1)
		{
			if(No == 0)
				return Index;
			No--;
		}
	}

	return -1;
}

void	MutatePPRate(RATES *Rates)
{
	int			Pos;
	double		NewR;
	PHYLOPLASTY *PP;
	int			Valid;
		
	PP = Rates->PhyloPlasty;

	Pos = GetPPRateMutPos(Rates);

	do
	{
		Valid = TRUE;
	//	NewR = PP->Rates[Pos] + (nrand(Rates->RandStates)*1.5);
		NewR = GenRandState(Rates->RandStates) * 10;

		if(GenRandState(Rates->RandStates) < 0.5)
			NewR = 1 + (GenRandState(Rates->RandStates) * 9);
		else
			NewR = GenRandState(Rates->RandStates);


		if(NewR < 0)
			Valid = FALSE;
		
		if(NewR > 10)
			Valid = FALSE;

	}while(Valid == FALSE);

	printf("Mutate\t%d\tOld\t%f\tNew\t%f\n", Pos, PP->Rates[Pos], NewR);
	fflush(stdout);
	PP->Rates[Pos] = NewR;
}

void	CollapseRate(RATES *Rates)
{
	PHYLOPLASTY *PP;
	int	Pos;

	PP = Rates->PhyloPlasty;


/*	Pos = rand() % (PP->NoRates - 1);
	Pos++;

	if(PP->Rates[Pos] == 1)
		return;
*/	
	Pos = GetPPRateMutPos(Rates);
	PP->Rates[Pos] = 1;

	printf("Colapes\t%d\n", Pos);
	fflush(stdout);


	PP->NoDiffRates--;
}

void	ExpandRate(RATES *Rates)
{
	PHYLOPLASTY *PP;
	int	Pos;

	PP = Rates->PhyloPlasty;

	do
	{
		Pos = rand() % (PP->NoRates - 1);
		Pos++;
	}while(PP->Rates[Pos] != 1);

	if(GenRandState(Rates->RandStates) < 0.5)
		PP->Rates[Pos] = 1 + (GenRandState(Rates->RandStates) * 9);
	else
		PP->Rates[Pos] = GenRandState(Rates->RandStates);

	
	printf("Expand\t%d\t%f\n", Pos, PP->Rates[Pos]);
	fflush(stdout);


	PP->NoDiffRates++;
}

void	PhyloPlasyMove(RATES *Rates)
{
	int Index;

	MutatePPRate(Rates);
	/*
	switch(PhyloPlasyMoveChoice(Rates))
	{
		case 0: MutatePPRate(Rates); break;
		case 1: CollapseRate(Rates); break;
		case 2: ExpandRate(Rates); break;
	}*/
}

/*
void	PhyloPlasyMove(RATES *Rates)
{
	PHYLOPLASTY *PP;
	int Pos;
	int	NC;

	PP = Rates->PhyloPlasty;
	Rates->Rates[0] = 3.720719;

	do
	{
		Pos = rand() % PP->NoCats;
	}while(Pos == 0);

	if(	(PP->Cats[Pos] != PP->OnePos) && 
		(GenRandState(Rates->RandStates) < 0.5))
	{
		PP->Cats[Pos] = PP->OnePos;
		return;
	}

	do
	{
		NC = rand() % PP->NoRates;
	} while(NC == PP->Cats[Pos]);


	PP->Cats[Pos] = NC;
}
*/

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

	Ret->NoRates = Trees->NoOfNodes;

	Ret->NoDiffRates = 0;
	Ret->Rates = (double*)malloc(sizeof(double) * Ret->NoRates);
	Ret->RealBL= (double*)malloc(sizeof(double) * Ret->NoRates);
	if((Ret->Rates == NULL) || (Ret->RealBL == NULL))
		MallocErr();
	
	for(Index=0;Index<Ret->NoRates;Index++)
	{
		Ret->Rates[Index] = 1;
		Ret->RealBL[Index] = Tree->NodeList[Index].Length;
	}
	Ret->InvV = TRUE;
/*
	for(Index=0;Index<Ret->NoRates;Index++)
	{
		printf("%d\t", Index);
		PrintNode(&Tree->NodeList[Index]);
		printf("\n");
	}
	exit(0); */
	return Ret;
}

void	MapRateToTree(TREES* Trees, RATES* Rates)
{
	TREE		*Tree;
	PHYLOPLASTY	*PP;
	int			Index;

	PP = Rates->PhyloPlasty;
	Tree = &Trees->Tree[Rates->TreeNo];

	for(Index=0;Index<PP->NoRates;Index++)
		Tree->NodeList[Index].Length = PP->RealBL[Index] * PP->Rates[Index];
}

void	SetPhyloPlastyV(TREES* Trees, RATES* Rates)
{
	TREE	*Tree;
	int		Index;

	Tree = &Trees->Tree[Rates->TreeNo];

	MapRateToTree(Trees, Rates);
//	CalcPVarCoVar(Trees, Tree); 
	MapPhyloPlastyToV(Trees, Tree);

	return;

	printf("Stat\t");
	PrintTime(stdout);
	printf("\n");
	fflush(stdout);
	for(Index=0;Index<10000;Index++)
	{
		MapRateToTree(Trees, Rates);
		MapPhyloPlastyToV(Trees, Tree);

		CopyMatrix(Tree->ConVars->TrueV, Tree->ConVars->V);
		if(Index%1000==0)
		{
			printf("\tdone\t%d\n", Index);
			fflush(stdout);
		}

//		FindInvV(Trees, Tree);
	}
	printf("End\t");
	PrintTime(stdout);
	exit(0);
}