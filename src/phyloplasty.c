#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "typedef.h"
#include "genlib.h"
#include "phyloplasty.h"
#include "matrix.h"
#include "rand.h"
#include "likelihood.h"

double**	MakeTrueBL(TREES *Trees)
{
	double **Ret;
	int		Index, NIndex;
	TREE	*T;
	Ret = (double**)malloc(sizeof(double*) * Trees->NoOfTrees);
	if(Ret == NULL)
		MallocErr();

	for(Index=0;Index<Trees->NoOfTrees;Index++)
	{
		Ret[Index] = (double*)malloc(sizeof(double) * Trees->NoOfNodes);
		if(Ret[Index] == NULL)
			MallocErr();
	}

	for(Index=0;Index<Trees->NoOfTrees;Index++)
	{
		T = &Trees->Tree[Index];
		for(NIndex=0;NIndex<Trees->NoOfNodes;NIndex++)
			Ret[Index][NIndex] = T->NodeList[NIndex].Length;
	}

	return Ret;
}

PLASTY*	CreatPlasty(RATES *Rates, TREES *Trees, OPTIONS *Opt)
{
	PLASTY* Ret;

	Ret = (PLASTY*)malloc(sizeof(PLASTY));
	if(Ret == NULL)
		MallocErr();

	Ret->NoTrees = Trees->NoOfTrees;
	Ret->NoNodes = 0;
	Ret->NodeList= NULL;

	Ret->TrueBL = MakeTrueBL(Trees);

	return Ret;
}

void	FreePlasty(PLASTY* Plasty)
{
	int Index;

	for(Index=0;Index<Plasty->NoTrees;Index++)
		free(Plasty->TrueBL[Index]);
	free(Plasty->TrueBL);

	if(Plasty->NodeList != NULL)
	{
		for(Index=0;Index<Plasty->NoNodes;Index++)
			free(Plasty->NodeList[Index]);
		
		free(Plasty->NodeList);
	}

	free(Plasty);
}

void	PlastySwap(RATES *Rates, TREES *Trees, OPTIONS *Opt)
{
	PLASTY		*Plasty;
	int			N1, N2;
	PLASTYNODE	*T;

	Plasty = Rates->Plasty;

	N1 = rand() % Plasty->NoNodes;
	do
	{
		N2 = rand() % Plasty->NoNodes;
	}while(N1 == N2);

	T = Plasty->NodeList[N1];
	Plasty->NodeList[N1] = Plasty->NodeList[N2];
	Plasty->NodeList[N2] = T;
}

void	PlastyAdd(RATES *Rates, TREES *Trees, OPTIONS *Opt)
{
	PLASTYNODE	*PNode;
	TREE		*T;
	NODE		N;
	PLASTY		*Plasty;

	Plasty = Rates->Plasty;

	PNode = (PLASTYNODE*)malloc(sizeof(PLASTYNODE));
	if(PNode == NULL)
		MallocErr();

	T = &Trees->Tree[0];

	N = &T->NodeList[rand() % Trees->NoOfNodes];

	PNode->Node = N;

	if(GenRandState(Rates->RandStates) < 0.20)
		PNode->Type = PPBRANCH;
	else
		PNode->Type = PPNODE;

	if(GenRandState(Rates->RandStates) < 0.5)
		PNode->Scale = GenRandState(Rates->RandStates);
	else
		PNode->Scale = 1 + (GenRandState(Rates->RandStates) * 9);

	Plasty->NodeList = (PLASTYNODE**) AddToList(&Plasty->NoNodes, Plasty->NodeList, (void*)PNode);

//	Rates->LhPrior = 300;
	Rates->LogHRatio= -PPCOST;
}

void	PlastyDel(RATES *Rates, TREES *Trees, OPTIONS *Opt)
{
	PLASTY		*Plasty;
	PLASTYNODE	**NList;
	int			Index, No;
	
	Plasty = Rates->Plasty;

	if(Plasty->NoNodes == 1)
	{
		free(Plasty->NodeList[0]);
		free(Plasty->NodeList);
		Plasty->NodeList = NULL;
		Plasty->NoNodes = 0;
		return;
	}

	No = rand() % Plasty->NoNodes;

	free(Plasty->NodeList[No]);
	Plasty->NodeList[No] = NULL;


	NList = (PLASTYNODE**)malloc(sizeof(PLASTYNODE*) * (Plasty->NoNodes - 1));
	if(NList == NULL)
		MallocErr();

	No = 0;
	for(Index=0;Index<Plasty->NoNodes;Index++)
	{
		if(Plasty->NodeList[Index] != NULL)
			NList[No++] = Plasty->NodeList[Index];
	}

	free(Plasty->NodeList);
	Plasty->NodeList = NList;
	Plasty->NoNodes--;

	Rates->LogHRatio= PPCOST;
}

double	ChangePlastyRate(double Rate, double dev, RANDSTATES *RS)
{
	double Ret;
	
	do
	{
		Ret = (GenRandState(RS) * dev) - (dev / 2.0); 
		Ret += Rate;
	} while(Ret <= 0);

	return Ret;
}

void	MutatePlasty(RATES *Rates, TREES *Trees, OPTIONS *Opt)
{
	PLASTY		*Plasty;
	PLASTYNODE	*Node;
	int			No;

	Plasty = Rates->Plasty;

	No = rand() % Plasty->NoNodes;
	Node = Plasty->NodeList[No];

	Node->Scale = ChangePlastyRate(Node->Scale, Opt->RateDev, Rates->RandStates);
}

void	PlastyMove(RATES *Rates, TREES *Trees, OPTIONS *Opt)
{
	PLASTY *Plasty;

	Plasty = Rates->Plasty;

	if(Plasty->NoNodes == 0)
	{
		PlastyAdd(Rates, Trees, Opt);
		return;
	}

	if((Plasty->NoNodes > 1) && (GenRandState(Rates->RandStates) < 0.05))
	{	
		PlastySwap(Rates, Trees, Opt);
		return;
	}

	if(GenRandState(Rates->RandStates) < 0.2)
	{
		PlastyDel(Rates, Trees, Opt);
		return;
	}

	if(GenRandState(Rates->RandStates) < 0.5)
		PlastyAdd(Rates, Trees, Opt);
	else
		MutatePlasty(Rates, Trees, Opt);

}

PLASTYNODE *ClonePlastyNode(PLASTYNODE *N2)
{
	PLASTYNODE *Ret;

	Ret = (PLASTYNODE*)malloc(sizeof(PLASTYNODE));
	if(Ret == NULL)
		MallocErr();

	Ret->Node = N2->Node;
	Ret->Scale= N2->Scale;
	Ret->Type = N2->Type;

	return Ret;
}

void	BlankPlasty(PLASTY *P)
{
	int Index;

	if(P->NoNodes != 0)
	{
		for(Index=0;Index<P->NoNodes;Index++)
			free(P->NodeList[Index]);
		free(P->NodeList);
	}

	P->NodeList = NULL;
	P->NoNodes = 0;
}

void	PlastyCopy(RATES *R1, RATES *R2)
{
	PLASTYNODE	**NList;
	PLASTY		*P1, *P2;
	int			Index;

	P1 = R1->Plasty;
	P2 = R2->Plasty;

	if(P2->NoNodes == 0)
	{
		BlankPlasty(P1);
		return;
	}

	NList = (PLASTYNODE**)malloc(sizeof(PLASTYNODE*) * P2->NoNodes);
	if(NList == NULL)
		MallocErr();

	for(Index=0;Index<P2->NoNodes;Index++)
		NList[Index] = ClonePlastyNode(P2->NodeList[Index]);

	BlankPlasty(P1);

	P1->NodeList = NList;
	P1->NoNodes = P2->NoNodes;
}

void	PlastyNode(NODE N, PLASTYNODE *P)
{
	N->Length = N->Length * P->Scale;
	if(P->Type == PPBRANCH)
		return;

	if(N->Tip == TRUE)
		return;

	PlastyNode(N->Left, P);
	PlastyNode(N->Right, P);
}

void	Plasty(OPTIONS *Opt, TREES *Trees, RATES *Rates)
{
	int Index;
	int	TNo;
	TREE *T;
	PLASTY *P;

	P = Rates->Plasty;
	TNo = Rates->TreeNo;
	T = &Trees->Tree[TNo];

	for(Index=0;Index<Trees->NoOfNodes;Index++)
		T->NodeList[Index].Length = P->TrueBL[TNo][Index];

	for(Index=0;Index<P->NoNodes;Index++)
		PlastyNode(P->NodeList[Index]->Node, P->NodeList[Index]);
}

void	InitPPTreeFile(OPTIONS *Opt, TREES *Trees)
{
	char	*Buffer;
	int		Index;

	Buffer = (char*)malloc(sizeof(char) * (strlen(Opt->LogFN) + 7));
	if(Buffer == NULL)
		MallocErr();
	sprintf(Buffer, "%s.trees", Opt->LogFN);
	Opt->PPTree = OpenWrite(Buffer);

	fprintf(Opt->PPTree, "#NEXUS\n");
	fprintf(Opt->PPTree, "\tBegin Trees;\n");
	fprintf(Opt->PPTree, "\t\tTranslate\n");

	for(Index=0;Index<Trees->NoOfTaxa-1;Index++)
		fprintf(Opt->PPTree, "\t\t%d\t%s,\n", Trees->Taxa[Index].No, Trees->Taxa[Index].Name);
	fprintf(Opt->PPTree, "\t\t%d\t%s\n\t\t;\n", Trees->Taxa[Index].No, Trees->Taxa[Index].Name);

	fflush(Opt->PPTree);
}

void	PrintPPNode(FILE *Out, NODE N)
{
	if(N->Tip == TRUE)
	{
		fprintf(Out, "%d:%f", N->Taxa->No, N->Length);
		return;
	}

	fprintf(Out, "(");

	PrintPPNode(Out, N->Left);
	fprintf(Out, ",");
	PrintPPNode(Out, N->Right);
	fprintf(Out, "):%f", N->Length);
}

void	PrintPPTree(OPTIONS *Opt, TREES *Trees, RATES *Rates)
{
	PLASTY *P;
	TREE	*T;

	P = Rates->Plasty;
	T = &Trees->Tree[Rates->TreeNo];

	Plasty(Opt, Trees, Rates);

	fprintf(Opt->PPTree, "\tTree T_%d = (", P->NoNodes);
	PrintPPNode(Opt->PPTree, T->Root->Left);
	fprintf(Opt->PPTree, ",");
	PrintPPNode(Opt->PPTree, T->Root->Right);
	fprintf(Opt->PPTree, ");\n");
	fflush(Opt->PPTree);
}