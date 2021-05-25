#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "typedef.h"
#include "genlib.h"
#include "phyloplasty.h"
#include "matrix.h"
#include "rand.h"
#include "likelihood.h"
#include "trees.h"

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

int		IsValidPPNode(NODE N)
{
	if(N->Length == 0)
		return FALSE;

	return TRUE;
}

int		FindNoValidNodes(TREES *Trees)
{
	TREE	*T;
	int		Index;
	int		Ret;

	Ret = 0;
	T = &Trees->Tree[0];
	for(Index=0;Index<Trees->NoOfNodes;Index++)
	{
		if(IsValidPPNode(&T->NodeList[Index]) == TRUE)
			Ret++;
	}

	return Ret;
}

void	MakeValidNodes(TREES *Trees, PLASTY* Plasty)
{
	TREE *T;
	int	Index;
	NODE N;
	Plasty->NoValidNode = FindNoValidNodes(Trees);

	Plasty->ValidNode = (NODE*)malloc(sizeof(NODE) * Plasty->NoValidNode);
	if(Plasty->ValidNode == NULL)
		MallocErr();

	Plasty->NoValidNode = 0;

	T = &Trees->Tree[0];
	for(Index=0;Index<Trees->NoOfNodes;Index++)
	{
		N = &T->NodeList[Index];
		if(IsValidPPNode(N) == TRUE)
			Plasty->ValidNode[Plasty->NoValidNode++] = N;
	}

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
	Ret->NoID = 0;

	MakeValidNodes(Trees, Ret);

	Ret->ScaleBL = (double*)malloc(sizeof(double) * Ret->NoValidNode);
	if(Ret->ScaleBL == NULL)
		MallocErr();

	Ret->TempList = (NODE*)malloc(sizeof(NODE) * Ret->NoValidNode);

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

	free(Plasty->ScaleBL);
	free(Plasty->ValidNode);
	free(Plasty);
}


void	PlastyAdd(RATES *Rates, TREES *Trees, OPTIONS *Opt, NODE N)
{
	PLASTYNODE	*PNode;
	PLASTY		*Plasty;

	Plasty = Rates->Plasty;

	PNode = (PLASTYNODE*)malloc(sizeof(PLASTYNODE));
	if(PNode == NULL)
		MallocErr();

	PNode->NodeID = Plasty->NoID++;

	PNode->Node = N;

	if(GenRandState(Rates->RandStates) < 0.1)
		PNode->Type = PPBRANCH;
	else
		PNode->Type = PPNODE;

	if(GenRandState(Rates->RandStates) < 0.5)
		PNode->Scale = GenRandState(Rates->RandStates);
	else
		PNode->Scale = 1 + (GenRandState(Rates->RandStates) * 9);

	Plasty->NodeList = (PLASTYNODE**) AddToList(&Plasty->NoNodes, Plasty->NodeList, (void*)PNode);

	Rates->LogHRatio += -PPCOST + log(PPJACOBIAN);
}

void	DelNNodes(PLASTY *Plasty, int NoDel)
{
	int No;
	int Pos;

	No = 0;
	do
	{
		do
		{
			Pos = rand() % Plasty->NoNodes;
		}while(Plasty->NodeList[Pos] == NULL);

		free(Plasty->NodeList[Pos]);
		Plasty->NodeList[Pos] = NULL;
		No++;
	} while(No != NoDel);
}

void	MassRemove(RATES *Rates, TREES *Trees, OPTIONS *Opt)
{
	PLASTY		*Plasty;
	PLASTYNODE	**NList;
	int			Index;
	int			No;

	Plasty = Rates->Plasty;

	if(Plasty->NoNodes <= PPMASSDEL)
		return;

	DelNNodes(Plasty, PPMASSDEL);

	NList = (PLASTYNODE**)malloc(sizeof(PLASTYNODE*) * (Plasty->NoNodes - PPMASSDEL));
	if(NList == NULL) MallocErr();

	No = 0;
	for(Index=0;Index<Plasty->NoNodes;Index++)
	{
		if(Plasty->NodeList[Index] != NULL)
			NList[No++] = Plasty->NodeList[Index];
	}

	free(Plasty->NodeList);
	Plasty->NodeList = NList;
	Plasty->NoNodes = Plasty->NoNodes - PPMASSDEL;

	Rates->LogHRatio= PPCOST * PPMASSDEL;
}

void	DelPlastyNode(RATES *Rates, TREES *Trees, OPTIONS *Opt, int No)
{
	PLASTY		*Plasty;
	PLASTYNODE	**NList;
	int			Index;
	
	Plasty = Rates->Plasty;

	if(Plasty->NoNodes == 1)
	{
		free(Plasty->NodeList[0]);
		free(Plasty->NodeList);
		Plasty->NodeList = NULL;
		Plasty->NoNodes = 0;
		Rates->LogHRatio= PPCOST;
		return;
	}

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


	Rates->LogHRatio = PPCOST + log(1.0 / PPJACOBIAN);
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

void	MakeTNodeList(NODE N, NODE* List, int *Size)
{
	if(IsValidPPNode(N) == TRUE)
	{
		List[*Size] = N;
		(*Size)++;
	}

	if(N->Tip == TRUE)
		return;

	MakeTNodeList(N->Left, List, Size);
	MakeTNodeList(N->Right, List, Size);
}

void	MoveNode(RATES *Rates, TREES *Trees, OPTIONS *Opt)
{
	PLASTY		*Plasty;
	PLASTYNODE	*PNode;
	NODE		N;
	int			No;

	Plasty = Rates->Plasty;

	No = rand() % Plasty->NoNodes;
	PNode = Plasty->NodeList[No];

	N = PNode->Node;

	Plasty->NoTempList = 0;
	if((N->Ans != NULL) && (IsValidPPNode(N->Ans) == TRUE))
		Plasty->TempList[Plasty->NoTempList++] = N->Ans;

	if(N->Tip == FALSE)
	{
		MakeTNodeList(N->Left, Plasty->TempList, &Plasty->NoTempList);
		MakeTNodeList(N->Right, Plasty->TempList, &Plasty->NoTempList);
	}

	if(Plasty->NoTempList > 0)
	{
		No = rand() % Plasty->NoTempList;
		PNode->Node = Plasty->TempList[No];
	}
	else
		PNode->Scale = ChangePlastyRate(PNode->Scale, Opt->RateDev, Rates->RandStates);
}

void	MutatePlasty(RATES *Rates, TREES *Trees, OPTIONS *Opt)
{
	PLASTY		*Plasty;
	PLASTYNODE	*Node;
	int			No;

	Plasty = Rates->Plasty;

	No = rand() % Plasty->NoNodes;
	Node = Plasty->NodeList[No];

	if((GenRandState(Rates->RandStates) < 0.05) && (Node->Node->Tip == FALSE))
	{
		if(Node->Type == PPNODE)
			Node->Type = PPBRANCH;
		else
			Node->Type = PPNODE;
		return;
	}

	if(GenRandState(Rates->RandStates) < 0.25)
		MoveNode(Rates, Trees, Opt);
	else
		Node->Scale = ChangePlastyRate(Node->Scale, Opt->RateDev, Rates->RandStates);
}

int GetPlastyNode(int ID, PLASTY *Plasty)
{
	int Index;

	for(Index=0;Index<Plasty->NoNodes;Index++)
		if(Plasty->NodeList[Index]->Node->ID == ID)
			return Index;

	return -1;
}

void	PlastyMove(RATES *Rates, TREES *Trees, OPTIONS *Opt)
{
	PLASTY *	Plasty;
	int			PNodeID;
	NODE		N;

	Plasty = Rates->Plasty;

	if((Plasty->NoNodes > 0) && (GenRandState(Rates->RandStates) < 0.3))
	{
		MutatePlasty(Rates, Trees, Opt);
		return;
	}
/*
	if((Plasty->NoNodes > PPMASSDEL) && (GenRandState(Rates->RandStates) < 0.01))
	{
		MassRemove(Rates, Trees, Opt);
		return;
	}
*/
	N = Plasty->ValidNode[rand() % Plasty->NoValidNode];
	PNodeID = GetPlastyNode(N->ID, Plasty);

	if(PNodeID == -1)
		PlastyAdd(Rates, Trees, Opt, N);		
	else
		DelPlastyNode(Rates, Trees, Opt, PNodeID);		
}

/*
void	PlastyMove(RATES *Rates, TREES *Trees, OPTIONS *Opt)
{
	PLASTY *	Plasty;
	int			PNodeID;
	NODE		N;

	Plasty = Rates->Plasty;

	if( (Plasty->NoNodes > 0) && 
		(GenRandState(Rates->RandStates) < 0.3))
	{
		MutatePlasty(Rates, Trees, Opt);
		return;
	}

	
	if((GenRandState(Rates->RandStates) < 0.2) && (Plasty->NoNodes > 0))
	{
		PNodeID = rand() % Plasty->NoNodes;		
		DelPlastyNode(Rates, Trees, Opt, PNodeID);
	}
	else
	{
		N = Plasty->ValidNode[rand() % Plasty->NoValidNode];
		PlastyAdd(Rates, Trees, Opt, N);	
	}
}
*/

PLASTYNODE *ClonePlastyNode(PLASTYNODE *N2)
{
	PLASTYNODE *Ret;

	Ret = (PLASTYNODE*)malloc(sizeof(PLASTYNODE));
	if(Ret == NULL)
		MallocErr();

	Ret->Node = N2->Node;
	Ret->Scale= N2->Scale;
	Ret->Type = N2->Type;
	Ret->NodeID = N2->NodeID;

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

	P1->NoID = P2->NoID;
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

	Buffer = (char*)malloc(sizeof(char) * (strlen(Opt->LogFN) + BUFFERSIZE));
	if(Buffer == NULL)
		MallocErr();
	sprintf(Buffer, "%s.PP.trees", Opt->DataFN);
	Opt->PPTree = OpenWrite(Buffer);
	free(Buffer);

	fprintf(Opt->PPTree, "#NEXUS\n");
	fprintf(Opt->PPTree, "\tBegin Trees;\n");
	fprintf(Opt->PPTree, "\t\tTranslate\n");

	for(Index=0;Index<Trees->NoOfTaxa-1;Index++)
		fprintf(Opt->PPTree, "\t\t%d\t%s,\n", Trees->Taxa[Index].No, Trees->Taxa[Index].Name);
	fprintf(Opt->PPTree, "\t\t%d\t%s\n\t\t;\n", Trees->Taxa[Index].No, Trees->Taxa[Index].Name);

	fflush(Opt->PPTree);
}


void	RecPrintPPNodes(FILE *Out, NODE N)
{
	if(N->Tip == TRUE)
	{
		fprintf(Out, "%d\t", N->Taxa->No);
		return;
	}

	RecPrintPPNodes(Out, N->Left);
	RecPrintPPNodes(Out, N->Right);
}


void	PPLogFileHeader(OPTIONS *Opt, TREES *Trees, RATES *Rates)
{
	TREE	*T;
	int		Index;
	NODE	N;
	PLASTY	*P;
	int		No;

	T = &Trees->Tree[0];
	P = Rates->Plasty;

	fprintf(Opt->PPLog, "%d\n", Trees->NoOfTaxa);
	for(Index=0;Index<Trees->NoOfTaxa;Index++)
		fprintf(Opt->PPLog, "%d\t%s\n", Trees->Taxa[Index].No, Trees->Taxa[Index].Name);

	fprintf(Opt->PPLog, "%d\n", P->NoValidNode);
	for(Index=0;Index<P->NoValidNode;Index++)
	{
		N = P->ValidNode[Index];
		No = 0;
		CTaxaBelow(N, &No);
		fprintf(Opt->PPLog, "%d\t%f\t%d\t", N->ID, N->Length, No);
		RecPrintPPNodes(Opt->PPLog, N);
		fprintf(Opt->PPLog, "\n");
	}
}

void	IntiPPLogFile(OPTIONS *Opt, TREES *Trees, RATES *Rates)
{
	char	*Buffer;

	Buffer = (char*)malloc(sizeof(char) * (strlen(Opt->LogFN) + BUFFERSIZE));
	if(Buffer == NULL)
		MallocErr();
	sprintf(Buffer, "%s.PP.txt", Opt->DataFN);
	Opt->PPLog = OpenWrite(Buffer);
	free(Buffer);

	PPLogFileHeader(Opt, Trees, Rates);

	fflush(Opt->PPLog);
}

void	GetNodeIDList(NODE N, int *Size, int *List)
{
	if(N->Tip == TRUE)
	{
		List[*Size] = N->Taxa->No;
		(*Size)++;
		return;
	}
	
	GetNodeIDList(N->Left, Size, List);
	GetNodeIDList(N->Right, Size, List);
}


void	InitPPFiles(OPTIONS *Opt, TREES *Trees, RATES* Rates)
{
	InitPPTreeFile(Opt, Trees);
	IntiPPLogFile(Opt, Trees, Rates);
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

void	PrintPPTree(OPTIONS *Opt, TREES *Trees, RATES *Rates, int It)
{
	PLASTY *P;
	TREE	*T;

	P = Rates->Plasty;
	T = &Trees->Tree[Rates->TreeNo];

	Plasty(Opt, Trees, Rates);

	fprintf(Opt->PPTree, "\tTree T_%010d_%d = (", It, P->NoNodes);
	PrintPPNode(Opt->PPTree, T->Root->Left);
	fprintf(Opt->PPTree, ",");
	PrintPPNode(Opt->PPTree, T->Root->Right);
	fprintf(Opt->PPTree, ");\n");
	fflush(Opt->PPTree);
}

void	LogPPResults(OPTIONS *Opt, TREES *Trees, RATES *Rates, int It)
{
	FILE		*Out;
	PLASTY		*P;
	int			Index;
	PLASTYNODE	*PNode;
	NODE		N;

	P = Rates->Plasty;

	Out = Opt->PPLog;

	fprintf(Out, "%d\t%f\t%d\t", It, Rates->Lh, P->NoNodes);
	fprintf(Out, "%f\t%f\t", Rates->Contrast->EstAlpha[0], Rates->Contrast->EstSigma[0]);

	for(Index=0;Index<P->NoNodes;Index++)
	{
		PNode = P->NodeList[Index];
		N = PNode->Node;

		fprintf(Out, "%d\t", N->ID);
		fprintf(Out, "%f\t", PNode->Scale);
		fprintf(Out, "%d\t", PNode->NodeID);

		if(PNode->Type == PPNODE)
			fprintf(Out, "Node\t");
		else
			fprintf(Out, "Branch\t");
	}
	fprintf(Out, "\n");

	fflush(Out);
}	

void	PrintPPOutput(OPTIONS *Opt, TREES *Trees, RATES *Rates, int It)
{
	PrintPPTree(Opt, Trees, Rates, It);
	LogPPResults(Opt, Trees, Rates, It);
}
