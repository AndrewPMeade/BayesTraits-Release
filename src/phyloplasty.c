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
#include "randdists.h"

void	TestGenRand(OPTIONS *Opt, TREES *Trees, RATES* Rates);

extern double gamma(double x);
extern double ndtr(a);

double	CalcNormalHasting(double x, double SD)
{
	return log(ndtr(x/SD));
}

double	PPSGammaPDF(double x, double Alpha, double Beta)
{
	double Ret, s, T1, T2;

	s = 1 / ((Alpha - 1) * Beta);

	T1 = exp(-(x/s) / Beta);
	T2 = pow(x/s, -1 + Alpha);
	T2 = T1 * T2 * pow(Beta, -Alpha);
	Ret = T2 / gamma(Alpha);
	
	return Ret / s;
}

double	RandGamma(double Shape, double Scale)
{
	double x;

	x = sgamma(Shape) * Scale;
	return x / (Scale* (Shape  - 1));
}


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
//	double		Cost;

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
/*
	if(GenRandState(Rates->RandStates) < 0.5)
		PNode->Scale = GenRandState(Rates->RandStates);
	else
		PNode->Scale = 1 + (GenRandState(Rates->RandStates) * (PPMAXSCALE-1));
*/
	PNode->Scale = GenRandState(Rates->RandStates) * PPMAXSCALE;

	Plasty->NodeList = (PLASTYNODE**) AddToList(&Plasty->NoNodes, Plasty->NodeList, (void*)PNode);

//	Cost = log(PPSGammaPDF(PNode->Scale, PPALPHA, PPBETA));
//	Rates->LnHastings = Cost;
//	Rates->LnHastings += -PPCOST;

	Rates->LnHastings = 0;
	Rates->LogJacobion = 0;
}

void	DelPlastyNode(RATES *Rates, TREES *Trees, OPTIONS *Opt, int No)
{
	PLASTY		*Plasty;
	PLASTYNODE	**NList;
	int			Index;
	double		Cost;
	

	Plasty = Rates->Plasty;
	
	if(Plasty->NoNodes == 1)
	{
		free(Plasty->NodeList[0]);
		free(Plasty->NodeList);
		Plasty->NodeList = NULL;
		Plasty->NoNodes = 0;
		Cost = -log(PPSGammaPDF(Plasty->NodeList[0]->Scale, PPALPHA, PPBETA));
		Rates->LnHastings= Cost;
		return;
	}

	Cost = -log(PPSGammaPDF(Plasty->NodeList[No]->Scale, PPALPHA, PPBETA));

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
	

//	Rates->LnHastings = PPCOST;
	Rates->LnHastings = 0;
	Rates->LogJacobion = 0;
}

double	ChangePlastyRate(RATES *Rates, double Scale, double SD)
{
	double		Ret;
	RANDSTATES	*RS;

	RS = Rates->RandStates;

	Rates->LnHastings = CalcNormalHasting(Scale, SD);
	
	do
	{
		Ret = (nrand(RS) * SD) + Scale;
	} while((Ret <= 0) || (Ret > PPMAXSCALE));

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

void 	PPMoveNode(RATES *Rates, TREES *Trees, OPTIONS *Opt)
{
	PLASTY		*Plasty;
	PLASTYNODE	*PNode;
	NODE		N;
	int			No;

	Plasty = Rates->Plasty;

	No = rand() % Plasty->NoNodes;
	PNode = Plasty->NodeList[No];

	N = PNode->Node;

	if((GenRandState(Rates->RandStates) < 0.05) && (PNode->Node->Tip == FALSE))
	{
		if(PNode->Type == PPNODE)
			PNode->Type = PPBRANCH;
		else
			PNode->Type = PPNODE;
		return;
	}

	Plasty->NoTempList = 0;
	if((N->Ans != NULL) && (IsValidPPNode(N->Ans) == TRUE))
		Plasty->TempList[Plasty->NoTempList++] = N->Ans;

	if(N->Tip == FALSE)
	{
		MakeTNodeList(N->Left, Plasty->TempList, &Plasty->NoTempList);
		MakeTNodeList(N->Right, Plasty->TempList, &Plasty->NoTempList);
	}


	if(Plasty->NoTempList == 0)
		return;

	No = rand() % Plasty->NoTempList;
	PNode->Node = Plasty->TempList[No];

	return;
}

void	PPChangeScale(RATES *Rates, TREES *Trees, OPTIONS *Opt)
{
	PLASTY		*Plasty;
	PLASTYNODE	*Node;
	int			No;
	double		POld, PNew;

	Plasty = Rates->Plasty;

	No = rand() % Plasty->NoNodes;
	Node = Plasty->NodeList[No];

	POld = PPSGammaPDF(Node->Scale, PPALPHA, PPBETA);
	Node->Scale = ChangePlastyRate(Rates, Node->Scale, Opt->RateDev);
	PNew = PPSGammaPDF(Node->Scale, PPALPHA, PPBETA);

	Rates->LhPrior = log(PNew / POld);
}

int GetPlastyNode(int ID, PLASTY *Plasty)
{
	int Index;

	for(Index=0;Index<Plasty->NoNodes;Index++)
		if(Plasty->NodeList[Index]->Node->ID == ID)
			return Index;

	return -1;
}

void	PPAddRemove(RATES *Rates, TREES *Trees, OPTIONS *Opt)
{
	PLASTY *	Plasty;
	int			PNodeID;
	NODE		N;

	Plasty = Rates->Plasty;

	N = Plasty->ValidNode[rand() % Plasty->NoValidNode];
	PNodeID = GetPlastyNode(N->ID, Plasty);

	if(PNodeID == -1)
		PlastyAdd(Rates, Trees, Opt, N);		
	else
		DelPlastyNode(Rates, Trees, Opt, PNodeID);		

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
	TestGenRand(Opt, Trees, Rates);
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
/*
 double gengam(double a,double r)
           GENerates random deviates from GAMma distribution
                              Function
     Generates random deviates from the gamma distribution whose
     density is
          (A**R)/Gamma(R) * X**(R-1) * Exp(-A*X)
                              Arguments
     a --> Location parameter of Gamma distribution
     JJV   (a > 0)
     r --> Shape parameter of Gamma distribution
     JJV   (r > 0)
                              Method
*/

double	GammPDF(double x, double Shape, double Scale)
{
	double T1, T2;

	T1 = 1.0 / (Scale * gamma(Shape));

	T2 = pow((x / Scale), Shape-1);
	T2 = T2 * exp((-x)/Scale);

	return T1 * T2;
}



void	NormalTest(void)
{
	double	x;

	printf("\n");

	for(x=0;x<10;x+=0.01)
		printf("Normal\t%f\t%f\n", x, CalcNormalHasting(x, 2));

	exit(0);
}

void	TestGenRand(OPTIONS *Opt, TREES *Trees, RATES* Rates)
{
	int Index;
	double	 dev;

	dev = 2;

//	NormalTest();

	printf("\n");
	for(Index=0;Index<100000;Index++)
//		printf("%f\n", RandGamma(8.85, 5.28));
		printf("%f\n", RandGamma(PPALPHA, PPBETA));
		
//		printf("%f\n", RandGamma(1.5, 6));
	exit(0);
}


double	CalcPPPriors(RATES *Rates, OPTIONS *Opt)
{
	double		Ret;
	int			Index;
	PLASTY		*Plasty;
	PLASTYNODE	*PNode;

	Plasty = Rates->Plasty;
	Ret = 0;

	Ret = Plasty->NoNodes * PPCOST;
	return Ret;

	for(Index=0;Index<Plasty->NoNodes;Index++)
	{
		PNode = Plasty->NodeList[Index];
		Ret += log(PPSGammaPDF(PNode->Scale, PPALPHA, PPBETA));
	}

	return Ret;
}