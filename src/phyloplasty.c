#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "typedef.h"
#include "genlib.h"
#include "phyloplasty.h"
#include "matrix.h"
#include "RandLib.h"
#include "likelihood.h"
#include "trees.h"
#include "randdists.h"

void	TestGenRand(OPTIONS *Opt, TREES *Trees, RATES* Rates);

extern double gamma(double x);
extern double ndtr(double a);

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

	Ret = Ret / s;
	
	Ret = Ret * PPPRIORSCALE;

	return Ret;
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
		Ret[Index] = (double*)malloc(sizeof(double) * Trees->Tree[Index]->NoNodes);
		if(Ret[Index] == NULL)
			MallocErr();
	}

	for(Index=0;Index<Trees->NoOfTrees;Index++)
	{
		T = Trees->Tree[Index];
		for(NIndex=0;NIndex<T->NoNodes;NIndex++)
			Ret[Index][NIndex] = T->NodeList[NIndex]->Length;
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
	T = Trees->Tree[0];
	for(Index=0;Index<T->NoNodes;Index++)
	{
		if(IsValidPPNode(T->NodeList[Index]) == TRUE)
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

	T = Trees->Tree[0];
	for(Index=0;Index<T->NoNodes;Index++)
	{
		N = T->NodeList[Index];
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

	MakeValidNodes(Trees, Ret);

	Ret->ScaleBL = (double*)malloc(sizeof(double) * Ret->NoValidNode);
	if(Ret->ScaleBL == NULL)
		MallocErr();

	Ret->TempList = (NODE*)malloc(sizeof(NODE) * Ret->NoValidNode);

#ifdef PPUNIFORM
	Ret->Alpha = -1;
#else
	Ret->Alpha = PPALPHA;
#endif

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


void	PlastyAdd(RATES *Rates, TREES *Trees, OPTIONS *Opt, NODE N, int It)
{
	PLASTYNODE	*PNode;
	PLASTY		*Plasty;

	Plasty = Rates->Plasty;

	PNode = (PLASTYNODE*)malloc(sizeof(PLASTYNODE));
	if(PNode == NULL)
		MallocErr();

	PNode->NodeID = It;

	PNode->Node = N;

#ifdef PPBLO
	PNode->Type = PPBRANCH;
#else
	if(RandDouble(Rates->RS) < 0.5)
		PNode->Type = PPBRANCH;
	else
		PNode->Type = PPNODE;
#endif

#ifdef PPUNIFORM
	PNode->Scale = GenRandState(Rates->RandStates) * PPMAXSCALE;
#else
	PNode->Scale = RandGamma(PPALPHA, PPBETA);
#endif

	Plasty->NodeList = (PLASTYNODE**) AddToList(&Plasty->NoNodes, (void**)Plasty->NodeList, (void*)PNode);

	Rates->LnHastings = 0;
	Rates->LnJacobion = 0;
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
	
	Rates->LnHastings = 0;
	Rates->LnJacobion = 0;
}

double	ChangePlastyRate(RATES *Rates, double Scale, double SD)
{
	double		Ret;
	RANDSTATES	*RS;

	RS = Rates->RS;

	do
	{
//		Ret = (nrand(RS) * SD) + Scale;
		Ret = ((RandDouble(RS) * SD) - (SD / 2.0)) + Scale; 
		
#ifdef PPUNIFORM
	} while((Ret <= 0) || (Ret > PPMAXSCALE));
#else
	} while(Ret <= 0);
#endif

	return Ret;
}

void	PPChangeScale(RATES *Rates, TREES *Trees, OPTIONS *Opt)
{
	PLASTY		*Plasty;
	PLASTYNODE	*Node;
	int			No;

	Plasty = Rates->Plasty;

	No = RandUSLong(Rates->RS) % Plasty->NoNodes;
	Node = Plasty->NodeList[No];

	Rates->LnHastings = CalcNormalHasting(Node->Scale, Opt->VarRatesScaleDev);
//	Rates->LnHastings = 0;
	Node->Scale = ChangePlastyRate(Rates, Node->Scale, Opt->VarRatesScaleDev);
}


int		NodeScaled(int NID, PLASTY *Plasty)
{
	int	Index;

	for(Index=0;Index<Plasty->NoNodes;Index++)
	{
		if(Plasty->NodeList[Index]->Node->ID == NID)
			return TRUE;
	}
	
	return FALSE;
}


void	MakeTNodeList(PLASTY *Plasty, NODE N, NODE* List, int *Size)
{
	int Index;
	if( (IsValidPPNode(N) == TRUE) && (NodeScaled(N->ID, Plasty) == FALSE))
	{
		List[*Size] = N;
		(*Size)++;
	}

	if(N->Tip == TRUE)
		return;

	for(Index=0;Index<N->NoNodes;Index++)
		MakeTNodeList(Plasty, N->NodeList[Index], List, Size);
}

void	SawpNodeType(PLASTYNODE	*PNode)
{
#ifdef PPBLO
	return;
#endif

	if(PNode->Type == PPNODE)
		PNode->Type = PPBRANCH;
	else
		PNode->Type = PPNODE;
}

void 	PPMoveNode(RATES *Rates, TREES *Trees, OPTIONS *Opt)
{
	PLASTY		*Plasty;
	PLASTYNODE	*PNode;
	NODE		N;
	int			No;
	int			Index;

	Plasty = Rates->Plasty;

	No = RandUSLong(Rates->RS) % Plasty->NoNodes;
	PNode = Plasty->NodeList[No];

	N = PNode->Node;

	if( (RandDouble(Rates->RS) < 0.05) && 
		(PNode->Node->Tip == FALSE))
	{
		SawpNodeType(PNode);
		return;
	}

	Plasty->NoTempList = 0;
	if((N->Ans != NULL) && (IsValidPPNode(N->Ans) == TRUE) && (N != Trees->Tree[0]->Root))
		if(NodeScaled(N->Ans->ID, Plasty) == FALSE)
			Plasty->TempList[Plasty->NoTempList++] = N->Ans;

	if(N->Tip == FALSE)
	{
		for(Index=0;Index<N->NoNodes;Index++)
			MakeTNodeList(Plasty, N->NodeList[Index], Plasty->TempList, &Plasty->NoTempList);
	}

	if(Plasty->NoTempList == 0)
	{
		SawpNodeType(PNode);
		return;
	}

	No = RandUSLong(Rates->RS) % Plasty->NoTempList;
	PNode->Node = Plasty->TempList[No];

	return;
}

int GetPlastyNode(int ID, PLASTY *Plasty)
{
	int Index;

	for(Index=0;Index<Plasty->NoNodes;Index++)
		if(Plasty->NodeList[Index]->Node->ID == ID)
			return Index;

	return -1;
}

void	PPAddRemove(RATES *Rates, TREES *Trees, OPTIONS *Opt, int It)
{
	PLASTY *	Plasty;
	int			PNodeID;
	NODE		N;

	Plasty = Rates->Plasty;

	do
	{
		N = Plasty->ValidNode[RandUSLong(Rates->RS) % Plasty->NoValidNode];
	}while(N == Trees->Tree[0]->Root);

	PNodeID = GetPlastyNode(N->ID, Plasty);

	if(PNodeID == -1)
		PlastyAdd(Rates, Trees, Opt, N, It);
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

	P1->Alpha = P2->Alpha;
}

void	PlastyNode(NODE N, PLASTYNODE *P)
{
	int Index;

	N->Length = N->Length * P->Scale;
	if(P->Type == PPBRANCH)
		return;

	if(N->Tip == TRUE)
		return;

	for(Index=0;Index<N->NoNodes;Index++)
		PlastyNode(N->NodeList[Index], P);
}

void	Plasty(OPTIONS *Opt, TREES *Trees, RATES *Rates, int Normalise)
{
	int Index;
	int	TNo;
	TREE *T;
	PLASTY *P;
	double SumBL, Scale;

	P = Rates->Plasty;
	TNo = Rates->TreeNo;
	T = Trees->Tree[TNo];

	for(Index=0;Index<T->NoNodes;Index++)
		T->NodeList[Index]->Length = P->TrueBL[TNo][Index];

	if(Normalise == TRUE)
		SumBL = SumNodeBL(T->Root);

	for(Index=0;Index<P->NoNodes;Index++)
		PlastyNode(P->NodeList[Index]->Node, P->NodeList[Index]);

	if(Normalise == FALSE)
		return;

	Scale = SumBL / SumNodeBL(T->Root);
	ScaleSubTree(T->Root, Scale);
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
		fprintf(Opt->PPTree, "\t\t%d\t%s,\n", Trees->Taxa[Index]->No, Trees->Taxa[Index]->Name);
	fprintf(Opt->PPTree, "\t\t%d\t%s\n\t\t;\n", Trees->Taxa[Index]->No, Trees->Taxa[Index]->Name);

	fflush(Opt->PPTree);
}


void	RecPrintPPNodes(FILE *Out, NODE N)
{
	int Index;

	if(N->Tip == TRUE)
	{
		fprintf(Out, "%d\t", N->Taxa->No);
		return;
	}

	for(Index=0;Index<N->NoNodes;Index++)
		RecPrintPPNodes(Out, N->NodeList[Index]);

}


void	PPLogFileHeader(OPTIONS *Opt, TREES *Trees, RATES *Rates)
{
	TREE	*T;
	int		Index;
	NODE	N;
	PLASTY	*P;
	int		No;

	T = Trees->Tree[0];
	P = Rates->Plasty;

	fprintf(Opt->PPLog, "%d\n", Trees->NoOfTaxa);
	for(Index=0;Index<Trees->NoOfTaxa;Index++)
		fprintf(Opt->PPLog, "%d\t%s\n", Trees->Taxa[Index]->No, Trees->Taxa[Index]->Name);

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

	fprintf(Opt->PPLog, "It\tLh\tLh + Prior\tNo Pram\tAlpha\tSigma\tAlpha Scale Prior\t");
	fprintf(Opt->PPLog, "Node ID\tScale\tCreat It\tNode / Branch\t");

	fprintf(Opt->PPLog, "\n");

	fflush(Opt->PPLog);
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
	int Index;

	if(N->Tip == TRUE)
	{
		List[*Size] = N->Taxa->No;
		(*Size)++;
		return;
	}
	
	for(Index=0;Index<N->NoNodes;Index++)
		GetNodeIDList(N->NodeList[Index], Size, List);

}

void	InitPPFiles(OPTIONS *Opt, TREES *Trees, RATES* Rates)
{
//	TestGenRand(Opt, Trees, Rates);
	InitPPTreeFile(Opt, Trees);
	IntiPPLogFile(Opt, Trees, Rates);
}

void	PrintPPNode(FILE *Out, NODE N)
{
	int Index;

	if(N->Tip == TRUE)
	{
		fprintf(Out, "%d:%f", N->Taxa->No, N->Length);
		return;
	}

	fprintf(Out, "(");
	for(Index=0;Index<N->NoNodes-1;Index++)
	{
		PrintPPNode(Out, N->NodeList[Index]);
		fprintf(Out, ",");
	}
	PrintPPNode(Out, N->NodeList[Index]);
	fprintf(Out, "):%f", N->Length);
}

void	PrintPPTree(OPTIONS *Opt, TREES *Trees, RATES *Rates, int It)
{
	PLASTY	*P;
	TREE	*T;
	int		Index;

	P = Rates->Plasty;
	T = Trees->Tree[Rates->TreeNo];

	Plasty(Opt, Trees, Rates, NORM_TRANSFORMS);

	fprintf(Opt->PPTree, "\tTree T_%010d_%d = (", It, P->NoNodes);

	for(Index=0;Index<T->Root->NoNodes-1;Index++)
	{
		PrintPPNode(Opt->PPTree, T->Root->NodeList[Index]);
		fprintf(Opt->PPTree, ",");
	}

	PrintPPNode(Opt->PPTree, T->Root->NodeList[Index]);
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
	double		Sigma;

	P = Rates->Plasty;

	Out = Opt->PPLog;

	fprintf(Out, "%d\t%f\t%f\t%d\t", It, Rates->Lh, Rates->Lh + Rates->LhPrior, P->NoNodes);

	if(Opt->Model == M_CONTRAST_STD)
	{
		Sigma = Rates->Contrast->SigmaMat->me[0][0];
		fprintf(Out, "%f\t%f\t%f\t", Rates->Contrast->Alpha[0], Sigma, P->Alpha);
	}

	
	if(Opt->Model == M_CONTRAST_FULL)
	{
		Sigma = Rates->Contrast->Sigma[0];
		fprintf(Out, "%f\t%f\t%f\t", Rates->Contrast->Alpha[0], Sigma, P->Alpha);
	}

	if(Opt->Model == M_CONTRAST_REG)
	{
		fprintf(Out, "%f\t%f\t%f\t", Rates->Contrast->RegAlpha, -1, P->Alpha);
	}



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

void	ChangePPHyperPrior(RATES *Rates, OPTIONS *Opt)
{
	double	NAlpha;
	PLASTY	*Plasty;

	Plasty = Rates->Plasty;

	Rates->LnHastings = CalcNormalHasting(Plasty->Alpha - 1, PPALPHASCLAE);
	
	do
	{
		NAlpha = RandNormal(Rates->RS, Plasty->Alpha, PPALPHASCLAE);
	} while(NAlpha <= 1);

	Plasty->Alpha = NAlpha;
}

double	CalcPPPriors(RATES *Rates, OPTIONS *Opt)
{
	double		Ret;
	int			Index;
	PLASTY		*Plasty;
	PLASTYNODE	*PNode;
	double		PVal;

	Plasty = Rates->Plasty;
	Ret = 0;

#ifdef PPUNIFORM
	Ret = Plasty->NoNodes * -PPUNICOST;
	return Ret;
#endif

	for(Index=0;Index<Plasty->NoNodes;Index++)
	{
		PNode = Plasty->NodeList[Index];
		PVal = PPSGammaPDF(PNode->Scale, Plasty->Alpha, PPBETA);
		Ret += log(PVal);
	}

	return Ret;
}


int		SeenNodeID(int NID, PLASTY *Plasty, int Size)
{
	int	Index;

	for(Index=0;Index<Size;Index++)
		if(Plasty->NodeList[Index]->Node->ID == NID)
			return TRUE;

	return FALSE;
}

void	CheckPlasty(RATES *Rates, TREES *Trees, OPTIONS *Opt)
{
	PLASTY	*Plasty;
	int		Index;
	
	Plasty = Rates->Plasty;

	for(Index=0;Index<Plasty->NoNodes;Index++)
	{
		if(SeenNodeID(Plasty->NodeList[Index]->Node->ID, Plasty, Index) == TRUE)
		{
			printf("Dup node in list %d.\n", Plasty->NodeList[Index]->Node->ID);
			exit(0);
		}
	}	
}

/*
PLASTY*	CreatPlasty(RATES *Rates, TREES *Trees, OPTIONS *Opt) { return NULL; }
void	FreePlasty(PLASTY* Plasty) { }

void	PPAddRemove(RATES *Rates, TREES *Trees, OPTIONS *Opt, int It) { }
void	PPChangeScale(RATES *Rates, TREES *Trees, OPTIONS *Opt) { }
void	PPMoveNode(RATES *Rates, TREES *Trees, OPTIONS *Opt) { }


void	PlastyCopy(RATES *R1, RATES *R2) { }
void	Plasty(OPTIONS *Opt, TREES *Trees, RATES *Rates) { }

void	InitPPFiles(OPTIONS *Opt, TREES *Trees, RATES *Rates) { }
void	PrintPPOutput(OPTIONS *Opt, TREES *Trees, RATES *Rates, int It) { }

double	CalcPPPriors(RATES *Rates, OPTIONS *Opt) { return 0;}
void	ChangePPHyperPrior(RATES *Rates, OPTIONS *Opt) { }

void	CheckPlasty(RATES *Rates, TREES *Trees, OPTIONS *Opt) { }
*/
