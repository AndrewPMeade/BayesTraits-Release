#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "typedef.h"
#include "genlib.h"
#include "VarRates.h"
#include "matrix.h"
#include "RandLib.h"
#include "likelihood.h"
#include "trees.h"
#include "randdists.h"
#include "RJLocalScalar.h"
#include "priors.h"
#include "ContrastsTrans.h"


extern double gamma(double x);
extern double ndtr(double a);

int		UseNonParametricMethods(OPTIONS *Opt)
{
	int Index;

	if(Opt->UseVarRates == TRUE)
		return TRUE;

	for(Index=0;Index<NO_RJ_LOCAL_SCALAR;Index++)
		if(Opt->UseRJLocalScalar[Index] == TRUE)
			return TRUE;

	return FALSE;
}

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

int		IsValidVarRatesNode(NODE N, RJ_VARRATE_TYPE	Type)
{
	PART *Part;
	
	if(N == NULL)
		return FALSE;

	Part = N->Part;

	if(N->Length == 0)
		return FALSE;

	if(Type == VR_BL || Type == VR_NODE)
	{
		if(N->Tip == TRUE && Type == VR_NODE)
			return FALSE;

		if(N->Ans == NULL)
			return FALSE;

		return TRUE;
	}

	if(Part->NoTaxa >= MIN_TAXA_VR_TRANS)
		return TRUE;

	return FALSE;
}
/*
int		FindNoValidNodes(OPTIONS *Opt, TREES *Trees)
{
	TREE	*T;
	int		Index;
	int		Ret;

	Ret = 0;
	T = Trees->Tree[0];
	for(Index=0;Index<T->NoNodes;Index++)
	{
		if(IsValidVarRatesNode(Opt, T->NodeList[Index]) == TRUE)
			Ret++;
	}

	return Ret;
}

void	MakeValidNodes(OPTIONS *Opt, TREES *Trees, PLASTY* Plasty)
{
	TREE *T;
	int	Index;
	NODE N;

	Plasty->NoValidNode = FindNoValidNodes(Opt, Trees);

	Plasty->ValidNode = (NODE*)malloc(sizeof(NODE) * Plasty->NoValidNode);
	if(Plasty->ValidNode == NULL)
		MallocErr();

	Plasty->NoValidNode = 0;

	T = Trees->Tree[0];
	for(Index=0;Index<T->NoNodes;Index++)
	{
		N = T->NodeList[Index];
		if(IsValidVarRatesNode(Opt, N) == TRUE)
			Plasty->ValidNode[Plasty->NoValidNode++] = N;
	}
}
*/

PLASTY*	CreatVarRates(RATES *Rates, TREES *Trees, OPTIONS *Opt)
{
	PLASTY* Ret;


	Ret = (PLASTY*)malloc(sizeof(PLASTY));
	if(Ret == NULL)
		MallocErr();

	Ret->NoTrees = Trees->NoOfTrees;
	Ret->NoNodes = 0;
	Ret->NodeList= NULL;

	Ret->TrueBL = MakeTrueBL(Trees);

//	MakeValidNodes(Opt, Trees, Ret);

//	Ret->ScaleBL = (double*)malloc(sizeof(double) * Ret->NoValidNode);
//	if(Ret->ScaleBL == NULL)
//		MallocErr();

	Ret->TempList = (NODE*)malloc(sizeof(NODE) * Trees->MaxNodes);
	if(Ret->TempList == NULL)
		MallocErr();

#ifdef PPUNIFORM
	Ret->Alpha = -1;
#else
	Ret->Alpha = PPALPHA;
#endif

	return Ret;
}

void	FreeVarRates(PLASTY* Plasty)
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

	free(Plasty->TempList);
//	free(Plasty->ScaleBL);
//	free(Plasty->ValidNode);
	free(Plasty);
}


void	ScaledGammaTest(RANDSTATES *RS)
{
	double C, N;
	double pC, pN, Heat;
	int Index;

	C = 1;
	for(Index=0;Index<10000000;Index++)
	{
			
		do
		{
			N = C;
			N += (RandDouble(RS) * 0.5) - 0.25;
//			N = C + (RandDouble(RS) - 0.5);
		} while(N < 0);

		pC = log(PPSGammaPDF(C, PPALPHA, PPBETA));
		pN = log(PPSGammaPDF(N, PPALPHA, PPBETA));

		Heat = pN - pC;

		if(log(RandDouble(RS)) <= Heat)
		{
			C = N;
			pC = pN;
		}
		
		if(Index % 100 == 0)
			printf("%d\t%f\t%f\n", Index, pC, C);
	}

	exit(0);

}

void	ScaleTest()
{
	double x;

	for(x=0.0000001;x<20;x+=0.0001)
		printf("%f\t%f\n", x, log(PPSGammaPDF(x, PPALPHA, PPBETA)));
	exit(0);
}

PLASTYNODE*	CreateVarRatesNode(long long It, NODE N)
{
	PLASTYNODE	*Ret;
	
	Ret = (PLASTYNODE*)malloc(sizeof(PLASTYNODE));
	if(Ret == NULL)
		MallocErr();

	Ret->NodeID = It;

	Ret->Node = N;

	Ret->Scale = -1;
	Ret->Type = VR_BL;

	return Ret;
}

void		FreeVarRatesNode(PLASTYNODE* VarRatesNode)
{
	free(VarRatesNode);
}

//RJ_VARRATE_TYPE	Get

RJ_VARRATE_TYPE	GetVarRatesType(RANDSTATES *RS, SCHEDULE *Shed)
{
	int Pos;

	Pos = RandInProportion(RS, Shed->FreqVarRatesOp, Shed->NoVarRatesOp);

	return Shed->VarRatesOp[Pos];
}

void	SetVRNodeBLRates(PLASTYNODE *PNode, RANDSTATES *RS)
{
#ifdef PPUNIFORM
	PNode->Scale = RandDouble(RS) * PPMAXSCALE;
#else
	PNode->Scale = RandGamma(PPALPHA, PPBETA);
#endif
}

void	SetVRScalar(OPTIONS *Opt, RATES *Rates, PLASTYNODE *PNode)
{
	PRIORS *Prior;

	Prior = GetPriorFromRJRatesScalar(Opt, PNode->Type);

	if(PNode->Type == VR_LAMBDA)
		PNode->Scale = RandDouble(Rates->RS);
	else
		PNode->Scale = RandFromPrior(Rates->RS, Prior);
//	PNode->Scale = 0;
}

int		CountPlasyID(long long ID, PLASTY *Plasty)
{
	int Ret, Index;

	Ret = 0;
	for(Index=0;Index<Plasty->NoNodes;Index++)
		if(Plasty->NodeList[Index]->Node->ID == ID)
			Ret++;

	return Ret;
}

void	CheckPlasyNodes(PLASTY *Plasty)
{
	int No, Index;
	PLASTYNODE	*PNode;

	for(Index=0;Index<Plasty->NoNodes;Index++)
	{
		No = CountPlasyID(Plasty->NodeList[Index]->Node->ID, Plasty);
		if(No != 1)
		{
			printf("Multiple hits of the same node.\n");
			exit(0);
		}

		printf("%d\n", No);
	}


	for(Index=0;Index<Plasty->NoNodes;Index++)
	{
		PNode = Plasty->NodeList[Index];
		if(PNode->Type == VR_LAMBDA && PNode->Scale > 1.0)
			printf("Err\n");
	}
}

void	VarRatesAddNode(RATES *Rates, TREES *Trees, OPTIONS *Opt, RJ_VARRATE_TYPE Type, NODE N, long long It)
{
	PLASTYNODE	*PNode;
	PLASTY		*Plasty;

	Plasty = Rates->Plasty;

	PNode = CreateVarRatesNode(It, N);

	PNode->Type = Type;

	if(PNode->Type == VR_BL || PNode->Type == VR_NODE)
		SetVRNodeBLRates(PNode, Rates->RS);
	else
		SetVRScalar(Opt, Rates, PNode);

	Plasty->NodeList = (PLASTYNODE**) AddToList(&Plasty->NoNodes, (void**)Plasty->NodeList, (void*)PNode);

	Rates->LnHastings = 0;
	Rates->LnJacobion = 0;
}

void	VarRatesDelNode(RATES *Rates, TREES *Trees, OPTIONS *Opt, int No)
{
	PLASTY		*Plasty;
	PLASTYNODE	**NList;
	int			Index;
	
	Plasty = Rates->Plasty;
	
	if(Plasty->NoNodes == 1)
	{
		FreeVarRatesNode(Plasty->NodeList[0]);
		free(Plasty->NodeList);
		Plasty->NodeList = NULL;
		Plasty->NoNodes = 0;
		return;
	}
	
	FreeVarRatesNode(Plasty->NodeList[No]);
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

double	ChangePlastyRate(RANDSTATES	*RS, double Scale, double SD)
{
	double		Ret;

#ifdef PPUNIFORM
		if(SD > PPMAXSCALE)
			SD = PPMAXSCALE;
#endif

	do
	{
		Ret = ((RandDouble(RS) * SD) - (SD / 2.0)) + Scale; 

#ifdef PPUNIFORM
	} while((Ret <= 0) || (Ret > PPMAXSCALE));
#else
	} while(Ret <= 0);
#endif

	return Ret;
}

double	ChangePlastyRateLambda(RANDSTATES	*RS, double Scale, double SD)
{
	double Ret;

	do
	{
		Ret = ChangePlastyRate(RS, Scale, SD);
	} while(Ret > 1.0);

	return Ret;
}

void	ChangeVarRatesScale(RATES *Rates, TREES *Trees, OPTIONS *Opt)
{
	PLASTY		*Plasty;
	PLASTYNODE	*Node;
	int			No;

	Plasty = Rates->Plasty;

	No = RandUSLong(Rates->RS) % Plasty->NoNodes;
	Node = Plasty->NodeList[No];

	Rates->LnHastings = CalcNormalHasting(Node->Scale, Opt->VarRatesScaleDev);
//	Rates->LnHastings = 0;

	if(Node->Type == VR_LAMBDA)
		Node->Scale = ChangePlastyRateLambda(Rates->RS, Node->Scale, Opt->VarRatesScaleDev);
	else
		Node->Scale = ChangePlastyRate(Rates->RS, Node->Scale, Opt->VarRatesScaleDev);
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

int		IsVarRateTypeRate(RJ_VARRATE_TYPE Type)
{
	if(Type == 	VR_NODE || Type == VR_BL)
		return TRUE;

	return FALSE;
}


int GetPlastyNode(int ID, PLASTY *Plasty, RJ_VARRATE_TYPE Type)
{
	int Index;

	for(Index=0;Index<Plasty->NoNodes;Index++)
	{

#ifdef VARRATES_ONE_OP_PER_NODE
		if(	Plasty->NodeList[Index]->Node->ID == ID)
			return Index;
#else
		if(	Plasty->NodeList[Index]->Node->ID == ID && 
			Plasty->NodeList[Index]->Type == Type)
		return Index;
#endif
	}

	return -1;
}

int		ValidMoveNode(PLASTY *Plasty, NODE N, RJ_VARRATE_TYPE Type)
{
	int ID;

	if(N == NULL)
		return FALSE;
	
	if(IsValidVarRatesNode(N, Type) == FALSE)
		return FALSE;
	
	ID = GetPlastyNode(N->ID, Plasty, Type);

	if(ID == -1)
		return TRUE;

	return FALSE;
}

void	MakeTNodeList(OPTIONS *Opt, PLASTY *Plasty, RJ_VARRATE_TYPE Type, NODE N, NODE* List, int *Size)
{
//	int Index;

	if(ValidMoveNode(Plasty, N, Type) == TRUE)
	{
		List[*Size] = N;
		(*Size)++;
	}

	if(N->Tip == TRUE)
		return;

//	for(Index=0;Index<N->NoNodes;Index++)
//		MakeTNodeList(Opt, Plasty,Type, N->NodeList[Index], List, Size);
}

void 	VarRatesMoveNode(RATES *Rates, TREES *Trees, OPTIONS *Opt)
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

	Plasty->NoTempList = 0;
	if(ValidMoveNode(Plasty, N->Ans, PNode->Type) == TRUE)
		Plasty->TempList[Plasty->NoTempList++] = N->Ans;

	if(N->Tip == FALSE)
	{
		for(Index=0;Index<N->NoNodes;Index++)
			MakeTNodeList(Opt, Plasty, PNode->Type, N->NodeList[Index], Plasty->TempList, &Plasty->NoTempList);
	}

	if(Plasty->NoTempList == 0)
		return;
	
	No = RandUSLong(Rates->RS) % Plasty->NoTempList;
	PNode->Node = Plasty->TempList[No];

	return;
}


NODE	GetVarRatesNode(RATES *Rates, TREES *Trees, RJ_VARRATE_TYPE	Type)
{
	PLASTY		*Plasty;
	TREE		*Tree;
	NODE		N;
	int			Pos;
	
	Tree = Trees->Tree[Rates->TreeNo];
	Plasty = Rates->Plasty;

	do
	{
		Pos = RandUSInt(Rates->RS) % Tree->NoNodes;
		N = Tree->NodeList[Pos];
	}while(IsValidVarRatesNode(N, Type) == FALSE); 
	
	return N;
}

void	VarRatesAddRemove(RATES *Rates, TREES *Trees, OPTIONS *Opt, SCHEDULE *Shed, long long It)
{
	PLASTY *	Plasty;
	int			PNodeID;
	NODE		N;
	RJ_VARRATE_TYPE		Type;
	


	Type = GetVarRatesType(Rates->RS, Shed);
		
	Plasty = Rates->Plasty;
	
	N = GetVarRatesNode(Rates, Trees, Type);

	PNodeID = GetPlastyNode(N->ID, Plasty, Type);

	if(PNodeID == -1)
		VarRatesAddNode(Rates, Trees, Opt, Type, N, It);
	else
		VarRatesDelNode(Rates, Trees, Opt, PNodeID);	

	CheckPlasyNodes(Plasty);
}


PLASTYNODE *CloneVarRatesNode(PLASTYNODE *N2)
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
			FreeVarRatesNode(P->NodeList[Index]);
		free(P->NodeList);
	}

	P->NodeList = NULL;
	P->NoNodes = 0;
}

void	VarRatesCopy(RATES *R1, RATES *R2)
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
		NList[Index] = CloneVarRatesNode(P2->NodeList[Index]);

	BlankPlasty(P1);

	P1->NodeList = NList;
	P1->NoNodes = P2->NoNodes;

	P1->Alpha = P2->Alpha;
}

void	ScaleNode(NODE N,  PLASTYNODE *P)
{
	int Index;

	N->Length = N->Length * P->Scale;

	if(N->Tip == TRUE)
		return;

	for(Index=0;Index<N->NoNodes;Index++)
		ScaleNode(N->NodeList[Index], P);
}

void	VarRatesNode(TREE *Tree, NODE N, PLASTYNODE *P)
{
//	int Index;
	int Norm;

//	Norm = FALSE;
	Norm = TRUE;

	if(P->Type == VR_BL)
		N->Length = N->Length * P->Scale;

	if(P->Type == VR_NODE)
		ScaleNode(N,  P);	

	if(P->Type == VR_KAPPA)
		TransContNodeKappa(N, P->Scale, Norm);

	if(P->Type == VR_LAMBDA)
	{
		SetTreeDistToRoot(Tree);
		TransContNodeLambda(N, P->Scale, Norm);
	}

	if(P->Type == VR_DELTA)
		TransContNodeDelta(N, P->Scale, Norm);

	if(P->Type == VR_OU)
	{
		SetTreeDistToRoot(Tree);
		TransContNodeOU(N, P->Scale, Norm);
	}
}

void	CheckVarRatesData(OPTIONS *Opt, TREES *Trees, RATES *Rates)	
{
	PLASTY *P;
	NODE N;
	int Index;

	P = Rates->Plasty;

	for(Index=0;Index<P->NoNodes;Index++)
	{
		N = P->NodeList[Index]->Node;

		if(N->Part->NoTaxa < MIN_TAXA_VR_TRANS)
		{
			printf("err.\n");
			exit(1);
		}

		if(IsValidVarRatesNode(P->NodeList[Index]->Node, P->NodeList[Index]->Type) == FALSE)
		{
			printf("err2.\n");
			exit(1);
		}
	
	}
}

void	VarRatesTree(OPTIONS *Opt, TREES *Trees, RATES *Rates, int Normalise)
{
	int Index;
	int	TNo;
	TREE *Tree;
	PLASTY *P;
	double SumBL, Scale;

//	CheckVarRatesData(Opt, Trees, Rates);

	P = Rates->Plasty;
	TNo = Rates->TreeNo;
	Tree = Trees->Tree[TNo];

//	This is done before. 
//	for(Index=0;Index<Tree->NoNodes;Index++)
//		Tree->NodeList[Index]->Length = P->TrueBL[TNo][Index];

	if(Normalise == TRUE)
		SumBL = SumNodeBL(Tree->Root);

	for(Index=0;Index<P->NoNodes;Index++)
		VarRatesNode(Tree, P->NodeList[Index]->Node, P->NodeList[Index]);

	if(Normalise == FALSE)
		return;

	Scale = SumBL / SumNodeBL(Tree->Root);
	ScaleSubTree(Tree->Root, Scale);
}

void	InitVarRatesTreeFile(OPTIONS *Opt, TREES *Trees)
{
	char	*Buffer;
	int		Index;
	
	Buffer = (char*)malloc(sizeof(char) * (strlen(Opt->LogFN) + BUFFERSIZE));
	if(Buffer == NULL)
		MallocErr();
		
	sprintf(Buffer, "%s.VarRates.trees", Opt->LogFN);
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

/*
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

	fprintf(Opt->PPLog, "It\tLh\tLh + Prior\tNo Pram\tAlpha\tSigma^2\tAlpha Scale Prior\t");
	fprintf(Opt->PPLog, "Node ID\tScale\tCreat It\tNode / Branch\t");

	fprintf(Opt->PPLog, "\n");

	fflush(Opt->PPLog);
}
*/

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
/*
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
	*/

	fprintf(Opt->PPLog, "%d\n", T->NoNodes);
	for(Index=0;Index<T->NoNodes;Index++)
	{
	//	N = P->ValidNode[Index];
		N = T->NodeList[Index];
		No = 0;
		CTaxaBelow(N, &No);
		fprintf(Opt->PPLog, "%d\t%f\t%d\t", N->ID, N->Length, No);
		RecPrintPPNodes(Opt->PPLog, N);
		fprintf(Opt->PPLog, "\n");
	}

	fprintf(Opt->PPLog, "It\tLh\tLh + Prior\tNo Pram\tAlpha\tSigma^2\tAlpha Scale Prior\t");
	fprintf(Opt->PPLog, "Node ID\tScale\tCreat It\tNode / Branch\t");

	fprintf(Opt->PPLog, "\n");

	fflush(Opt->PPLog);
}

void	IntiVarRatesLogFile(OPTIONS *Opt, TREES *Trees, RATES *Rates)
{
	char	*Buffer;

	Buffer = (char*)malloc(sizeof(char) * (strlen(Opt->LogFN) + BUFFERSIZE));
	if(Buffer == NULL)
		MallocErr();

	sprintf(Buffer, "%s.VarRates.txt", Opt->LogFN);
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

void	InitVarRatesFiles(OPTIONS *Opt, TREES *Trees, RATES* Rates)
{
	InitVarRatesTreeFile(Opt, Trees);
	IntiVarRatesLogFile(Opt, Trees, Rates);
}

void	FinishVarRatesFiles(OPTIONS *Opt)
{
	fprintf(Opt->PPTree, "end;");
	fclose(Opt->PPTree);
	fclose(Opt->PPLog);
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

void	PrintPPTree(OPTIONS *Opt, TREES *Trees, RATES *Rates, long long It)
{
	PLASTY	*P;
	TREE	*Tree;
	int		Index;

	P = Rates->Plasty;
	Tree = Trees->Tree[Rates->TreeNo];

	ReSetBranchLength(Tree);
	VarRatesTree(Opt, Trees, Rates, NORMALISE_TREE_CON_SCALING);

	fprintf(Opt->PPTree, "\tTree VarRates_%llu_%d = (", It, P->NoNodes);

	for(Index=0;Index<Tree->Root->NoNodes-1;Index++)
	{
		PrintPPNode(Opt->PPTree, Tree->Root->NodeList[Index]);
		fprintf(Opt->PPTree, ",");
	}

	PrintPPNode(Opt->PPTree, Tree->Root->NodeList[Index]);
	fprintf(Opt->PPTree, ");\n");

	fflush(Opt->PPTree);
}

void	OutputVarRatesType(FILE *Out, RJ_VARRATE_TYPE Type)
{
	if(Type == VR_NODE)
		fprintf(Out, "Node\t");
	
	if(Type == VR_BL)
		fprintf(Out, "Branch\t");

	if(Type == VR_KAPPA)
		fprintf(Out, "Kappa\t");

	if(Type == VR_LAMBDA)
		fprintf(Out, "Lambda\t");

	if(Type == VR_DELTA)
		fprintf(Out, "Delta\t");

	if(Type == VR_OU)
		fprintf(Out, "OU\t");
}

void	LogPPResults(OPTIONS *Opt, TREES *Trees, RATES *Rates, long long It)
{
	FILE		*Out;
	PLASTY		*P;
	int			Index;
	PLASTYNODE	*PNode;
	NODE		N;
	double		Sigma;

	P = Rates->Plasty;
//	CheckPlasyNodes(P);

	Out = Opt->PPLog;

	fprintf(Out, "%lld\t%f\t%f\t%d\t", It, Rates->Lh, Rates->Lh + Rates->LhPrior, P->NoNodes);

	if(Opt->Model == M_CONTRAST_CORREL)
	{
		Sigma = Rates->Contrast->SigmaMat->me[0][0];
		fprintf(Out, "%f\t%f\t%f\t", Rates->Contrast->Alpha[0], Sigma, P->Alpha);
	}

	if(Opt->Model == M_CONTRAST)
	{
		Sigma = Rates->Contrast->Sigma[0];
		fprintf(Out, "%f\t%f\t%f\t", Rates->Contrast->Alpha[0], Sigma, P->Alpha);
	}

	if(Opt->Model == M_CONTRAST_REG)
	{
		Sigma = Rates->Contrast->GlobalVar;
		fprintf(Out, "%f\t%f\t%f\t", Rates->Contrast->RegAlpha, Sigma, P->Alpha);
	}
	
//	fprintf(Out, "\n");
	for(Index=0;Index<P->NoNodes;Index++)
	{
		PNode = P->NodeList[Index];
		N = PNode->Node;

		fprintf(Out, "%d\t", N->ID);
		fprintf(Out, "%f\t", PNode->Scale);
		fprintf(Out, "%llu\t", PNode->NodeID);

		OutputVarRatesType(Out, PNode->Type);
//		fprintf(Out, "\n");
	}
	fprintf(Out, "\n");
//	exit(0);
	fflush(Out);
}	

void	PrintVarRatesOutput(OPTIONS *Opt, TREES *Trees, RATES *Rates, long long It)
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

void	ChangeVarRatesHyperPrior(RATES *Rates, OPTIONS *Opt)
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

double	CalcVarRatesPriors(RATES *Rates, OPTIONS *Opt)
{
	double		Ret;
	int			Index;
	PLASTY		*Plasty;
	PLASTYNODE	*PNode;
	double		PVal;

//	return Rates->Plasty->NoNodes * -2.302585093;

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

*/
