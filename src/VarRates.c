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
#include "TransformTree.h"

#include <gsl/gsl_cdf.h>

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

TRANSFORM_TYPE	StrToVarRatesType(char *Str)
{
	MakeLower(Str);

	if(strcmp("node", Str) == 0)
		return VR_NODE;

	if(strcmp("branch", Str) == 0)
		return VR_BL;

	if(strcmp("kappa", Str) == 0)
		return VR_KAPPA;

	if(strcmp("lambda", Str) == 0)
		return VR_LAMBDA;

	if(strcmp("delta", Str) == 0)
		return VR_DELTA;

	if(strcmp("ou", Str) == 0)
		return VR_OU;


	printf("uknown varaible rate type %s\n", Str); 
	exit(0);
}

char* VarRatesTypeToStr(TRANSFORM_TYPE Type)
{
	if(Type == VR_NODE)
		return "Node";

	if(Type == VR_BL)
		return "Branch";

	if(Type == VR_KAPPA)
		return "Kappa";

	if(Type == VR_LAMBDA)
		return "Lambda";

	if(Type == VR_DELTA)
		return "Delta";

	if(Type == VR_OU)
		return "OU";

	printf("%s::%d unkonwn RJ Variable type\n", __FILE__, __LINE__);
	exit(0);
	return NULL;
}

double	CalcNormalHasting(double x, double SD)
{
	double Ret;

	Ret = gsl_cdf_gaussian_P(x, SD);
//	Ret = ndtr(x/SD);

	return log(Ret);
}

double	RandGamma(double Shape, double Scale)
{
	double x;

	x = sgamma(Shape) * Scale;
	return x / (Scale* (Shape  - 1));
}


int		IsValidVarRatesNode(NODE N, TRANSFORM_TYPE	Type)
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

VARRATES*	CreatVarRates(RATES *Rates, TREES *Trees, OPTIONS *Opt)
{
	VARRATES* Ret;
	
	Ret = (VARRATES*)malloc(sizeof(VARRATES));
	if(Ret == NULL)
		MallocErr();

	Ret->NoNodes = 0;
	Ret->NodeList= NULL;

	Ret->TempList = (NODE*)malloc(sizeof(NODE) * Trees->MaxNodes);
	if(Ret->TempList == NULL)
		MallocErr();

#ifdef PPUNIFORM
	Ret->Alpha = -1;
#else
	Ret->Alpha = VARRATES_ALPHA;
#endif

	return Ret;
}

void	FreeVarRates(VARRATES* Plasty)
{
	int Index;

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


VAR_RATES_NODE*	CreateVarRatesNode(long long It, NODE N)
{
	VAR_RATES_NODE	*Ret;
	
	Ret = (VAR_RATES_NODE*)malloc(sizeof(VAR_RATES_NODE));
	if(Ret == NULL)
		MallocErr();

	Ret->NodeID = It;

	Ret->Node = N;

	Ret->Scale = -1;
	Ret->Type = VR_BL;

	return Ret;
}

void		FreeVarRatesNode(VAR_RATES_NODE* VarRatesNode)
{
	free(VarRatesNode);
}

//RJ_VARRATE_TYPE	Get

TRANSFORM_TYPE	GetVarRatesType(RANDSTATES *RS, SCHEDULE *Shed)
{
	int Pos;

	Pos = RandInProportion(RS, Shed->FreqVarRatesOp, Shed->NoVarRatesOp);

	return Shed->VarRatesOp[Pos];
}

void	SetVRNodeBLRates(VAR_RATES_NODE *PNode, RANDSTATES *RS)
{
#ifdef PPUNIFORM
	PNode->Scale = RandDouble(RS) * PPMAXSCALE;
#else
	PNode->Scale = RandGamma(VARRATES_ALPHA, VARRATES_BETA);
#endif
}

void	SetVRScalar(OPTIONS *Opt, RATES *Rates, VAR_RATES_NODE *PNode)
{
	PRIOR *Prior;

	Prior = GetPriorFromRJRatesScalar(Opt, PNode->Type);
	
	PNode->Scale = RandFromPrior(Rates->RNG, Prior);
}

int		CountPlasyID(long long ID, VARRATES *Plasty)
{
	int Ret, Index;

	Ret = 0;
	for(Index=0;Index<Plasty->NoNodes;Index++)
		if(Plasty->NodeList[Index]->Node->ID == ID)
			Ret++;

	return Ret;
}

void	CheckPlasyNodes(VARRATES *Plasty)
{
	int No, Index;
	VAR_RATES_NODE	*PNode;

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



void	VarRatesAddNode(RATES *Rates, TREES *Trees, OPTIONS *Opt, TRANSFORM_TYPE Type, NODE N, long long It)
{
	VAR_RATES_NODE	*PNode;
	VARRATES	*VarRates;

	VarRates = Rates->VarRates;

	PNode = CreateVarRatesNode(It, N);

	PNode->Type = Type;

	if(PNode->Type == VR_BL || PNode->Type == VR_NODE)
		SetVRNodeBLRates(PNode, Rates->RS);
	else
		SetVRScalar(Opt, Rates, PNode);

	VarRates->NodeList = (VAR_RATES_NODE**)AddToList(&VarRates->NoNodes, (void**)VarRates->NodeList, (void*)PNode);

	Rates->LnHastings = 0;
	Rates->LnJacobion = 0;
}

void	VarRatesDelNode(RATES *Rates, TREES *Trees, OPTIONS *Opt, int No)
{
	VARRATES	*VarRates;
	VAR_RATES_NODE	**NList;
	int			Index;
	
	VarRates = Rates->VarRates;
	
	if(VarRates->NoNodes == 1)
	{
		FreeVarRatesNode(VarRates->NodeList[0]);
		free(VarRates->NodeList);
		VarRates->NodeList = NULL;
		VarRates->NoNodes = 0;
		return;
	}
	
	FreeVarRatesNode(VarRates->NodeList[No]);
	VarRates->NodeList[No] = NULL;

	NList = (VAR_RATES_NODE**)malloc(sizeof(VAR_RATES_NODE*) * (VarRates->NoNodes - 1));
	if(NList == NULL)
		MallocErr();

	No = 0;
	for(Index=0;Index<VarRates->NoNodes;Index++)
	{
		if(VarRates->NodeList[Index] != NULL)
			NList[No++] = VarRates->NodeList[Index];
	}

	free(VarRates->NodeList);
	VarRates->NodeList = NList;
	VarRates->NoNodes--;
	
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

void	TestNormHasting(void)
{
	double X, Dev;
	double LhP;
	double LhG;

	Dev = 2.0;
	X = 0.1;

	LhP = CalcNormalHasting(X, Dev);
	LhP = exp(LhP);

	LhG = gsl_cdf_gaussian_P(X, Dev);
	LhG = gsl_cdf_gaussian_Q(X, Dev);

	exit(0);
}

void	ChangeVarRatesScale(RATES *Rates, TREES *Trees, OPTIONS *Opt, SCHEDULE* Shed)
{
	VARRATES	*VarRates;
	VAR_RATES_NODE	*Node;
	int			No;
	double		Dev;
	
	Shed->CurrentAT = Shed->VarRateAT;
	Dev = Shed->CurrentAT->CDev;

	VarRates = Rates->VarRates;

	No = RandUSLong(Rates->RS) % VarRates->NoNodes;
	Node = VarRates->NodeList[No];

//	TestNormHasting();


	Rates->LnHastings = CalcNormalHasting(Node->Scale, Dev);
//	Rates->LnHastings = 0;

	if(Node->Type == VR_LAMBDA)
		Node->Scale = ChangePlastyRateLambda(Rates->RS, Node->Scale, Dev);
	else
		Node->Scale = ChangePlastyRate(Rates->RS, Node->Scale, Dev);
}

int		NodeScaled(int NID, VARRATES *VarRates)
{
	int	Index;

	for(Index=0;Index<VarRates->NoNodes;Index++)
	{
		if(VarRates->NodeList[Index]->Node->ID == NID)
			return TRUE;
	}
	
	return FALSE;
}

int		IsVarRateTypeRate(TRANSFORM_TYPE Type)
{
	if(Type == 	VR_NODE || Type == VR_BL)
		return TRUE;

	return FALSE;
}


int GetPlastyNode(int ID, VARRATES *VarRates, TRANSFORM_TYPE Type)
{
	int Index;

	for(Index=0;Index<VarRates->NoNodes;Index++)
	{

#ifdef VARRATES_ONE_OP_PER_NODE
		if(	VarRates->NodeList[Index]->Node->ID == ID)
			return Index;
#else
		if(	Plasty->NodeList[Index]->Node->ID == ID && 
			Plasty->NodeList[Index]->Type == Type)
		return Index;
#endif
	}

	return -1;
}

int		ValidMoveNode(VARRATES *VarRates, NODE N, TRANSFORM_TYPE Type)
{
	int ID;

	if(N == NULL)
		return FALSE;
	
	if(IsValidVarRatesNode(N, Type) == FALSE)
		return FALSE;
	
	ID = GetPlastyNode(N->ID, VarRates, Type);

	if(ID == -1)
		return TRUE;

	return FALSE;
}

void	MakeTNodeList(OPTIONS *Opt, VARRATES *Plasty, TRANSFORM_TYPE Type, NODE N, NODE* List, int *Size)
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

int		CompVarRatesNode(const void *Vr1, const void *Vr2)
{
	VAR_RATES_NODE **VR1, **VR2;

	VR1 = (VAR_RATES_NODE**)Vr1;
	VR2 = (VAR_RATES_NODE**)Vr2;


	if((*VR1)->Node->Part->NoTaxa >= (*VR2)->Node->Part->NoTaxa)
		return -1;

	if((*VR1)->Node->Part->NoTaxa < (*VR2)->Node->Part->NoTaxa)
		return 1;

	return 0;
}

void 	VarRatesMoveNode(RATES *Rates, TREES *Trees, OPTIONS *Opt)
{
	VARRATES	*VarRates;
	VAR_RATES_NODE	*PNode;
	NODE		N;
	int			No;
	int			Index;

	VarRates = Rates->VarRates;

	No = RandUSLong(Rates->RS) % VarRates->NoNodes;
	PNode = VarRates->NodeList[No];

	N = PNode->Node;

	VarRates->NoTempList = 0;
	if(ValidMoveNode(VarRates, N->Ans, PNode->Type) == TRUE)
		VarRates->TempList[VarRates->NoTempList++] = N->Ans;

	if(N->Tip == FALSE)
	{
		for(Index=0;Index<N->NoNodes;Index++)
			MakeTNodeList(Opt, VarRates, PNode->Type, N->NodeList[Index], VarRates->TempList, &VarRates->NoTempList);
	}

	if(VarRates->NoTempList == 0)
		return;
	
	No = RandUSLong(Rates->RS) % VarRates->NoTempList;
	PNode->Node = VarRates->TempList[No];

	qsort(VarRates->NodeList, VarRates->NoNodes, sizeof(VAR_RATES_NODE*), CompVarRatesNode);

	return;
}


NODE	GetVarRatesNode(RATES *Rates, TREES *Trees, TRANSFORM_TYPE	Type)
{
	VARRATES	*VarRates;
	TREE		*Tree;
	NODE		N;
	int			Pos;
	
	Tree = Trees->Tree[Rates->TreeNo];
	VarRates = Rates->VarRates;

	do
	{
		Pos = RandUSInt(Rates->RS) % Tree->NoNodes;
		N = Tree->NodeList[Pos];
	}while(IsValidVarRatesNode(N, Type) == FALSE); 
	
	return N;
}



void	DumpVRNodes(VARRATES	*VarRates)
{
	int Index;

	for(Index=0;Index<VarRates->NoNodes;Index++)
		printf("%d\t%d\n", Index, VarRates->NodeList[Index]->Node->Part->NoTaxa);

	exit(0);
}

void	VarRatesAddRemove(RATES *Rates, TREES *Trees, OPTIONS *Opt, SCHEDULE *Shed, long long It)
{
	VARRATES	*VarRates;
	int			PNodeID;
	NODE		N;
	TRANSFORM_TYPE		Type;
	
	Type = GetVarRatesType(Rates->RS, Shed);
		
	VarRates = Rates->VarRates;
	
	N = GetVarRatesNode(Rates, Trees, Type);

	PNodeID = GetPlastyNode(N->ID, VarRates, Type);

	if(PNodeID == -1)
		VarRatesAddNode(Rates, Trees, Opt, Type, N, It);
	else
		VarRatesDelNode(Rates, Trees, Opt, PNodeID);	
	
	qsort(VarRates->NodeList, VarRates->NoNodes, sizeof(VAR_RATES_NODE*), CompVarRatesNode);
	
//	CheckPlasyNodes(Plasty);
}


VAR_RATES_NODE *CloneVarRatesNode(VAR_RATES_NODE *N2)
{
	VAR_RATES_NODE *Ret;

	Ret = (VAR_RATES_NODE*)malloc(sizeof(VAR_RATES_NODE));
	if(Ret == NULL)
		MallocErr();

	Ret->Node = N2->Node;
	Ret->Scale= N2->Scale;
	Ret->Type = N2->Type;
	Ret->NodeID = N2->NodeID;

	return Ret;
}

void	BlankPlasty(VARRATES *P)
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
	VAR_RATES_NODE	**NList;
	VARRATES	*P1, *P2;
	int			Index;

	P1 = R1->VarRates;
	P2 = R2->VarRates;

	if(P2->NoNodes == 0)
	{
		BlankPlasty(P1);
		return;
	}

	NList = (VAR_RATES_NODE**)malloc(sizeof(VAR_RATES_NODE*) * P2->NoNodes);
	if(NList == NULL)
		MallocErr();

	for(Index=0;Index<P2->NoNodes;Index++)
		NList[Index] = CloneVarRatesNode(P2->NodeList[Index]);

	BlankPlasty(P1);

	P1->NodeList = NList;
	P1->NoNodes = P2->NoNodes;

	P1->Alpha = P2->Alpha;
}

void	ScaleNode(NODE N,  double Scale)
{
	int Index;

	N->Length = N->Length * Scale;

	if(N->Tip == TRUE)
		return;

	for(Index=0;Index<N->NoNodes;Index++)
		ScaleNode(N->NodeList[Index], Scale);
}

void	VarRatesNode(TREES *Trees, TREE *Tree, NODE N, double Scale, TRANSFORM_TYPE Type)
{
	int Norm;

	Norm = FALSE;

	if(Type == VR_BL)
		N->Length = N->Length * Scale;

	if(Type == VR_NODE)
		ScaleNode(N,  Scale);

	if(Type == VR_KAPPA)
		TransformTreeKappa(N, Scale, Norm);

	if(Type == VR_LAMBDA)
	{
		SetTreeDistToRoot(Tree);
		TransformTreeLambda(N, Scale, Norm);
	}

	if(Type == VR_DELTA)
		TransformTreeDelta(N, Scale, Norm);

	if(Type == VR_OU)
	{
		SetTreeDistToRoot(Tree);
		TransformTreeOU(Trees, N, Scale, Norm);
	}
}

void	CheckVarRatesData(OPTIONS *Opt, TREES *Trees, RATES *Rates)	
{
	VARRATES *VarRates;
	NODE N;
	int Index;

	VarRates = Rates->VarRates;

	for(Index=0;Index<VarRates->NoNodes;Index++)
	{
		N = VarRates->NodeList[Index]->Node;

		if(N->Part->NoTaxa < MIN_TAXA_VR_TRANS)
		{
			printf("err.\n");
			exit(1);
		}

		if(IsValidVarRatesNode(VarRates->NodeList[Index]->Node, VarRates->NodeList[Index]->Type) == FALSE)
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
	VARRATES *VarRates;
	double SumBL, Scale;

//	CheckVarRatesData(Opt, Trees, Rates);

	VarRates = Rates->VarRates;
	TNo = Rates->TreeNo;
	Tree = Trees->Tree[TNo];

	if(Normalise == TRUE)
		SumBL = SumNodeBL(Tree->Root);

	for(Index=0;Index<VarRates->NoNodes;Index++)
		VarRatesNode(Trees, Tree, VarRates->NodeList[Index]->Node, VarRates->NodeList[Index]->Scale, VarRates->NodeList[Index]->Type);

	if(Normalise == FALSE)
		return;

	Scale = SumBL / SumNodeBL(Tree->Root);
	ScaleSubTree(Tree->Root, Scale);
}

void	InitVarRatesTreeFile(OPTIONS *Opt, TREES *Trees)
{
	char	*Buffer;
	int		Index;
	
	Buffer = (char*)SMalloc(sizeof(char) * (strlen(Opt->LogFN) + BUFFERSIZE));
			
	sprintf(Buffer, "%s.VarRates.trees", Opt->LogFN);
	Opt->PPTree = OpenWrite(Buffer);
	free(Buffer);

	fprintf(Opt->PPTree, "#NEXUS\n");
	fprintf(Opt->PPTree, "\tBegin Trees;\n");
	fprintf(Opt->PPTree, "\t\tTranslate\n");

	for(Index=0;Index<Trees->NoOfTaxa-1;Index++)
		fprintf(Opt->PPTree, "\t\t%d\t%s,\n", Trees->Taxa[Index]->No+1, Trees->Taxa[Index]->Name);
	fprintf(Opt->PPTree, "\t\t%d\t%s\n\t\t;\n", Trees->Taxa[Index]->No+1, Trees->Taxa[Index]->Name);

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
	VARRATES	*VarRates;
	int		No;

	T = Trees->Tree[0];
	VarRates = Rates->VarRates;

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
		fprintf(Out, "%d:%f", N->Taxa->No+1, N->Length);
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
	VARRATES	*VarRates;
	TREE	*Tree;
	int		Index;

	VarRates = Rates->VarRates;
	Tree = Trees->Tree[Rates->TreeNo];

	ReSetBranchLength(Tree);
	VarRatesTree(Opt, Trees, Rates, NORMALISE_TREE_CON_SCALING);

	fprintf(Opt->PPTree, "\tTree VarRates_%llu_%d = (", It, VarRates->NoNodes);

	for(Index=0;Index<Tree->Root->NoNodes-1;Index++)
	{
		PrintPPNode(Opt->PPTree, Tree->Root->NodeList[Index]);
		fprintf(Opt->PPTree, ",");
	}

	PrintPPNode(Opt->PPTree, Tree->Root->NodeList[Index]);
	fprintf(Opt->PPTree, ");\n");

	fflush(Opt->PPTree);
}

void	OutputVarRatesType(FILE *Out, TRANSFORM_TYPE Type)
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
	VARRATES		*VarRates;
	int			Index;
	VAR_RATES_NODE	*PNode;
	NODE		N;
	double		Sigma;

	VarRates = Rates->VarRates;
//	CheckPlasyNodes(P);

	Out = Opt->PPLog;

	fprintf(Out, "%lld\t%f\t%f\t%d\t", It, Rates->Lh, Rates->Lh + Rates->LhPrior, VarRates->NoNodes);

	if(Opt->Model == M_CONTRAST_CORREL)
	{
		Sigma = Rates->Contrast->SigmaMat->me[0][0];
		fprintf(Out, "%f\t%f\t%f\t", Rates->Contrast->Alpha[0], Sigma, VarRates->Alpha);
	}

	if(Opt->Model == M_CONTRAST) 
	{
		Sigma = Rates->Contrast->Sigma[0];
		fprintf(Out, "%f\t%f\t%f\t", Rates->Contrast->Alpha[0], Sigma, VarRates->Alpha);
	}

	if(Opt->Model == M_CONTRAST_REG)
	{
		Sigma = Rates->Contrast->GlobalVar;
		fprintf(Out, "%f\t%f\t%f\t", Rates->Contrast->RegAlpha, Sigma, VarRates->Alpha);
	}

	if(Opt->ModelType == MT_FATTAIL || Opt->ModelType == MT_DISCRETE)
		fprintf(Out, "NA\tNA\tNA\t");
		
//	fprintf(Out, "\n");
	for(Index=0;Index<VarRates->NoNodes;Index++)
	{
		PNode = VarRates->NodeList[Index];
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
	VARRATES	*VarRates;

	VarRates = Rates->VarRates;

	Rates->LnHastings = CalcNormalHasting(VarRates->Alpha - 1, VARRATES_HP_ALPHA_SCLAE);
	
	do
	{
		NAlpha = RandNormal(Rates->RS, VarRates->Alpha, VARRATES_HP_ALPHA_SCLAE);
	} while(NAlpha <= 1);

	VarRates->Alpha = NAlpha;
}

PRIOR*	GetVRPrior(TRANSFORM_TYPE Type, RATES *Rates)
{
	if(Type == VR_BL)
		return GetPriorFromName("VRBL", Rates->Priors, Rates->NoPriors);

	if(Type == VR_NODE)
		return GetPriorFromName("VRNode", Rates->Priors, Rates->NoPriors);

	if(Type == VR_KAPPA)
		return GetPriorFromName("Kappa", Rates->Priors, Rates->NoPriors);

	if(Type == VR_DELTA)
		return GetPriorFromName("Delta", Rates->Priors, Rates->NoPriors);

	if(Type == VR_LAMBDA)
		return GetPriorFromName("Lambda", Rates->Priors, Rates->NoPriors);

	if(Type == VR_OU)
		return GetPriorFromName("OU", Rates->Priors, Rates->NoPriors);

	return NULL;
}

double	CaclVRPrior(VAR_RATES_NODE *PNode, RATES *Rates)
{
	double Ret;
	PRIOR *Prior;

	if(PNode->Type == VR_KAPPA && PNode->Scale >= MAX_VR_KAPPA)
		return ERRLH;

	if(PNode->Type == VR_DELTA && PNode->Scale >= MAX_VR_DELTA)
		return ERRLH;

	if(PNode->Type == VR_LAMBDA && PNode->Scale >= MAX_VR_LAMBDA)
		return ERRLH;

	if(PNode->Type == VR_OU && PNode->Scale >= MAX_VR_OU)
		return ERRLH;
	
	Prior = GetVRPrior(PNode->Type, Rates);
	
	Ret = CalcLhPriorP(PNode->Scale, Prior);
	
	return Ret;
}

double	CalcVarRatesPriors(RATES *Rates, OPTIONS *Opt)
{
	double		Ret;
	int			Index;
	VARRATES	*VarRates;
	VAR_RATES_NODE	*PNode;
	double		PVal;

	VarRates = Rates->VarRates;
	Ret = 0;

	for(Index=0;Index<VarRates->NoNodes;Index++)
	{
		PNode = VarRates->NodeList[Index];

		PVal = CaclVRPrior(PNode, Rates);

		if(PVal == ERRLH)
			return ERRLH;

		Ret += PVal;
	}

	return Ret;
}


int		SeenNodeID(int NID, VARRATES *Plasty, int Size)
{
	int	Index;

	for(Index=0;Index<Size;Index++)
		if(Plasty->NodeList[Index]->Node->ID == NID)
			return TRUE;

	return FALSE;
}

VAR_RATES_NODE*	CreateTextPNode(NODE N, double Scale, long long CIt, TRANSFORM_TYPE Type)
{
	VAR_RATES_NODE*	Ret;
	
	Ret = CreateVarRatesNode(CIt, N);

	Ret->Scale = Scale;
	Ret->Type = Type;

	return Ret;
}

TRANSFORM_TYPE	StrToRJVarRatesType(char *Str)
{
	MakeLower(Str);

	if(strcmp("node", Str) == 0)
		return VR_NODE;
	
	if(strcmp("branch", Str) == 0)
		return VR_BL;
	
	if(strcmp("kappa", Str) == 0)
		return VR_KAPPA;

	if(strcmp("lambda", Str) == 0)
		return VR_LAMBDA;

	if(strcmp("delta", Str) == 0)
		return VR_DELTA;

	if(strcmp("ou", Str) == 0)
		return VR_OU;

	printf("Unkown string (%s) in %s::%d.\n", Str, __FILE__, __LINE__);
	exit(0);
}

void	AddTextVarRate(TREE *Tree, RATES *Rates, int Tokes, char **Passed)
{
	VAR_RATES_NODE* PNode;
	NODE N;
	int NodeID;
	double Scale;
	long long Itter;
	TRANSFORM_TYPE Type;
	VARRATES		*Plasty;

	Plasty = Rates->VarRates;

	NodeID = atoi(Passed[0]);
	Scale = atof(Passed[1]);
	sscanf(Passed[2], "%lld", &Itter);
	
	Type = StrToRJVarRatesType(Passed[3]);

	N = Tree->NodeList[NodeID];
		
	PNode = CreateTextPNode(N, Scale, Itter, Type);

	Plasty->NodeList = (VAR_RATES_NODE**) AddToList(&Plasty->NoNodes, (void**)Plasty->NodeList, (void*)PNode);
}

VAR_RATES_NODE**	MakeNewList(VARRATES *Plasty, int Pos)
{
	VAR_RATES_NODE** Ret;
	int CPos, Index;

	Ret = (VAR_RATES_NODE**)malloc(sizeof(VAR_RATES_NODE*) * (Plasty->NoNodes - 1));


	CPos = 0;
	for(Index=0;Index<Plasty->NoNodes;Index++)
	{
		if(Index != Pos)
			Ret[CPos++] = Plasty->NodeList[Index];
	}
	
	return Ret;
}

void	SetVarRatesList(VARRATES *Plasty, VAR_RATES_NODE** VRateList, int NoVRates)
{
	free(Plasty->NodeList);

	Plasty->NodeList = VRateList;
	Plasty->NoNodes = NoVRates;
}


double	LhVarRatesList(OPTIONS *Opt, RATES *Rates, TREES *Trees, VAR_RATES_NODE** VRateList, int NoVRates)
{
	VAR_RATES_NODE** OList;
	int NoOld;
	VARRATES		*VarRates;
	double		Ret;

	VarRates = Rates->VarRates;
	
	NoOld = VarRates->NoNodes;
	OList = VarRates->NodeList;

	VarRates->NodeList = VRateList;
	VarRates->NoNodes = NoVRates;
	
	Ret = Likelihood(Rates, Trees, Opt);

	VarRates->NodeList = OList;
	VarRates->NoNodes = NoOld;


	return Ret;
}

void	RemoveEachVarRate(OPTIONS *Opt, RATES *Rates)
{
	VARRATES		*VarRates;
	VAR_RATES_NODE	**RList;
	int Index;
	double Lh;

	VarRates = Rates->VarRates;
	
	printf("\n\n\n");

	printf("Removing nodes one at a time\n");

	for(Index=0;Index<VarRates->NoNodes;Index++)
	{

		RList = MakeNewList(VarRates, Index);
		Lh = LhVarRatesList(Opt, Rates, Opt->Trees, RList, VarRates->NoNodes - 1);

		printf("%d\tLh\t%f\n", Index, Lh);
				
		free(RList);
	}

//	exit(0);
}

VAR_RATES_NODE** CreateVRateList1(VARRATES *Plasty, int No1)
{
	VAR_RATES_NODE** Ret;

	Ret = (VAR_RATES_NODE**)malloc(sizeof(VAR_RATES_NODE*) * 1);
	
	Ret[0] = Plasty->NodeList[No1];

	return Ret;
}

VAR_RATES_NODE** CreateVRateList2(VARRATES *Plasty, int No1, int No2)
{
	VAR_RATES_NODE** Ret;

	Ret = (VAR_RATES_NODE**)malloc(sizeof(VAR_RATES_NODE*) * 2);
	
	Ret[0] = Plasty->NodeList[No1];
	Ret[1] = Plasty->NodeList[No2];


	return Ret;
}

void	OneRateAtOnce(OPTIONS *Opt, RATES *Rates)
{
	int Index;
	VARRATES		*VarRates;
	VAR_RATES_NODE**	List;
	double Lh;

	VarRates = Rates->VarRates;

	printf("Testing each node, one at a time\n");

//	for(X=0;X<Plasty->NoNodes;X++)
	for(Index=0;Index<VarRates->NoNodes;Index++)
	{
		List = CreateVRateList1(VarRates, Index);
		Lh = LhVarRatesList(Opt, Rates, Opt->Trees, List, 1);

		printf("Using:\t%d\t%f\n", Index, Lh);

		free(List);
	}
}

void	SpecifcPairTest(OPTIONS *Opt, RATES *Rates)
{
	double		Lh;
	VARRATES		*VarRates;
	VAR_RATES_NODE**	List;

	VarRates = Rates->VarRates;

	List = CreateVRateList2(VarRates, 1, 6);

	SetVarRatesList(VarRates, List, 2);

	Lh = Likelihood(Rates, Opt->Trees, Opt);
	printf("Lh:\t%f\n", Lh);

	CalcPriors(Rates, Opt);

	PrintPPTree(Opt, Opt->Trees, Rates, 1);
		
	exit(0);
}

void	TestEachVarRatePair(OPTIONS *Opt, RATES *Rates)
{
	int X, Y;
	VARRATES		*VarRates;
	VAR_RATES_NODE**	List;
	double Lh;

	VarRates = Rates->VarRates;

	SpecifcPairTest(Opt, Rates);
/*
	List = CreateVRateList(Plasty, 1, 3);


	SetVarRatesList(Plasty, List, 2);

	Lh = Likelihood(Rates, Opt->Trees, Opt);
	printf("Lh:\t%f\n", Lh);

	CalcPriors(Rates, Opt);

	PrintPPTree(Opt, Opt->Trees, Rates, 1);
		
	exit(0);
*/	

	printf("\n\n\n");
	printf("Testing each pair of node\n");

	for(X=0;X<VarRates->NoNodes;X++)
	{
		for(Y=0;Y<VarRates->NoNodes;Y++)
		{
			if(X != Y)
			{
				List = CreateVRateList2(VarRates, X, Y);
				Lh = LhVarRatesList(Opt, Rates, Opt->Trees, List, 2);

				printf("%d\t%d\t%f\n", X, Y, Lh);

				free(List);
			}
		}
	}

	exit(0);
}

void	TestVarRates(OPTIONS *Opt, RATES *Rates)
{
	TREES *Trees;
	TREE *Tree;
	double Lh;

	Trees = Opt->Trees;
	Tree = Trees->Tree[0];

	ReSetBranchLength(Tree);

	Lh = Likelihood(Rates, Opt->Trees, Opt);

	printf("Lh = %f\n", Lh);
//	exit(0);

	OneRateAtOnce(Opt, Rates);
	
	RemoveEachVarRate(Opt, Rates);



	TestEachVarRatePair(Opt, Rates);
	
//	PrintPPTree(Opt, Opt->Trees, Rates, 1);

	exit(0);
}

void	DumpVarRates(RATES *Rates)
{
	int Index;
	NODE N;
	VAR_RATES_NODE	*PNode;
	VARRATES		*VarRates;
	
	VarRates = Rates->VarRates;

	for(Index=0;Index<VarRates->NoNodes;Index++)
	{
		PNode = VarRates->NodeList[Index];
		N = PNode->Node;
		
		printf("%d\t", Index);
		printf("%d\t", N->ID);
		printf("%f\t", PNode->Scale);
		printf("%llu\t", PNode->NodeID);

		OutputVarRatesType(stdout, PNode->Type);

		printf("\n");
	}

	printf("\n\n\n\n\n\n");

	fflush(stdout);
//	exit(0);
}


void	SetVarRatesFromStr(char *Str, RATES *Rates, OPTIONS *Opt)
{
	char **Passed;
	int Index, Tokes;
	char *S;

	PrintPPTree(Opt, Opt->Trees, Rates, 1);
	

	S = StrMake(Str);

	Passed = (char**)malloc(sizeof(char*) * strlen(Str));
	if(Passed == NULL)
		MallocErr();

	Tokes = MakeArgv(S, Passed, (int)strlen(S));

	Rates->Contrast->Alpha[0] = atof(Passed[4]);
	Rates->Contrast->Sigma[0] = atof(Passed[5]);

	Index = 7;
	for(;Index<Tokes;Index+=4)
		AddTextVarRate(Opt->Trees->Tree[0], Rates, 4, &Passed[Index]);
//		printf("%s\n", Passed[Index]);
	
	DumpVarRates(Rates);

	TestVarRates(Opt, Rates);
	
	free(Passed);
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


