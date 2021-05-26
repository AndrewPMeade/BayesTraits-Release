/*
*  BayesTriats 3.0
*
*  copyright 2017
*
*  Andrew Meade
*  School of Biological Sciences
*  University of Reading
*  Reading
*  Berkshire
*  RG6 6BX
*
* BayesTriats is free software: you can redistribute it and/or modify
* it under the terms of the GNU General Public License as published by
* the Free Software Foundation, either version 3 of the License, or
* (at your option) any later version.
* 
* This program is distributed in the hope that it will be useful,
* but WITHOUT ANY WARRANTY; without even the implied warranty of
* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
* GNU General Public License for more details.
* 
* You should have received a copy of the GNU General Public License
* along with this program.  If not, see <http://www.gnu.org/licenses/>
*
*/

#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "TypeDef.h"
#include "GenLib.h"
#include "VarRates.h"
#include "Matrix.h"
#include "RandLib.h"
#include "Likelihood.h"
#include "Trees.h"
#include "RandDists.h"
#include "RJLocalScalar.h"
#include "Priors.h"
#include "TransformTree.h"
#include "Part.h"
#include "StableDist.h"
#include "StochasticBeta.h"

#include <gsl/gsl_cdf.h>

void	OutputVarRatesType(FILE *Out, TRANSFORM_TYPE Type);

int		UseRJLocalScalar(OPTIONS* Opt)
{
	int Index;

	for (Index = 0; Index < NO_RJ_LOCAL_SCALAR; Index++)
		if (Opt->UseRJLocalScalar[Index] == TRUE)
			return TRUE;

	return FALSE;
}

int		UseNonParametricMethods(OPTIONS *Opt)
{
	if(UseRJLocalScalar(Opt) == TRUE)
		return TRUE;

	if(Opt->UseStochasticBeta == TRUE)
		return TRUE;

	return FALSE;
}

TRANSFORM_TYPE	StrToVarRatesType(char *Str)
{
	if(StrICmp("Node", Str) == 0)
		return VR_NODE;

	if(StrICmp("Branch", Str) == 0)
		return VR_BL;

	if(StrICmp("Kappa", Str) == 0)
		return VR_KAPPA;

	if(StrICmp("Lambda", Str) == 0)
		return VR_LAMBDA;

	if(StrICmp("Delta", Str) == 0)
		return VR_DELTA;

	if(StrICmp("OU", Str) == 0)
		return VR_OU;

	if(StrICmp("LandscapeBL", Str) == 0)
		return VR_LS_BL;

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

	if(Type == VR_LS_BL)
		return "LandscapeBL";

	printf("%s::%d unkonwn RJ Variable type\n", __FILE__, __LINE__);
	exit(0);
	return NULL;
}

NODE	GetVRNode(TREES *Trees, int TreeNo, VAR_RATES_NODE *VR_Node)
{
	if(VR_Node->NodeList[TreeNo] == NULL)
		VR_Node->NodeList[TreeNo] = PartGetMRCA(Trees->Tree[TreeNo], VR_Node->Part);

	return VR_Node->NodeList[TreeNo];
}

double	RandGamma(double Shape, double Scale)
{
	double x;

	x = sgamma(Shape) * Scale;
	return x / (Scale* (Shape  - 1));
}


int		IsValidVarRatesNode(NODE N, TRANSFORM_TYPE	Type, OPTIONS *Opt)
{
	PART *Part;
	
	if(N == NULL)
		return FALSE;
	
	if(N->Length == 0)
		return FALSE;

	Part = N->Part;

	/* Don't scale the root */
	if(Type == VR_BL || Type == VR_NODE || Type == VR_LS_BL)
	{
		if(N->Ans == NULL)
			return FALSE;
	}

//	if(Type == VR_LS_BL && N->Tip == TRUE)
//		return FALSE;

	if(Type == VR_LS_BL)
		return TRUE;

	if(Type == VR_BL || Type == VR_NODE)
	{
		if(N->Tip == TRUE && Type == VR_NODE)
			return FALSE;
		
		if(Type == VR_NODE && Part->NoTaxa < MIN_TAXA_VR_NODE)
			return FALSE;

		return TRUE;
	}

	if(Part->NoTaxa >= Opt->MinTransTaxaNo)
		return TRUE;

	return FALSE;
}

VARRATES*	CreatVarRates(RATES *Rates, TREES *Trees, OPTIONS *Opt)
{
	VARRATES* Ret;
	
	Ret = (VARRATES*)SMalloc(sizeof(VARRATES));

	Ret->NoNodes = 0;
	Ret->NodeList= NULL;

	Ret->TempList = (NODE*)SMalloc(sizeof(NODE) * Trees->MaxNodes);

#ifdef PPUNIFORM
	Ret->Alpha = -1;
#else
	Ret->Alpha = VARRATES_ALPHA;
#endif

	return Ret;
}

void		FreeVarRatesNode(VAR_RATES_NODE* VarRatesNode)
{
	free(VarRatesNode->NodeList);
	free(VarRatesNode);
}

void	FreeVarRates(VARRATES* Plasty)
{
	int Index;

	if(Plasty->NodeList != NULL)
	{
		for(Index=0;Index<Plasty->NoNodes;Index++)
			FreeVarRatesNode(Plasty->NodeList[Index]);
		
		free(Plasty->NodeList);
	}

	free(Plasty->TempList);
	free(Plasty);
}


void			BlankNodeList(NODE *NList, int NoTrees)
{
	int TIndex;

	for(TIndex=0;TIndex<NoTrees;TIndex++)
		NList[TIndex] = NULL;
}

VAR_RATES_NODE*	AllocVarRatesNode(int NoTrees)
{
	VAR_RATES_NODE	*Ret;

	Ret = (VAR_RATES_NODE*)SMalloc(sizeof(VAR_RATES_NODE));
	Ret->NodeList = (NODE*)SMalloc(sizeof(NODE) * NoTrees);

	
	BlankNodeList(Ret->NodeList, NoTrees);

	return Ret;
}


VAR_RATES_NODE*	CreateVarRatesNode(long long It, PART *Part, int NoTrees)
{
	VAR_RATES_NODE	*Ret;
	
	Ret = AllocVarRatesNode(NoTrees);

	Ret->NodeID = It;
		
	Ret->Part = Part;

	Ret->Scale = -1;
	Ret->Type = VR_BL;

	return Ret;
}

TRANSFORM_TYPE	GetVarRatesType(RANDSTATES *RS, SCHEDULE *Shed)
{
	int Pos;

	Pos = RandInProportion(RS, Shed->FreqVarRatesOp, Shed->NoVarRatesOp);

//	if(Shed->VarRatesOp[Pos] == VR_BL)
//		return VR_NODE;

//	if(Shed->VarRatesOp[Pos] == VR_NODE)
//		return VR_BL;

	
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

	if(PNode->Type == VR_LS_BL && RandDouble(Rates->RS) < 0.5)
		PNode->Scale = PNode->Scale * -1;
}

int		CountPlasyID(VAR_RATES_NODE *VRNode, VARRATES *VarRates)
{
	int Ret, Index;

	Ret = 0;
	for(Index=0;Index<VarRates->NoNodes;Index++)
#ifdef VARRATES_ONE_OP_PER_NODE
		if(VRNode->Part->PartID == VarRates->NodeList[Index]->Part->PartID )
#else
		if( VRNode->Part->PartID == VarRates->NodeList[Index]->Part->PartID && 
			VRNode->Type == VarRates->NodeList[Index]->Type)
#endif
			Ret++;

	return Ret;
}

void	CheckPlasyNodes(VARRATES *VarRates)
{
	int No, Index;
	VAR_RATES_NODE	*PNode;

	for(Index=0;Index<VarRates->NoNodes;Index++)
	{
		PNode = VarRates->NodeList[Index];
		No = CountPlasyID(PNode, VarRates);
		if(No != 1)
		{
			printf("Multiple hits of the same node %d.\n", No);
			exit(1);
		}

		if(PNode->Type == VR_NODE && PNode->Part->NoTaxa < MIN_TAXA_VR_NODE)
		{
			printf("Node no taxa < minimum.\n");
			exit(1);
		}

		/*if(PNode->Type == VR_LS_BL && PNode->Part->NoTaxa == 1)
		{
			printf("VR LS BL placed on a Tip\n");
			exit(1);
		}*/
	}

	for(Index=0;Index<VarRates->NoNodes;Index++)
	{
		PNode = VarRates->NodeList[Index];
		if(PNode->Type == VR_LAMBDA && PNode->Scale > 1.0)
			printf("Err\n");
	}
}

void	SetScalar(OPTIONS *Opt, RATES *Rates,  VAR_RATES_NODE *PNode)
{
	if(PNode->Type == VR_BL || PNode->Type == VR_NODE)
	{
		SetVRNodeBLRates(PNode, Rates->RS);
		return;
	}

	SetVRScalar(Opt, Rates, PNode);
}

double	SetLandscapeBeta(RATES *Rates, double t)
{
	double Z, Sig2, Beta;

	Z = RandNormal(Rates->RS, 0.0, 1.0);
	Sig2 = Rates->Rates[1];

	Beta = (Z * sqrt(Sig2 * t)) / t;

	return Beta;
}

void	VarRatesAddNode(RATES *Rates, TREES *Trees, OPTIONS *Opt, TRANSFORM_TYPE Type, NODE N, long long It)
{
	VAR_RATES_NODE	*PNode;
	VARRATES	*VarRates;

	VarRates = Rates->VarRates;

	PNode = CreateVarRatesNode(It, N->Part, Trees->NoTrees);

	PNode->Type = Type;

	SetScalar(Opt, Rates, PNode);
	
	VarRates->NodeList = (VAR_RATES_NODE**)AddToList(&VarRates->NoNodes, (void**)VarRates->NodeList, (void*)PNode);

	Rates->LnHastings = 0;
	Rates->LnJacobion = 0;
}

int GetVRPosFromNode(VARRATES* VarRates, VAR_RATES_NODE *VRNode)
{
	int Index;

	for(Index=0;Index<VarRates->NoNodes;Index++)
	{
		if(VarRates->NodeList[Index] == VRNode)
			return Index;
	}

	return -1;
}

void	VarRatesDelNode(RATES *Rates, TREES *Trees, OPTIONS *Opt,  VAR_RATES_NODE *VRNode)
{
	VARRATES	*VarRates;
	VAR_RATES_NODE	**NList;
	int			Index;
	int			No;
	
	VarRates = Rates->VarRates;
	
	if(VarRates->NoNodes == 1)
	{
		FreeVarRatesNode(VarRates->NodeList[0]);
		free(VarRates->NodeList);
		VarRates->NodeList = NULL;
		VarRates->NoNodes = 0;
		return;
	}
	
	No = GetVRPosFromNode(VarRates, VRNode);
	FreeVarRatesNode(VRNode);
	VarRates->NodeList[No] = NULL;

	NList = (VAR_RATES_NODE**)SMalloc(sizeof(VAR_RATES_NODE*) * (VarRates->NoNodes - 1));

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
	{
		Node->Scale = ChangePlastyRateLambda(Rates->RS, Node->Scale, Dev);
		return;
	}

	if(Node->Type == VR_LS_BL )
	{
		Node->Scale = RandNormal(Rates->RS, Node->Scale, Dev);
		return;
	}
	
	Node->Scale = ChangePlastyRate(Rates->RS, Node->Scale, Dev);
}

int		IsVarRateTypeRate(TRANSFORM_TYPE Type)
{
	if(Type == 	VR_NODE || Type == VR_BL)
		return TRUE;

	return FALSE;
}

/*
int GetPlastyNode(int ID, VARRATES *VarRates, TRANSFORM_TYPE Type)
{
	int Index;

	for(Index=0;Index<VarRates->NoNodes;Index++)
	{

#ifdef VARRATES_ONE_OP_PER_NODE
		if(	VarRates->NodeList[Index]->Node->ID == ID)
			return Index;
#else
		if(	VarRates->NodeList[Index]->Node->ID == ID && 
			VarRates->NodeList[Index]->Type == Type)
		return Index;
#endif
	}

	return -1;
}*/

VAR_RATES_NODE* NodeHasVRSclar(TREES *Trees, NODE N, VARRATES *VarRates, TRANSFORM_TYPE Type, int TreeNo)
{
	NODE VNode;
	int Index;

	for(Index=0;Index<VarRates->NoNodes;Index++)
	{
		VNode = GetVRNode(Trees, TreeNo, VarRates->NodeList[Index]);
		if(VNode == N && Type == VarRates->NodeList[Index]->Type)
			return VarRates->NodeList[Index];
	}

	return NULL;
}


int		ValidMoveNode(VARRATES *VarRates, NODE N, TRANSFORM_TYPE Type, OPTIONS *Opt, int TreeNo)
{
	if(N == NULL)
		return FALSE;
	
	if(IsValidVarRatesNode(N, Type, Opt) == FALSE)
		return FALSE;
	
	if(NodeHasVRSclar(Opt->Trees, N, VarRates, Type, TreeNo) == NULL)
		return TRUE;

	return FALSE;
}

void	MakeTNodeList(OPTIONS *Opt, int TreeNo, VARRATES *Plasty, TRANSFORM_TYPE Type, NODE N, NODE* List, int *Size)
{
	if(ValidMoveNode(Plasty, N, Type, Opt, TreeNo) == TRUE)
	{
		List[*Size] = N;
		(*Size)++;
	}

	if(N->Tip == TRUE)
		return;

}

int		CompVarRatesNodeSize(int S1, int S2)
{
	return S2 - S1;
}

int		CompVarRatesNodeType(VAR_RATES_NODE *VR1, VAR_RATES_NODE *VR2)
{
	int Size1, Size2;
	
	Size1 = VR1->Part->NoTaxa;
	Size2 = VR2->Part->NoTaxa;

	if((VR1->Type == VR_BL || VR1->Type == VR_NODE) && (VR2->Type == VR_BL || VR2->Type == VR_NODE))
		return CompVarRatesNodeSize(Size1, Size2);

	if(VR1->Type == VR_BL || VR1->Type == VR_NODE)
		return -1;

	if(VR2->Type == VR_BL || VR2->Type == VR_NODE)
		return 1;
	
	return CompVarRatesNodeSize(Size1, Size2);
}

int		CompVarRatesNode(const void *Vr1, const void *Vr2)
{
	VAR_RATES_NODE **VR1, **VR2;

	VR1 = (VAR_RATES_NODE**)Vr1;
	VR2 = (VAR_RATES_NODE**)Vr2;

	if((*VR1)->Type != (*VR2)->Type)
		return CompVarRatesNodeType((*VR1), (*VR2));

	return CompVarRatesNodeSize((*VR1)->Part->NoTaxa, (*VR2)->Part->NoTaxa);
}

void 	VarRatesMoveNode(RATES *Rates, TREES *Trees, OPTIONS *Opt)
{
	VARRATES	*VarRates;
	VAR_RATES_NODE	*PNode;
	NODE		N;
	int			No;
	int			Index;
	int			TreeNo;

	VarRates = Rates->VarRates;
	TreeNo = Rates->TreeNo;

	No = RandUSLong(Rates->RS) % VarRates->NoNodes;
	PNode = VarRates->NodeList[No];


	N = PNode->NodeList[TreeNo];

	VarRates->NoTempList = 0;
	if(ValidMoveNode(VarRates, N->Ans, PNode->Type, Opt, TreeNo) == TRUE)
		VarRates->TempList[VarRates->NoTempList++] = N->Ans;

	if(N->Tip == FALSE)
	{
		for(Index=0;Index<N->NoNodes;Index++)
			MakeTNodeList(Opt, TreeNo, VarRates, PNode->Type, N->NodeList[Index], VarRates->TempList, &VarRates->NoTempList);
	}
	
	if(VarRates->NoTempList == 0)
		return;
	
	No = RandUSLong(Rates->RS) % VarRates->NoTempList;

	PNode->Part = VarRates->TempList[No]->Part;
	BlankNodeList(PNode->NodeList, Trees->NoTrees);

	PNode->NodeList[TreeNo] = VarRates->TempList[No];
	
	qsort(VarRates->NodeList, VarRates->NoNodes, sizeof(VAR_RATES_NODE*), CompVarRatesNode);
	
	if(N->Tip == TRUE)
		Rates->LnHastings = log(1.0 / 3.0);

	if(N->Ans != NULL)
		if(N->Ans->Ans == NULL)
			Rates->LnHastings = log(2.0 / 3.0);
	
	return;
}

NODE	GetVarRatesNode(OPTIONS *Opt, RATES *Rates, TREES *Trees, TRANSFORM_TYPE	Type)
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
	}while(IsValidVarRatesNode(N, Type, Opt) == FALSE); 
	
	return N;
}

void	VarRatesAddRemove(RATES *Rates, TREES *Trees, OPTIONS *Opt, SCHEDULE *Shed, long long It)
{
	VARRATES	*VarRates;
	VAR_RATES_NODE	*VRNode;
	NODE		N;
	TRANSFORM_TYPE		Type;
	
	Type = GetVarRatesType(Rates->RS, Shed);

	VarRates = Rates->VarRates;
	
	N = GetVarRatesNode(Opt, Rates, Trees, Type);
	
	VRNode = NodeHasVRSclar(Trees, N, VarRates, Type, Rates->TreeNo);

	if(VRNode == NULL)
		VarRatesAddNode(Rates, Trees, Opt, Type, N, It);
	else
		VarRatesDelNode(Rates, Trees, Opt, VRNode);	

	qsort(VarRates->NodeList, VarRates->NoNodes, sizeof(VAR_RATES_NODE*), CompVarRatesNode);
	
	CheckPlasyNodes(VarRates);

}

VAR_RATES_NODE *CloneVarRatesNode(VAR_RATES_NODE *VR_Node, int NoTrees)
{
	VAR_RATES_NODE *Ret;

	Ret = AllocVarRatesNode(NoTrees);

	memcpy(Ret->NodeList, VR_Node->NodeList, sizeof(NODE) * NoTrees);
	

	Ret->Scale= VR_Node->Scale;
	Ret->Type = VR_Node->Type;
	Ret->NodeID = VR_Node->NodeID;
	Ret->Part = VR_Node->Part;

	return Ret;
}

void	EmptyPlasty(VARRATES *P)
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

void	VarRatesCopy(TREES *Trees, RATES *R1, RATES *R2)
{
	VAR_RATES_NODE	**NList;
	VARRATES	*VR1, *VR2;
	int			Index;

	VR1 = R1->VarRates;
	VR2 = R2->VarRates;

	if(VR2->NoNodes == 0)
	{
		EmptyPlasty(VR1);
		return;
	}

	NList = (VAR_RATES_NODE**)SMalloc(sizeof(VAR_RATES_NODE*) * VR2->NoNodes);

	for(Index=0;Index<VR2->NoNodes;Index++)
		NList[Index] = CloneVarRatesNode(VR2->NodeList[Index], Trees->NoTrees);

	EmptyPlasty(VR1);

	VR1->NodeList = NList;
	VR1->NoNodes = VR2->NoNodes;

	VR1->Alpha = VR2->Alpha;
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

	// this will have been applyed. 
	if(Type == VR_LS_BL)
		return;

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
	PART *Part;
	int Index;

	VarRates = Rates->VarRates;

	for(Index=0;Index<VarRates->NoNodes;Index++)
	{
		Part = VarRates->NodeList[Index]->Part;

		if(Part->NoTaxa < Opt->MinTransTaxaNo)
		{
			printf("err.\n");
			exit(1);
		}

		if(IsValidVarRatesNode(VarRates->NodeList[Index]->NodeList[Rates->TreeNo], VarRates->NodeList[Index]->Type, Opt) == FALSE)
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
	VAR_RATES_NODE *VR_Node;
	NODE Node;
	
	VarRates = Rates->VarRates;
	TNo = Rates->TreeNo;
	Tree = Trees->Tree[TNo];

	if(Normalise == TRUE)
		SumBL = SumNodeBL(Tree->Root);

	for(Index=0;Index<VarRates->NoNodes;Index++)
	{
		VR_Node = VarRates->NodeList[Index];

		Node = GetVRNode(Trees, TNo, VR_Node);
	
		VarRatesNode(Trees, Tree, Node, VR_Node->Scale, VR_Node->Type);
	}

	if(Normalise == FALSE)
		return;

	Scale = SumBL / SumNodeBL(Tree->Root);
	ScaleSubTree(Tree->Root, Scale);
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


void	VarRatesLogFileHeaderSingleTree(OPTIONS *Opt, TREES *Trees, RATES *Rates)
{
	TREE	*T;
	int		Index;
	NODE	N;
	VARRATES	*VarRates;
	int		No;

	T = Trees->Tree[0];
	VarRates = Rates->VarRates;

	fprintf(Opt->VarRatesLog, "%d\n", Trees->NoTaxa);
	for(Index=0;Index<Trees->NoTaxa;Index++)
		fprintf(Opt->VarRatesLog, "%d\t%s\n", Trees->Taxa[Index]->No, Trees->Taxa[Index]->Name);

	fprintf(Opt->VarRatesLog, "%d\n", T->NoNodes);
	for(Index=0;Index<T->NoNodes;Index++)
	{
		N = T->NodeList[Index];
		No = 0;
		CTaxaBelow(N, &No);
		fprintf(Opt->VarRatesLog, "%zu\t%f\t%d\t", N->Part->PartID, N->Length, No);
		RecPrintPPNodes(Opt->VarRatesLog, N);
		fprintf(Opt->VarRatesLog, "\n");
	}

	fprintf(Opt->VarRatesLog, "It\tLh\tLh + Prior\tNo Pram\tAlpha\tSigma^2\tAlpha Scale Prior\t");
	fprintf(Opt->VarRatesLog, "Node ID\tScaler\tCreat It\tNode / Branch\t");

	fprintf(Opt->VarRatesLog, "\n");

	fflush(Opt->VarRatesLog);
}

void	VarRatesLogFileHeader(OPTIONS *Opt, TREES *Trees, RATES *Rates)
{
	PART		*Part;	
	int			Index;
	VARRATES	*VarRates;
	FILE		*Str;


	Str = Opt->VarRatesLog;
	VarRates = Rates->VarRates;


	fprintf(Str, "No Parts:\t%zu\n", Trees->NoParts);
	fprintf(Str, "PartID\tFrequency\tProbability\tNoTaxa\tTaxa\n");

	for(Index=0;Index<Trees->NoParts;Index++)
	{
		Part = Trees->PartList[Index];
		PrintPart(Str, Trees, Part);
		fprintf(Str, "\n");
	}

	fprintf(Str, "It\tLh\tLh + Prior\tTree No\tNo Pram\tAlpha\tSigma^2\tAlpha Scale Prior\t");
	fprintf(Str, "Part ID\tScaler\tCreat It\tNode / Branch\t");
	fprintf(Str, "\n");
	
/*
	fprintf(Opt->VarRatesLog, "%d\n", Trees->NoTaxa);
	for(Index=0;Index<Trees->NoTaxa;Index++)
		fprintf(Opt->VarRatesLog, "%d\t%s\n", Trees->Taxa[Index]->No, Trees->Taxa[Index]->Name);

	fprintf(Opt->VarRatesLog, "%d\n", T->NoNodes);
	for(Index=0;Index<T->NoNodes;Index++)
	{
		N = T->NodeList[Index];
		No = 0;
		CTaxaBelow(N, &No);
		fprintf(Opt->VarRatesLog, "%d\t%f\t%d\t", N->ID, N->Length, No);
		RecPrintPPNodes(Opt->VarRatesLog, N);
		fprintf(Opt->VarRatesLog, "\n");
	}


	fprintf(Opt->VarRatesLog, "\n");
	*/
}

void	IntiVarRatesLogFile(OPTIONS *Opt, TREES *Trees, RATES *Rates)
{
	Opt->VarRatesLog = OpenWriteWithExt(Opt->BaseOutputFN, OUTPUT_EXT_VAR_RATES);

	if(Trees->NoTrees == 1)
		VarRatesLogFileHeaderSingleTree(Opt, Trees, Rates);
	else
		VarRatesLogFileHeader(Opt, Trees, Rates);


	fflush(Opt->VarRatesLog);
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
	IntiVarRatesLogFile(Opt, Trees, Rates);
}

void	FinishVarRatesFiles(OPTIONS *Opt)
{
	fclose(Opt->VarRatesLog);
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
	
	if(Type == VR_LS_BL)
		fprintf(Out, "LS_Beta\t");
}

void	LogVarRatesResults(OPTIONS *Opt, TREES *Trees, RATES *Rates, long long It)
{
	FILE		*Out;
	VARRATES		*VarRates;
	int			Index;
	VAR_RATES_NODE	*PNode;
	NODE		N;
	double		Sigma;

	if(Opt->UseStochasticBeta == TRUE && UseRJLocalScalar(Opt) == FALSE)
	{
		fprintf(Opt->VarRatesLog, "%lld\t%f\t%f\t%d\t", It, Rates->Lh, Rates->Lh + Rates->LhPrior, GetNoStochasticBetaType(Rates, SB_RJ));
		LogStochasticBetaResults(Opt, Trees, Rates, It);
		fprintf(Opt->VarRatesLog,"\n");
		return;
	}

	VarRates = Rates->VarRates;
//	CheckPlasyNodes(P);

	Out = Opt->VarRatesLog;
	
	if(Trees->NoTrees == 1)
		fprintf(Out, "%lld\t%f\t%f\t%d\t", It, Rates->Lh, Rates->Lh + Rates->LhPrior, VarRates->NoNodes);
	else
		fprintf(Out, "%lld\t%f\t%f\t%d\t%d\t", It, Rates->Lh, Rates->Lh + Rates->LhPrior, Rates->TreeNo, VarRates->NoNodes);

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
		N = PNode->NodeList[Rates->TreeNo];

		fprintf(Out, "%zu\t", PNode->Part->PartID);
		fprintf(Out, "%f\t", PNode->Scale);
		fprintf(Out, "%llu\t", PNode->NodeID);

		OutputVarRatesType(Out, PNode->Type);
//		fprintf(Out, "\n");
	}

	if(Opt->UseStochasticBeta == TRUE)
		LogStochasticBetaResults(Opt,Trees,Rates,It);

	fprintf(Out, "\n");
//	exit(0);
	fflush(Out);
}	

void	PrintVarRatesOutput(OPTIONS *Opt, TREES *Trees, RATES *Rates, long long It)
{
	LogVarRatesResults(Opt, Trees, Rates, It);
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

	if(Type == VR_LS_BL)
		return GetPriorFromName("VR_LS_BL", Rates->Priors, Rates->NoPriors);

	return NULL;
}

double	CaclVRPrior(double X, TRANSFORM_TYPE Type, RATES *Rates)
{
	double Ret;
	PRIOR *Prior;

//	return -2.0;

	if(Type == VR_KAPPA && X >= MAX_VR_KAPPA)
		return ERRLH;

	if(Type == VR_DELTA && X >= MAX_VR_DELTA)
		return ERRLH;

	if(Type == VR_LAMBDA && X >= MAX_VR_LAMBDA)
		return ERRLH;

	if(Type == VR_OU && X >= MAX_VR_OU)
		return ERRLH;

	Prior = GetVRPrior(Type, Rates);
	Ret = CalcLhPriorP(X, Prior);

	return Ret;
}

double	CalcVRLandPriorLh(double Beta, double T, double Sig2)
{
	double Lh, Z;

	Z = (Beta * T) / sqrt(Sig2 * T);
	Z = fabs(Z);

	Z = Z - 3;
	if(Z < 0)
		Z = 0;

	Lh = CaclNormalLogLh(Z, 1.0, 1.0);

	return Lh;
}

double	CaclVRLandPrior(RATES *Rates, VAR_RATES_NODE *VRNode)
{
	double Ret, X;
	PRIOR *Prior;

	X = VRNode->Scale;
	X = fabs(X);

	Prior = GetVRPrior(VRNode->Type, Rates);
	Ret = CalcLhPriorP(X, Prior);

	return Ret;
/*
	double T, Beta, Sig2, Z, LhP;
	NODE N;

//	return 0;

	N = VRNode->NodeList[Rates->TreeNo];
	T = N->Length;

	Beta = VRNode->Scale;

	Sig2 = Rates->Rates[1];

	LhP = CalcVRLandPriorLh(Beta, T, Sig2);

	return LhP;*/
}

void	TestR(RATES *Rates)
{
	TRANSFORM_TYPE Type;
	PRIOR *Prior;
	double X, Lh;

	Type = VR_LS_BL;

	Prior = GetVRPrior(Type, Rates);
	
	for(X=0.0001;X<3;X+=0.01)
	{
		Lh = CalcLhPriorP(X, Prior);

		printf("%f\t%f\t%f\n", X, Lh, exp(Lh));
	}
	
	exit(0);
}

double CalcDiffTHoldCost(RATES *Rates, OPTIONS* Opt)
{
	double Ret;
	int Index;
	VARRATES* VarRates;
	VAR_RATES_NODE* PNode;

	VarRates = Rates->VarRates;

	Ret = 0;
	for(Index=0;Index<VarRates->NoNodes;Index++)
	{
		PNode = VarRates->NodeList[Index];

		if(PNode->Type == VR_NODE)
			Ret += -1.0;
		else
			Ret += Opt->RJThreshold;

	}
	
	return Ret;
}

double	CalcVarRatesPriors(RATES *Rates,OPTIONS *Opt)
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

		if(PNode->Type != VR_LS_BL)
			PVal = CaclVRPrior(PNode->Scale,PNode->Type,Rates);
		else
			PVal = CaclVRLandPrior(Rates,PNode);

		if(PVal == ERRLH)
			return ERRLH;

		Ret += PVal;
	}

	// Add a theshold cost for all VR paramiters
	//Ret += VarRates->NoNodes * Opt->RJThreshold;
	Ret += CalcDiffTHoldCost(Rates, Opt);

	return Ret;
}


/*
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

		if(PNode->Type != VR_LS_BL)
			PVal = CaclVRPrior(PNode->Scale, PNode->Type, Rates);
		else
			PVal = CaclVRLandPrior(Rates, PNode);

		if(PVal == ERRLH)
			return ERRLH;

		Ret += PVal;
	}
	
	// Add a theshold cost for all VR paramiters
	Ret += VarRates->NoNodes * Opt->RJThreshold;
//	Ret += CalcDiffTHoldCost(Rates, Opt);
	
	return Ret;
}
*/
VAR_RATES_NODE*	CreateTextPNode(NODE N, double Scale, long long CIt, TRANSFORM_TYPE Type)
{
	VAR_RATES_NODE*	Ret;
	
	Ret = CreateVarRatesNode(CIt, N->Part, 1);

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

	if(strcmp("ls_beta", Str) == 0)
		return VR_LS_BL;

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

	Plasty->NodeList = (VAR_RATES_NODE**)AddToList(&Plasty->NoNodes, (void**)Plasty->NodeList, (void*)PNode);
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

	PrintVarRatesTree(Opt, Opt->Trees, Rates, 1);
		
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

	SetUserBranchLength(Tree);

	Lh = Likelihood(Rates, Opt->Trees, Opt);

	printf("Lh = %f\n", Lh);
//	exit(0);

	OneRateAtOnce(Opt, Rates);
	
	RemoveEachVarRate(Opt, Rates);



	TestEachVarRatePair(Opt, Rates);
	
//	PrintVarRatesTree(Opt, Opt->Trees, Rates, 1);

	exit(0);
}


void	SetVarRatesFromStr(RATES *Rates, OPTIONS *Opt, char *Str)
{
	char **Passed;
	int Index, Tokes;
	char *S;
	double Lh;
	
	S = StrMake(Str);

	Passed = (char**)SMalloc(sizeof(char*) * strlen(Str));
	
	Tokes = MakeArgv(S, Passed, (int)strlen(S));

	Rates->Rates[0] = atof(Passed[4]);
	Rates->Rates[1] = atof(Passed[5]);

	Index = 7;
	for(;Index<Tokes;Index+=4)
		AddTextVarRate(Opt->Trees->Tree[0], Rates, 4, &Passed[Index]);


	Lh = Likelihood(Rates, Opt->Trees, Opt);

//	printf("%f\n", Lh);
//	exit(0);

	
//	PrintVarRatesOutput(Opt, Opt->Trees, Rates, 1); exit(0);

//	TestVarRates(Opt, Rates);

	free(S);
	free(Passed);
}