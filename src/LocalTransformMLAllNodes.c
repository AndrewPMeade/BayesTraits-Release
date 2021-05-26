#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "LocalTransformMLAllNodes.h"
#include "LocalTransform.h"
#include "VarRates.h"
#include "ML.h"
#include "GenLib.h"
#include "Tag.h"
#include "Likelihood.h"
#include "Trees.h"


LT_ALL_NODES*	AllocLTMLAllNodes(int NoNodes)
{
	int Index;

	LT_ALL_NODES* Ret;

	Ret = (LT_ALL_NODES*)SMalloc(sizeof(LT_ALL_NODES));

	Ret->NoNodes = NoNodes;
	Ret->TagList = (TAG**)SMalloc(sizeof(TAG*) * NoNodes);

	Ret->ValidBetaNodes = (int*)SMalloc(sizeof(int) * NoNodes);
	Ret->ValidBLNodes = (int*)SMalloc(sizeof(int) * NoNodes);
	Ret->ValidNodeNodes = (int*)SMalloc(sizeof(int) * NoNodes);

	for(Index=0;Index<NoNodes;Index++)
	{
		Ret->TagList[Index] = NULL;

		Ret->ValidBetaNodes[Index] = FALSE;
		Ret->ValidBLNodes[Index] = FALSE;
		Ret->ValidNodeNodes[Index] = FALSE;
	}


	Ret->EstBeta = FALSE;
	Ret->EstBL = FALSE;
	Ret->EstNodes = FALSE;

	return Ret;
}

void	FreeLTMLAllNodes(LT_ALL_NODES *Data)
{
	int Index;

	for(Index=0;Index<Data->NoNodes;Index++)
		FreeTag(Data->TagList[Index]);
	free(Data->TagList);

	free(Data->ValidBetaNodes);
	free(Data->ValidBLNodes);
	free(Data->ValidNodeNodes);

	free(Data);
}

int		ValidNode(NODE Node, TRANSFORM_TYPE TType)
{
	if(Node->Ans == NULL)
		return FALSE;

	if(Node->Tip == FALSE || TType == VR_BL)
		return TRUE;

	return FALSE;
}

void	SetValidNodes(LT_ALL_NODES *Data, TREE *Tree)
{
	int Index;

	for(Index=0;Index<Data->NoNodes;Index++)
	{
		Data->ValidNodeNodes[Index] = ValidNode(Tree->NodeList[Index], VR_NODE);
		Data->ValidBLNodes[Index] = ValidNode(Tree->NodeList[Index], VR_BL);
		Data->ValidBetaNodes[Index] = ValidNode(Tree->NodeList[Index], VR_LS_BL);
	}
}

void	GetTagTaxaNames(NODE Node, int *NoTaxa, char **TagNames)
{
	int Index;
	NODE N;

	if(Node->Tip == TRUE)
	{
		TagNames[*NoTaxa] = Node->Taxa->Name;
		(*NoTaxa)++;
		return;
	}

	for(Index=0;Index<Node->NoNodes;Index++)
	{
		N = Node->NodeList[Index];
		GetTagTaxaNames(N, NoTaxa, TagNames);
	}
}

void	SetAllTags(LT_ALL_NODES *Data, TREES *Trees)
{
	TREE *Tree;
	NODE Node;
	int Index, NoTaxa;
	char **TaxaNames;
	char *TagName;

	Tree = Trees->Tree[0];
	TagName = (char*)SMalloc(sizeof(char) * MAX_LT_NAME_SIZE);
	TaxaNames = (char**)SMalloc(sizeof(char*) * Tree->NoNodes);


	for(Index=0;Index<Tree->NoNodes;Index++)
	{
		Node = Tree->NodeList[Index];

		NoTaxa = 0;
		GetTagTaxaNames(Node, &NoTaxa, TaxaNames);

		sprintf(TagName, "Tag-%010d", Index);

		Data->TagList[Index] = CreateTag(Trees, TagName, NoTaxa, TaxaNames);

	}
	
	free(TagName);
	free(TaxaNames);
}

LT_ALL_NODES*	CreatLTMLAllNodes(TREES *Trees, int EstNodes, int EstBL, int EstBeta)
{
	LT_ALL_NODES *Ret;
	TREE *Tree;

	Tree = Trees->Tree[0];

	Ret = AllocLTMLAllNodes(Tree->NoNodes);

	if(EstNodes == TRUE)
		Ret->EstNodes = TRUE;

	if(EstBL == TRUE)
		Ret->EstBL = TRUE;

	if(EstBeta == TRUE)
		Ret->EstBeta = TRUE;

	SetValidNodes(Ret, Tree);
	SetAllTags(Ret, Trees);

	return Ret;
}

LOCAL_TRANSFORM*	CreateNewLocalTransforms(int No, TRANSFORM_TYPE TType)
{
	char *Name;
	LOCAL_TRANSFORM *LT;
	TAG **TList;

	Name = (char*)SMalloc(sizeof(char) * MAX_LT_NAME_SIZE);
	sprintf(Name, "LT-%10d", No);

	TList = (TAG**)SMalloc(sizeof(TAG*));
	TList[0] = NULL;

	LT = CreateLocalTransforms(Name, TList, 1, TType, TRUE, 1.0);

	free(TList);

	free(Name);
	return(LT);
}

void	AddScalar(LT_ALL_NODES* Data, TRANSFORM_TYPE TType, OPTIONS *Opt, TREES *Trees, RATES *Rates)
{
	LOCAL_TRANSFORM *LT;
	int NoLT;

	NoLT = Rates->NoLocalTransforms;

	LT = CreateNewLocalTransforms(NoLT, TType);

	Rates->LocalTransforms[NoLT] = LT;
	Rates->NoLocalTransforms++;
}

void	RemoveScalar(RATES *Rates)
{
	FreeLocalTransforms(Rates->LocalTransforms[Rates->NoLocalTransforms-1]);
	Rates->NoLocalTransforms--;
}

double	MLLTScalar(NODE N, LOCAL_TRANSFORM *LT, OPTIONS *Opt, TREES *Trees, RATES *Rates)
{
	double InitLh, OptLh;

	InitLh = Likelihood(Rates, Trees, Opt);

	MLTree(Opt, Trees, Rates);

	OptLh = Likelihood(Rates, Trees, Opt);

	return OptLh - InitLh;
}

int		GetMaxLhGainPos(LT_ALL_NODES* Data, double *LhList, int *ValidNodes)
{
	int Index, Ret;

	Ret = -1;

	for(Index=0;Index<Data->NoNodes;Index++)
	{
		if(ValidNodes[Index] == TRUE)
		{
			if(Ret == -1)
				Ret = Index;
			else
			{
				if(LhList[Index] > LhList[Ret])
					Ret = Index;
			}
		}
	}

	return Ret;
}

void	PostProcNewScalar(LT_ALL_NODES* Data, int *ValidNodes, LOCAL_TRANSFORM *LT, RATES *Rates, double *LhList, double *ScaleList)
{
	int MaxPos;

	MaxPos = GetMaxLhGainPos(Data, LhList, ValidNodes);

	LT->TagList[0] = Data->TagList[MaxPos];
	LT->Scale = ScaleList[MaxPos];
	LT->Est = FALSE;
}

void	LhTestSurce(OPTIONS *Opt, TREES *Trees, RATES *Rates)
{
	double Val, Lh;
	LOCAL_TRANSFORM *LT;

	LT = Rates->LocalTransforms[Rates->NoLocalTransforms-1];

	for(Val=5;Val<10;Val+=1.0)
	{
		LT->Scale = Val;

		Lh = Likelihood(Rates, Trees, Opt);

		printf("%f\t%f\n",Val, Lh);
	}

	exit(0);

}

void	IncludeNewScalarType(LT_ALL_NODES* Data, TRANSFORM_TYPE TType, int *ValidNodes, OPTIONS *Opt, TREES *Trees, RATES *Rates)
{
	TREE *Tree;
	int NIndex;
	NODE Node;
	TAG	*Tag;
	LOCAL_TRANSFORM *LT;
	double *LhGain, *Scale;
	
	Tree = Trees->Tree[0];
	AddScalar(Data, TType, Opt, Trees, Rates);
	LT = Rates->LocalTransforms[Rates->NoLocalTransforms-1];

	LhGain = (double*)SMalloc(sizeof(double)*Data->NoNodes);
	Scale = (double*)SMalloc(sizeof(double)*Data->NoNodes);

	for(NIndex=0;NIndex<Data->NoNodes;NIndex++)
	{
		if(ValidNodes[NIndex] == TRUE)
		{
			Node = Tree->NodeList[NIndex];
			Tag = Data->TagList[NIndex];
			LT->TagList[0] = Tag;

			LT->Scale = 1.0;
			if(LT->Type == VR_LS_BL)
				LT->Scale = 0.0;			

			LhGain[NIndex] = MLLTScalar(Node, LT, Opt, Trees, Rates);
			Scale[NIndex] = LT->Scale;
		}
	}

	PostProcNewScalar(Data, ValidNodes, LT, Rates, LhGain, Scale);
	
	free(LhGain);
	free(Scale);
}

void	AddLTNode(RATES *Rates, LOCAL_TRANSFORM *LTNode)
{
	Rates->LocalTransforms[Rates->NoLocalTransforms++] = LTNode;
}

void	SetNodeAsInValid(LT_ALL_NODES* Data, LOCAL_TRANSFORM *LTNode, int *NodeList)
{
	int Index;

	for(Index=0;Index<Data->NoNodes;Index++)
	{
		if(LTNode->TagList[0] == Data->TagList[Index])
		{
			NodeList[Index] = FALSE;
			return;
		}
	}
	
}

void	IncludeNewScalar(LT_ALL_NODES* Data, OPTIONS *Opt, TREES *Trees, RATES *Rates)
{
	double InitLh, LhGainNode, LhGainBeta, LhGainBL;
	LOCAL_TRANSFORM *LTNode, *LTBeta, *LTBL;

	InitLh = Likelihood(Rates, Trees, Opt);

	LhGainNode = LhGainBeta = LhGainBL = -1;
	LTNode = LTBeta = LTBL = NULL;

	if(Data->EstNodes == TRUE)
	{
		IncludeNewScalarType(Data, VR_NODE, Data->ValidNodeNodes, Opt, Trees, Rates);
		LTNode = Rates->LocalTransforms[Rates->NoLocalTransforms-1];
		LhGainNode = Likelihood(Rates, Trees, Opt) - InitLh;
		Rates->NoLocalTransforms--;
	}

	if(Data->EstBeta == TRUE)
	{
		IncludeNewScalarType(Data, VR_LS_BL, Data->ValidBetaNodes, Opt, Trees, Rates);
		LTBeta = Rates->LocalTransforms[Rates->NoLocalTransforms-1];
		LhGainBeta = Likelihood(Rates, Trees, Opt) - InitLh;
		Rates->NoLocalTransforms--;
	}

	if(Data->EstBL == TRUE)
	{
		IncludeNewScalarType(Data, VR_BL, Data->ValidBLNodes, Opt, Trees, Rates);
		LTBL = Rates->LocalTransforms[Rates->NoLocalTransforms-1];
		LhGainBL = Likelihood(Rates, Trees, Opt) - InitLh;
		Rates->NoLocalTransforms--;
	}


	if(LhGainNode > LhGainBeta && LhGainNode > LhGainBL)
	{		
		AddLTNode(Rates, LTNode);
		SetNodeAsInValid(Data, LTNode, Data->ValidNodeNodes);
		FreeLocalTransforms(LTBeta);
		FreeLocalTransforms(LTBL);
		return;
	}

	if(LhGainBeta > LhGainBL)
	{		
		AddLTNode(Rates, LTBeta);
		SetNodeAsInValid(Data, LTBeta, Data->ValidBetaNodes);
		FreeLocalTransforms(LTNode);
		FreeLocalTransforms(LTBL);
		return;
	}

	if(LhGainBL != -1)
	{
		AddLTNode(Rates, LTBL);
		SetNodeAsInValid(Data, LTBL, Data->ValidBetaNodes);
		FreeLocalTransforms(LTNode);
		FreeLocalTransforms(LTBeta);
		return;
	}
}

void	IntiLocalTransforms(TREES *Trees, RATES *Rates)
{
	TREE *Tree;

	Tree = Trees->Tree[0];
	Rates->NoLocalTransforms = 0;

	Rates->LocalTransforms = (LOCAL_TRANSFORM**)SMalloc(sizeof(LOCAL_TRANSFORM*) * Tree->NoNodes * 10);
}

void	OutputNewScalar(double BaseLh, double OptLh, double GOptLh, RATES *Rates)
{
	LOCAL_TRANSFORM *LT;

	LT = Rates->LocalTransforms[Rates->NoLocalTransforms-1];

	printf("Adding\t%d\t%f\t%f\t%f\t%f\t", Rates->NoLocalTransforms, BaseLh, OptLh, GOptLh, LT->Scale);
	OutputVarRatesType(stdout, LT->Type);

	PrintTag(stdout, LT->TagList[0]);

	fflush(stdout);
}

void	OutputAllScalar(RATES *Rates)
{
	int Index;
	LOCAL_TRANSFORM *LT;

	for(Index=0;Index<Rates->NoLocalTransforms;Index++)
	{
		LT = Rates->LocalTransforms[Index];
		printf("Final\t%d\t%f\t", Index, LT->Scale);
		OutputVarRatesType(stdout, LT->Type);
		PrintTag(stdout, LT->TagList[0]);
	}
	fflush(stdout);
}

void	OutputTrees(OPTIONS *Opt, TREES *Trees, RATES *Rates)
{
	char *FName;

	FName = (char*)SMalloc(sizeof(char) * (strlen(OUTPUT_EXT_TREES) + strlen(Opt->BaseOutputFN) + 1));
	sprintf(FName, "%s%s", Opt->BaseOutputFN, OUTPUT_EXT_TREES);

	SaveTrees(FName, Trees);

	free(FName);
}

void	SetAllEstLT(RATES *Rates, int Est)
{
	int Index;

	for(Index=0;Index<Rates->NoLocalTransforms;Index++)
		Rates->LocalTransforms[Index]->Est = Est;
}

void	LTMLTreeGlobal(OPTIONS *Opt, TREES *Trees, RATES *Rates)
{
	ML_MAP*	CMap, *BMap;
	double CLh, BLh;
	int Index;
	

	BMap = AllocMLMap();

	BuildMLMap(BMap, Opt, Trees, Rates);

	for(Index=0;Index<Rates->NoLocalTransforms;Index++)
		BMap->PVal[Index] = Rates->LocalTransforms[Index]->Scale;

	BLh = LikelihoodML(BMap, Opt, Trees, Rates);

	for(Index=0;Index<Opt->MLTries;Index++)
	{
		CMap = MLMapTreeTry(Opt, Trees, Rates, BMap);
		CLh = LikelihoodML(CMap, Opt, Trees, Rates);
			
		if(CLh > BLh)
		{
			CopyMLMap(BMap, CMap);
			BLh = CLh;
		}

		FreeMLMap(CMap);
	}
	
	Rates->Lh = LikelihoodML(BMap, Opt, Trees, Rates);
	FreeMLMap(BMap);
}



double	GlobalOpt(OPTIONS *Opt, TREES *Trees, RATES *Rates)
{
	double Pre, Post;

	Pre = Likelihood(Rates, Trees, Opt);
	
	SetAllEstLT(Rates, TRUE);

	LTMLTreeGlobal(Opt, Trees, Rates);

	SetAllEstLT(Rates, FALSE);

	Post = Likelihood(Rates, Trees, Opt);

	if(Pre > Post)
	{
		printf("eee\n");
		exit(0);
	}

	return Likelihood(Rates, Trees, Opt);
}

void	LocalTransformMLAllNodes(OPTIONS *Opt, TREES *Trees, RATES *Rates, int EstNodes, int EstBL, int EstBeta)
{
	LT_ALL_NODES* Data;
	int Valid;
	double InitLh, OptLh, GOptLh;

	IntiLocalTransforms(Trees, Rates);

	Data = CreatLTMLAllNodes(Trees, EstNodes, EstBL, EstBeta);
	
	InitLh = Likelihood(Rates, Trees, Opt);
	printf("IntLH:\t%f\n", InitLh);

	do
	{
		InitLh = Likelihood(Rates, Trees, Opt);
		IncludeNewScalar(Data, Opt, Trees, Rates);
		OptLh = Likelihood(Rates, Trees, Opt);
		Rates->Lh = OptLh;

		
		if(OptLh - InitLh > ML_CUT_POINT)
			Valid = TRUE;
		else
		{
			RemoveScalar(Rates);
			Valid= FALSE;
		}

		if(Valid == TRUE)
		{
			GOptLh = GlobalOpt(Opt, Trees, Rates);
			OutputNewScalar(InitLh, OptLh, GOptLh, Rates);
		}
	}while(Valid == TRUE);

	OptLh = Likelihood(Rates, Trees, Opt);
	
	printf("OptLh:\t%f\n", OptLh);

	OutputTrees(Opt, Trees, Rates);

	OutputAllScalar(Rates);	
	
	FreeLTMLAllNodes(Data);

	exit(0);
}