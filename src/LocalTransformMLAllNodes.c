#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "LocalTransformMLAllNodes.h"
#include "LocalTransform.h"
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
	Ret->ValidNodes = (int*)SMalloc(sizeof(int) * NoNodes);

	for(Index=0;Index<NoNodes;Index++)
	{
		Ret->TagList[Index] = NULL;
		Ret->ValidNodes[Index] = FALSE;
	}

	return Ret;
}

int		ValidNode(NODE Node, TRANSFORM_TYPE TType)
{
	if(Node->Ans == NULL)
		return FALSE;

	if(Node->Tip == FALSE || TType == VR_BL)
		return TRUE;

	return FALSE;
}

void	SetValidNodes(LT_ALL_NODES *Data, TREE *Tree, TRANSFORM_TYPE TType)
{
	int Index;

	for(Index=0;Index<Data->NoNodes;Index++)
		Data->ValidNodes[Index] = ValidNode(Tree->NodeList[Index], TType);
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

LT_ALL_NODES*	CreatLTMLAllNodes(TREES *Trees, TRANSFORM_TYPE TType)
{
	LT_ALL_NODES *Ret;
	TREE *Tree;

	Tree = Trees->Tree[0];

	Ret = AllocLTMLAllNodes(Tree->NoNodes);

	SetValidNodes(Ret, Tree, TType);
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

int		GetMaxLhGainPos(LT_ALL_NODES* Data, double *LhList)
{
	int Index, Ret;

	Ret = -1;

	for(Index=0;Index<Data->NoNodes;Index++)
	{
		if(Data->ValidNodes[Index] == TRUE)
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

int		PostProcNewScalar(LT_ALL_NODES* Data, LOCAL_TRANSFORM *LT, RATES *Rates, double *LhList, double *ScaleList)
{
	int MaxPos;

	MaxPos = GetMaxLhGainPos(Data, LhList);

	if(LhList[MaxPos] > ML_CUT_POINT)
	{
		LT->TagList[0] = Data->TagList[MaxPos];
		LT->Scale = ScaleList[MaxPos];
		LT->Est = FALSE;
		Data->ValidNodes[MaxPos] = FALSE;

		return TRUE;
	}

	return FALSE;
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

int		IncludeNewScalar(LT_ALL_NODES* Data, TRANSFORM_TYPE TType, OPTIONS *Opt, TREES *Trees, RATES *Rates)
{
	TREE *Tree;
	int NIndex, Ret;
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
		if(Data->ValidNodes[NIndex] == TRUE)
		{
			Node = Tree->NodeList[NIndex];
			Tag = Data->TagList[NIndex];
			LT->TagList[0] = Tag;

			LT->Scale = 1.0;
			if(LT->Type == VR_LS_BL)
				LT->Scale = 0.0;			

//			LhTestSurce(Opt, Trees, Rates);
			LhGain[NIndex] = MLLTScalar(Node, LT, Opt, Trees, Rates);
			Scale[NIndex] = LT->Scale;

		}
	}

	Ret = PostProcNewScalar(Data, LT, Rates, LhGain, Scale);
	
	if(Ret == FALSE)
		RemoveScalar(Rates);
	
	free(LhGain);
	free(Scale);

	return Ret;
}

void	IntiLocalTransforms(TREES *Trees, RATES *Rates)
{
	TREE *Tree;

	Tree = Trees->Tree[0];
	Rates->NoLocalTransforms = 0;

	Rates->LocalTransforms = (LOCAL_TRANSFORM**)SMalloc(sizeof(LOCAL_TRANSFORM*) * Tree->NoNodes);
}

void	OutputNewScalar(double BaseLh, double OptLh, double GOptLh, RATES *Rates)
{
	LOCAL_TRANSFORM *LT;

	LT = Rates->LocalTransforms[Rates->NoLocalTransforms-1];

	printf("%d\t%f\t%f\t%f\t%f\t", Rates->NoLocalTransforms, BaseLh, OptLh, GOptLh, LT->Scale);
	PrintTag(stdout, LT->TagList[0]);

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

double	GlobalOpt(OPTIONS *Opt, TREES *Trees, RATES *Rates)
{
	SetAllEstLT(Rates, TRUE);

	MLTree(Opt, Trees, Rates);

	SetAllEstLT(Rates, FALSE);

	return Likelihood(Rates, Trees, Opt);
}

void	LocalTransformMLAllNodes(OPTIONS *Opt, TREES *Trees, RATES *Rates)
{
	LT_ALL_NODES* Data;
	TRANSFORM_TYPE TType;
	int Valid;
	double InitLh, OptLh, GOptLh;

//	TType = VR_NODE;

	// you need to make sure a node is allready being scaled.. 
	TType = VR_LS_BL;

	IntiLocalTransforms(Trees, Rates);

	Data = CreatLTMLAllNodes(Trees, TType);
	
	do
	{
		InitLh = Likelihood(Rates, Trees, Opt);
		Valid = IncludeNewScalar(Data, TType, Opt, Trees, Rates);
		OptLh = Likelihood(Rates, Trees, Opt);
		if(Valid == TRUE)
		{
			GOptLh = GlobalOpt(Opt, Trees, Rates);
			OutputNewScalar(InitLh, OptLh, GOptLh, Rates);
		}
	}while(Valid == TRUE);

	Likelihood(Rates, Trees, Opt);
	OutputTrees(Opt, Trees, Rates);

	exit(0);
}