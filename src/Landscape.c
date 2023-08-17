#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

#include "TypeDef.h"
#include "Trees.h"
#include "GenLib.h"
#include "Part.h"
#include "VarRates.h"
#include "Likelihood.h"
#include "ML.h"
#include "NLOptBT.h"
#include "Contrasts.h"

int			UseLandscapeBeta(OPTIONS* Opt, RATES *Rates)
{
	int Index;

	if(Opt->UseRJLocalScalar[VR_LS_BL] == TRUE)
		return TRUE;

	for(Index=0;Index<Rates->NoLocalTransforms;Index++)
		if(Rates->LocalTransforms[Index]->Type == VR_LS_BL)
			return TRUE;

	return FALSE;
}

void		ResetTreeLandscape(TREE *Tree)
{
	int Index;
	for(Index=0;Index<Tree->NoNodes;Index++)
		Tree->NodeList[Index]->LandscapeBeta = 0.0;
}

void		SetLandscapeBetaTrait(NODE Node, int NoSites, double Change)
{
	CONDATA *Con;
	int SIndex;

	Con = Node->ConData;

	for(SIndex=0;SIndex<NoSites;SIndex++)
		Con->Contrast[0]->Data[SIndex] = Node->Taxa->ConData[SIndex] - Change;
}

void		PropLandscapeBeta(TREES *Trees, NODE Node, double Change)
{
	int Index;

	Change += Node->LandscapeBeta;

	if(Node->Tip == TRUE)
	{
		SetLandscapeBetaTrait(Node, Trees->NoSites, Change);
		return;
	}
	
	for(Index=0;Index<Node->NoNodes;Index++)
		PropLandscapeBeta(Trees, Node->NodeList[Index], Change);
}

void		SetBetaNode(NODE N, double Beta)
{
	int Index;

	N->LandscapeBeta = Beta;

	for(Index=0;Index<N->NoNodes;Index++)
		SetBetaNode(N->NodeList[Index], Beta);
}

int			NodeSubSet(OPTIONS *Opt, TREES *Trees, RATES *Rates)
{
	int i,j;
	VARRATES* VarRates;
	VAR_RATES_NODE *VRN_i, *VRN_j;

	VarRates = Rates->VarRates;

	for(i=0;i<VarRates->NoNodes;i++)
	{
		VRN_i = VarRates->NodeList[i];

		if(VRN_i->Type == VR_LS_BL)
		{
			for(j=0;j<VarRates->NoNodes;j++)
			{
				VRN_j = VarRates->NodeList[j];

				if(i != j)
				{
					if(PartSubSet(VRN_i->Part, VRN_j->Part) == TRUE)
						return TRUE;
				}
			}

		}
	}

	return FALSE;
}

void		MapRJLandscape(OPTIONS *Opt, TREES *Trees, RATES *Rates)
{
	NODE Node;
	
	VARRATES* VarRates;
	VAR_RATES_NODE *VR_Node;
	int TreeNo, Index;
	
	TreeNo = Rates->TreeNo;

	if(Rates->VarRates == NULL)
		return;

	VarRates = Rates->VarRates;

	for(Index=0;Index<VarRates->NoNodes;Index++)
	{
		VR_Node = VarRates->NodeList[Index];

		if(VR_Node->Type == VR_LS_BL)
		{
			Node = GetVRNode(Trees, TreeNo, VR_Node);

			Node->LandscapeBeta = VR_Node->Scale;
		}
	}
}

void		MapLocalTranfromsBeta(OPTIONS *Opt, TREES *Trees, RATES *Rates)
{
	LOCAL_TRANSFORM *LRate;
	int TIndex, Index;
	TREE *Tree;
	NODE N;

	Tree = Trees->Tree[Rates->TreeNo];

	for(Index=0;Index<Rates->NoLocalTransforms;Index++)
	{
		LRate = Rates->LocalTransforms[Index];

		if(LRate->Type == VR_LS_BL)
		{
			for(TIndex=0;TIndex<LRate->NoTags;TIndex++)
			{
				N = LRate->TagList[TIndex]->NodeList[Rates->TreeNo];
				N->LandscapeBeta = LRate->Scale;
			}
		}
	}	
}
void		MapLandscape(OPTIONS *Opt, TREES *Trees, RATES *Rates)
{
	TREE *Tree;
	int TreeNo;


	TreeNo = Rates->TreeNo;
	Tree = Trees->Tree[TreeNo];
	
	ResetTreeLandscape(Tree);

	MapLocalTranfromsBeta(Opt, Trees, Rates);
	
	MapRJLandscape(Opt, Trees, Rates);

	PropLandscapeBeta(Trees, Tree->Root, 0.0);
}



void	AncRecCalcContrast(NODE Node, int NoSites, double *TempAnc, double GlobalBeta)
{
	int		Index;
	CONDATA *Con;

	if(Node->Tip == TRUE)
		return;

	for(Index=0;Index<Node->NoNodes;Index++)
		AncRecCalcContrast(Node->NodeList[Index], NoSites, TempAnc, GlobalBeta);

	CalcNodeContrast(Node, NoSites);

	Con = Node->ConData;
	TempAnc[Node->ID] = Con->Contrast[0]->Data[0];
	Con->Contrast[0]->Data[0] = Con->Contrast[0]->Data[0] - Node->LandscapeBeta - (Node->Length * GlobalBeta);
}

void	SetFabricAncStatesTipData(OPTIONS *Opt, TREE *Tree, RATES *Rates)
{
	int Index;
	NODE Node;


	for(Index=0;Index<Tree->NoNodes;Index++)
	{
		Node = Tree->NodeList[Index];
		if(Node->Tip == TRUE)
			Node->ConData->Contrast[0]->Data[0] = Node->ConData->Contrast[0]->Data[0] - Node->LandscapeBeta - (Node->Length * Rates->GlobalTrend);
	}

}

// use this to test for it. 
// 	if(UseLandscapeBeta(Opt, Rates) == TRUE)
void	SetFabricAncStates(OPTIONS *Opt, TREES *Trees, RATES *Rates)
{
	TREE *Tree;
	int Index;
	NODE Node;
	double *TempAnc;

//	return;

	Tree = Trees->Tree[Rates->TreeNo];

	TempAnc = (double*)SMalloc(sizeof(double) * Tree->NoNodes);

	LhTransformTree(Rates, Trees, Opt);

	RetSetConTraitData(Trees->Tree[Rates->TreeNo], Trees->NoSites);

	SetFabricAncStatesTipData(Opt, Tree, Rates);

	AncRecCalcContrast(Tree->Root, Trees->NoSites, TempAnc, Rates->GlobalTrend);

	for(Index=0;Index<Tree->NoNodes;Index++)
	{
		Node = Tree->NodeList[Index];
		if(Node->Tip == FALSE)
			Node->ConData->Contrast[0]->Data[0] = TempAnc[Node->ID];
	}

	free(TempAnc);

//	This will need to be done, but not at this point. 
//	Likelihood(Rates, Trees, Opt);
}