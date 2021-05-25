#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

#include "TypeDef.h"
#include "Trees.h"
#include "GenLib.h"
#include "Part.h"
#include "VarRates.h"

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

	Change += Node->LandscapeBeta * Node->Length;

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
	PART Pi, Pj;
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

void		MapLandscape(OPTIONS *Opt, TREES *Trees, RATES *Rates)
{
	VARRATES* VarRates;
	VAR_RATES_NODE *VR_Node;
	TREE *Tree;
	int TreeNo, Index, NIndex;
	NODE Node;
	
	TreeNo = Rates->TreeNo;
	Tree = Trees->Tree[TreeNo];
	
	VarRates = Rates->VarRates;

	ResetTreeLandscape(Tree);

	for(Index=0;Index<VarRates->NoNodes;Index++)
	{
		VR_Node = VarRates->NodeList[Index];

		if(VR_Node->Type == VR_LS_BL)
		{
			Node = GetVRNode(Trees, TreeNo, VR_Node);

	//		for(NIndex=0;NIndex<Node->NoNodes;NIndex++)
	//			SetBetaNode(Node->NodeList[NIndex], Landscape->NodeList[Index]->Beta);

			Node->LandscapeBeta = VR_Node->Scale;
		}
	}
	
//	exit(0);

//	for(Index=0;Index<Tree->Root->NoNodes;Index++)
//		PropLandscapeBeta(Trees, Tree->NodeList[Index], 0.0);

	PropLandscapeBeta(Trees, Tree->Root, 0.0);
}
