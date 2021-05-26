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

int			UseLandscapeBeta(OPTIONS* Opt)
{
	int Index;

	if(Opt->UseRJLocalScalar[VR_LS_BL] == TRUE)
		return TRUE;

	for(Index=0;Index<Opt->NoLocalTransforms;Index++)
		if(Opt->LocalTransforms[Index]->Type == VR_LS_BL)
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

/*
LANDSCAPE_NODE*	CreatMLLandscapeNode(NODE N)
{
	LANDSCAPE_NODE*	 Ret;

	Ret =  (LANDSCAPE_NODE*)SMalloc(sizeof(LANDSCAPE_NODE));

	Ret->Beta = 0;
	Ret->NodeID = 0;
	Ret->Part = NULL;

	Ret->NodeList = (NODE*) SMalloc(sizeof(NODE*));
	Ret->NodeList[0] = N;
	Ret->MLEst = TRUE;
	
	return Ret;
}

void AddLandscapeNode(LANDSCAPE *Landscape, LANDSCAPE_NODE *LandNode)
{
	Landscape->NodeList[Landscape->NoNodes] = LandNode;
	Landscape->NoNodes++;
}

LANDSCAPE_NODE* GetLandscapeNode(LANDSCAPE *Landscape, NODE N)
{
	int Index;

	for(Index=0;Index<Landscape->NoNodes;Index++)
		if(N == Landscape->NodeList[Index]->NodeList[0])
			return Landscape->NodeList[Index];

	return NULL;
}

NODE*	MakeNodeList(TREE *Tree, LANDSCAPE *Landscape, int *NoNodes)
{
	NODE*	Ret;
	int Index;
	NODE N;

	Ret = (NODE*)SMalloc(sizeof(NODE) * Tree->NoNodes);

	*NoNodes = 0;

	for(Index=0;Index<Tree->NoNodes;Index++)
	{
		N = Tree->NodeList[Index];

		if(N->Tip == FALSE && N != Tree->Root && GetLandscapeNode(Landscape, N) == NULL)
//		if(N != Tree->Root && GetLandscapeNode(Landscape, N) == NULL)
		{
			Ret[*NoNodes] = N;
			(*NoNodes)++;
		}
	}

	return Ret;
}

int	GetMaxLhPos(double *Lh, int NoNodes)
{
	int Index, Ret;

	Ret = 0;
	for(Index=0;Index<NoNodes;Index++)
		if(Lh[Index] > Lh[Ret])
			Ret = Index;

	return Ret;
}

int	SaveNewNode(RATES *Rates, double *Beta, double *Lh, NODE* NodeList, int NoNodes)
{
	int MLPos; 
	LANDSCAPE *Landscape;
	LANDSCAPE_NODE	*LandNode;

	Landscape = Rates->Landscape;
	LandNode = Landscape->NodeList[Landscape->NoNodes-1];

	MLPos = GetMaxLhPos(Lh, NoNodes);
	if(Lh[MLPos] < 1.98)
	{
		FreeLandScapeNode(LandNode);
		Landscape->NoNodes--;
		return FALSE;
	}

	LandNode->NodeList[0] = NodeList[MLPos];
	LandNode->Beta = Beta[MLPos];
	LandNode->MLEst = FALSE;

	return TRUE;
}



int	AddMLLandscapeNNode(OPTIONS *Opt, TREES *Trees, RATES *Rates)
{
	LANDSCAPE *Landscape;
	LANDSCAPE_NODE *LandNode;
	TREE *Tree;
	NODE	*NodeList, Node;
	int		NoNodes, Index, Valid;
	ML_MAP *ML_Map;
	double	Lh, InitLh;
	double	*Beta, *LhList;
	
	// Find the current Lh
	InitLh = Likelihood(Rates, Trees, Opt);
		
	Landscape = Rates->Landscape;
	Tree = Trees->Tree[Rates->TreeNo];

	// Get All valid nodes
	NodeList = MakeNodeList(Tree, Landscape, &NoNodes);

	// Space to save Beta and Lh
	Beta = (double*)SMalloc(sizeof(double) * NoNodes);
	LhList  = (double*)SMalloc(sizeof(double) * NoNodes);

	// Create a New land node and add it to the model
	LandNode  = CreatMLLandscapeNode(NULL);
	AddLandscapeNode(Landscape, LandNode);

	// Setup the ML map
	ML_Map = AllocMLMap();
	BuildMLMap(ML_Map, Opt, Trees, Rates);
	
	// For all the valid odes
	for(Index=0;Index<NoNodes;Index++)
	{
//		Index = 1;

		Node = 	NodeList[Index];
	
		LandNode->NodeList[0] = Node;
		LandNode->Beta = 0.0;
		
		LandNode->Beta = 1.0;
		Lh = Likelihood(Rates,Trees,Opt);

		// get the ML value
		NLOptBT(Rates, Opt, Trees, ML_Map);

		// map it back to the model 
		MLMapToRates(ML_Map, Opt, Rates);

		// calck the new lh
		Lh = Likelihood(Rates, Trees, Opt);

		// Save the Lh and Beta
		LhList[Index] = Lh - InitLh;
		Beta[Index] = LandNode->Beta;

		printf("%d\t%f\t%f\n", Index, LhList[Index], Beta[Index]);
	}


	Valid = SaveNewNode(Rates, Beta, LhList, NodeList, NoNodes);

	free(NodeList);
	free(LhList);
	free(Beta);
	FreeMLMap(ML_Map);

	return Valid;
}

void	SetLandRateEst(LANDSCAPE *Land, int Val)
{
	int Index;

	for(Index=0;Index<Land->NoNodes;Index++)
		Land->NodeList[Index]->MLEst = Val;
}

void	GlobalMLLandscapeOpt(OPTIONS *Opt, TREES *Trees, RATES *Rates)
{
	LANDSCAPE *Land;
	ML_MAP *ML_Map;
	double Lh;

	Land = Rates->Landscape;

	SetLandRateEst(Land, TRUE);

	ML_Map = AllocMLMap();
	BuildMLMap(ML_Map, Opt, Trees, Rates);

	NLOptBT(Rates, Opt, Trees, ML_Map);

	MLMapToRates(ML_Map, Opt, Rates);

	// calck the new lh
	Lh = Likelihood(Rates, Trees, Opt);

	FreeMLMap(ML_Map);

	SetLandRateEst(Land, FALSE);
}

void	PrintMLLandscape(TREES *Trees, RATES *Rates)
{
	LANDSCAPE *Landscape;
	LANDSCAPE_NODE *LandNode;
	int Index;

	Landscape = Rates->Landscape;

	printf("No Betas\t%d\n", Landscape->NoNodes);
	
	printf("Index\tBeta X BL\tBL\n");
	for(Index=0;Index<Landscape->NoNodes;Index++)
	{
		LandNode = Landscape->NodeList[Index];
		printf("%d\t%f\t%f\t", Index, LandNode->Beta, LandNode->NodeList[0]->Length);
		PrintPart(stdout, Trees, LandNode->NodeList[0]->Part);
		printf("\n");
	}

}

void	OutputLandscapeInfo(OPTIONS *Opt, TREES *Trees, RATES *Rates)
{
	double LH, Alpha, Sig;
	

	LH = Likelihood(Rates, Trees, Opt);

	Alpha = Rates->Contrast->Alpha[0];
	Sig = Rates->Contrast->Sigma[0];


	printf("Lh\tAlpha\tSig2\n");
	printf("%f\t%f\t%f\n", LH, Alpha, Sig);
}



void	MLLandscape(OPTIONS *Opt, TREES *Trees, RATES *Rates)
{
	int Valid;
	double InitLh, OptLh;

	OutputLandscapeInfo(Opt, Trees, Rates);

	printf("Start Lh\t1 D Opt\tGlobal Opt\n");
	do
	{
		InitLh = Likelihood(Rates, Trees, Opt);
		Valid = AddMLLandscapeNNode(Opt, Trees, Rates);
				
		if(Valid == TRUE)
		{
			OptLh = Likelihood(Rates, Trees, Opt);
			printf("%f\t%f\t", InitLh, OptLh);
			
			GlobalMLLandscapeOpt(Opt, Trees, Rates);
			
			OptLh = Likelihood(Rates, Trees, Opt);
			printf("%f\t", OptLh);

			printf("\n");
			fflush(stdout); 
		}


		fflush(stdout);
	}while(Valid == TRUE);
	
	PrintMLLandscape(Trees, Rates);
	
	OutputLandscapeInfo(Opt, Trees, Rates);

	exit(0);
}
*/



void	MLLandscape(OPTIONS *Opt,TREES *Trees,RATES *Rates)
{

//	LOCAL_TRANSFORM*	CreateLocalTransforms(char *Name,TAG **TagList,int NoTags,TRANSFORM_TYPE Type,int Est,double Scale);
//	void				FreeLocalTransforms(LOCAL_TRANSFORM* LTrans);

}
