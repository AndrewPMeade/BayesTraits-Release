#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

#include "TypeDef.h"
#include "Trees.h"
#include "GenLib.h"
#include "Part.h"

void	FreeLandscapeNode(LANDSCAPE_NODE* SNode);

LANDSCAPE*	CreateLandscape(OPTIONS *Opt)
{
	LANDSCAPE *Ret;

	Ret = (LANDSCAPE*)SMalloc(sizeof(LANDSCAPE));


	Ret->NoNodes = 0;
	Ret->NodeList = NULL;
	
	return Ret;
}

void		FreeLandscape(LANDSCAPE *Landscape)
{
	int Index;

	for(Index=0;Index<Landscape->NoNodes;Index++)
		FreeLandscapeNode(Landscape->NodeList[Index]);
	
	if(Landscape->NodeList != NULL)
		free(Landscape->NodeList);
	free(Landscape);
}

void	FreeLandscapeNode(LANDSCAPE_NODE* SNode)
{
	free(SNode->NodeList);
	free(SNode);
}

LANDSCAPE_NODE*	AllocLandscapNode(int NoTrees)
{
	LANDSCAPE_NODE*	Ret;

	Ret = (LANDSCAPE_NODE*)SMalloc(sizeof(LANDSCAPE_NODE));
	Ret->NodeList = (NODE*)SMalloc(sizeof(NODE) * NoTrees);

	return Ret;
}

LANDSCAPE_NODE*	CreateLandscapNode(PART *Part, int NoTrees)
{
	LANDSCAPE_NODE*	Ret;
	int Index;

	Ret = AllocLandscapNode(NoTrees);

	for(Index=0;Index<NoTrees;Index++)
		Ret->NodeList[Index] = NULL;

	Ret->Beta = 0.0;
	Ret->Part = Part;
	Ret->NodeID = 0;

	return Ret;
}

LANDSCAPE_NODE*	CloneLandscapNode(LANDSCAPE_NODE *LNode , int NoTrees)
{
	LANDSCAPE_NODE* Ret;

	Ret = AllocLandscapNode(NoTrees);

	Ret->Beta	= LNode->Beta;
	Ret->Part	= LNode->Part;
	Ret->NodeID = LNode->NodeID;

	memcpy(Ret->NodeList, LNode->NodeList, sizeof(NODE) * NoTrees);

	return Ret;
}

void	BlankLandscape(LANDSCAPE *Land)
{
	int Index;

	for(Index=0;Index<Land->NoNodes;Index++)
		FreeLandscapeNode(Land->NodeList[Index]);

	if(Land->NodeList != NULL)
		free(Land->NodeList);

	Land->NodeList = NULL;
	Land->NoNodes = 0;
}

void	CopyLandscape(LANDSCAPE *A, LANDSCAPE *B, int NoTrees)
{
	int Index;

	BlankLandscape(A);

	A->NodeList = (LANDSCAPE_NODE**)SMalloc(sizeof(LANDSCAPE_NODE*) * B->NoNodes);
	for(Index=0;Index<B->NoNodes;Index++)
		A->NodeList[Index] = CloneLandscapNode(B->NodeList[Index], NoTrees);

	A->NoNodes = B->NoNodes;
}

void		ResetTreeLandscape(TREE *Tree)
{
	int Index;
	for(Index=0;Index<Tree->NoNodes;Index++)
		Tree->NodeList[Index]->LandscapeBeta = 0.0;
}

NODE		GetLandscapeNode(TREES *Trees, int TreeNo, LANDSCAPE_NODE* LNode)
{
	if(LNode->NodeList[TreeNo] == NULL)
		LNode->NodeList[TreeNo] = PartGetMRCA(Trees->Tree[TreeNo], LNode->Part);

	return LNode->NodeList[TreeNo];
}

void		SetLandscapeBetaTrait(NODE Node, int NoSites, double Change)
{
	CONDATA *Con;
	int SIndex;

	Con = Node->ConData;

	for(SIndex=0;SIndex<NoSites;SIndex++)
		Con->Contrast[0]->Data[SIndex] = Node->Taxa->ConData[SIndex] + Change;
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

void		MapLandscape(OPTIONS *Opt, TREES *Trees, RATES *Rates)
{
	LANDSCAPE *Landscape;
	TREE *Tree;
	int TreeNo, Index, NIndex;
	NODE Node;
	
	
	TreeNo = Rates->TreeNo;
	Tree = Trees->Tree[TreeNo];
	Landscape = Rates->Landscape;

	ResetTreeLandscape(Tree);
	
	for(Index=0;Index<Landscape->NoNodes;Index++)
	{
		Node = GetLandscapeNode(Trees, TreeNo, Landscape->NodeList[Index]);

//		SetBetaNode(Node, Landscape->NodeList[Index]->Beta);

//		for(NIndex=0;NIndex<Node->NoNodes;NIndex++)
//			SetBetaNode(Node->NodeList[NIndex], Landscape->NodeList[Index]->Beta);

		Node->LandscapeBeta = Landscape->NodeList[Index]->Beta;
	}

	PropLandscapeBeta(Trees, Tree->Root, 0.0);
}

void		AddNodeToLandscape(LANDSCAPE *Land, LANDSCAPE_NODE *NewNode)
{
	LANDSCAPE_NODE	**NList;

	NList = (LANDSCAPE_NODE**)SMalloc(sizeof(LANDSCAPE_NODE*) * (Land->NoNodes + 1));

	memcpy(NList, Land->NodeList, sizeof(LANDSCAPE_NODE*) * Land->NoNodes);

	NList[Land->NoNodes] = NewNode;
	if(Land->NodeList != NULL)
		free(Land->NodeList);
	Land->NodeList = NList;
	Land->NoNodes++;
}

void		AddLandscapeFromPart(LANDSCAPE *Land, PART *Part, TREES *Trees, double Beta)
{
	LANDSCAPE_NODE *NewNode;

	NewNode = CreateLandscapNode(Part, Trees->NoTrees);

	NewNode->Beta = Beta;
	NewNode->NodeID = 0;
	NewNode->Part = Part;

	AddNodeToLandscape(Land, NewNode);
}