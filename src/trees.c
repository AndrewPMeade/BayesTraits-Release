#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "typedef.h"
#include "trees.h"
#include "genlib.h"
#include "data.h"
#include "likelihood.h"
#include "matrix.h"
#include "linalg.h"
#include "continuous.h"
#include "treepasser.h"
#include "BigLh.h"
#include "part.h"
#include "RandLib.h"
#include "QuadDouble.h"
#include "contrasts.h"
#include "FatTail.h"

#ifdef BTOCL
#include "btocl_discrete.h"
#endif

void	WriteTreeToFile(NODE n, NODE Root, FILE *F);

void	GetNoNodes(NODE N, int *No)
{
	int Index;

	(*No)++;
	for(Index=0;Index<N->NoNodes;Index++)
		GetNoNodes(N->NodeList[Index], No);
}

void	SetNodes(NODE *List, NODE N, int *Pos)
{
	int Index;

	List[*Pos] = N;
	(*Pos)++;

	for(Index=0;Index<N->NoNodes;Index++)
		SetNodes(List, N->NodeList[Index], Pos);
}

void	SetNodeList(TREE *Tree)
{
	if(Tree->NodeList != NULL)
		free(Tree->NodeList);

	Tree->NoNodes = 0;
	GetNoNodes(Tree->Root, &Tree->NoNodes);
	
	Tree->NodeList = (NODE*)malloc(sizeof(NODE) * Tree->NoNodes);
	if(Tree->NodeList == NULL)
		MallocErr();

	Tree->NoNodes = 0;
	SetNodes(Tree->NodeList, Tree->Root, &Tree->NoNodes);
}

NODE	AllocNode(void)
{
	NODE	Ret;

	Ret = (NODE)SMalloc(sizeof(struct INODE));

	Ret->ID				=	-1;
	Ret->Ans			=	NULL;
	Ret->Tip			=	FALSE;
	Ret->TipID			=	-1;
	Ret->Length			=	-1;
	Ret->DistToRoot		=	-1;
	Ret->Taxa			=	NULL;
	Ret->Partial		=	NULL;
	Ret->FossilMask		=	NULL;
	Ret->Part			=	NULL;
	Ret->GammaPartial	=	NULL;

	Ret->NodeList		=	NULL;
	Ret->NoNodes		=	-1;

	Ret->VPosX			=	-1;
	Ret->VPosY			=	-1;

	return Ret;
}

TAXA*	GetTaxaFromID(int ID, TAXA **Taxa, int NoOfTaxa)
{
	int	Index;

	for(Index=0;Index<NoOfTaxa;Index++)
	{
		if(Taxa[Index]->No == ID)
			return Taxa[Index];
	}
	return NULL;
}

TAXA*	GetTaxaFromName(char *Name, TAXA **Taxa, int NoOfTaxa)
{
	int	Index;

	for(Index=0;Index<NoOfTaxa;Index++)
	{
		if(strcmp(Taxa[Index]->Name, Name) == 0)
			return Taxa[Index];
	}
	return NULL;
}


void	LinkTipsToTaxa(NODE N, TAXA **Taxa, int NoOfTaxa)
{
	int NIndex;

	if(N->Tip == TRUE)
	{
		N->Taxa = GetTaxaFromID(N->TipID, Taxa, NoOfTaxa);
	}
	else
	{
		for(NIndex=0;NIndex<N->NoNodes;NIndex++)
			LinkTipsToTaxa(N->NodeList[NIndex], Taxa, NoOfTaxa);
	}

}

void	FreeNode(NODE N)
{
	if(N->Partial != NULL)
		free(N->Partial);

	if(N->FatTailNode != NULL)
		free(N->FatTailNode);

	if(N->GammaPartial != NULL)
		free(N->GammaPartial);

	if(N->NodeList != NULL)
		free(N->NodeList);

	if(N->FossilMask != NULL)
		free(N->FossilMask);

	free(N);
}

void	FreeTree(TREE* Tree, int NoOfSites, int NoOfTaxa)
{
	int		NIndex;
	NODE	N;

	N = Tree->Root;
	if(N->Partial != NULL)
		free(N->Partial[0]);

	if(N->GammaPartial != NULL)
		free(N->GammaPartial[0]);
	
	for(NIndex=0;NIndex<Tree->NoNodes;NIndex++)
		FreeNode(Tree->NodeList[NIndex]);

	free(Tree->NodeList);

	if(Tree->FatTailTree != NULL)
		FreeFatTailTree(Tree->FatTailTree);

	if(Tree->ConVars!= NULL)
		FreeConVar(Tree->ConVars, NoOfTaxa);

	for(NIndex=0;NIndex<Tree->NoFGroups;NIndex++)
		free(Tree->FNodes[NIndex]);
	free(Tree->FNodes);
	free(Tree->NoFNodes);

	if(Tree->PNodes != NULL)
		free(Tree->PNodes);

	free(Tree);
}

void	FreeTrees(TREES* Trees, OPTIONS *Opt)
{
	int	Index;

	if(Opt->ModelType == MT_CONTRAST)
		FreeAllContrast(Opt, Trees);

	FreeTreeBigLh(Opt, Trees);

#ifdef QUAD_DOUBLE
	if(Opt->ModelType == MT_DISCRETE)
		FreeQuadLh(Opt, Trees);
#endif

#ifdef BTOCL
	// Do this before individual trees are deallocated
	if(Opt->ModelType == MT_DISCRETE) {
		// not used, but keeping code
		//btocl_FreeLhInfo(Trees);
	}
#endif

	FreeData(Opt);
	free(Trees->Taxa);

	FreeParts(Trees);

	for(Index=0;Index<Trees->NoOfTrees;Index++)
		FreeTree(Trees->Tree[Index], Trees->NoOfSites, Trees->NoOfTaxa);
	free(Trees->Tree);

	if(Trees->PList != NULL)
	{
		FreeMultiMatrixLinMem(Trees->PList, Trees->MaxNodes);
		free(Trees->PList);
	}

	if(Trees->InvInfo != NULL)
		FreeInvInfo(Trees->InvInfo);

	if(Trees->RemovedTaxa != NULL)
	{
		for(Index=0;Index<Trees->NoOfRemovedTaxa;Index++)
			free(Trees->RemovedTaxa[Index]);
		free(Trees->RemovedTaxa);
	}

	if(Trees->SymbolList != NULL)
		free(Trees->SymbolList);

	if(Trees->TempConVars != NULL)
		FreeTempConVars(Trees->TempConVars);

	if(Trees->PMean != NULL)
		free(Trees->PMean);

	if(Trees->PSD != NULL)
		free(Trees->PSD);

	free(Trees);
}

void	SetTipHit(NTREES *Trees, NNODE N, char *THits)
{
	int	Index;

	for(Index=0;Index<Trees->NoOfTaxa;Index++)
		if(&Trees->Taxa[Index] == N->Taxa)
		{
			THits[Index] = TRUE;
			return;
		}

	return;
}

void	CheckTaxaPresent(NTREES *Trees)
{
	char	*THits;
	int		TIndex;
	int		Index;
	NTREE	*Tree;
	int		Err;

	THits = (char*)malloc(sizeof(char) * Trees->NoOfTaxa);
	if(THits== NULL)
		MallocErr();

	for(TIndex=0;TIndex<Trees->NoOfTrees;TIndex++)
	{
		Tree = &Trees->Trees[TIndex];

		for(Index=0;Index<Trees->NoOfTaxa;Index++)
			THits[Index] = FALSE;

		for(Index=0;Index<Tree->NoOfNodes;Index++)
		{
			if(Tree->NodeList[Index]->Tip == TRUE)
				SetTipHit(Trees, Tree->NodeList[Index], THits);
		}

		Err = FALSE;

		for(Index=0;Index<Trees->NoOfTaxa;Index++)
			if(THits[Index] == FALSE)
				Err = TRUE;

		if(Err == TRUE)
		{
			for(Index=0;Index<Trees->NoOfTaxa;Index++)
				if(THits[Index] == FALSE)
					printf("Error: Taxa %s is not present in tree %d\n",  Trees->Taxa[Index].Name, TIndex+1);
		
			exit(0);
		}		
	}

	free(THits);
}

void	CheckBaseTirs(NTREES *Trees)
{
	int		TIndex;
	NTREE	*Tree;

	for(TIndex=0;TIndex<Trees->NoOfTrees;TIndex++)
	{
		Tree = &Trees->Trees[TIndex];

		if(Tree->Root->NoOfNodes != 2)
		{
			printf("Error: Tree %d has a %d-way basal polytomy, only rooted binary trees are allowed.\n", TIndex+1, Tree->Root->NoOfNodes);                                
			exit(0);
		}
	}
}

void	CheckPoly(NTREES *Trees)
{
	int		TIndex;
	int		NIndex;
	int		Index;
	NTREE	*Tree;
	NNODE	Node;

	for(TIndex=0;TIndex<Trees->NoOfTrees;TIndex++)
	{
		Tree = &Trees->Trees[TIndex];
		for(NIndex=0;NIndex<Tree->NoOfNodes;NIndex++)
		{
			Node = Tree->NodeList[NIndex];
			if((Node->Tip == FALSE) && (Node->NoOfNodes != 2))
			{
				printf("Error: Tree %d has a %d-way polymotme.\n", TIndex+1, Node->NoOfNodes);
				for(Index=0;Index<Node->NoOfNodes;Index++)
				{
					printf("\tNode %d\t", Index+1);
					PrintNTree(stdout, Node->NodeList[Index]);
					printf("\n");
				}
				exit(0);
			}
		}
	}
}

void	CheckPresentBL(NTREES *Trees)
{
	int		TIndex;
	int		NIndex;
	NNODE	N;
	NTREE	*Tree;

	for(TIndex=0;TIndex<Trees->NoOfTrees;TIndex++)
	{
		Tree = &Trees->Trees[TIndex];

		for(NIndex=0;NIndex<Tree->NoOfNodes;NIndex++)
		{
			N = Tree->NodeList[NIndex];

			if((N->Length == -1) && (N != Tree->Root))
			{
				printf("Error: Node in tree %d does not have a valid branch length.\n", TIndex++);
				printf("Tree:\t");
				PrintNTree(stdout, N);
				printf("\n");
				exit(0);
			}
		}
	}
}

void	ValidateTrees(NTREES *Trees)
{
	/* Check each tree has all taxa */
	CheckTaxaPresent(Trees);

	/* No basis trification */
/*	CheckBaseTirs(Trees); */
	
	/* Only binary trees */
/*	CheckPoly(Trees); */

	/* More than one tree */
	if(Trees->NoOfTrees < 1)
	{
		printf("Error: BayesTraits requires one or more valid trees\n");
		exit(0);
	}

	/* More than one Taxa */
	if(Trees->NoOfTaxa < 1)
	{
		printf("Error: BayesTraits requires one or more valid taxa\n");
		exit(0);
	}
	
	/* All nodes barch lengths */
	CheckPresentBL(Trees);
}

int		GetArrPos(NNODE *Start, NNODE Pos)
{
	int	Index;

	Index = 0;
	while(1)
	{
		if(Start[Index] == Pos)
			return Index;
		Index++;
	}

	return -1;
}

void	MakeNewTree(TREE *Tree, NTREE *PTree)
{
	int		Index, NIndex;
	NODE	Node;
	NNODE	NNode;
	int		Pos;

	Tree->NoNodes = PTree->NoOfNodes;

	Tree->NodeList = (NODE*)malloc(sizeof(NODE) * Tree->NoNodes);
	if(Tree->NodeList == NULL)
		MallocErr();

	for(Index=0;Index<PTree->NoOfNodes;Index++)
		Tree->NodeList[Index] = AllocNode();

	for(Index=0;Index<PTree->NoOfNodes;Index++)
	{
		Node = Tree->NodeList[Index];
		NNode= PTree->NodeList[Index];

		Node->Length		= NNode->Length;
		Node->UserLength	= NNode->Length;
		Node->Tip			= NNode->Tip;
		Node->ConData		= NULL;
		Node->FatTailNode	= NULL;
	

		Node->Part = NULL;
	//	Node->PSize= 0;

		if(Node->Tip == TRUE)
		{
			Node->TipID = NNode->TaxaID;
		}
		else
		{
			Node->NodeList = (NODE*)malloc(sizeof(NODE) * NNode->NoOfNodes);
			if(Node->NodeList == NULL)
				MallocErr();
			Node->NoNodes = NNode->NoOfNodes;

			for(NIndex=0;NIndex<Node->NoNodes;NIndex++)
			{
				Pos = GetArrPos(&PTree->NodeList[0], NNode->NodeList[NIndex]);
				Node->NodeList[NIndex] = Tree->NodeList[Pos];
			} 
		}

		if(NNode->Ans == NULL)
			Node->Ans = NULL;
		else
		{
			Pos = GetArrPos(&PTree->NodeList[0], NNode->Ans);
			Node->Ans = Tree->NodeList[Pos];
		}
	}

	Pos = GetArrPos(&PTree->NodeList[0], PTree->Root);
	Tree->Root = Tree->NodeList[Pos];
}

int		FindMaxPoly(TREES *Trees)
{
	int TIndex,NIndex;
	int	Ret;
	TREE *Tree;
	NODE Node;

	Ret = 0;

	for(TIndex=0;TIndex<Trees->NoOfTrees;TIndex++)
	{
		Tree = Trees->Tree[TIndex];
		for(NIndex=0;NIndex<Tree->NoNodes;NIndex++)
		{
			Node = Tree->NodeList[NIndex];
			if(Node->NoNodes > Ret)
				Ret = Node->NoNodes;
		}
	}

	return Ret;
}

 int	FindMaxNodes(TREES *Trees)
 {
	int TIndex;
	int	Ret;
	TREE *Tree;
	
	Ret = 0;


	for(TIndex=0;TIndex<Trees->NoOfTrees;TIndex++)
	{
		Tree = Trees->Tree[TIndex];
		if(Tree->NoNodes > Ret)
			Ret = Tree->NoNodes;
	}

	return Ret;
 }


 TAXA*	InitTaxa(int No, char *Name)
 {
	TAXA *Ret;

	Ret = (TAXA*)malloc(sizeof(TAXA));
	if(Ret == NULL)
		MallocErr();

	Ret->No				=	No;
	Ret->Name			=	StrMake(Name);

	Ret->DesDataChar	=	NULL;
	Ret->ConData		=	NULL;
	Ret->Exclude		=	FALSE;

	Ret->EstData		=	FALSE;
	Ret->EstDepData		=	FALSE;
	Ret->EstDataP		=	NULL;
	Ret->RealData		=	NULL;
	
	 return Ret;
 }


 double	CaclAveBL(TREE *Tree)
 {
	 double Ret;
	 int Index;

	 Ret = 0.0;
	 for(Index=1;Index<Tree->NoNodes;Index++)
		 Ret += Tree->NodeList[Index]->Length;
	 
	 Ret = Ret / (Tree->NoNodes - 1);
	 return Ret;
 }

 TREE*	InitTree(TREES *Trees, NTREES *PTrees, int TNo)
 {
	TREE	*Ret;
	
	Ret = (TREE*)malloc(sizeof(TREE));
	if(Ret == NULL)
		MallocErr();

	Ret->ConVars	= NULL;
	Ret->NodeList	= NULL;
	Ret->Root		= NULL;

	Ret->FNodes		= NULL;
	Ret->NoFNodes	= NULL;
	Ret->NoFGroups	= -1;

	Ret->NoPNodes	= -1;
	Ret->PNodes		= NULL;

	Ret->FatTailTree= NULL;

	Ret->NoContrast	= -1;
	
	MakeNewTree(Ret, &PTrees->Trees[TNo]);

	LinkTipsToTaxa(Ret->Root, Trees->Taxa, Trees->NoOfTaxa);
		
	SetNodeIDs(Ret);

	Ret->AveBL = CaclAveBL(Ret);

//	printf("AveBL\t%f\n", Ret->AveBL);

	return Ret;
 }


 void	ReSetTreeTaxaID(TREE *Tree)
 {
	int Index;
	NODE N;

	for(Index=0;Index<Tree->NoNodes;Index++)
	{
		N = Tree->NodeList[Index];
		if(N->Tip == TRUE)
			N->TipID = N->Taxa->No;
	}
 }

 void	ReSetTaxaID(TREES *Trees)
 {
	int Index;
	
	for(Index=0;Index<Trees->NoOfTaxa;Index++)
		Trees->Taxa[Index]->No = Index;
	
	for(Index=0;Index<Trees->NoOfTrees;Index++)
		ReSetTreeTaxaID(Trees->Tree[Index]);
 }

void	InitialTrees(TREES *Trees, NTREES *PTrees)
{
	int		Index;

	Trees->NoOfTaxa		= PTrees->NoOfTaxa;
	Trees->NoOfTrees	= PTrees->NoOfTrees;
	Trees->NOSPerSite	= FALSE;
	
	Trees->Taxa = (TAXA**) malloc(sizeof(TAXA*) * Trees->NoOfTaxa);
	if(Trees->Taxa == NULL)
		MallocErr();
	
	for(Index=0;Index<Trees->NoOfTaxa;Index++)
		Trees->Taxa[Index] = InitTaxa(PTrees->Taxa[Index].No, PTrees->Taxa[Index].Name);

	Trees->Tree = (TREE**)malloc(sizeof(TREE*) * PTrees->NoOfTrees);
	if(Trees->Tree == NULL)
		MallocErr();

	for(Index=0;Index<Trees->NoOfTrees;Index++)
		Trees->Tree[Index] = InitTree(Trees, PTrees, Index);
	
	Trees->MaxNodes = FindMaxNodes(Trees);

	ReSetTaxaID(Trees);

	SetParts(Trees);
}

void	TestTreeLoad(TREES *Trees)
{
	int Index;
	TAXA *Taxa;
	TREE *T;
	NODE N;

	for(Index=0;Index<Trees->NoOfTaxa;Index++)
	{
		Taxa = Trees->Taxa[Index];
		printf("%d\t%s\n", Taxa->No, Taxa->Name);
	}

	T = Trees->Tree[0];
	for(Index=0;Index<T->NoNodes;Index++)
	{
		N = T->NodeList[Index];
		if(N->Tip == TRUE)
		{
			printf("%d\t%d\t%s\n", N->TipID, N->Taxa->No, N->Taxa->Name);
		}
	}

	exit(0);
}

TREES*	LoadTrees(char* FileName)
{
	TREES*		Ret=NULL;
	NTREES*		PTrees;
	char		*Err;

	Ret = (TREES*)malloc(sizeof(TREES));
	if(Ret==NULL)
		MallocErr();

	Ret->MaxNodes			= -1;
	Ret->PMem				= NULL;
	Ret->PList				= NULL;

	Ret->Taxa				= NULL;
	Ret->Tree				= NULL;
	Ret->InvInfo			= NULL;
	Ret->SymbolList			= NULL;
	Ret->RemovedTaxa		= NULL;
	Ret->NoOfRemovedTaxa	= 0;
	Ret->ValidCData			= TRUE;
	Ret->ValidDData			= TRUE;

	Ret->PMean				= NULL;
	Ret->PSD				= NULL;
	Ret->TempConVars		= NULL;

	PTrees = LoadNTrees(FileName, &Err);

	if(PTrees == NULL)
	{
		printf("Err: %s\n", Err);
		free(Err);
		exit(0);
	}

	ValidateTrees(PTrees);

	InitialTrees(Ret, PTrees);

	FreeNTrees(PTrees);
	
	Ret->NormConst = FindTreeNormalise(Ret);

	Ret->JStop	=	FALSE;

//	TestTreeLoad(Ret);


	return Ret;
}

void	PrintTreesInfo(FILE*	Str, TREES *Trees, DATATYPE DataType)
{
	fprintf(Str, "Tree Information\n");
	fprintf(Str, "     Trees:                      %d\n", Trees->NoOfTrees);
	fprintf(Str, "     Taxa:                       %d\n", Trees->NoOfTaxa);
	fprintf(Str, "     Sites:                      %d\n", Trees->NoOfSites);

	if(DataType == DISCRETE)
		fprintf(Str, "     States:                     %d\n", Trees->NoOfStates);
}

int		SymbolToPos(char Symbol, char *List)
{
	int	Index;

	for(Index=0;Index<(int)strlen(List);Index++)
	{
		if(List[Index] == Symbol)
			return Index;
	}

	return -1;
}

int		SiteHadUnKnownState(char *StatList)
{
	int	Index;

	for(Index=0;Index<(int)strlen(StatList);Index++)
	{
		if(StatList[Index] == UNKNOWNSTATE)
			return TRUE;
	}

	return FALSE;
}

void	SetDescCVTipData(int SiteNo, NODE N)
{
	N->Partial[SiteNo][4] = N->Partial[SiteNo][0];
	N->Partial[SiteNo][5] = N->Partial[SiteNo][1];
	N->Partial[SiteNo][6] = N->Partial[SiteNo][2];
	N->Partial[SiteNo][7] = N->Partial[SiteNo][3];
}

void	SetNodeTipData(OPTIONS *Opt, NODE N, TREE* Tree, TREES *Trees)
{
	int		SiteIndex;
	int		StateIndex;
	int		Pos;
	int		StrLen;
	int		NOS;

	NOS = Trees->NoOfStates;
	if(Trees->UseCovarion == TRUE)
		NOS = NOS / 2;

	for(SiteIndex=0;SiteIndex<Trees->NoOfSites;SiteIndex++)
	{
		if(SiteHadUnKnownState(N->Taxa->DesDataChar[SiteIndex]) == FALSE)
		{
			/* Set all boxes to 0  */
			for(StateIndex=0;StateIndex<NOS;StateIndex++)
				N->Partial[SiteIndex][StateIndex] = 0;

			/* Set the corrispoding boxes to 1  */
			StrLen = (int)strlen(N->Taxa->DesDataChar[SiteIndex]);
			for(StateIndex=0;StateIndex<StrLen;StateIndex++)
			{
				Pos = SymbolToPos(N->Taxa->DesDataChar[SiteIndex][StateIndex], Trees->SymbolList);
				N->Partial[SiteIndex][Pos] = 1;
			}
		}
		else
		{
			for(StateIndex=0;StateIndex<NOS;StateIndex++)
				N->Partial[SiteIndex][StateIndex] = 1;
		}

		/* Copy for the Dep CV model */
		if(Opt->Model == M_DESCCV)
			SetDescCVTipData(SiteIndex, N);

		/* Copy the sites for the covarion mode */
		if(Trees->UseCovarion == TRUE)
		{
			for(StateIndex=0;StateIndex<NOS;StateIndex++)
				N->Partial[SiteIndex][StateIndex+NOS] = N->Partial[SiteIndex][StateIndex];
		}
	}
}
void	SetTipData(OPTIONS *Opt, TREE *Tree, TREES *Trees)
{
	int		NIndex;
	NODE	N;
	
	for(NIndex=0;NIndex<Tree->NoNodes;NIndex++)
	{
		N = Tree->NodeList[NIndex];
		if(N->Tip == TRUE)
			SetNodeTipData(Opt, N, Tree, Trees);
	}
}

/*
void	SetTipData(TREE *Tree, TREES *Trees)
{
	int		NIndex;
	int		SiteIndex;
	int		StateIndex;
	int		Pos;
	NODE	N;
	int		StrLen;
	int		NOS;

	NOS = Trees->NoOfStates;
	if(Trees->UseCovarion == TRUE)
		NOS = NOS / 2;

	for(NIndex=0;NIndex<Trees->NoOfNodes;NIndex++)
	{
		N = &Tree->NodeList[NIndex];
		if(N->Tip == TRUE)
		{
			for(SiteIndex=0;SiteIndex<Trees->NoOfSites;SiteIndex++)
			{
				if(SiteHadUnKnownState(N->Taxa->DesDataChar[SiteIndex]) == FALSE)
				{
					
					for(StateIndex=0;StateIndex<NOS;StateIndex++)
						N->Partial[SiteIndex][StateIndex] = 0;
					
					StrLen = (int)strlen(N->Taxa->DesDataChar[SiteIndex]);

					for(StateIndex=0;StateIndex<StrLen;StateIndex++)
					{
						Pos = SymbolToPos(N->Taxa->DesDataChar[SiteIndex][StateIndex], Trees->SymbolList);
						N->Partial[SiteIndex][Pos] = 1;
					}
				}
				else
				{
					for(StateIndex=0;StateIndex<NOS;StateIndex++)
						N->Partial[SiteIndex][StateIndex] = 1;
				}

				
				if(Trees->UseCovarion == TRUE)
				{
					for(StateIndex=0;StateIndex<NOS;StateIndex++)
						N->Partial[SiteIndex][StateIndex+NOS] = N->Partial[SiteIndex][StateIndex];
				}
			}
		}
	}
}
*/
void	SetNOSPerSiteTipData(TREE *Tree, TREES *Trees)
{
	int		NIndex;
	int		SiteIndex;
	int		StateIndex;
	int		Pos;
	NODE	N;
	int		StrLen;
	int		NOS;
	double	**Part;

	for(NIndex=0;NIndex<Tree->NoNodes;NIndex++)
	{
		N = Tree->NodeList[NIndex];
		if(N->Tip == TRUE)
		{
			Part = N->Partial;
			for(SiteIndex=0;SiteIndex<Trees->NoOfSites;SiteIndex++)
			{
				NOS = Trees->NOSList[SiteIndex];
				if(Trees->UseCovarion == TRUE)
					NOS = NOS / 2;

				if(SiteHadUnKnownState(N->Taxa->DesDataChar[SiteIndex]) == FALSE)
				{
					/* Set all boxes to 0 */
					for(StateIndex=0;StateIndex<NOS;StateIndex++)
						Part[SiteIndex][StateIndex] = 0;
					/* Set the corrispoding boxes to 1 */
					StrLen = (int)strlen(N->Taxa->DesDataChar[SiteIndex]);
					for(StateIndex=0;StateIndex<StrLen;StateIndex++)
					{
						Pos = SymbolToPos(N->Taxa->DesDataChar[SiteIndex][StateIndex], Trees->SiteSymbols[SiteIndex]);
						Part[SiteIndex][Pos] = 1;
					}
				}
				else
				{
					for(StateIndex=0;StateIndex<NOS;StateIndex++)
						Part[SiteIndex][StateIndex] = 1;
				}

				/* Copy the sites for the covarion mode */
				if(Trees->UseCovarion == TRUE)
				{
					for(StateIndex=0;StateIndex<NOS;StateIndex++)
						Part[SiteIndex][StateIndex+NOS] = Part[SiteIndex][StateIndex];
				}
			}
		}
	}
}

/*
void	AllocNodePartial(NODE N, TREES *Trees, int Gamma)
{
	int		Outter;

	N->GammaPartial = NULL;

	N->Partial = (double**)malloc(sizeof(double*)*Trees->NoOfSites);
	if(N->Partial == NULL)
		MallocErr();

	for(Outter=0;Outter<Trees->NoOfSites;Outter++)
	{
		N->Partial[Outter] = (double*)malloc(sizeof(double) * Trees->NoOfStates);
		if(N->Partial[Outter] == NULL)
			MallocErr();
	}

	if(Gamma == FALSE)
		return;

	if(N->Tip == TRUE)
		return;
	
	N->GammaPartial = (double**)malloc(sizeof(double*)*Trees->NoOfSites);
	if(N->GammaPartial == NULL)
		MallocErr();

	for(Outter=0;Outter<Trees->NoOfSites;Outter++)
	{
		N->GammaPartial[Outter] = (double*)malloc(sizeof(double) * Trees->NoOfStates);
		if(N->GammaPartial[Outter] == NULL)
			MallocErr();
	}
}

void	AllocPartial(OPTIONS *Opt, TREES* Trees, int Gamma)
{
	int		TIndex;
	int		NIndex;
	NODE	N;

	
	for(TIndex=0;TIndex<Trees->NoOfTrees;TIndex++)
	{
		for(NIndex=0;NIndex<Trees->Tree[TIndex]->NoNodes;NIndex++)
		{
			N = Trees->Tree[TIndex]->NodeList[NIndex];
			AllocNodePartial(N, Trees, Gamma);
		}

		SetTipData(Opt, Trees->Tree[TIndex], Trees);
	}
}
*/

void	AllocNodePartial(NODE N, TREES *Trees, int Gamma)
{
	int		Outter;

	N->GammaPartial = NULL;

	N->Partial = (double**)malloc(sizeof(double*)*Trees->NoOfSites);
	if(N->Partial == NULL)
		MallocErr();

	for(Outter=0;Outter<Trees->NoOfSites;Outter++)
	{
		N->Partial[Outter] = (double*)malloc(sizeof(double) * Trees->NoOfStates);
		if(N->Partial[Outter] == NULL)
			MallocErr();
	}

	if(Gamma == FALSE)
		return;

	if(N->Tip == TRUE)
		return;
	
	N->GammaPartial = (double**)malloc(sizeof(double*)*Trees->NoOfSites);
	if(N->GammaPartial == NULL)
		MallocErr();

	for(Outter=0;Outter<Trees->NoOfSites;Outter++)
	{
		N->GammaPartial[Outter] = (double*)malloc(sizeof(double) * Trees->NoOfStates);
		if(N->GammaPartial[Outter] == NULL)
			MallocErr();
	}
}

void	AllocTreePartail(TREES* Trees, TREE *Tree)
{
	double *Mem;
	NODE	N;
	int		Index, SIndex;
		
	Mem = (double*)malloc(sizeof(double) * Trees->NoOfStates * Trees->NoOfSites * Tree->NoNodes);
	if(Mem == NULL)
		MallocErr();

	for(Index=0;Index<Tree->NoNodes;Index++)
	{
		N = Tree->NodeList[Index];

		N->Partial = (double**)malloc(sizeof(double*) * Trees->NoOfSites);
		if(N->Partial == NULL)
			MallocErr();

		N->Partial[0] = &Mem[Index * Trees->NoOfSites * Trees->NoOfStates];
		for(SIndex=1;SIndex<Trees->NoOfSites;SIndex++)
			N->Partial[SIndex] = N->Partial[0] + (SIndex * Trees->NoOfStates);
	}
}

void	SetGammaPartials(TREES* Trees, TREE *Tree)
{
	int Index;
	NODE N;

	for(Index=0;Index<Tree->NoNodes;Index++)
	{
		N = Tree->NodeList[Index];
		N->GammaPartial = N->Partial;
	}
}

void	AllocPartial(OPTIONS *Opt, TREES* Trees, int Gamma)
{
	int		TIndex;
	
	for(TIndex=0;TIndex<Trees->NoOfTrees;TIndex++)
	{
		AllocTreePartail(Trees, Trees->Tree[TIndex]);

		if(Gamma == TRUE)
		{
			SetGammaPartials(Trees, Trees->Tree[TIndex]);
			AllocTreePartail(Trees, Trees->Tree[TIndex]);
		}

		SetTipData(Opt, Trees->Tree[TIndex], Trees);
	}
}


double	GetStateProbPct(int State, int NoOfStates, double *Part)
{
	double	Tot=0;
	int		Index;

	for(Index=0;Index<NoOfStates;Index++)
		Tot +=	Part[Index];
	 
	return Part[State] / Tot;
}

void	PrintTaxaNames(NODE N)
{
	int	Index;

	if(N->Tip == TRUE)
		printf("%s\t%f\n", N->Taxa->Name, N->Length);
	else
	{
		for(Index=0;Index<N->NoNodes;Index++)
			PrintTaxaNames(N->NodeList[Index]);
	}
}

void	ReBuildTree(NODE N, NODE NewNodeList)
{
	NODE	NewNode;
	int		NIndex;

	if(N->ID == -1)
	{
		if(N->Taxa == FALSE)
			ReBuildTree(N->NodeList[0], NewNodeList);
		return;
	}

	NewNode = &NewNodeList[N->ID];

	NewNode->Tip		=	N->Tip;
	NewNode->TipID		=	N->TipID;
	NewNode->Length		=	N->Length;
	NewNode->Visited	=	N->Visited;
	NewNode->Taxa		=	NULL;
	NewNode->Partial	=	N->Partial;

	if(N->Ans != NULL)
		NewNode->Ans		=	&NewNodeList[N->Ans->ID];
	else
		NewNode->Ans		=	NULL;

	if(N->Tip == FALSE)
	{
		NewNode->NodeList = (NODE*)malloc(sizeof(NODE) * N->NoNodes);
		if(NewNode->NodeList == NULL)
			MallocErr();
		NewNode->NoNodes = N->NoNodes;

		for(NIndex=0;NIndex<N->NoNodes;NIndex++)
			NewNode->NodeList[NIndex] = &NewNodeList[N->NodeList[NIndex]->ID];

		for(NIndex=0;NIndex<N->NoNodes;NIndex++)
			ReBuildTree(N->NodeList[NIndex], NewNodeList);
	}
}


void	SetDelNodeIDs(NODE N, int *ID)
{
	int Index;

	if(N->ID == -1)
	{
		N->ID = *ID;
		*ID = (*ID) + 1;
	}

	if(N->Tip == FALSE)
	{
		for(Index=0;Index<N->NoNodes;Index++)
			SetDelNodeIDs(N->NodeList[Index], ID);
	}
}

void	SwapNode(NODE Node, NODE With)
{
	int NIndex;
	NODE Ans;

	Ans = Node->Ans;
	for(NIndex=0;NIndex<Ans->NoNodes;NIndex++)
		if(Ans->NodeList[NIndex] == Node)
			Ans->NodeList[NIndex] = With;
}




int	GetTaxaPos(TREES *Trees, char *TName)
{
	int TIndex;

	for(TIndex=0;TIndex<Trees->NoOfTaxa;TIndex++)
	{
		if(strcmp(TName, Trees->Taxa[TIndex]->Name)==0)
			return TIndex;
	}

	return -1;
}

NODE	GetTaxaNode(TREE *Tree, char *TName)
{
	int Index;
	NODE N;

	for(Index=0;Index<Tree->NoNodes;Index++)
	{
		N = Tree->NodeList[Index];
		if(N->Tip == TRUE)
		{
			if(strcmp(N->Taxa->Name, TName) == 0)
				return N;
		}
	}

	return NULL;
}


void	RemoveTipNode(NODE TipNode)
{
	int		NIndex, Pos;
	NODE	*NList, Ans;


	Ans = TipNode->Ans;
	NList = (NODE*)malloc(sizeof(NODE) * (Ans->NoNodes - 1));
	if(NList == NULL)
		MallocErr();
	
	Pos = 0;
	for(NIndex=0;NIndex<Ans->NoNodes;NIndex++)
		if(Ans->NodeList[NIndex] != TipNode)
			NList[Pos++] = Ans->NodeList[NIndex];

	free(Ans->NodeList);
	Ans->NodeList = NList;

	Ans->NoNodes--;
}

void	AddOrphanNode(NODE Node)
{
	NODE Des, Ans;
	int	Index;
		
	Ans = Node->Ans;

	assert(Ans != NULL);

	Des = Node->NodeList[0];
	Des->Length += Node->Length;

	for(Index=0;Index<Ans->NoNodes;Index++)
		if(Ans->NodeList[Index] == Node)
			Ans->NodeList[Index] = Des;

	Des->Ans = Ans;
}

void	ReSetOrphanNodes(TREE *Tree, NODE Node)
{
	int Index;
	
	if(Node->Tip == TRUE)
		return;

	for(Index=0;Index<Node->NoNodes;Index++)
		ReSetOrphanNodes(Tree, Node->NodeList[Index]);

	if(Node->NoNodes > 1)
		return;

	assert(Node->NoNodes == 1);
		
	if(Node->Ans != NULL)
		AddOrphanNode(Node);
	else
	{
		Tree->Root = Node->NodeList[0];
		Tree->Root->Ans = NULL;
	}

	FreeNode(Node);
}

void	RemoveTaxaFromTree(TREES *Trees, TREE *Tree, char *TName)
{
	NODE	DNode, DNodeAns;

	DNode = GetTaxaNode(Tree, TName);
	if(DNode == NULL)
	{
		printf("%s::%d Cound not find taxa %s\n", __FILE__, __LINE__, TName);
		exit(1);
	}
	
	DNodeAns = DNode->Ans;

	RemoveTipNode(DNode);
	FreeNode(DNode);

	ReSetOrphanNodes(Tree, Tree->Root);

	SetNodeList(Tree);
}

void	RemoveTaxaRec(TREES *Trees, int TaxaPos)
{
	TAXA	**NewTaxaList;
	char	**NewUsedTaxa;
	char	*TaxaName;
	int		NewTIndex, OldTIndex, TIndex;

	TaxaName = StrMake(Trees->Taxa[TaxaPos]->Name);

	NewTaxaList = (TAXA**)malloc(sizeof(TAXA*) * (Trees->NoOfTaxa - 1));
	if(NewTaxaList == NULL)
		MallocErr();

	NewTIndex = 0;
	for(OldTIndex=0;OldTIndex<Trees->NoOfTaxa;OldTIndex++)
	{
		if(OldTIndex != TaxaPos)
		{
			NewTaxaList[NewTIndex] = Trees->Taxa[OldTIndex];
			NewTIndex++;
		}
	}

	FreeTaxa(Trees->Taxa[TaxaPos], Trees->NoOfSites);
	
	free(Trees->Taxa);
	Trees->Taxa	= NewTaxaList;
	Trees->NoOfTaxa--;

	for(TIndex=0;TIndex<Trees->NoOfTrees;TIndex++)
		LinkTipsToTaxa(Trees->Tree[TIndex]->Root, Trees->Taxa, Trees->NoOfTaxa);
	
	NewUsedTaxa = (char**)malloc(sizeof(char*) * (Trees->NoOfRemovedTaxa + 1));
	if(NewUsedTaxa == NULL)
		MallocErr();

	for(TIndex=0;TIndex<Trees->NoOfRemovedTaxa;TIndex++)
		NewUsedTaxa[TIndex] = Trees->RemovedTaxa[TIndex];

	NewUsedTaxa[TIndex] = TaxaName;

	free(Trees->RemovedTaxa);
	Trees->RemovedTaxa = NewUsedTaxa;
	Trees->NoOfRemovedTaxa++;
}

void	CheckDelTaxa(OPTIONS *Opt, TREES *Trees, char *TName)
{
	int TaxaPos;

	if(Trees->NoOfTaxa <= 2)
	{
		printf("There must be two or more taxa in a tree\n");
		exit(0);
	}

	TaxaPos = GetTaxaPos(Trees, TName);

	if(TaxaPos == -1)
	{
		printf("Could not find taxa %s for deleting\n", TName);
		exit(0);
	}
	
	if(Opt == NULL)
		return;

	if(Opt->NoOfRecNodes != 0)
	{
		printf("There must be no node to reconstruct before a taxa can be removed for a tree \n");
		exit(0);
	}


	if(Opt->NoLocalTransforms != 0)
	{
		printf("There must be no local transforms before a taxa can be removed.");
		exit(0);
	}
}

void	RemoveTaxa(TREES *Trees, char *TName)
{
	int TaxaPos, TIndex;

	TaxaPos = GetTaxaPos(Trees, TName);

	for(TIndex=0;TIndex<Trees->NoOfTrees;TIndex++)
		RemoveTaxaFromTree(Trees, Trees->Tree[TIndex], TName);

	RemoveTaxaRec(Trees, TaxaPos);
		
	ReSetTaxaID(Trees);
}

void WriteTreeToFile(NODE n, NODE Root, FILE *F)
{
	double	BL;
	int		Index;

	BL = n->Length;

	if(n->Tip == TRUE)
		fprintf(F, "%d:%10.10f", n->Taxa->No+1, BL);
	else
	{
		fprintf(F, "(");

		for(Index=0;Index<n->NoNodes-1;Index++)
		{
			WriteTreeToFile(n->NodeList[Index], Root, F);
			fprintf(F, ",");
		}
		WriteTreeToFile(n->NodeList[Index], Root, F);
 
		
		if(n != Root)
			fprintf(F, "):%10.10f", BL);	
		else
			fprintf(F, ")");	
	}
}

void	SaveTreesHeader(FILE *Out, TREES* Trees)
{
	int TIndex;

	fprintf(Out, "#NEXUS\n");
	fprintf(Out, "begin trees;\n");
	fprintf(Out, "\ttranslate\n");

	for(TIndex=0;TIndex<Trees->NoOfTaxa;TIndex++)
	{
		fprintf(Out, "\t\t%d %s", Trees->Taxa[TIndex]->No+1, Trees->Taxa[TIndex]->Name);

		if(TIndex<Trees->NoOfTaxa-1)
			fprintf(Out, ",\n");
		else
			fprintf(Out, ";\n");
	}
}

void	SaveTreeFotter(FILE *Out)
{
	fprintf(Out, "end;");
	fclose(Out);
}

void	SaveTrees(char	*FileName, TREES* Trees)
{

	FILE	*TreeFile=NULL;
	int		TIndex;
	char	*Name;
	char	Buffer[128];
	size_t	NoOfChar;


	TreeFile = OpenWrite(FileName);
	
	SaveTreesHeader(TreeFile, Trees);

	sprintf(&Buffer[0], "%d", Trees->NoOfTrees);
	NoOfChar = strlen(&Buffer[0]);

	for(TIndex=0;TIndex<Trees->NoOfTrees;TIndex++)
	{
		Name = FormatInt(TIndex+1, (int)NoOfChar);

		fprintf(TreeFile, "\t\ttree No_%s = ", Name);
		WriteTreeToFile(Trees->Tree[TIndex]->Root, Trees->Tree[TIndex]->Root, TreeFile);
		fprintf(TreeFile, ";\n");
		
		free(Name);
	}

	SaveTreeFotter(TreeFile);
}	

void	SetNodeIDs(TREE *Tree)
{
	int Index;

	for(Index=0;Index<Tree->NoNodes;Index++)
		Tree->NodeList[Index]->ID = Index;
}

void SetMinBL(TREES *Trees)
{
	int		TIndex;
	int		NIndex;
	NODE	Node;
	TREE	*Tree;
	int		NoErr;

	NoErr=0;
	for(TIndex=0;TIndex<Trees->NoOfTrees;TIndex++)
	{
		Tree = Trees->Tree[TIndex];
		for(NIndex=0;NIndex<Tree->NoNodes;NIndex++)
		{
			Node = Tree->NodeList[NIndex];
			if(Node != Tree->Root)
			{
				if((Node->Length == -1) || (Node->Length < MINBL))
				{
					if(Node->Length == -1)
					{
						printf("Error: Tree %d has a node without a branch length.\n", TIndex);
						exit(0);
					}

					if(NoErr < 5)
						printf("Tree %d has a Branch length less then %f.\nSetting branch length to %f\n", TIndex, MINBL, MINBL);
					
					
					if(NoErr == 5)
						printf("5 or more branch length error have occurred. Suppressing warnings.\n");

					NoErr++;

					Node->Length = MINBL;
				}
			}
		}
	}
}

void	AllocNOSPerSite(OPTIONS *Opt)
{
	TREES	*Trees;
	int		Index;

	Trees = Opt->Trees;
	
	Trees->NOSPerSite	= TRUE;
	Trees->MaxNOS		= Trees->NoOfStates;
	Trees->NOSList		= (int*)malloc(sizeof(int) * Trees->NoOfSites);
	Trees->SiteSymbols	= (char**)malloc(sizeof(char*) * Trees->NoOfSites);

	if((Trees->NOSList == NULL) || (Trees->SiteSymbols == NULL))
		MallocErr();

	for(Index=0;Index<Trees->NoOfSites;Index++)
		FindSiteSymbols(Trees, Index);
}

void	SetNOSPerSite(OPTIONS *Opt)
{
	TREES	*Trees;
	int		Index;

	if(Opt->NOSPerSite == FALSE)
		return;
	
	Trees = Opt->Trees;
	AllocNOSPerSite(Opt);

	for(Index=0;Index<Trees->NoOfTrees;Index++)
		SetNOSPerSiteTipData(Trees->Tree[Index], Trees);
}

void	ReLinkNodes(NODE NewList, NODE OldList, int Size)
{
	int		Index, NIndex;
	NODE	N;
	NODE	C;

	for(Index=0;Index<Size;Index++)
		OldList[Index].ID = Index;

	for(Index=0;Index<Size;Index++)
	{
		C = &OldList[Index];
		N = &NewList[Index];

		if(C->Tip == FALSE)
		{
			for(NIndex=0;NIndex<N->NoNodes;NIndex++)
				N->NodeList[NIndex] = &NewList[C->NodeList[NIndex]->ID];

		}
		N->Ans = NULL;
	}
}

void	ReLinkAns(NODE N)
{
	int NIndex;

	if(N->Tip == TRUE)
		return;

	for(NIndex=0;NIndex<N->NoNodes;NIndex++)
		N->NodeList[NIndex]->Ans = N;

	for(NIndex=0;NIndex<N->NoNodes;NIndex++)
		ReLinkAns(N->NodeList[NIndex]);

}

void	AddNewRecNodeTree(TREES *Trees, TREE *Tree, RECNODE *RecNode)
{
	NODE	NewTaxa, Node;
	NODE	*NList;
	int		NoNodes, Index;
	
	NoNodes = 0;
	Node = FindNode(RecNode, Tree, &NoNodes);
	
	NewTaxa = AllocNode();

	NList = (NODE*)malloc(sizeof(NODE) * (Node->NoNodes + 1));
	if(NList == NULL)
		MallocErr();

	memcpy(NList, Node->NodeList, sizeof(NODE) * Node->NoNodes);
	NList[Node->NoNodes] = NewTaxa;
	free(Node->NodeList);
	Node->NodeList = NList;
	Node->NoNodes++;

	NewTaxa->Ans	= Node;
	NewTaxa->Length = MINBL;
	NewTaxa->Tip	= TRUE;
	NewTaxa->Taxa	= GetTaxaFromName(RecNode->Name, Trees->Taxa, Trees->NoOfTaxa);
	NewTaxa->TipID	= NewTaxa->Taxa->No;
	NewTaxa->ID		= Tree->NoNodes;
	
	SetNodeList(Tree);

	for(Index=0;Index<Tree->NoNodes;Index++)
	{
		Node = Tree->NodeList[Index];
		if(Node->Tip == TRUE)
			Node->Taxa = GetTaxaFromID(Node->TipID, Trees->Taxa, Trees->NoOfTaxa);
	}

}

void	AddNewRecNode(TREES* Trees, RECNODE *RecNode)
{
	int	TIndex;
	TREE	*Tree;

	for(TIndex=0;TIndex<Trees->NoOfTrees;TIndex++)
	{
		Tree = Trees->Tree[TIndex];
		AddNewRecNodeTree(Trees, Tree, RecNode);
	}
} 

double	GetRootToTip(NODE N)
{
	double Ret;
	
	Ret = 0;
	do
	{
		N = N->NodeList[0];
		Ret += N->Length;
	}while(N->Tip == FALSE);

	return Ret;
}

double	GetTipPathLength(NODE Node)
{
	double	Ret;

	Ret = 0;

	do
	{
		Ret += Node->Length;
		Node = Node->Ans;
	} while(Node->Ans != NULL);
	
	return Ret;
}

void	MakeTreeUM(TREES *Trees, int TNo, double RootTip)
{
	int		Index;
	double	Dist;
	NODE	N;
	TREE	*Tree;

	Tree = Trees->Tree[TNo];

	for(Index=0;Index<Tree->NoNodes;Index++)
	{
		N = Tree->NodeList[Index];
		if(N->Tip == TRUE)
		{
			Dist = RootTip - GetTipPathLength(N);
			N->Length += Dist;
			if(N->Length < 0)
			{
				printf("Tree\t%d\tTip\t%s\tcannot be made into an ultrametric tree.\n", TNo, N->Taxa->Name);
				exit(0);
			}
		}
	}
}

void	MakeUM(TREES* Trees)
{
	int	Index;
	double RootTip;

	for(Index=0;Index<Trees->NoOfTrees;Index++)
	{
		RootTip = GetRootToTip(Trees->Tree[Index]->Root);

//		printf("Tree\t%d\t%f\n", Index, RootTip); 

		MakeTreeUM(Trees, Index, RootTip);
	}

//	exit(0);
}


NODE	GetNodeFromTID(TREE *Tree, int ID)
{
	int Index;
	NODE	N;

	
	for(Index=0;Index<Tree->NoNodes;Index++)
	{
		if(Tree->NodeList[Index]->Tip == TRUE)
		{
			N = Tree->NodeList[Index];
			if(N->Taxa->No == ID)
				return N;
		}
	}

	return NULL;
}

void	ListOddPPTaxa(TREES *Trees)
{
	TREE	*Tree;
	int		TIndex, Index;
	TAXA	*CT, *T;
	NODE	CN, N;
	int		No;
	int		GID;

	Tree = Trees->Tree[0];


	for(TIndex=0;TIndex<Tree->NoNodes;TIndex++)
	{
		printf("%d\t%f\t%f\t", TIndex, Tree->NodeList[TIndex]->Length, log(Trees->Tree[0]->NodeList[TIndex]->Length));

		if(Tree->NodeList[TIndex]->Length == 0)
			printf("Zeor\n");
		else
			printf("NonZero\n");

	}

	exit(0);

	printf("\n\n");
	for(TIndex=0;TIndex<Trees->NoOfTaxa;TIndex++)
	{
		CT = Trees->Taxa[TIndex];
		CN = GetNodeFromTID(Tree, CT->No);

		printf("%d\t%s\t%f\t%f\n", TIndex, CT->Name, CT->ConData[0], N->Length);
				
	}
	exit(0);

	GID = 0;
	for(TIndex=0;TIndex<Trees->NoOfTaxa;TIndex++)
	{
		CT = Trees->Taxa[TIndex];
		CN = GetNodeFromTID(Tree, CT->No);
		if(CT->ConData[0] != -1)
		{
			No = 0;
			for(Index=TIndex+1;Index<Trees->NoOfTaxa;Index++)
			{
				T = Trees->Taxa[Index];
				if(T->ConData[0] == CT->ConData[0])
				{
					N = GetNodeFromTID(Tree, T->No);
					if(N->Length == CN->Length)
					{

						printf("%d\t%d\t%s\t%f\t%f\t0\n", GID, No+1, T->Name, T->ConData[0], N->Length);
						T->ConData[0] = -1;
						No++;
					}
				}					
			}

			if(No > 0)
			{
				printf("%d\t%d\t%s\t%f\t%f\t1\n\n\n", GID, No+1, CT->Name, CT->ConData[0], CN->Length);
				GID++;
			}
		}
	}

	exit(0);
}

void	CTaxaBelow(NODE N, int *No)
{
	int	Index;

	if(N->Tip == TRUE)
	{
		(*No)++;
		return;
	}

	for(Index=0;Index<N->NoNodes;Index++)
		CTaxaBelow(N->NodeList[Index], No);
}

char*	GetTName(FILE *Str, TREES *Trees, int TNo)
{
	int Index;

	for(Index=0;Index<Trees->NoOfTaxa;Index++)
	{
		if(Trees->Taxa[Index]->No == TNo)
			return Trees->Taxa[Index]->Name;
	}

	return NULL;
}

void	NormaliseTrees(double NormC, TREES *Trees)
{
	int Index, TIndex;
	NODE N;
	TREE *T;

	for(TIndex=0;TIndex<Trees->NoOfTrees;TIndex++)
	{
		T = Trees->Tree[TIndex];
		for(Index=0;Index<T->NoNodes;Index++)
		{
			N = T->NodeList[Index];
			if(N != T->Root)
				N->Length = NormC * N->Length;
		}
	}
}

double	FindTreeNormalise(TREES *Trees)
{
	int NoBL, Index, TIndex;
	double SumBL, Ret;
	NODE N;
	TREE *T;
	
	Ret = 0;
	
	
	for(TIndex=0;TIndex<Trees->NoOfTrees;TIndex++)
	{
		T = Trees->Tree[TIndex];
		SumBL = 0;
		NoBL = 0;
		for(Index=0;Index<T->NoNodes;Index++)
		{
			N = T->NodeList[Index];
			if(N != T->Root)
			{
				NoBL++;
				SumBL += N->Length;
			}
		}

		SumBL = SumBL / NoBL;
		Ret += SumBL * (1.0 / Trees->NoOfTrees); 
	}

//	Ret = SumBL / NoBL;
	
	Ret = NORM_MEAN_BL / Ret;

	return Ret;
}

void	RecSetDistToRoot(NODE N)
{
	int Index;

	if(N->Ans == NULL)
		N->DistToRoot = 0;
	else
		N->DistToRoot = N->Ans->DistToRoot + N->Length;

	if(N->Tip == TRUE)
		return;

	for(Index=0;Index<N->NoNodes;Index++)
		RecSetDistToRoot(N->NodeList[Index]);	
}

void	SetTreeDistToRoot(TREE *Tree)
{
	RecSetDistToRoot(Tree->Root);
}

void	SetTreesDistToRoot(TREES *Trees)
{
	int Index;

	for(Index=0;Index<Trees->NoOfTrees;Index++)
		SetTreeDistToRoot(Trees->Tree[Index]);
}

void	SetVisitedNode(NODE N, int Val)
{
	int Index;

	N->Visited = Val;
	if(N->Tip == TRUE)
		return;

	for(Index=0;Index<N->NoNodes;Index++)
		SetVisitedNode(N->NodeList[Index], Val);
}

void	SetVisitedTree(TREE *Tree, int Val)
{
	int Index;

	for(Index=0;Index<Tree->NoNodes;Index++)
		Tree->NodeList[Index]->Visited = Val;
}

void	SimRandChange(NODE N, double RateDev, RANDSTATES *RS)
{

}

double*	GetPhyChanges(TREES *Trees, TREE *Tree, double RateDev, RANDSTATES *RS)
{
	double	*Ret;

	Ret = NULL;

	return Ret;
}

void	AddTaxaErr(TREES *Trees, int TaxaID, double Err)
{
	int TIndex;
	TREE *Tree;
	NODE N;

	for(TIndex=0;TIndex<Trees->NoOfTrees;TIndex++)
	{
		Tree = Trees->Tree[TIndex];
		N = GetNodeFromTID(Tree, TaxaID);

		if(N == NULL)
		{
			printf("%s::%d cound not find taxa.\n", __FILE__, __LINE__);
			exit(0);
		}
		N->Length += Err;
	}
}

int		TaxaIndexToNo(TREES *Trees, int Index)
{
	return Trees->Taxa[Index]->No;
}

int		TaxaNoToIndex(TREES *Trees, int No)
{
	int Index;

	for(Index=0;Index<Trees->NoOfTaxa;Index++)
		if(Trees->Taxa[Index]->No == No)
			return Index;

	printf("%s::%d Taxa No %d not valid.\n", __FILE__, __LINE__, No);
	exit(0);
	return -1;
}


void	SetUserBranchLength(TREE *Tree)
{
	int Index;
	NODE N;

	for(Index=0;Index<Tree->NoNodes;Index++)
	{
		N = Tree->NodeList[Index];
		N->Length = N->UserLength;
	}
}

void	ReSetBranchLength(TREE *Tree)
{
	SetUserBranchLength(Tree);
	SetTreeDistToRoot(Tree);
}


void	RecGetSumBL(NODE N, double *SumBL)
{
	int Index;

	if(N->Length >= 0)
		*SumBL += N->Length;

	if(N->Tip == TRUE)
		return;

	for(Index=0;Index<N->NoNodes;Index++)
		RecGetSumBL(N->NodeList[Index], SumBL);
}

double	SumNodeBL(NODE N)
{
	double Ret;

	Ret = 0;
	RecGetSumBL(N, &Ret);

//	Ret = Ret - N->Length;

	return Ret;
}

void	RecScaleSubTree(NODE N, double Scale)
{
	int Index;

	N->Length = N->Length * Scale;
	
	if(N->Tip == TRUE)
		return;

	for(Index=0;Index<N->NoNodes;Index++)
		RecScaleSubTree(N->NodeList[Index], Scale);
}

void	ScaleSubTree(NODE N, double Scale)
{
	if(N->Tip == TRUE)
	{
		printf("Sub Tree can only scale internal nodes.\n");
		exit(0);
	}


	RecScaleSubTree(N, Scale);

//	for(Index=0;Index<N->NoNodes;Index++)
//		RecScaleSubTree(N->NodeList[Index], Scale);
}

void	ScaleTrees(TREES *Trees, double Scale)
{
	int Index;

	for(Index=0;Index<Trees->NoOfTrees;Index++)
		RecScaleSubTree(Trees->Tree[Index]->Root, Scale);

}

void	RecScaleUserSubTree(NODE N, double Scale)
{
	int Index;

	N->Length = N->Length * Scale;
	N->UserLength = N->UserLength * Scale;

	if(N->Tip == TRUE)
		return;

	for(Index=0;Index<N->NoNodes;Index++)
		RecScaleUserSubTree(N->NodeList[Index], Scale);
}


void	ScaleUserTrees(TREES *Trees, double Scale)
{
	int Index;

	for(Index=0;Index<Trees->NoOfTrees;Index++)
		RecScaleUserSubTree(Trees->Tree[Index]->Root, Scale);
}	

NODE	GetTreeTaxaNode(TREE *Tree, int TaxaNo)
{
	int Index;

	for(Index=0;Index<Tree->NoNodes;Index++)
		if(Tree->NodeList[Index]->Tip == TRUE)
			if(Tree->NodeList[Index]->Taxa->No == TaxaNo)
				return Tree->NodeList[Index];

	return NULL;
}

void	InitialiseOutputTrees(OPTIONS *Opt, TREES *Trees)
{
	char	*Buffer;

	Buffer = (char*)SMalloc(sizeof(char) * (strlen(Opt->LogFN) + BUFFERSIZE));

	sprintf(Buffer, "%s.Output.trees", Opt->LogFN);
	Opt->OutTrees = OpenWrite(Buffer);
	free(Buffer);

	SaveTreesHeader(Opt->OutTrees, Trees);

	fflush(Opt->OutTrees);
}

void	OutputTree(OPTIONS *Opt, TREES *Trees, RATES *Rates, long long No, FILE *Out)
{
	TREE *Tree;
	NODE Root;
	 
	fprintf(Out, "\t\ttree No_%lld = ", No);
	Tree = Trees->Tree[Rates->TreeNo];

	Root = Tree->Root;

	WriteTreeToFile(Root, Root, Out);
	fprintf(Out, ";\n");

	fflush(Out);
}