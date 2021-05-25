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
#include "treenode.h"
#include "treepasser.h"
#include "treenode.h"


void	WriteTreeToFile(NODE n, NODE Root, FILE *F);

void	BlankNode(NODE N)
{
		N->ID		=	-1;
		N->Ans		=	NULL;
		N->Left		=	NULL;
		N->Right	=	NULL;
		N->Tip		=	FALSE;
		N->TipID	=	-1;
		N->Length	=	-1;
		N->Taxa		=	NULL;
		N->Partial	=	NULL;
		N->FossilState=	-1;
		N->Part		=	NULL;
		N->PSize	=	-1;
		N->GammaPartial = NULL;
}


TAXA*	GetTaxaFromID(int ID, TAXA *Taxa, int NoOfTaxa)
{
	int	Index;

	for(Index=0;Index<NoOfTaxa;Index++)
	{
		if(Taxa[Index].No == ID)
			return &Taxa[Index];
	}
	return NULL;
}

TAXA*	GetTaxaFromName(char *Name, TAXA *Taxa, int NoOfTaxa)
{
	int	Index;

	for(Index=0;Index<NoOfTaxa;Index++)
	{
		if(strcmp(Taxa[Index].Name, Name) == 0)
			return &Taxa[Index];
	}
	return NULL;
}


void	LinkTipsToTaxa(NODE N, TAXA *Taxa, int NoOfTaxa)
{
	if(N->Tip == TRUE)
	{
		N->Taxa = GetTaxaFromID(N->TipID, Taxa, NoOfTaxa);
	}
	else
	{
		LinkTipsToTaxa(N->Left, Taxa, NoOfTaxa);
		LinkTipsToTaxa(N->Right, Taxa, NoOfTaxa);
	}

}

void	FreeTree(TREE* Tree, int NoOfNodes, int NoOfSites, int NoOfTaxa)
{
	int		NIndex;
	int		DIndex;
	NODE	N=NULL;
	int		TIndex=0;

	for(NIndex=0;NIndex<NoOfNodes;NIndex++)
	{
		N = &Tree->NodeList[NIndex];

		if(N->Partial != NULL)
		{
			for(DIndex=0;DIndex<NoOfSites;DIndex++)
			{
				free(N->Partial[DIndex]);
				if(N->GammaPartial != NULL)
					free(N->GammaPartial[DIndex]);
			}
		}

		free(N->Partial);

		if(N->GammaPartial != NULL)
			free(N->GammaPartial);
		TIndex++;
	}

	free(Tree->NodeList);

	if(Tree->ConVars!= NULL)
		FreeConVar(Tree->ConVars, NoOfTaxa);
}

void	FreeTrees(TREES* Trees, OPTIONS *Opt)
{
	int	Index;

	FreeData(Opt);
	free(Trees->Taxa);

	FreePartitions(Trees);

	for(Index=0;Index<Trees->NoOfTrees;Index++)
		FreeTree(&Trees->Tree[Index], Trees->NoOfNodes, Trees->NoOfSites, Trees->NoOfTaxa);

	if(Trees->PLeft != NULL)
		FreeMatrix(Trees->PLeft);

	if(Trees->PRight != NULL)
		FreeMatrix(Trees->PRight);
	
	if(Trees->InvInfo != NULL)
		FreeInvInfo(Trees->InvInfo);

	if(Trees->RemovedTaxa != NULL)
	{
		for(Index=0;Index<Trees->NoOfRemovedTaxa;Index++)
			free(Trees->RemovedTaxa[Index]);
		free(Trees->RemovedTaxa);
	}

	free(Trees->Tree);
	if(Trees->SymbolList != NULL)
		free(Trees->SymbolList);
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
	CheckBaseTirs(Trees);
	
	/* Only binary trees */
	CheckPoly(Trees);

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
	int		Index;
	NODE	Node;
	NNODE	NNode;
	int		Pos;

	Tree->NodeList = (NODE)malloc(sizeof(struct INODE) * PTree->NoOfNodes);
	if(Tree->NodeList == NULL)
		MallocErr();
	memset(Tree->NodeList , '\0', sizeof(struct INODE) * PTree->NoOfNodes);

	for(Index=0;Index<PTree->NoOfNodes;Index++)
		BlankNode(&Tree->NodeList[Index]);

	for(Index=0;Index<PTree->NoOfNodes;Index++)
	{
		Node = &Tree->NodeList[Index];
		NNode= PTree->NodeList[Index];

		Node->Length	= NNode->Length;
		Node->Tip		= NNode->Tip;
		Node->Contrast	= NULL;

		Node->Part = NULL;
		Node->PSize= 0;

		if(Node->Tip == TRUE)
		{
			Node->TipID = NNode->TaxaID;
		}
		else
		{
			Pos = GetArrPos(&PTree->NodeList[0], NNode->Left);
			Node->Left = &Tree->NodeList[Pos];

			Pos = GetArrPos(&PTree->NodeList[0], NNode->Right);
			Node->Right = &Tree->NodeList[Pos];
		}

		if(NNode->Ans == NULL)
			Node->Ans = NULL;
		else
		{
			Pos = GetArrPos(&PTree->NodeList[0], NNode->Ans);
			Node->Ans = &Tree->NodeList[Pos];
		}
	}

	Pos = GetArrPos(&PTree->NodeList[0], PTree->Root);
	Tree->Root = &Tree->NodeList[Pos];
}

void	InitialTrees(TREES *Trees, NTREES *PTrees)
{
	int		Index;
	TAXA	*Taxa;
	TREE	*Tree;
	int		NodeID;

	Trees->NoOfTaxa		= PTrees->NoOfTaxa;
	Trees->NoOfTrees	= PTrees->NoOfTrees;
	Trees->NoOfNodes	= PTrees->Trees[0].NoOfNodes;
	Trees->NOSPerSite	= FALSE;
	 

	Trees->Taxa = (TAXA*) malloc(sizeof(TAXA) * Trees->NoOfTaxa);
	if(Trees->Taxa == NULL)
		MallocErr();
	/* Keep Purify Happy, because of padding */
	memset(Trees->Taxa , '\0', sizeof(TAXA) * Trees->NoOfTaxa);

	for(Index=0;Index<Trees->NoOfTaxa;Index++)
	{
		Taxa = &Trees->Taxa[Index];

		Taxa->Name			= StrMake(PTrees->Taxa[Index].Name);
		Taxa->No			= PTrees->Taxa[Index].No;
		Taxa->DesDataChar	= NULL;
		Taxa->ConData		= NULL;
		Taxa->Exclude		= FALSE;


		Taxa->EstData		=	FALSE;
		Taxa->EstDepData	=	FALSE;
		Taxa->EstDataP		=	NULL;
		Taxa->RealData		=	NULL;
	}

	Trees->Tree = (TREE*)malloc(sizeof(TREE) * PTrees->NoOfTrees);
	if(Trees->Tree == NULL)
		MallocErr();

	for(Index=0;Index<Trees->NoOfTrees;Index++)
	{
		Tree = &Trees->Tree[Index];

		Tree->ConVars	= NULL;
		Tree->NodeList	= NULL;
		Tree->Root		= NULL;

		MakeNewTree(Tree, &PTrees->Trees[Index]);

		LinkTipsToTaxa(Trees->Tree[Index].Root, Trees->Taxa, Trees->NoOfTaxa);

		NodeID = 0;
		SetNodeIDs(Trees->Tree[Index].Root, &NodeID);
	}

	SetPartitions(Trees);
}

TREES*	LoadTrees(char* FileName)
{
	TREES*		Ret=NULL;
	NTREES*		PTrees;
	char		*Err;

	Ret = (TREES*)malloc(sizeof(TREES));
	if(Ret==NULL)
		MallocErr();

	Ret->PLeft				= NULL;
	Ret->PRight				= NULL;
	Ret->Taxa				= NULL;
	Ret->Tree				= NULL;
	Ret->InvInfo			= NULL;
	Ret->SymbolList			= NULL;
	Ret->RemovedTaxa		= NULL;
	Ret->NoOfRemovedTaxa	= 0;
	Ret->ValidCData			= TRUE;
	Ret->ValidDData			= TRUE;

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

	Ret->JStop	=	FALSE;

	return Ret;
}

void	PrintTrees(FILE*	Str, TREES *Trees, DATATYPE DataType)
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

void	SetNodeTipData(NODE N, TREE* Tree, TREES *Trees)
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

		/* Copy the sites for the covarion mode */
		if(Trees->UseCovarion == TRUE)
		{
			for(StateIndex=0;StateIndex<NOS;StateIndex++)
				N->Partial[SiteIndex][StateIndex+NOS] = N->Partial[SiteIndex][StateIndex];
		}
	}
}
void	SetTipData(TREE *Tree, TREES *Trees)
{
	int		NIndex;
	NODE	N;
	

	for(NIndex=0;NIndex<Trees->NoOfNodes;NIndex++)
	{
		N = &Tree->NodeList[NIndex];
		if(N->Tip == TRUE)
			SetNodeTipData(N, Tree, Trees);
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

	for(NIndex=0;NIndex<Trees->NoOfNodes;NIndex++)
	{
		N = &Tree->NodeList[NIndex];
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

void	AllocPartial(TREES* Trees, int Gamma)
{
	int		TIndex;
	int		NIndex;
	NODE	N;

	
	for(TIndex=0;TIndex<Trees->NoOfTrees;TIndex++)
	{
		for(NIndex=0;NIndex<Trees->NoOfNodes;NIndex++)
		{
			N = &Trees->Tree[TIndex].NodeList[NIndex];
			AllocNodePartial(N, Trees, Gamma);
		}

		SetTipData(&Trees->Tree[TIndex], Trees);
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

void	CopyNode(NODE A, NODE B)
{
	A->Tip		=	B->Tip;
	A->TipID	=	B->TipID;

	A->Length	=	B->Length;

	A->Left		=	B->Left;
	A->Right	=	B->Right;
	A->Ans		=	B->Ans;

	A->Visited	=	B->Visited;

	A->Partial	=	B->Partial;

	A->Taxa		=	B->Taxa;
}

void	PrintTaxaNames(NODE N)
{
	if(N->Tip == TRUE)
		printf("%s\t%f\n", N->Taxa->Name, N->Length);
	else
	{
		PrintTaxaNames(N->Left);
		PrintTaxaNames(N->Right);
	}
}

void	ReBuildTree(NODE N, NODE NewNodeList)
{
	NODE	NewNode;

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
		NewNode->Left		=	&NewNodeList[N->Left->ID];
		NewNode->Right		=	&NewNodeList[N->Right->ID];

		ReBuildTree(N->Left, NewNodeList);
		ReBuildTree(N->Right, NewNodeList);
	}
	else
	{
		NewNode->Left = NULL;
		NewNode->Right= NULL;
	}
}



void	RemoveTaxaFromTree(TREES *Trees, TREE *Tree, char *TName)
{
	NODE	DelNode=NULL;
	NODE	KeepNode=NULL;
	NODE	DelTaxaNode=NULL;
	NODE	NewNodeList=NULL;
	NODE	N;
	int		NIndex;
	int		OnLeft;

	for(NIndex=0;NIndex<Trees->NoOfNodes;NIndex++)
	{
		N = &Tree->NodeList[NIndex];
		if(N->Tip == TRUE)
		{
			if(strcmp(N->Taxa->Name, TName) == 0)
			{
				DelTaxaNode		= N;
				DelNode			= N->Ans;
				if(DelNode->Left == N)
					OnLeft = TRUE;
				else
					OnLeft = FALSE;
			}
		}
	}

	if(DelNode == NULL)
	{
		printf("Could not find node to del in %s::RemoveTaxaFromTree %d\n", __FILE__, __LINE__);
		exit(1);
	}

	if(DelNode != Tree->Root)
	{
		if(OnLeft == TRUE)
		KeepNode = DelNode->Right;
		else
			KeepNode = DelNode->Left;

		KeepNode->Length += DelNode->Length;
		KeepNode->Ans = DelNode->Ans;
	
		if(DelNode == DelNode->Ans->Left)
			DelNode->Ans->Left = KeepNode;
		else
			DelNode->Ans->Right= KeepNode;
	}
	else
	{
		if(DelNode->Left == DelTaxaNode)
			Tree->Root = DelNode->Right;
		else
			Tree->Root = DelNode->Left;

		Tree->Root->Ans = NULL;
	}

	BlankNode(DelNode);
	BlankNode(DelTaxaNode);

	NewNodeList = (NODE)malloc(sizeof(struct INODE) * (Trees->NoOfNodes - 2));
	if(NewNodeList == NULL)
		MallocErr();
	
	for(NIndex=0;NIndex<Trees->NoOfNodes-2;NIndex++)
	{
		N = &NewNodeList[NIndex];
		BlankNode(N);
	}

	NIndex=0;
	SetNodeIDs(Tree->Root, &NIndex);

	ReBuildTree(Tree->Root, NewNodeList);

	free(Tree->NodeList);
	Tree->NodeList = NewNodeList;
	Tree->Root = &Tree->NodeList[0];
}

int	RemoveTaxa(OPTIONS *Opt, TREES *Trees, char *TName)
{
	int		TIndex;
	int		OldTIndex;
	int		NewTIndex;
	char	**NewUsedTaxa=NULL;
	TAXA	*NewTaxaList=NULL;
	int		FoundTaxa=FALSE;
	int		TPos;
	char	*TaxaName;

	if(Trees->NoOfTaxa <= 2)
	{
		printf("There must be two or more taxa in a tree\n");
		return FALSE;
	}

	if(Opt != NULL)
	{
		if(Opt->NoOfRecNodes != 0)
		{
			printf("There must be no node to reconsutct befor a taxa can be removed for a tree\n");
			return FALSE;
		}
	}

	for(TIndex=0;TIndex<Trees->NoOfTaxa;TIndex++)
	{
		if(strcmp(TName, Trees->Taxa[TIndex].Name)==0)
		{
			FoundTaxa = TRUE;
			TPos = TIndex;
		}
	}

	
	if(FoundTaxa == FALSE)
	{
		printf("Could not find taxa %s\n", TName);
		return FALSE;
	}

	
	TaxaName = StrMake(TName);
	
	for(TIndex=0;TIndex<Trees->NoOfTrees;TIndex++)
		RemoveTaxaFromTree(Trees, &Trees->Tree[TIndex], TName);

	Trees->NoOfNodes = Trees->NoOfNodes - 2;

	/*
		Must delete the taxa form the List 
	*/

	NewTaxaList = (TAXA*)malloc(sizeof(TAXA) * (Trees->NoOfTaxa - 1));
	if(NewTaxaList == NULL)
		MallocErr();

	NewTIndex = 0;
	for(OldTIndex=0;OldTIndex<Trees->NoOfTaxa;OldTIndex++)
	{
		if(OldTIndex != TPos)
		{
			NewTaxaList[NewTIndex].Name			= Trees->Taxa[OldTIndex].Name;
			NewTaxaList[NewTIndex].ConData		= Trees->Taxa[OldTIndex].ConData;
			NewTaxaList[NewTIndex].DesDataChar	= Trees->Taxa[OldTIndex].DesDataChar;
			NewTaxaList[NewTIndex].No			= Trees->Taxa[OldTIndex].No;
			NewTaxaList[NewTIndex].Exclude		= Trees->Taxa[OldTIndex].Exclude;
			NewTaxaList[NewTIndex].EstData		= Trees->Taxa[OldTIndex].EstData;
			NewTaxaList[NewTIndex].EstDataP		= Trees->Taxa[OldTIndex].EstDataP;

			NewTIndex++;
		}
	}

	FreeTaxa(&Trees->Taxa[TPos], Trees->NoOfSites);

	free(Trees->Taxa);
	Trees->Taxa	= NewTaxaList;
	Trees->NoOfTaxa--;

	for(TIndex=0;TIndex<Trees->NoOfTrees;TIndex++)
		LinkTipsToTaxa(Trees->Tree[TIndex].Root, Trees->Taxa, Trees->NoOfTaxa);
	
	NewUsedTaxa = (char**)malloc(sizeof(char*) * (Trees->NoOfRemovedTaxa + 1));
	if(NewUsedTaxa == NULL)
		MallocErr();

	for(TIndex=0;TIndex<Trees->NoOfRemovedTaxa;TIndex++)
		NewUsedTaxa[TIndex] = Trees->RemovedTaxa[TIndex];

	NewUsedTaxa[TIndex] = TaxaName;

	free(Trees->RemovedTaxa);
	Trees->RemovedTaxa = NewUsedTaxa;
	Trees->NoOfRemovedTaxa++;

	return TRUE;
}

void WriteTreeToFile(NODE n, NODE Root, FILE *F)
{
	double	BL;

	BL = n->Length;

	if(n->Tip == TRUE)
		fprintf(F, "%d:%10.10f", n->TipID, BL);
	else
	{
		fprintf(F, "(");

		WriteTreeToFile(n->Left, Root, F);

		fprintf(F, ",");
 
		WriteTreeToFile(n->Right, Root, F);
		
		if(n != Root)
			fprintf(F, "):%10.10f", BL);	
		else
			fprintf(F, ")");	

	}
}

void	PrintTree(char	*FileName, TREES* Trees, OPTIONS *Opt)
{

	FILE	*TreeFile=NULL;
	int		TIndex;
	char	*Name;
	char	Buffer[128];
	int		NoOfChar;


	TreeFile = fopen(FileName, "w");
	if(TreeFile == NULL)
	{
		printf("Could not open file %s for writting\n", FileName);
		return;
	}

	sprintf(&Buffer[0], "%d", Trees->NoOfTrees);
	NoOfChar = strlen(&Buffer[0]);

	fprintf(TreeFile, "#NEXUS\n");

	fprintf(TreeFile, "begin trees;\n");
	fprintf(TreeFile, "\ttranslate\n");

	for(TIndex=0;TIndex<Trees->NoOfTaxa;TIndex++)
	{

		fprintf(TreeFile, "\t\t%d %s", Trees->Taxa[TIndex].No, Trees->Taxa[TIndex].Name);

		if(TIndex<Trees->NoOfTaxa-1)
			fprintf(TreeFile, ",\n");
		else
			fprintf(TreeFile, ";\n");
	}


	for(TIndex=0;TIndex<Trees->NoOfTrees;TIndex++)
	{
		Name = FormatInt(TIndex+1, NoOfChar);

		fprintf(TreeFile, "\t\ttree No_%s = ", Name);
		WriteTreeToFile(Trees->Tree[TIndex].Root, Trees->Tree[TIndex].Root, TreeFile);
		fprintf(TreeFile, ";\n");
		
		free(Name);
	}

	fprintf(TreeFile, "end;\n");
	fclose(TreeFile);

}	

void	SetNodeIDs(NODE N, int *No)
{
	N->ID = *No;

	*No = (*No) + 1;

	if(N->Tip == FALSE)
	{
		SetNodeIDs(N->Left, No);
		SetNodeIDs(N->Right, No);
	}
}

void	NormaliseTrees(TREE *Tree, int NoOfNodes)
{
	double	Total=0;
	int		Index=0;
	
	for(Index=1;Index<NoOfNodes;Index++)
		Total += Tree->NodeList[Index].Length;
		
	for(Index=1;Index<NoOfNodes;Index++)
		Tree->NodeList[Index].Length = Tree->NodeList[Index].Length / Total;

}

void	SetFossiles(TREES *Trees, OPTIONS *Opt)
{
	RECNODE	RNode;
	int		TIndex;

	RNode = Opt->RecNode;

	while(RNode != NULL)
	{
		if(RNode->NodeType == FOSSIL)
		{
			for(TIndex=0;TIndex<Trees->NoOfTrees;TIndex++)
				RNode->TreeNodes[TIndex]->FossilState = RNode->FossilState;
		}

		RNode = RNode->Next;
	}
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
		Tree = &Trees->Tree[TIndex];
		for(NIndex=0;NIndex<Trees->NoOfNodes;NIndex++)
		{
			Node = &Tree->NodeList[NIndex];
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
		SetNOSPerSiteTipData(&Trees->Tree[Index], Trees);
}

void	ReLinkNodes(NODE NewList, NODE OldList, int Size)
{
	int	Index;
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
			N->Left = &NewList[C->Left->ID];
			N->Right= &NewList[C->Right->ID];
		}
		N->Ans = NULL;
	}
}

void	ReLinkAns(NODE N)
{
	if(N->Tip == TRUE)
		return;

	N->Left->Ans = N;
	N->Right->Ans= N;

	ReLinkAns(N->Left);
	ReLinkAns(N->Right);	
}

void	AddNewRecNodeTree(TREES *Trees, TREE *Tree, RECNODE RecNode)
{
	NODE	New;
	NODE	NewTaxa;
	NODE	NodeList;
	NODE	Node; 
	int		NoNodes;
	int		Index;

	NoNodes = 0;
	Node = FindNode(RecNode, Tree, &NoNodes, Trees->NoOfNodes);

	NoNodes = Trees->NoOfNodes;

	NodeList = (NODE)malloc(sizeof(struct INODE) * (NoNodes + 2));
	if(NodeList == NULL)
		MallocErr();
	memset(NodeList , '\0', sizeof(struct INODE) * (NoNodes + 2));

	memcpy(NodeList, Tree->NodeList, sizeof(struct INODE) * NoNodes);
	ReLinkNodes(NodeList, Tree->NodeList, Trees->NoOfNodes);

	Tree->Root = &NodeList[Tree->Root->ID];
	ReLinkAns(Tree->Root);

	Node = &NodeList[Node->ID];

	New = &NodeList[NoNodes];
	NewTaxa = &NodeList[NoNodes+1];

	BlankNode(New);
	BlankNode(NewTaxa);

	New->Length = MINBL;
	NewTaxa->Length = MINBL;

	New->Left = Node->Left;
	New->Right= Node->Right;

	New->Left->Ans = New;
	New->Right->Ans= New;

	New->Ans = Node;
	NewTaxa->Ans = Node;

	Node->Left = New;
	Node->Right = NewTaxa;


	New->Tip = FALSE;

	NewTaxa->Tip = TRUE;
	NewTaxa->Taxa= GetTaxaFromName(RecNode->Name, Trees->Taxa, Trees->NoOfTaxa);
	NewTaxa->TipID= NewTaxa->Taxa->No;

	New->ID = Trees->NoOfNodes;
	NewTaxa->ID = Trees->NoOfNodes+1;

	free(Tree->NodeList);
	Tree->NodeList = NodeList;

	for(Index=0;Index<NoNodes;Index++)
	{
		Node = &NodeList[Index];
		if(Node->Tip == TRUE)
			Node->Taxa = GetTaxaFromID(Node->TipID, Trees->Taxa, Trees->NoOfTaxa);
	}
}

/*
void	AddNewRecNodeTree(TREES *Trees, TREE *Tree, RECNODE RecNode)
{
	NODE	New;
	NODE	NewTaxa;
	NODE	Ans;
	NODE	NodeList;
	NODE	Node; 
	int		NoNodes;
	int		Index;

	NoNodes = 0;
	Node = FindNode(RecNode, Tree, &NoNodes, Trees->NoOfNodes);

	NoNodes = Trees->NoOfNodes;

	NodeList = (NODE)malloc(sizeof(struct INODE) * (NoNodes + 2));
	if(NodeList == NULL)
		MallocErr();

	memcpy(NodeList, Tree->NodeList, sizeof(struct INODE) * NoNodes);
	ReLinkNodes(NodeList, Tree->NodeList, Trees->NoOfNodes);

	Tree->Root = &NodeList[Tree->Root->ID];
	ReLinkAns(Tree->Root);

	Node = &NodeList[Node->ID];
	Ans = Node->Ans;

	New = &NodeList[NoNodes];
	NewTaxa = &NodeList[NoNodes+1];

	BlankNode(New);
	BlankNode(NewTaxa);

	New->Length = MINBL;
	NewTaxa->Length = MINBL;


	New->Left = Node;
	New->Right= NewTaxa;
	New->Ans = Ans;

	if(Node == Ans->Left)
		Ans->Left = New;
	else
		Ans->Right= New;
	
	New->Left->Ans = New;
	New->Right->Ans= New;

	New->Tip = FALSE;

	NewTaxa->Ans = Ans;
	NewTaxa->Tip = TRUE;
	NewTaxa->Taxa= GetTaxaFromName(RecNode->Name, Trees->Taxa, Trees->NoOfTaxa);
	NewTaxa->TipID= NewTaxa->Taxa->No;

	New->ID = Trees->NoOfNodes;
	NewTaxa->ID = Trees->NoOfNodes+1;

	free(Tree->NodeList);
	Tree->NodeList = NodeList;

	for(Index=0;Index<NoNodes;Index++)
	{
		Node = &NodeList[Index];
		if(Node->Tip == TRUE)
			Node->Taxa = GetTaxaFromID(Node->TipID, Trees->Taxa, Trees->NoOfTaxa);
	}
}
*/

void	AddNewRecNode(TREES* Trees, RECNODE RecNode)
{
	int	TIndex;
	TREE	*Tree;

	for(TIndex=0;TIndex<Trees->NoOfTrees;TIndex++)
	{
		Tree = &Trees->Tree[TIndex];
		AddNewRecNodeTree(Trees, Tree, RecNode);
	}

	Trees->NoOfNodes += 2;
} 

double	GetRootToTip(NODE N)
{
	double Ret;
	
	Ret = 0;
	do
	{
		N = N->Left;
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

	Tree = &Trees->Tree[TNo];

	for(Index=0;Index<Trees->NoOfNodes;Index++)
	{
		N = &Tree->NodeList[Index];
		if(N->Tip == TRUE)
		{
			Dist = RootTip - GetTipPathLength(N);
			N->Length += Dist;
			if(N->Length < 0)
			{
				printf("Tree\t%d\tTip\t%s\tcannot be made into an ultrametric tree.\n", TNo, N->Taxa->Name);
				exit(0);
			}
	//		printf("%s\t%30.30f\n", N->Taxa->Name, Dist);
		}
	}

//	exit(0);
}

void	MakeUM(TREES* Trees)
{
	int	Index;
	double RootTip;

	for(Index=0;Index<Trees->NoOfTrees;Index++)
	{
		RootTip = GetRootToTip(Trees->Tree[Index].Root);

//		printf("Tree\t%d\t%f\n", Index, RootTip); 

		MakeTreeUM(Trees, Index, RootTip);
	}

//	exit(0);
}


NODE	GetNodeFromTID(TREE *Tree, int NoNodes, int ID)
{
	int Index;
	NODE	N;

	for(Index=0;Index<NoNodes;Index++)
	{
		if(Tree->NodeList[Index].Tip == TRUE)
		{
			N = &Tree->NodeList[Index];
			if(N->Taxa->No == ID)
				return N;
		}
	}

	return NULL;
}

void	ListOddPPTaxa(TREES *Trees)
{
	int		TIndex, Index;
	TAXA	*CT, *T;
	NODE	CN, N;
	int		No;
	int		GID;


	for(TIndex=0;TIndex<Trees->NoOfNodes;TIndex++)
	{
		printf("%d\t%f\t%f\t", TIndex, Trees->Tree[0].NodeList[TIndex].Length, log(Trees->Tree[0].NodeList[TIndex].Length));

		if(Trees->Tree[0].NodeList[TIndex].Length == 0)
			printf("Zeor\n");
		else
			printf("NonZero\n");

	}

	exit(0);

	printf("\n\n");
	for(TIndex=0;TIndex<Trees->NoOfTaxa;TIndex++)
	{
		CT = &Trees->Taxa[TIndex];
		CN = GetNodeFromTID(&Trees->Tree[0], Trees->NoOfNodes, CT->No);

		printf("%d\t%s\t%f\n", TIndex, CT->Name, CT->ConData[0], N->Length);


		
	}
	exit(0);

	GID = 0;
	for(TIndex=0;TIndex<Trees->NoOfTaxa;TIndex++)
	{
		CT = &Trees->Taxa[TIndex];
		CN = GetNodeFromTID(&Trees->Tree[0], Trees->NoOfNodes, CT->No);
		if(CT->ConData[0] != -1)
		{
			No = 0;
			for(Index=TIndex+1;Index<Trees->NoOfTaxa;Index++)
			{
				T = &Trees->Taxa[Index];
				if(T->ConData[0] == CT->ConData[0])
				{
					N = GetNodeFromTID(&Trees->Tree[0], Trees->NoOfNodes, T->No);
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
	
	if(N->Tip == TRUE)
	{
		(*No)++;
		return;
	}
	CTaxaBelow(N->Left, No);
	CTaxaBelow(N->Right, No);
}