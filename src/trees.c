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

void	SetNodeIDs(NODE N, int *No);
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

	FreeData(Trees, Opt->Model);
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
			printf("Error: Tree %d has a %d-way bases polytomey, only rooted binary trees are allowed.\n", TIndex+1, Tree->Root->NoOfNodes);
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

	for(Index=0;Index<PTree->NoOfNodes;Index++)
		BlankNode(&Tree->NodeList[Index]);

	for(Index=0;Index<PTree->NoOfNodes;Index++)
	{
		Node = &Tree->NodeList[Index];
		NNode= PTree->NodeList[Index];

		Node->Length	= NNode->Length;
		Node->Tip		= NNode->Tip;

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

	Trees->NoOfTaxa = PTrees->NoOfTaxa;
	Trees->NoOfTrees = PTrees->NoOfTrees;
	Trees->NoOfNodes= PTrees->Trees[0].NoOfNodes;

	Trees->Taxa = (TAXA*) malloc(sizeof(TAXA) * Trees->NoOfTaxa);
	if(Trees->Taxa == NULL)
		MallocErr();

	for(Index=0;Index<Trees->NoOfTaxa;Index++)
	{
		Taxa = &Trees->Taxa[Index];

		Taxa->Name			= StrMake(PTrees->Taxa[Index].Name);
		Taxa->No			= PTrees->Taxa[Index].No;
		Taxa->DesDataChar	= NULL;
		Taxa->ConData		= NULL;
		Taxa->EstData		= FALSE;
		Taxa->EstDataP		= NULL;
		Taxa->Exclude		= FALSE;
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
					/* Set all boxes to 0 */
					for(StateIndex=0;StateIndex<NOS;StateIndex++)
						N->Partial[SiteIndex][StateIndex] = 0;
					/* Set the corrispoding boxes to 1 */
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
			FoundTaxa = TRUE;
	}

	if(FoundTaxa == FALSE)
	{
		printf("Could not find taxa %s\n", TName);
		return FALSE;
	}
	
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
		if(strcmp(Trees->Taxa[OldTIndex].Name, TName) != 0)
		{
			NewTaxaList[NewTIndex].Name			= Trees->Taxa[OldTIndex].Name;
			NewTaxaList[NewTIndex].ConData		= Trees->Taxa[OldTIndex].ConData;
			NewTaxaList[NewTIndex].DesDataChar	= Trees->Taxa[OldTIndex].DesDataChar;
			NewTaxaList[NewTIndex].No			= Trees->Taxa[OldTIndex].No;
			NewTaxaList[NewTIndex].Exclude		= Trees->Taxa[OldTIndex].Exclude;

			NewTIndex++;
		}
	}

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

	NewUsedTaxa[TIndex] = (char*)malloc(sizeof(char) * strlen(TName) + 1);
	if(NewUsedTaxa[TIndex] == NULL)
		MallocErr();

	strcpy(NewUsedTaxa[TIndex], TName);

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