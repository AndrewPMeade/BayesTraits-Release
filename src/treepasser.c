#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#include "genlib.h"
#include "typedef.h"
#include "treepasser.h"

#define	NEXUSTAG	"#NEXUS"
#define	TREESTAG	"TREES"
#define BEGINTAG	"BEGIN"
#define	TRANSTAG	"TRANSLATE"


void	GetPartitons(NTREES *Trees);
void	FreePartitons(NTREES* Trees);
int		IsPartitionEqual(PARTITION *PartA, PARTITION *PartB);
void	SetBinaryTree(NTREE *Tree);

int		FindNexusStart(char** TreeFile, int NoOfLines, int MaxLine)
{
	char	*Buffer;
	int		Index;
	
	Buffer = (char*)malloc(sizeof(char) * (MaxLine+1));
	if(Buffer == NULL)
		MallocErr();

	for(Index=0;Index<NoOfLines;Index++)
	{
		strcpy(Buffer, TreeFile[Index]);

		RemoveChar(';', Buffer);
		RemoveChar(' ', Buffer);
		RemoveChar('\t', Buffer);
		
		MakeUpper(Buffer);

		if(strcmp(Buffer, NEXUSTAG) == 0)
		{
			free(Buffer);
			return Index;
		}
	}

	free(Buffer);
	return -1;
}

void	FindBeginTrees(char** TreeFile, int NoOfLines, int MaxLine, int NexusStart, int *TBStart, int *TRStart)
{
	int		Index;
	char	*Buffer;
	char	**Passed;
	int		Tokes;

	Buffer	= (char*)malloc(sizeof(char) * (MaxLine+1));
	Passed	= (char**)malloc(sizeof(char*) * (MaxLine+1));
	if((Buffer == NULL) || (Passed == NULL))
		MallocErr();

	*TBStart = -1;
	*TRStart = -1;


	for(Index=NexusStart;Index<NoOfLines;Index++)
	{
		strcpy(Buffer, TreeFile[Index]);
		RemoveChar(';', Buffer);
		RemoveChar(',', Buffer);
		Tokes = MakeArgv(Buffer, Passed, MaxLine);

		if(Tokes == 2)
		{
			MakeUpper(Passed[0]);
			MakeUpper(Passed[1]);
			if((strcmp(Passed[0], BEGINTAG) == 0) && (strcmp(Passed[1], TREESTAG) == 0))
				*TBStart = Index;
		}

		if((Tokes == 1) && (*TBStart != -1))
		{
			MakeUpper(Passed[0]);
			if(strcmp(Passed[0], TRANSTAG) == 0)
			{
				free(Buffer);
				free(Passed);
				*TRStart = Index;
				return;
			}
		}
	}

	free(Buffer);
	free(Passed);
}

char*	CopyStr(char* Str)
{
	char*	Ret;

	Ret = (char*)malloc(sizeof(char) * (strlen(Str) + 1));
	if(Ret == NULL)
		MallocErr();

	strcpy(Ret, Str);

	return Ret;
}	

void	AddTaxa(NTREES* Trees, int TaxaNo, int Tokes, char **Passed, int LineSize)
{
	NTAXA*	Taxa;
	char*	Buffer;
	int		Index;

	Buffer = (char*)malloc(sizeof(char) * (LineSize + 1));
	if(Buffer == NULL)
		MallocErr();

	Taxa = &Trees->Taxa[TaxaNo];

	Taxa->No = atoi(Passed[0]);

	strcpy(Buffer, Passed[1]);

	for(Index=2;Index<Tokes;Index++)
	{
		strcat(Buffer, " ");
		strcat(Buffer, Passed[Index]);
	}

	Index = strlen(Buffer);
	if((Buffer[Index-1] == ';') || (Buffer[Index-1] == ','))
		Buffer[Index-1] = '\0';

	Taxa->Name = (char*)malloc(sizeof(char) * (strlen(Buffer) + 1));
	if(Taxa->Name == NULL)
		MallocErr();
	
	strcpy(Taxa->Name, Buffer);

	free(Buffer);
}

int		ProcessTaxa(NTREES* Trees, char** TreeFile, int NoOfLines, int MaxLine, int TransStart, int *TreeStart)
{
	int		Ret=0;
	char*	Buffer;
	char**	Passed;
	int		Index;
	int		Tokes;
	int		NoTaxa;

	Buffer = (char*)malloc(sizeof(char) * (MaxLine+1));
	Passed	= (char**)malloc(sizeof(char*) * (MaxLine+1));
	if((Buffer == NULL) || (Passed == NULL))
		MallocErr();

	NoTaxa = 0;
	for(Index=TransStart;Index<NoOfLines;Index++)
	{
		strcpy(Buffer, TreeFile[Index]);

		Tokes = MakeArgv(Buffer, Passed, MaxLine);

		if(Tokes >= 2)
		{
			if(IsValidInt(Passed[0]) == TRUE)
			{
				if(Trees != NULL)
					AddTaxa(Trees, NoTaxa, Tokes, Passed, MaxLine);

				NoTaxa++;
			}

			if(Passed[Tokes-1][strlen(Passed[Tokes-1])-1] == ';')
			{
				free(Buffer);
				free(Passed);
				(*TreeStart) = Index + 1;
				return NoTaxa;
			}
		}

		if(Tokes == 1)
		{
			if(Passed[0][0] == ';')
			{
				free(Buffer);
				free(Passed);
				(*TreeStart) = Index + 1;
				return NoTaxa;
			}
		}
	}

	free(Buffer);
	free(Passed);
	return -1;
}

void	GetTaxa(NTREES* Trees, char** TreeFile, int NoOfLines, int MaxLine, int TransStart, int *TreeStart)
{
	Trees->NoOfTaxa = ProcessTaxa(NULL, TreeFile, NoOfLines, MaxLine, TransStart, TreeStart);

	Trees->Taxa = (NTAXA*)malloc(sizeof(NTAXA) * Trees->NoOfTaxa);
	if(Trees->Taxa == NULL)
		MallocErr();

	Trees->NoOfTaxa = ProcessTaxa(Trees, TreeFile, NoOfLines, MaxLine, TransStart, TreeStart);	
}

char*	GetTaxaName(NTREES* Trees, int No)
{
	int	Index;

	for(Index=0;Index<Trees->NoOfTaxa;Index++)
	{
		if(Trees->Taxa[No].No == No)
			return Trees->Taxa[No].Name;
	}

	return NULL;
}

int		GetTaxaNo(NTREES* Trees, char* Name)
{
	int	Index;

	for(Index=0;Index<Trees->NoOfTaxa;Index++)
	{
		if(strcmp(Trees->Taxa[Index].Name, Name) == 0)
			return Trees->Taxa[Index].No;
	}

	return -1;
}


int		TaxaOK(NTREES* Trees, char** Err)
{
	int		x,y;
	
	if(Trees->NoOfTaxa < 2)
	{
		
		*Err = CopyStr("Must have more than 1 taxa");
		return FALSE;
	}

	for(x=0;x<Trees->NoOfTaxa;x++)
	{
		for(y=x+1;y<Trees->NoOfTaxa;y++)
		{
			if(Trees->Taxa[x].No == Trees->Taxa[y].No)
			{
				*Err = (char*) malloc(sizeof(char) * (BUFFERSIZE + 1));
				if(*Err == NULL)
					MallocErr();

				sprintf(*Err, "Taxa %s and %s have the same taxa number %d\n", Trees->Taxa[y].Name, Trees->Taxa[x].Name, Trees->Taxa[y].No);
				return FALSE;
			}
		}
	}

	for(x=0;x<Trees->NoOfTaxa;x++)
	{
		if(Trees->Taxa[x].No < 0)
		{
			*Err = (char*) malloc(sizeof(char) * (BUFFERSIZE + 1));
			if(*Err == NULL)
				MallocErr();

			sprintf(*Err, "Taxa %s has a taxa number (%d) less than 0", Trees->Taxa[x].Name, Trees->Taxa[x].No);
			return FALSE;
		}
	}

	return TRUE;
}

int	IsValidTreeTag(char* TreeTag)
{
	int	Index=0;

	MakeLower(TreeTag);
	
	while(strlen(TREETAGS[Index]) != 0)
	{
		if(strcmp(TreeTag, TREETAGS[Index]) == 0)
			return TRUE;
		Index++;
	}	

	return FALSE;
}

void	GetTreeName(NTREE *Tree, char** Passed)
{
	int	Index;

	Index = 1;
	
	if(strcmp(Passed[Index], "*") == 0)
		Index++;

	Tree->Tag = (char*)malloc(sizeof(char) * (strlen(Passed[Index]) + 1));
	if(Tree->Tag == NULL)
		MallocErr();

	strcpy(Tree->Tag, Passed[Index]);
}

NNODE	AllocNodeI()
{
	NNODE	Ret;
	
	Ret = (NNODE) malloc(sizeof(struct NINODE));
	if(Ret == NULL)
		MallocErr();

	Ret->Left		= NULL;
	Ret->Right		= NULL;
	Ret->Ans		= NULL;

	Ret->Tip		= FALSE;
	Ret->TaxaID		= -1;
	Ret->Taxa		= NULL;
	Ret->Length		= -1;

	Ret->NodeList	= NULL;
	Ret->NoOfNodes	= 0;

	Ret->Part		= NULL;

	return Ret;
}

NNODE	AddDescendant(NNODE N)
{
	NNODE	Ret;
	int		Index;
	NNODE*	NewList;

	Ret = AllocNodeI();
	Ret->Ans = N;

	if(N->NoOfNodes == 0)
	{
		N->NodeList		= (NNODE*)malloc(sizeof(NNODE));
		if(N->NodeList == NULL)
			MallocErr();
		N->NodeList[0]	= Ret;
		N->NoOfNodes++;
		return Ret;
	}

	NewList = (NNODE*)malloc(sizeof(struct NINODE*) * (N->NoOfNodes + 1));
	if(NewList == NULL)
		MallocErr();

	for(Index=0;Index<N->NoOfNodes;Index++)
		NewList[Index] = N->NodeList[Index];
	
	NewList[Index] = Ret;
	free(N->NodeList);
	N->NodeList = NewList;
	N->NoOfNodes++;

	return Ret;
}

double	GetLength(char *TreeString, int *Pos, char* Buffer)
{
	int		Index;
	double	Ret;
	
	Index=0;
	do
	{
		Buffer[Index] = TreeString[*Pos];
		(*Pos)++;
		Index++;
	} while((TreeString[*Pos] != ',') &&
			(TreeString[*Pos] != ')') && 
			(TreeString[*Pos] != ';'));

	Buffer[Index] = '\0';

	Ret = atof(Buffer);

	return Ret;
}

void	GetTip(NNODE Node, char *TreeString, int* Pos, char *Buffer)
{
	int		Index;

	Node->Tip = TRUE;

	Index=0;
	do
	{
		Buffer[Index] = TreeString[*Pos];
		(*Pos)++;
		Index++;
	} while((TreeString[*Pos] != ',') &&
			(TreeString[*Pos] != ')') &&
			(TreeString[*Pos] != ':'));

	Buffer[Index] = '\0';

	Node->TaxaID = atoi(Buffer);

	if(TreeString[*Pos] == ':')
	{
		(*Pos)++;
		Node->Length = GetLength(TreeString, Pos, Buffer);
	}
}

void	PassTree(NNODE CNode, char *TreeString, int* Pos, char *Buffer)
{
	NNODE	NNode;

	do 
	{
		(*Pos)++;

		NNode = AddDescendant(CNode);
		
		if(TreeString[(*Pos)] == '(')
		{
			PassTree(NNode, TreeString, Pos, Buffer);
		}
		else
		{
			GetTip(NNode, TreeString , Pos, Buffer);
		}
	} while (TreeString[*Pos] == ',') ;

	(*Pos)++;
	if(TreeString[(*Pos)] == ':')
	{
		(*Pos)++;
		CNode->Length = GetLength(TreeString, Pos, Buffer);
	}
}

void	CheckBrackets(char* TreeString, int TreeNo)
{
	int	Open;
	int Close;
	int	Index;
	int	Len;

	Open = 0;
	Close = 0;

	Len = strlen(TreeString);
	for(Index=0;Index<Len;Index++)
	{
		if(TreeString[Index] == '(')
			Open++;

		if(TreeString[Index] == ')')
			Close++;
	}

	if(Open != Close)
	{
		printf("Tree no (%d) has %d opening brackets and %d closing brackets.\n", TreeNo, Open, Close);
		printf("Brackets must be balanced\n");
		exit(0);
	}
}

int	GetTrees(NTREES *Trees, char** TreeFile, int NoOfLines, int MaxLine, int TreeStart)
{
	int		Index;
	char	*Buffer, *PBuffer;
	char**	Passed;
	int		Tokes;
	int		NoOfTrees;
	int		Pos;


	NoOfTrees = 0;

	Buffer = (char*)malloc(sizeof(char) * (MaxLine+1));
	Passed = (char**)malloc(sizeof(char*) * (MaxLine));
	PBuffer = (char*)malloc(sizeof(char) * BUFFERSIZE);

	if((Buffer == NULL) || (Passed == NULL) || (PBuffer == NULL))
		MallocErr();


	for(Index=TreeStart;Index<NoOfLines;Index++)
	{		
		strcpy(Buffer, TreeFile[Index]);
		Tokes = MakeArgv(Buffer, Passed, MaxLine);

		if(Tokes >= 3)
		{
			if(IsValidTreeTag(Passed[0]) == TRUE)
			{
				if(Trees != NULL)
				{
					GetTreeName(&Trees->Trees[NoOfTrees], Passed);

					Pos = 0;

					Trees->Trees[NoOfTrees].Root = AllocNodeI();

					CheckBrackets(Passed[Tokes-1], NoOfTrees);
					PassTree(Trees->Trees[NoOfTrees].Root, Passed[Tokes-1], &Pos, PBuffer);

					Trees->Trees[NoOfTrees].NoOfNodes = 0;
					FindNoOfNodes(Trees->Trees[NoOfTrees].Root, &Trees->Trees[NoOfTrees].NoOfNodes);
				}

				NoOfTrees++;
			}
		}
	}

	free(Buffer);
	free(Passed);
	free(PBuffer);

	return NoOfTrees;
}


NTREES*	AllocTrees()
{
	NTREES* Ret;

	Ret = (NTREES*)malloc(sizeof(NTREES));
	if(Ret == NULL)
		MallocErr();

	Ret->NoOfParts	= 0;
	Ret->PartList	= NULL;
	
	Ret->NoOfTrees	= 0;
	Ret->Taxa		= NULL;
	Ret->Trees		= NULL;

	return Ret;
}

NTAXA*	GetTaxaPtr(NTREES* Trees, int TaxaID)
{
	int	TIndex;

	for(TIndex=0;TIndex<Trees->NoOfTaxa;TIndex++)
		if(Trees->Taxa[TIndex].No == TaxaID)
			return &Trees->Taxa[TIndex];

	return NULL;
}

void	SetTaxaPtr(NTREES*	Trees)
{
	int		TIndex;
	int		NIndex;
	NTREE	*Tree;
	NNODE	N;

	for(TIndex=0;TIndex<Trees->NoOfTrees;TIndex++)
	{
		Tree = &Trees->Trees[TIndex];
		for(NIndex=0;NIndex<Tree->NoOfNodes;NIndex++)
		{
			N = Tree->NodeList[NIndex];

			if(N->Tip == TRUE)
				N->Taxa = GetTaxaPtr(Trees, N->TaxaID);
		}
	}
}

NTREES*	LoadNTrees(char* FileName, char** Err)
{
	TEXTFILE*	TreeFile;
	NTREES*		Ret;
	int			MaxLine;
	int			NexusStart;
	int			TreeBlockStart;
	int			TranslateStart;
	int			TreeStart;	
	int			Index;
	NTREE*		Tree;

	TreeFile = LoadTextFile(FileName, TRUE);

	MaxLine = TreeFile->MaxLine;

	NexusStart = FindNexusStart(TreeFile->Data, TreeFile->NoOfLines, MaxLine);
	if(NexusStart == -1)		
	{
		FreeTextFile(TreeFile);
		*Err = CopyStr("Tree file does not have a valid nexus tag.\n");
		return NULL;
	}

	FindBeginTrees(TreeFile->Data, TreeFile->NoOfLines, MaxLine, NexusStart, &TreeBlockStart, &TranslateStart);

	if(TreeBlockStart == -1)
	{
		FreeTextFile(TreeFile);
		*Err = CopyStr("Does have a Tree Block.");
		return NULL;
	}

	if(TranslateStart == -1)
	{
		FreeTextFile(TreeFile);
		*Err = CopyStr("Does have a translate line.");
		return NULL;
	}

	Ret = (NTREES*)malloc(sizeof(NTREES));
	if(Ret==NULL)
		MallocErr();

	Ret->PartList = NULL;

	GetTaxa(Ret, TreeFile->Data, TreeFile->NoOfLines, MaxLine, TranslateStart, &TreeStart);

	if(TaxaOK(Ret, Err) == FALSE)
	{
		FreeTextFile(TreeFile);
		free(Ret);
		return NULL;
	}

	Ret->NoOfTrees = GetTrees(NULL, TreeFile->Data, TreeFile->NoOfLines, MaxLine, TreeStart);

	Ret->Trees = (NTREE*) malloc(sizeof(NTREE) * Ret->NoOfTrees);
	if(Ret->Trees == NULL)
		MallocErr();

	for(Index=0;Index<Ret->NoOfTrees;Index++)
	{
		Tree = &Ret->Trees[Index];
		
		Tree->Binary	= FALSE;
		Tree->NodeList	= NULL;
		Tree->NoOfNodes	= 0;
		Tree->Tag		= NULL;
		Tree->PartList	= NULL;
	}

	GetTrees(Ret, TreeFile->Data, TreeFile->NoOfLines, MaxLine, TreeStart);

	for(Index=0;Index<Ret->NoOfTrees;Index++)
		FlattenNodes(&Ret->Trees[Index]);

	for(Index=0;Index<Ret->NoOfTrees;Index++)
	{
		Ret->Trees[Index].Binary = BinaryTree(&Ret->Trees[Index]);
		if(Ret->Trees[Index].Binary == TRUE)
			SetBinaryTree(&Ret->Trees[Index]);
	}

	SetTaxaPtr(Ret);

	FreeTextFile(TreeFile);
/*	
	for(Index=0;Index<Ret->NoOfTrees;Index++)
	{
		printf(" tree %s = ", Ret->Trees[Index].Tag);
		PrintTree(stdout, Ret->Trees[Index].Root);
		printf("\n");
	}
*/
	GetPartitons(Ret);
	return Ret;
}

void	FreeNTree(NTREE* Tree)
{
	int	Index;
	
	for(Index=0;Index<Tree->NoOfNodes;Index++)
	{
		if(Tree->NodeList[Index]->NodeList != NULL)
			free(Tree->NodeList[Index]->NodeList);
		free(Tree->NodeList[Index]);
	}

	free(Tree->NodeList);
}

void	FreeNTrees(NTREES* Trees)
{
	int	Index;

	if(Trees == NULL)
		return;

	FreePartitons(Trees);

	for(Index=0;Index<Trees->NoOfTaxa;Index++)
		free(Trees->Taxa[Index].Name);

	for(Index=0;Index<Trees->NoOfTrees;Index++)
	{
		free(Trees->Trees[Index].Tag);
		FreeNTree(&Trees->Trees[Index]);
	}

	free(Trees->Trees);
	free(Trees->Taxa);
	
	free(Trees);
}


void	FindNoOfNodes(NNODE N, int *No)
{
	int	Index;

	(*No) += N->NoOfNodes;

	for(Index=0;Index<N->NoOfNodes;Index++)
	{
		FindNoOfNodes(N->NodeList[Index], No);
	}
}

void	SetFlatNodes(NNODE N, NNODE* List, int *Pos)
{
	int	Index;

	for(Index=0;Index<N->NoOfNodes;Index++)
	{
		List[*Pos] = N->NodeList[Index];
		(*Pos)++;
		SetFlatNodes(N->NodeList[Index], List, Pos);
	}
}

void	FlattenNodes(NTREE*Tree)
{
	int	Index;

	if(Tree->NodeList != NULL)
		free(Tree->NodeList);

	Tree->NoOfNodes = 1;

	FindNoOfNodes(Tree->Root, &Tree->NoOfNodes);

	Tree->NodeList = (NNODE*)malloc(sizeof(struct NINODE*) * Tree->NoOfNodes);
	if(Tree->NodeList == NULL)
		MallocErr();

	Tree->NodeList[0] = Tree->Root;
	Index=1;

	SetFlatNodes(Tree->Root, Tree->NodeList, &Index);
}

int	FindNoOfTaxa(NNODE N)
{
	int	Index;
	int	Ret;

	if(N->Tip == TRUE)
		return 1;
	
	Ret = 0;
	for(Index=0;Index<N->NoOfNodes;Index++)
		Ret += FindNoOfTaxa(N->NodeList[Index]);
	
	return Ret;
}

void	SetPartitionTaxa(NNODE N, int* List, int *No)
{
	int	Index;

	if(N->Tip == TRUE)
	{
		List[(*No)] = N->TaxaID;
		(*No)++;
		return;
	}

	for(Index=0;Index<N->NoOfNodes;Index++)
		SetPartitionTaxa(N->NodeList[Index], List, No);
}

int	TaxaIDComp(const void *Pv1, const void *Pv2)
{
	int *P1, *P2;

	P1 = (int*)Pv1;
	P2 = (int*)Pv2;

	if(*P1 > *P2)
		return 1;

	if(*P2 > *P1)
		return -1;

	return 0;
}

void	PrintPartition(NTREES* Trees, PARTITION *Part)
{
	int	Index;

	printf("%d\t%f\t%d\t", Part->Freq, Part->Prob, Part->No);
	for(Index=0;Index<Part->No;Index++)
		printf("\t%s", GetTaxaName(Trees, Part->TaxaNo[Index])); 
	printf("\n");

	fflush(stdout);
}

void	PrintTreePartition(NTREES* Trees, NTREE* Tree)
{
	int	Index;

	printf("%s\n", Tree->Tag);

	for(Index=0;Index<Tree->NoOfParts;Index++)
		PrintPartition(Trees, Tree->PartList[Index]);
}

PARTITION* GetPartition(NNODE N)
{
	PARTITION*	Ret;
	int			No;
	
	Ret = (PARTITION*)malloc(sizeof(PARTITION));
	if(Ret == NULL)
		MallocErr();

	Ret->Freq	= 0;
	Ret->Prob	= 0;

	Ret->No		= FindNoOfTaxa(N);

	Ret->TaxaNo	= (int*)malloc(sizeof(int) * Ret->No);
	if(Ret->TaxaNo == NULL)
		MallocErr();
	memset(Ret->TaxaNo, 0, sizeof(int)*Ret->No);

	No = 0;
	SetPartitionTaxa(N, Ret->TaxaNo, &No);

	qsort(Ret->TaxaNo, Ret->No, sizeof(int), TaxaIDComp);

	return Ret;
}

void	GetTreePartitions(NTREE* Tree)
{
	int		Index;
	int		PIndex;
	NNODE	N;

	Tree->NoOfParts = 0;
	for(Index=0;Index<Tree->NoOfNodes;Index++)
	{
		N = Tree->NodeList[Index];

		if(N->Tip == FALSE)
		{
			N->Part = GetPartition(N);
			Tree->NoOfParts++;
		}
	}

	Tree->PartList = (PARTITION**)malloc(sizeof(PARTITION*) * Tree->NoOfParts);
	if(Tree->PartList == NULL)
		MallocErr();

	PIndex=0;
	for(Index=0;Index<Tree->NoOfNodes;Index++)
	{
		N = Tree->NodeList[Index];
		if(N->Tip == FALSE)
		{
			Tree->PartList[PIndex] = N->Part;
			PIndex++;
		}
	}

	fflush(stdout);
}


void	AddPartitonsToList(PARTITION**	PList, int *No, PARTITION* Part)
{
	int	Index;

	for(Index=0;Index<*No;Index++)
	{
		if(IsPartitionEqual(PList[Index], Part) == TRUE)
		{
			PList[Index]->Freq++;
			return;
		}
	}

	PList[*No] = Part;
	PList[*No]->Freq = 1;

	(*No)++;
}

//int		CompPartFreq(PARTITION **Pv1, PARTITION **P2)
int		CompPartFreq(const void *Pv1, const void *Pv2)
{
	PARTITION **P1, **P2;

	P1 = (PARTITION**)Pv1;
	P2 = (PARTITION**)Pv2;


	if((*P1)->Prob>= (*P2)->Prob)
		return -1;

	if((*P1)->Prob < (*P2)->Prob)
		return 1;

	return 0;
}

void	GetPartitons(NTREES *Trees)
{
	int			TIndex;
	int			MaxParts;
	PARTITION**	TempList;
	int			PIndex;
	int			Index;
	

	MaxParts = 0;
	for(TIndex=0;TIndex<Trees->NoOfTrees;TIndex++)
	{
		GetTreePartitions(&Trees->Trees[TIndex]);
		MaxParts += Trees->Trees[TIndex].NoOfParts;
	}

	TempList = (PARTITION**)malloc(sizeof(PARTITION*) * MaxParts);
	if(TempList == NULL)
		MallocErr();

	PIndex = 0;
	for(TIndex=0;TIndex<Trees->NoOfTrees;TIndex++)
	{
		for(Index=0;Index<Trees->Trees[TIndex].NoOfParts;Index++)
			AddPartitonsToList(TempList, &PIndex, Trees->Trees[TIndex].PartList[Index]);
	}

	Trees->NoOfParts = PIndex;

	Trees->PartList = (PARTITION**)malloc(sizeof(PARTITION*) * Trees->NoOfParts);
	if(Trees->PartList == NULL)
		MallocErr();

	for(PIndex=0;PIndex<Trees->NoOfParts;PIndex++)
	{
		Trees->PartList[PIndex] = TempList[PIndex];
		Trees->PartList[PIndex]->Prob = (double)Trees->PartList[PIndex]->Freq / (double)Trees->NoOfTrees;
	}

	qsort(Trees->PartList, Trees->NoOfParts, sizeof(PARTITION*), CompPartFreq);

	free(TempList);
}

void	FreePartitons(NTREES* Trees)
{
	int		TIndex;
	int		NIndex;
	NTREE*	Tree;
	NNODE	N;

	if(Trees->PartList != NULL)
		free(Trees->PartList);
	Trees->PartList = NULL;

	for(TIndex=0;TIndex<Trees->NoOfTrees;TIndex++)
	{
		Tree = &Trees->Trees[TIndex];

		if(Tree->PartList != NULL)
		{
			free(Tree->PartList);
			Tree->PartList = NULL;
		}

		for(NIndex=0;NIndex<Tree->NoOfNodes;NIndex++)
		{
			N = Tree->NodeList[NIndex];
			if(N->Part != NULL)
			{
				free(N->Part->TaxaNo);
				free(N->Part);
				N->Part = NULL;
			}
		}
	}
}

void	ResolveNode(NNODE N, double Len)
{
	NNODE	NewNode;
	int		Index;
	NNODE*	TempList;

	NewNode = (NNODE)malloc(sizeof(struct NINODE));
	if(NewNode == NULL)
		MallocErr();

	NewNode->Tip  = FALSE;
	NewNode->Part = NULL;

	NewNode->NodeList = (NNODE*)malloc(sizeof(struct NINODE*) * (N->NoOfNodes - 1));
	if(NewNode->NodeList == NULL)
		MallocErr();
	NewNode->NoOfNodes = N->NoOfNodes - 1;

	for(Index=1;Index<N->NoOfNodes;Index++)
		NewNode->NodeList[Index - 1] = N->NodeList[Index];

	TempList = (NNODE*)malloc(sizeof(struct NINODE*) * 2);
	if(TempList == NULL)
		MallocErr();

	TempList[0] = N->NodeList[0];
	TempList[1] = NewNode;

	free(N->NodeList);
	N->NodeList = TempList;
	N->NoOfNodes= 2;

	NewNode->Ans = N;
	NewNode->Length = Len;

	NewNode->Part = GetPartition(NewNode);
}

void	ResolveTreeDet(NNODE N, double Len)
{
	int	I;

	if(N->NoOfNodes > 2)
		ResolveNode(N, Len);

	if(N->Tip == FALSE)
	{
		for(I=0;I<N->NoOfNodes;I++)
			ResolveTreeDet(N->NodeList[I], Len);

		N->Right = N->NodeList[0];
		N->Left	 = N->NodeList[1];


		ResolveTreeDet(N->Right, Len);
		ResolveTreeDet(N->Left, Len);
	}
	else
	{
		N->Right= NULL;
		N->Left = NULL;
	}
}

void	ResolveTreesDet(NTREES *Trees, double Len)
{
	int TIndex;	

	for(TIndex=0;TIndex<Trees->NoOfTrees;TIndex++)
	{
		ResolveTreeDet(Trees->Trees[TIndex].Root, Len);
		FlattenNodes(&Trees->Trees[TIndex]);
	}
}

NNODE	AllocPartNode(PARTITION* Part)
{
	NNODE	Ret;

	Ret = (NNODE)malloc(sizeof(struct NINODE));
	if(Ret == NULL)
		MallocErr();

	Ret->Part = Part;


	Ret->Ans		= NULL;
	Ret->Left		= NULL;
	Ret->Right		= NULL;

	Ret->Length		= -1;

	Ret->NodeList	= NULL;
	Ret->NoOfNodes	= 0;

	Ret->Tip		= FALSE;

	return Ret;
}

void	AddNewPartNode(NNODE Node, NNODE New)
{
	NNODE*	NewList;
	int		Index;

	NewList = (NNODE*)malloc(sizeof(struct NINODE*) * (Node->NoOfNodes + 1));
	if(NewList == NULL)
		MallocErr();

	New->Ans = Node;
	NewList[0] = New;

	for(Index=1;Index<Node->NoOfNodes+1;Index++)
		NewList[Index] = Node->NodeList[Index-1];

	free(Node->NodeList);
	Node->NodeList = NewList;

	
	Node->NoOfNodes++;
}

int		IsTaxaInPartition(PARTITION *Part, int TaxaID)
{
	int	Index;

	for(Index=0;Index<Part->No;Index++)
	{
		if(Part->TaxaNo[Index] == TaxaID)
			return TRUE;
	}
	
	return FALSE;
}

int		IsPartitionSubSet(PARTITION *Part, PARTITION *SubSet)
{
	int	Index;

	if(SubSet->No > Part->No)
		return FALSE;

	for(Index=0;Index<SubSet->No;Index++)
	{
		if(IsTaxaInPartition(Part, SubSet->TaxaNo[Index]) == FALSE)
			return FALSE;
	}

	return TRUE;
}

int		IsPartitionEqual(PARTITION *PartA, PARTITION *PartB)
{
	int	Index;

	if(PartA->No != PartB->No)
		return FALSE;

	for(Index=0;Index<PartA->No;Index++)
	{
		if(PartA->TaxaNo[Index] != PartB->TaxaNo[Index])
			return FALSE;
	}

	return TRUE;
}

void	PrintNTree(FILE* Str, NNODE N)
{
	int	Index;

	if(N->Tip == TRUE)
	{
		fprintf(Str, "%d", N->TaxaID);
		if(N->Length != -1)
			fprintf(Str, ":%f", N->Length);
		return;
	}

	fprintf(Str, "(");

	for(Index=0;Index<N->NoOfNodes;Index++)
	{
		if(N->NodeList[Index]->Tip == FALSE)
		{
			PrintNTree(Str, N->NodeList[Index]);
		}
		else
		{
			fprintf(Str, "%d", N->NodeList[Index]->TaxaID);
			if(N->NodeList[Index]->Length != -1)
				fprintf(Str, ":%f", N->NodeList[Index]->Length);
		}
		if(Index != N->NoOfNodes - 1)
			fprintf(Str, ",");
	}

	if((N->NoOfNodes == 0) && (N->Tip == FALSE))
	{
		for(Index=0;Index<N->Part->No;Index++)
		{
			printf("%d", N->Part->TaxaNo[Index]);
			if(Index!=N->Part->No-1)
				printf(",");
		}
	}

	fprintf(Str, ")");

	if(N->Length != -1)
		fprintf(Str, ":%f", N->Length);
}

void	PrintNTreeHeadder(FILE* Str, NTREES *Trees)
{
	int	Index;

	fprintf(Str, "#NEXUS\n");
	fprintf(Str, "begin trees;\n");
	fprintf(Str, "\ttranslate\n");

	for(Index=0;Index<Trees->NoOfTaxa-1;Index++)
		fprintf(Str, "\t\t%d %s,\n", Trees->Taxa[Index].No, Trees->Taxa[Index].Name);
	fprintf(Str, "\t\t%d %s;\n", Trees->Taxa[Index].No, Trees->Taxa[Index].Name);
}

void	PrintNTrees(FILE* Str, NTREES *Trees)
{
	int		Index;
	char	*Buffer;
	int		MaxDig;
	char	*TNo;

	Buffer = (char*)malloc(sizeof(char) * BUFFERSIZE);
	if(Buffer == NULL)
		MallocErr();

	sprintf(Buffer, "%d", Trees->NoOfTrees);
	MaxDig = strlen(Buffer);
	free(Buffer);

	PrintNTreeHeadder(Str, Trees);

	for(Index=0;Index<Trees->NoOfTrees;Index++)
	{
		if(Trees->Trees[Index].Tag != NULL)
			fprintf(Str, "\t\ttree %s = ", Trees->Trees[Index].Tag);
		else
		{
			TNo = FormatInt(Index+1, MaxDig);
			fprintf(Str, "\t\ttree tree.%s = ", TNo);
			free(TNo);
		}

		PrintNTree(Str, Trees->Trees[Index].Root);

		fprintf(Str, ";\n");
	}

	fprintf(Str, "end;\n");
}

void	SaveNTrees(char* FileName, NTREES *Trees)
{
	FILE*	OutFile;

	OutFile = OpenWrite(FileName);
	PrintNTrees(OutFile, Trees);
	fclose(OutFile);
}

void	SetBinaryTree(NTREE *Tree)
{
	int		Index;
	NNODE	N;

	for(Index=0;Index<Tree->NoOfNodes;Index++)
	{
		N = Tree->NodeList[Index]; 
		if(N->Tip == FALSE)
		{
			N->Right = N->NodeList[0];
			N->Left  = N->NodeList[1];
		}
	}
}

int		BinaryTree(NTREE *Tree)
{
	int Index;

	for(Index=0;Index<Tree->NoOfNodes;Index++)
	{
		if(Tree->NodeList[Index]->NoOfNodes > 2)
			return FALSE;
	}

	return TRUE;
}

double	GetTreeLength(NTREES* Trees, int TreeNo)
{
	double	Ret=0;
	int		Index;
	NTREE*	Tree;
	NNODE	Node;

	Tree = &Trees->Trees[TreeNo];
	for(Index=0;Index<Tree->NoOfNodes;Index++)
	{
		Node = Tree->NodeList[Index];
		if(Node->Length > 0)
			Ret += Node->Length;
	}

	return Ret;
}

int		GetTaxaID(NTREES* Trees, char *TaxaName)
{
	int	Index;

	for(Index=0;Index<Trees->NoOfTaxa;Index++)
	{
		if(strcmp(Trees->Taxa[Index].Name, TaxaName) == 0)
			return Trees->Taxa[Index].No;
	}

	return -1;
}


int		NodesBelow(NNODE N)
{
	if(N->Tip == TRUE)
		return 1;

	return (NodesBelow(N->Left) + NodesBelow(N->Right));

}

void	MClockSetBL(NNODE N)
{
	N->Length = NodesBelow(N->Ans) - NodesBelow(N);

	if(N->Tip == FALSE)
	{
		MClockSetBL(N->Left);
		MClockSetBL(N->Right);
	}
}

void	MakeUltrametric(NTREE* Tree, double Length)
{
	int		i;
	double	PathLen;

	MClockSetBL(Tree->Root->Left);
	MClockSetBL(Tree->Root->Right);

	i = 0;
	while(Tree->NodeList[i]->Tip == FALSE)
		i++;

	PathLen = GetPathLength(Tree->NodeList[i]);

	PathLen = Length / PathLen;

	for(i=0;i<Tree->NoOfNodes;i++)
		Tree->NodeList[i]->Length *= PathLen;
}


double	GetPathLength(NNODE Node)
{
	double	Ret;

	Ret = 0;

	while(Node != NULL)
	{
		Ret += Node->Length;
		Node = Node->Ans;
	}

	return Ret;
}


double	AveBL(NTREE* Tree)
{
	double	Sum;
	int		Index;

	Sum = 0;
	for(Index=0;Index<Tree->NoOfNodes;Index++)
		Sum += Tree->NodeList[Index]->Length;
	
	Sum = Sum / Tree->NoOfNodes;

	return Sum;
}
