#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>

#include "FatTail.h"

#include "typedef.h"

#include "genlib.h"
#include "StableDist.h"
#include "likelihood.h"
#include "part.h"
#include "trees.h"
#include "praxis.h"
#include "Geo.h"
#include "SliceSampler.h"

#define NO_SLICE_STEPS 100
#define	MAX_STEP_DIFF	5.0

double	NodeSliceSampler(NODE N, RANDSTATES *RS);
void	SetInitAnsStates(OPTIONS *Opt, TREES *Trees, TREE *Tree);

//void MapRatesToTree(TREE *Tree, int NoSites, FATTAILRATES *FTR)
void	FatTailSetAnsSates(TREE *Tree, int NoSites, FATTAILRATES *FTR)
{
	size_t Size;
	
	Size = sizeof(double) * Tree->NoNodes * NoSites;

	memcpy(Tree->FatTailTree->AnsVect, FTR->AnsVect, Size);
}

void FatTailGetAnsSates(TREE *Tree, int NoSites, FATTAILRATES *FTR)
{
	size_t Size;
	
	Size = sizeof(double) * Tree->NoNodes * NoSites;

	memcpy(FTR->AnsVect, Tree->FatTailTree->AnsVect, Size);
}

void MapRatesToFatTailRate(RATES *Rates, FATTAILRATES *FatTailRates)
{
	int Index, Pos;

	Pos = 0;
	for(Index=0;Index<FatTailRates->NoSD;Index++)
	{
		FatTailRates->Alpha[Index] = Rates->Rates[Pos++];
		FatTailRates->Scale[Index] = Rates->Rates[Pos++];
	}
}


void MapFatTailRateToRates(RATES *Rates, FATTAILRATES *FatTailRates)
{
	int Index, Pos;

	Pos = 0;
	for(Index=0;Index<FatTailRates->NoSD;Index++)
	{
		Rates->Rates[Pos++] = FatTailRates->Alpha[Index];
		Rates->Rates[Pos++] = FatTailRates->Scale[Index];
	}
}

FATTAILRATES*	AllocFatTailRates(OPTIONS *Opt, TREES *Trees)
{
	FATTAILRATES* Ret;
	int Index, NoSites;

	NoSites = Opt->Trees->NoOfSites;

	Ret = (FATTAILRATES*)malloc(sizeof(FATTAILRATES));
	if(Ret == NULL)
		MallocErr();

	Ret->Alpha = (double*)malloc(sizeof(double) * NoSites);
	Ret->Scale = (double*)malloc(sizeof(double) * NoSites);
	Ret->SiteLh = (double*)malloc(sizeof(double) * NoSites);

	if( (Ret->Alpha == NULL) ||
		(Ret->Scale == NULL) ||
		(Ret->SiteLh == NULL))
		MallocErr();

	Ret->SiteMin = (double*)malloc(sizeof(double) * NoSites);
	Ret->SiteMax = (double*)malloc(sizeof(double) * NoSites);
	Ret->SiteSD  = (double*)malloc(sizeof(double) * NoSites);

	if( (Ret->SiteMin == NULL) ||
		(Ret->SiteMax == NULL) ||
		(Ret->SiteSD  == NULL))
		MallocErr();

	Ret->SliceSampler = CrateSliceSampler(NO_SLICE_STEPS);


	Ret->AnsVect = (double*)malloc(sizeof(double) * NoSites * Trees->MaxNodes);
	if(Ret->AnsVect == NULL)
		MallocErr();	

	Ret->NoSD = NoSites;
	if(Opt->UseGeoData == TRUE)
		Ret->NoSD = 1;

	Ret->SDList = (STABLEDIST**)malloc(sizeof(double) * Ret->NoSD);
	if(Ret->SDList == NULL)
		MallocErr();
	for(Index=0;Index<Ret->NoSD;Index++)
		Ret->SDList[Index] = CreatStableDist();
	
	return Ret;
}

void			GetSiteInfo(int SiteNo, TREES *Trees, FATTAILRATES* FTR)
{
	int Index;
	double Data, Mean, SD;


	FTR->SiteMax[SiteNo] = Trees->Taxa[0]->ConData[SiteNo];
	FTR->SiteMin[SiteNo] = Trees->Taxa[0]->ConData[SiteNo];

	Mean = 0;

	for(Index=0;Index<Trees->NoOfTaxa;Index++)
	{
		Data = Trees->Taxa[Index]->ConData[SiteNo];
				
		Mean += Data;

		if(Data > FTR->SiteMax[SiteNo])
			FTR->SiteMax[SiteNo] = Data;

		if(Data < FTR->SiteMin[SiteNo])
			FTR->SiteMin[SiteNo] = Data;
	}

	Mean = Mean / Trees->NoOfTaxa;

	SD = 0;
	for(Index=0;Index<Trees->NoOfTaxa;Index++)
	{
		Data = Trees->Taxa[Index]->ConData[SiteNo];
		SD += (Data - Mean) * (Data - Mean);
	}
	SD = SD / Trees->NoOfTaxa;
	SD = sqrt(SD);

	FTR->SiteSD[SiteNo] = SD;
}

FATTAILRATES*	CreateFatTailRates(OPTIONS *Opt, TREES *Trees)
{
	FATTAILRATES* Ret;
	int	Index;

	Ret = AllocFatTailRates(Opt, Trees);

	for(Index=0;Index<Ret->NoSD;Index++)
	{
		if(Opt->FatTailNormal == FALSE)
			Ret->Alpha[Index] = 0.5;
		else
			Ret->Alpha[Index] = 2.0;

		Ret->Scale[Index] = 0.5;
	}

	for(Index=0;Index<Opt->Trees->NoOfSites;Index++)
		GetSiteInfo(Index, Trees, Ret);

	SetInitAnsStates(Opt, Trees, Trees->Tree[0]);
	FatTailGetAnsSates(Trees->Tree[0], Trees->NoOfSites, Ret);
			
	return Ret;
}

void	FreeFatTailRates(FATTAILRATES* FTR, int NoSites)
{
	int Index;

	free(FTR->Alpha);
	free(FTR->Scale);
	free(FTR->SiteLh);

	free(FTR->SiteMin);
	free(FTR->SiteMax);
	free(FTR->SiteSD);

//	FreeSliceSampler(FTR->SliceSampler);

	for(Index=0;Index<FTR->NoSD;Index++)
		FreeStableDist(FTR->SDList[Index]);
	free(FTR->SDList);
		
	free(FTR->AnsVect);
	
	free(FTR);
}

void			CopyFatTailRates(TREES *Trees, FATTAILRATES *A, FATTAILRATES *B)
{
	// Alpha and Scale should be mapped from globabal rates
	
	memcpy(A->AnsVect, B->AnsVect, sizeof(double) * Trees->MaxNodes * Trees->NoOfSites);
}

FATTAILNODE*	InitFatTailNode(int NoSites, NODE N, double *AnsVect)
{
	FATTAILNODE*	Ret;
	int			Index;
	size_t		Pos;

	Ret = (FATTAILNODE*)malloc(sizeof(FATTAILNODE));
	if(Ret == NULL)
		MallocErr();
	
	Pos = N->ID * NoSites;
	Ret->Ans = &AnsVect[Pos];

	for(Index=0;Index<NoSites;Index++)
	{
		Ret->Ans[Index] = 0;
		if(N->Tip == TRUE)
			Ret->Ans[Index] = N->Taxa->ConData[Index];
	}

	return Ret;
}

void	InitFatTailTree(OPTIONS *Opt, TREE *Tree)
{
	int Index;
	NODE N;

	for(Index=0;Index<Tree->NoNodes;Index++)
	{
		N = Tree->NodeList[Index];
		N->FatTailNode = InitFatTailNode(Opt->Trees->NoOfSites, N, Tree->FatTailTree->AnsVect);
	}
}

FATTAILTREE*	AllocFatTailTree(TREE *Tree, int NoOfSites)
{
	FATTAILTREE* Ret;

	Ret = (FATTAILTREE*)malloc(sizeof(FATTAILTREE));
	if(Ret == NULL)
		MallocErr();

	Ret->AnsVect = (double*)malloc(sizeof(double) * Tree->NoNodes * NoOfSites);
	if(Ret->AnsVect == NULL)
		MallocErr();

	return Ret;
}

void			FreeFatTailTree(FATTAILTREE *FatTailTree)
{
	free(FatTailTree->AnsVect);
	free(FatTailTree);
}

void	SetInitAnsStateNodes(int SiteNo, TREE *Tree)
{
	NODE N;
	int NIndex;
	FATTAILNODE *FTN;

	for(NIndex=0;NIndex<Tree->NoNodes;NIndex++)
	{
		N = Tree->NodeList[NIndex];
		FTN = N->FatTailNode;

		if(N->Tip == TRUE)
			FTN->Data = FTN->Ans[SiteNo];

		FTN->Cont = FTN->Err = FTN->Var = FTN->v = 0.0;
	}
}

void	GetInitAnsStateNodes(int SiteNo, TREE *Tree)
{
	NODE N;
	int NIndex;
	FATTAILNODE *FTN;

	for(NIndex=0;NIndex<Tree->NoNodes;NIndex++)
	{
		N = Tree->NodeList[NIndex];
		FTN = N->FatTailNode;
		FTN->Ans[SiteNo] = FTN->Data;
	}
}

void	SetContrastAnsStates(NODE N)
{
	int Index;
	NODE	N0, N1;
	FATTAILNODE *C, *C0, *C1;
	double	l0, l1, t;

	C = N->FatTailNode;


	if(N->Tip == TRUE)
		return;

	for(Index=0;Index<N->NoNodes;Index++)
		SetContrastAnsStates(N->NodeList[Index]);
/*
	if(N->NoNodes != 2)
	{
		printf("Code is only for biforcating tree.\n");
		exit(0);
	}
*/
	N0 = N->NodeList[0];
	N1 = N->NodeList[1];

	C0 = N0->FatTailNode;
	C1 = N1->FatTailNode;

	l0 = N0->Length + C0->Err;
	l1 = N1->Length + C1->Err;

	t = (l0 * C1->Data) +  (l1 * C0->Data);
	t = t / (l0 + l1);
		
	C->Data = t;
	C->Cont = C0->Data - C1->Data;

	C->Err = (l0 * l1) / (l0 + l1);
	C->Var = l0 + l1;

	C->v = C->Err;
	if(N->Length > 0)
		C->v += N->Length;
}


void	SetInitAnsStates(OPTIONS *Opt, TREES *Trees, TREE *Tree)
{
	int SIndex;

	for(SIndex=0;SIndex<Trees->NoOfSites;SIndex++)
	{
		SetInitAnsStateNodes(SIndex, Tree);
		SetContrastAnsStates(Tree->Root);
		GetInitAnsStateNodes(SIndex, Tree);
	}

	if(Opt->UseGeoData == TRUE)
		CorrectIntGeoNodes(Tree);
}

void	InitFatTailTrees(OPTIONS *Opt, TREES *Trees)
{
	int Index;
	TREE *Tree;

	for(Index=0;Index<Trees->NoOfTrees;Index++)
	{
		Tree = Trees->Tree[Index];
		Tree->FatTailTree = AllocFatTailTree(Tree, Trees->NoOfSites);
		InitFatTailTree(Opt, Tree);
	}
}


double	CalcNodeStableLh(NODE N, int NoSites, STABLEDIST **SDList, int UseGeoData)
{
	int Index, SIndex;
	double Ret;
	double L, x;
	STABLEDIST *SD;
	
	if(N->Tip == TRUE)
		return 0;

	Ret = 0;
	SD = SDList[0];
	for(SIndex=0;SIndex<NoSites;SIndex++)
	{
		if(UseGeoData == FALSE)
			SD = SDList[SIndex];

		for(Index=0;Index<N->NoNodes;Index++)
		{
			x = N->FatTailNode->Ans[SIndex]- N->NodeList[Index]->FatTailNode->Ans[SIndex];
			
			L = StableDistTPDF(SD, x , N->NodeList[Index]->Length);
			Ret += L;
		}
	}

	return Ret;
}

double	CalcTreeStableLh(OPTIONS *Opt, TREES *Trees, RATES *Rates)
{
	int NoSites, Index;
	double Ret;
	FATTAILRATES *FTR;
	TREE *Tree;

	Tree = Trees->Tree[Rates->TreeNo];
	NoSites = Trees->NoOfSites;
	FTR = Rates->FatTailRates;

	FatTailSetAnsSates(Tree, NoSites, FTR);

	for(Index=0;Index<FTR->NoSD;Index++)
		SetStableDist(FTR->SDList[Index], FTR->Alpha[Index], FTR->Scale[Index]);
	
	Ret = 0;

	#pragma omp parallel for num_threads(Opt->Cores) reduction(+:Ret)
	for(Index=0;Index<Tree->NoNodes;Index++)
	{
		if(Tree->NodeList[Index]->Tip == FALSE)
			Ret += CalcNodeStableLh(Tree->NodeList[Index], NoSites, FTR->SDList, Opt->UseGeoData);
	}

	if(Ret != Ret || Ret == Ret + 1.0)
		return ERRLH;

	return Ret;
}

NODE	GetSliceSampleNode(TREE *Tree, RANDSTATES *RS)
{
	NODE Ret;
	int Pos;

	do
	{
		Pos = RandUSInt(RS) % Tree->NoNodes;
		Ret = Tree->NodeList[Pos];
	}while(Ret->Tip == TRUE);

	return Ret;
}

double	AnsStateLh(double X, int SiteNo, NODE N, STABLEDIST *SD)
{
	double Ret, Val;
	int Index;

	Ret = 0;

	for(Index=0;Index<N->NoNodes;Index++)
	{
		Val = X - N->NodeList[Index]->FatTailNode->Ans[SiteNo];
		Ret += StableDistTPDF(SD, Val, N->NodeList[Index]->Length);
	}

	if(N->Ans != NULL)
	{
		Val = X - N->Ans->FatTailNode->Ans[SiteNo];
		Ret += StableDistTPDF(SD, Val, N->Length);
	}

	return Ret;
}

void	SetTestVect(NODE N, int SiteNo, STABLEDIST *SD, int NoSteps, double *XVect, double *YVect)
{
	int Index;

	for(Index=0;Index<NoSteps;Index++)
		YVect[Index] = 0.0;

	for(Index=100;Index<200;Index++)
		YVect[Index] = 10.0;

	for(Index=500;Index<600;Index++)
		YVect[Index] = 20;
}

int		FatTailSetYPosVect(SLICESAMPLER *SS, OPTIONS *Opt, NODE N, int SiteNo, STABLEDIST *SD)
{
	int Index; 
	
	#pragma omp parallel for num_threads(Opt->Cores) 
	for(Index=0;Index<SS->NoSteps;Index++)
		SS->SliceY[Index] = AnsStateLh(SS->SliceX[Index], SiteNo, N, SD);

	for(Index=0;Index<SS->NoSteps;Index++)
	{
		if(Index > 1)
		{
			if(SS->SliceY[Index] - SS->SliceY[Index-1] > MAX_STEP_DIFF)
				return FALSE;

			if(SS->SliceY[Index] == SS->SliceY[0])
				return FALSE;
		}

		if(SS->SliceY[Index] != SS->SliceY[Index] || SS->SliceY[Index] == SS->SliceY[Index] + 1.0)
			return FALSE;
	}

	return TRUE;
}

void	GetSiteMinMax(FATTAILRATES *FTR, int SiteNo, double *Min, double *Max)
{
	*Min = FTR->SiteMin[SiteNo] - (FTR->SiteMin[SiteNo] * 0.1);
	*Max = FTR->SiteMax[SiteNo] + (FTR->SiteMax[SiteNo] * 0.1);
}



void	PrintAnsVect(TREE *Tree)
{
	int Index;
	NODE N;

	for(Index=0;Index<Tree->NoNodes;Index++)
	{
		N = Tree->NodeList[Index];
		printf("Node:\t%d\t%d\t%f\t%f\n", Index,N->Tip, N->FatTailNode->Ans[0], N->FatTailNode->Ans[1]);
	}

	exit(0);
}

void	PrintAnsStates(TREES *Trees, NODE N)
{
	int Index;

	if(N->Tip == FALSE)
		for(Index=0;Index<N->NoNodes;Index++)
			PrintAnsStates(Trees, N->NodeList[Index]);


	printf("Ans =\t");
	for(Index=0;Index<Trees->NoOfSites;Index++)
		printf("%f\t", N->FatTailNode->Ans[Index]);

	PrintPart(stdout, Trees, N->Part);


	printf("\n");
}

void	SetNodeAns(TREES *Trees, NODE N, double Val)
{
	int Index;

	if(N->Tip == TRUE)
		return;

	for(Index=0;Index<N->NoNodes;Index++)
		SetNodeAns(Trees, N->NodeList[Index], Val);

	N->FatTailNode->Ans[0] = Val;
	N->FatTailNode->Ans[1] = Val;
}

void	TestMapping(OPTIONS *Opt, TREES *Trees, RATES *Rates)
{
	TREE *Tree;
	FATTAILRATES *FTR;
	double Lh;
	
	FTR = Rates->FatTailRates;

	Tree = Trees->Tree[Rates->TreeNo];

	PrintAnsStates(Trees, Tree->Root);

	

	Lh = Likelihood(Rates, Trees, Opt);
	printf("Lh:\t%f\n", Lh);

	FatTailSetAnsSates(Tree, Trees->NoOfSites, FTR);
	
	SetNodeAns(Trees, Tree->Root, 0.0);
	
	FatTailGetAnsSates(Tree, Trees->NoOfSites, FTR);
	Lh = Likelihood(Rates, Trees, Opt);
	printf("Lh:\t%f\n", Lh);
		PrintAnsStates(Trees, Tree->Root);

	exit(0);
}


void	SliceSampleFatTail(OPTIONS *Opt, TREES *Trees, RATES *Rates)
{
	TREE *Tree;
	int TSite, Valid;
	NODE TNode;
	double NAns, CAns, PCAns, Min, Max;
	FATTAILRATES *FTR;
	
//	TestMapping(Opt, Trees, Rates);

//	return;

	FTR = Rates->FatTailRates;

	Tree = Trees->Tree[Rates->TreeNo];

	FatTailSetAnsSates(Tree, Trees->NoOfSites, FTR);

//	PrintAnsVect(Tree);
	
	TSite = RandUSInt(Rates->RS) % Trees->NoOfSites;
	TNode = GetSliceSampleNode(Tree, Rates->RS);
	
	SetStableDist(FTR->SDList[TSite], FTR->Alpha[TSite], FTR->Scale[TSite]);
	
	CAns = TNode->FatTailNode->Ans[TSite];
	PCAns = AnsStateLh(CAns, TSite, TNode, FTR->SDList[TSite]);
	
	GetSiteMinMax(FTR, TSite, &Min, &Max);

//	TNode->FatTailNode->Ans[TSite] = Min + (RandDouble(Rates->RS) * (Max - Min));
//	MapTreeToRates(Tree, Trees->NoOfSites, FTR);

	SSSetXPosVect(FTR->SliceSampler, Min, Max); 
	
	Valid = FatTailSetYPosVect(FTR->SliceSampler, Opt, TNode, TSite, FTR->SDList[TSite]);

	if(Valid == FALSE)
		NAns = Min + (RandDouble(Rates->RS) * (Max - Min));
	else
		NAns = SSGetNewPoint(FTR->SliceSampler, Rates->RS, PCAns);

	TNode->FatTailNode->Ans[TSite] = NAns;
	
	FatTailGetAnsSates(Tree, Trees->NoOfSites, FTR);
	
	return;
}

void	NodeSliceSampleFatTail(NODE N, int SiteNo, OPTIONS *Opt, TREES *Trees, RATES *Rates) 
{
	int Changed, Valid;
	double CLh, CAns, NAns, NLh, Min, Max;
	FATTAILRATES *FTR;

	FTR = Rates->FatTailRates;
	
	CAns = N->FatTailNode->Ans[SiteNo];
	CLh = AnsStateLh(CAns, SiteNo, N, FTR->SDList[SiteNo]);

	FTR = Rates->FatTailRates;
	
	GetSiteMinMax(FTR, SiteNo, &Min, &Max);

	SSSetXPosVect(FTR->SliceSampler, Min, Max); 
	
	Valid = FatTailSetYPosVect(FTR->SliceSampler, Opt, N, SiteNo, FTR->SDList[SiteNo]);
	
	Changed = FALSE;
	do
	{
		if(Valid == FALSE)
			NAns = Min + (RandDouble(Rates->RS) * (Max - Min));
		else
			NAns = SSGetNewPoint(FTR->SliceSampler, Rates->RS, CLh);

		NLh = AnsStateLh(NAns, SiteNo, N, FTR->SDList[SiteNo]);

		if(log(RandDouble(Rates->RS)) < (NLh - CLh))
		{
		//	printf("New:\t%f\tOld\t%f\tDiffer\t%f\t%f\t%f\n", NLh, CLh, NLh - CLh, NAns, CAns);fflush(stdout);
			Changed = TRUE;
		}
	} while(Changed == FALSE);

	N->FatTailNode->Ans[SiteNo] = NAns;
}


void	AllSliceSampleFatTail(OPTIONS *Opt, TREES *Trees, RATES *Rates)
{
	int SIndex, NIndex;
	TREE *Tree;
	NODE N;
	FATTAILRATES *FTR;

	Tree = Trees->Tree[Rates->TreeNo];

	FTR = Rates->FatTailRates;


	FatTailSetAnsSates(Tree, Trees->NoOfSites, FTR);
	MapRatesToFatTailRate(Rates, FTR);

	for(SIndex=0;SIndex<Trees->NoOfSites;SIndex++)
	{
		SetStableDist(FTR->SDList[SIndex], FTR->Alpha[SIndex], FTR->Scale[SIndex]);

		for(NIndex=0;NIndex<Tree->NoNodes;NIndex++)
		{
			N = Tree->NodeList[NIndex];

			if(N->Tip == FALSE)
				NodeSliceSampleFatTail(N, SIndex, Opt, Trees, Rates);
		}
	}

	FatTailGetAnsSates(Tree, Trees->NoOfSites, FTR);

	Rates->Lh = Likelihood(Rates, Trees, Opt);
}

int	GetMutateFatTailRatesPos(OPTIONS *Opt, TREES* Trees, RATES* Rates, SCHEDULE* Shed)
{
	int Pos;



	if(Opt->FatTailNormal == FALSE)
		return RandUSInt(Rates->RS) % Shed->NoParm;

	do
	{
		Pos = RandUSInt(Rates->RS) % Shed->NoParm;
	}while(Rates->Rates[Pos] == FAT_TAIL_NORMAL_VAL);

	return Pos;
}
/*
void MutateFatTailRates(OPTIONS *Opt, TREES* Trees, RATES* Rates, SCHEDULE*	Shed)
{
	int Pos;
	double NewR, OldR, Dev;

//	Shed->PNo = RandUSInt(Rates->RS) % Shed->NoParm;
	Shed->PNo = GetMutateFatTailRatesPos(Opt, Trees, Rates, Shed);

	Pos = Shed->PNo;
	Dev = Opt->RateDevList[Shed->PNo];
//	Dev = 0.5;
	OldR = Rates->Rates[Pos];

	do
	{
		NewR = OldR + (RandDouble(Rates->RS) * Dev) - (Dev / 2.0);
	} while(NewR < 0.0);
	
	Rates->Rates[Pos] = NewR;

//	MapRatesToFatTailRate(Rates, Rates->FatTailRates);
}
*/

void MutateFatTailRates(OPTIONS *Opt, TREES* Trees, RATES* Rates, SCHEDULE*	Shed)
{
	int Index, Pos;
	double NewR, OldR, Dev;

//	Shed->PNo = RandUSInt(Rates->RS) % Shed->NoParm;
	Shed->PNo = GetMutateFatTailRatesPos(Opt, Trees, Rates, Shed);

	Pos = Shed->PNo;
	Dev = Opt->RateDevList[Shed->PNo];
	OldR = Rates->Rates[Pos];

	do
	{
		NewR = OldR + (RandDouble(Rates->RS) * Dev) - (Dev / 2.0);
	} while(NewR < 0.0);
	
	Rates->Rates[Pos] = NewR;
	
//	MapRatesToFatTailRate(Rates, Rates->FatTailRates);
}

void	InitFattailFile(OPTIONS *Opt, TREES *Trees)
{
	char	*Buffer;
	TREE	*Tree;
	NODE	N;
	int		Index, SIndex, TIndex, NID;
	PART	*Part;
	TAXA	*Taxa;

	Buffer = (char*)malloc(sizeof(char) * (strlen(Opt->LogFN) + BUFFERSIZE));
	if(Buffer == NULL)
		MallocErr();
	sprintf(Buffer, "%s.AnsStates.txt", Opt->LogFN);

	Opt->LogFatTail = OpenWrite(Buffer);
	
	free(Buffer);

	Tree = Trees->Tree[0];
	NID = 0;
	for(Index=0;Index<Tree->NoNodes;Index++)
	{
		N = Tree->NodeList[Index];
		if(N->Tip == FALSE)
		{
			Part = N->Part;

			fprintf(Opt->LogFatTail, "Node-%05d\t", NID++);
			
			for(TIndex=0;TIndex<Part->NoTaxa;TIndex++)
			{
				Taxa = Trees->Taxa[Part->Taxa[TIndex]];
				fprintf(Opt->LogFatTail, "%s\t", Taxa->Name);
			}
			fprintf(Opt->LogFatTail, "\n");
		}
	}

	fprintf(Opt->LogFatTail, "Itter\tLh\t");

	if(Opt->UseGeoData == TRUE)
		fprintf(Opt->LogFatTail, "Alpha\tScale\t");
	else
	{
		for(Index=0;Index<Trees->NoOfSites;Index++)
			fprintf(Opt->LogFatTail, "Alpha %d\tScale %d\t", Index+1, Index+1);
	}

	for(Index=0;Index<NID;Index++)
	{
		if(Opt->UseGeoData == TRUE)
			fprintf(Opt->LogFatTail, "Node-%05d - Long\tNode-%05d - Lat\t", Index, Index);
		//	fprintf(Opt->LogFatTail, "Node-%05d - X\tNode-%05d - Y\tNode-%05d - Z\t", Index, Index, Index);
		else
		{
			for(SIndex=0;SIndex<Trees->NoOfSites;SIndex++)
				fprintf(Opt->LogFatTail, "Node-%05d - %d\t",Index, SIndex+1);
		}
	}
	
	fprintf(Opt->LogFatTail, "\n");

	fflush(Opt->LogFatTail);
}

void	OutputFatTail(long long Itter, OPTIONS *Opt, TREES *Trees, RATES *Rates)
{
	int Index, SIndex;
	NODE N;
	TREE *Tree;

	double Long, Lat;
	
	Tree = Trees->Tree[Rates->TreeNo];

	FatTailSetAnsSates(Tree, Trees->NoOfSites, Rates->FatTailRates);

	fprintf(Opt->LogFatTail, "%lld\t%f\t", Itter, Rates->Lh);

	for(Index=0;Index<Rates->NoOfRates;Index++)
		fprintf(Opt->LogFatTail, "%f\t", Rates->Rates[Index]);

	for(Index=0;Index<Tree->NoNodes;Index++)
	{
		N = Tree->NodeList[Index];
		if(N->Tip == FALSE)
		{
			if(Opt->UseGeoData == TRUE)
			{
				NodeToLongLat(N, &Long, &Lat);
				fprintf(Opt->LogFatTail, "%f\t%f\t", Long, Lat);
			}
			else
			{
				NodeToLongLat(N, &Long, &Lat);
				for(SIndex=0;SIndex<Trees->NoOfSites;SIndex++)
					fprintf(Opt->LogFatTail, "%f\t", N->FatTailNode->Ans[SIndex]);
			}
		}
	}

	fprintf(Opt->LogFatTail, "\n");
	fflush(Opt->LogFatTail);
}

double	FatTailLhPraxis(void* P, double *List)
{
	PRAXSTATE	*PState;
	double		Ret;

	PState = (PRAXSTATE*)P;

	memcpy(PState->Rates->Rates, List, sizeof(double) * PState->n);

	Ret = Likelihood(PState->Rates, PState->Trees, PState->Opt);

	printf("Lh:\t%f\n", Ret);fflush(stdout);

	return Ret;
}

void	TestLhNorm2(RATES *Rates)
{
	double p, x, a, c;
	STABLEDIST *SDist;
	
	x = 0.1;
	a = 0.2;
	c = 1.0;

	SDist = Rates->FatTailRates->SDList[0];

	SetStableDist(SDist, a, c);
	p = StableDistPDF(SDist, x);
//	p = exp(p);

	printf("p:\t%f\n", p);
	exit(0);
}

void	TestLhNorm(RATES *Rates)
{
	double X, P, SD;
	STABLEDIST *SDist;

//	TestLhNorm2(Rates);
	
	SDist = Rates->FatTailRates->SDList[0];
	
	SD = 5.0;
	SD = sqrt(SD);

	SetStableDist(SDist, FAT_TAIL_NORMAL_VAL, SD/sqrt(2.0));
	
	for(X = -5;X < 5; X+=0.001)
	{
		P = StableDistPDF(SDist, X);
		printf("%f\t%f\n", X, exp(P));
	}
	exit(0);
}


void	SetRandFatTail(OPTIONS *Opt, RATES *Rates, int SiteNo)
{
	int Pos;
	PRIORS *P;

//	TestLhNorm(Rates);

	Pos = SiteNo * 2;
		
	P = Rates->Prios[Pos];
	if(Opt->FatTailNormal == FALSE)
		Rates->Rates[Pos] = RandUniDouble(Rates->RS, P->DistVals[0], P->DistVals[1]);
	else
		Rates->Rates[Pos] = FAT_TAIL_NORMAL_VAL;
	Pos++;

	P = Rates->Prios[Pos];
//	Rates->Rates[Pos] = RandUniDouble(Rates->RS, P->DistVals[0], P->DistVals[1]);
	Rates->Rates[Pos] = RandUniDouble(Rates->RS, 0, 10);
	Pos++;

}

void	InitFatTailRates(OPTIONS *Opt, TREES *Trees, RATES *Rates)
{
	double Lh;
	int Index;
	
	do
	{
		for(Index=0;Index<Rates->FatTailRates->NoSD;Index++)
			SetRandFatTail(Opt, Rates, Index);

		Lh = Likelihood(Rates, Trees, Opt);
	} while(Lh == ERRLH);

	return;
/*
	PRAXSTATE* PS;
	double		*TempVect;
	double		Lh;

	TempVect = (double*)malloc(sizeof(double) * Rates->NoOfRates);
	if(TempVect == NULL)
		MallocErr();
	memcpy(TempVect, Rates->Rates, sizeof(double) * Rates->NoOfRates);
	
	PS = IntiPraxis(FatTailLhPraxis, TempVect, Rates->NoOfRates, 0, 0, 4, 10000);
	
	PS->Opt		= Opt;
	PS->Trees	= Trees;
	PS->Rates	= Rates;
	
	Lh = praxis(PS);

	memcpy(Rates->Rates, TempVect, sizeof(double) * Rates->NoOfRates);
		
	FreePracxStates(PS);

	printf("%f\t%f\t%f\n", Lh, Rates->Rates[0], Rates->Rates[1]);
	exit(0);*/
}