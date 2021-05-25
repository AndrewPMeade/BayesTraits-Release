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

#define NO_SLICE_STEPS 1000
#define	MAX_STEP_DIFF	5.0

double	NodeSliceSampler(NODE N, RANDSTATES *RS);
void	SetInitAnsStates(OPTIONS *Opt, TREES *Trees, TREE *Tree);

void MapRatesToTree(TREE *Tree, int NoSites, FATTAILRATES *FTR)
{
	size_t Size;
	
	Size = sizeof(double) * Tree->NoNodes * NoSites;

	memcpy(Tree->FatTailTree->AnsVect, FTR->AnsVect, Size);
}

void MapTreeToRates(TREE *Tree, int NoSites, FATTAILRATES *FTR)
{
	size_t Size;
	
	Size = sizeof(double) * Tree->NoNodes * NoSites;

	memcpy(FTR->AnsVect, Tree->FatTailTree->AnsVect, Size);
}

void MapRatesToFatTailRate(int NoSites, RATES *Rates, FATTAILRATES *FatTailRates)
{
	int Index, Pos;

	Pos = 0;
	for(Index=0;Index<NoSites;Index++)
	{
		FatTailRates->Alpha[Index] = Rates->Rates[Pos++];
		FatTailRates->Scale[Index] = Rates->Rates[Pos++];
	}
}


void MapFatTailRateToRates(int NoSites, RATES *Rates, FATTAILRATES *FatTailRates)
{
	int Index, Pos;

	Pos = 0;
	for(Index=0;Index<NoSites;Index++)
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

	Ret->SliceX = (double*)malloc(sizeof(double) * NO_SLICE_STEPS);
	Ret->SliceY = (double*)malloc(sizeof(double) * NO_SLICE_STEPS);

	if(Ret->SliceX == NULL || Ret->SliceY == NULL)
		MallocErr();

	Ret->SliceMin = (double*)malloc(sizeof(double) * NO_SLICE_STEPS);
	Ret->SliceMax = (double*)malloc(sizeof(double) * NO_SLICE_STEPS);
	if(Ret->SliceMin == NULL || Ret->SiteMax == NULL)
		MallocErr();

	Ret->AnsVect = (double*)malloc(sizeof(double) * NoSites * Trees->MaxNodes);
	if(Ret->AnsVect == NULL)
		MallocErr();	

	Ret->SDList = (STABLEDIST**)malloc(sizeof(double) * NoSites);
	if(Ret->SDList == NULL)
		MallocErr();
	for(Index=0;Index<NoSites;Index++)
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

	for(Index=0;Index<Opt->Trees->NoOfSites;Index++)
	{
		Ret->Alpha[Index] = 0.5;
		Ret->Scale[Index] = 0.5;
	}

	for(Index=0;Index<Opt->Trees->NoOfSites;Index++)
		GetSiteInfo(Index, Trees, Ret);

	SetInitAnsStates(Opt, Trees, Trees->Tree[0]);
	MapTreeToRates(Trees->Tree[0], Trees->NoOfSites, Ret);
			
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

	free(FTR->SliceX);
	free(FTR->SliceY);

	for(Index=0;Index<NoSites;Index++)
		FreeStableDist(FTR->SDList[Index]);
	free(FTR->SDList);
	
	free(FTR->SliceMin);
	free(FTR->SliceMax);

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


double	CalcNodeStableLh(NODE N, int NoSites, STABLEDIST **SDList)
{
	int Index, SIndex;
	double Ret;
	double L, x;
	
	if(N->Tip == TRUE)
		return 0;

	Ret = 0;
	
	for(Index=0;Index<N->NoNodes;Index++)
	{
		for(SIndex=0;SIndex<NoSites;SIndex++)
		{
			x = N->FatTailNode->Ans[SIndex]- N->NodeList[Index]->FatTailNode->Ans[SIndex];
			L = StableDistTPDF(SDList[SIndex], x , N->NodeList[Index]->Length);
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

	MapRatesToTree(Tree, NoSites, FTR);

	for(Index=0;Index<NoSites;Index++)
		SetStableDist(FTR->SDList[Index], FTR->Alpha[Index], FTR->Scale[Index]);
	
	Ret = 0;

	#pragma omp parallel for num_threads(Opt->Cores) reduction(+:Ret)
	for(Index=0;Index<Tree->NoNodes;Index++)
	{
		if(Tree->NodeList[Index]->Tip == FALSE)
			Ret += CalcNodeStableLh(Tree->NodeList[Index], NoSites, FTR->SDList);
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

void	SetXPosVect(double Min, double Max, int NoSteps, double *XVect)
{
	double StepSize;
	int Index;

//	for(Index=0;Index<NoSteps;Index++)
//		XVect[Index] = Index;
//	return;

	StepSize = (Max - Min) / (NoSteps - 1);

	for(Index=0;Index<NoSteps;Index++)
		XVect[Index] = Min + (StepSize * Index);
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

int		SetYPosVect(OPTIONS *Opt, NODE N, int SiteNo, STABLEDIST *SD, int NoSteps, double *XVect, double *YVect)
{
	int Index; 
	
	#pragma omp parallel for num_threads(Opt->Cores) 
	for(Index=0;Index<NoSteps;Index++)
		YVect[Index] = AnsStateLh(XVect[Index], SiteNo, N, SD);

	for(Index=0;Index<NoSteps;Index++)
	{
		if(Index > 1)
		{
			if(YVect[Index] - YVect[Index-1] > MAX_STEP_DIFF)
				return FALSE;

			if(YVect[Index] == YVect[0])
				return FALSE;
		}

		if(YVect[Index] != YVect[Index] || YVect[Index] == YVect[Index] + 1.0)
			return FALSE;
	}

	return TRUE;
}

void	GetSiteMinMax(FATTAILRATES *FTR, int SiteNo, double *Min, double *Max)
{
	*Min = FTR->SiteMin[SiteNo] - (FTR->SiteMin[SiteNo] * 0.1);
	*Max = FTR->SiteMax[SiteNo] + (FTR->SiteMax[SiteNo] * 0.1);
}

void	GetMinMaxY(double *List, int Size, double *Min, double *Max)
{
	int Index;

	*Min = List[0];
	*Max = List[0];

	for(Index=1;Index<Size;Index++)
	{
		if(List[Index] > *Max)
			*Max = List[Index];

		if(List[Index] < *Min)
			*Min = List[Index];
	}
}

double	GetXCrossPoint(double YPoint, int Index, double *YList, double *XList)
{
	double Ret;
	
	double Diff1, Diff2;

//	Slope = (YList[Index-1] - YList[Index]) / (XList[Index-1] - XList[Index]);
//	YDiff = YPoint - YList[Index-1];

	Diff1 = fabs(YPoint - YList[Index-1]);
	Diff2 = fabs(YPoint - YList[Index]);


	Ret = Diff1 / (Diff1 + Diff2);
	Ret = Ret * (XList[Index] - XList[Index-1]);
	Ret = Ret + XList[Index-1];

	return Ret;
}

double FindSliceStart(double YPoint, FATTAILRATES *FTR, int *Pos)
{
	double	*YList, *XList;

	YList = FTR->SliceY;
	XList = FTR->SliceX;

	if(*Pos == 0 && YList[0] > YPoint)
		return XList[0];

	while(*Pos < NO_SLICE_STEPS)
	{
		if(YList[*Pos] > YPoint)
			return GetXCrossPoint(YPoint, *Pos, YList, XList);

		(*Pos)++;
	}

	return 0.0;
}

double FindSliceEnd(double YPoint, FATTAILRATES *FTR, int *Pos)
{
	double	*YList, *XList;

	YList = FTR->SliceY;
	XList = FTR->SliceX;


	while(*Pos < NO_SLICE_STEPS)
	{
		if(*Pos == NO_SLICE_STEPS - 1)
			return XList[*Pos];
		
		if(YList[*Pos] < YPoint)
			return GetXCrossPoint(YPoint, *Pos, YList, XList);

		(*Pos)++;
	}

	return 0.0;
}

double	FindNewPoint(FATTAILRATES *FTR, RANDSTATES *RS)
{
	double Sum, Point;
	int Index;

	Sum = 0.0;

	for(Index=0;Index<FTR->NoSlices;Index++)
		Sum += FTR->SliceMax[Index] - FTR->SliceMin[Index];

	Point = Sum * RandDouble(RS);

	Sum = 0.0;

	for(Index=0;Index<FTR->NoSlices;Index++)
	{
		if(Point < Sum + (FTR->SliceMax[Index] - FTR->SliceMin[Index]))
			return FTR->SliceMin[Index] + (Point - Sum);

		Sum += FTR->SliceMax[Index] - FTR->SliceMin[Index];
	}


	for(Index=0;Index<NO_SLICE_STEPS;Index++)
		printf("%d\t%f\t%f\n", Index, FTR->SliceX[Index], FTR->SliceY[Index]);
	fflush(stdout);

	printf("Error in (%s::%d)\n", __FILE__, __LINE__);
	exit(0);

	return 1.0;
}


double	GetSlices(FATTAILRATES *FTR, RANDSTATES *RS, double PCAns)
{
	int CPos, Exit, NoSlices;
	double	SPos, EPos, Ret, YPoint, MinY, MaxY;
	double	*YList, *XList;

	YList = FTR->SliceY;
	XList = FTR->SliceX;

	Ret = 0;
	
	GetMinMaxY(FTR->SliceY, NO_SLICE_STEPS, &MinY, &MaxY);

	do
	{
		YPoint = MinY + ((PCAns - MinY) * RandDouble(RS));
	} while(YPoint > MaxY);
	
	Exit = FALSE;

	CPos = 0;
	NoSlices = 0;

	do
	{
		SPos = FindSliceStart(YPoint, FTR, &CPos);

		if(CPos >= NO_SLICE_STEPS)
			Exit = TRUE;
		else
		{
			EPos = FindSliceEnd(YPoint, FTR, &CPos);

			FTR->SliceMin[NoSlices] = SPos;
			FTR->SliceMax[NoSlices] = EPos;
			NoSlices++;
			CPos+=1;
		}

	}while(Exit == FALSE);

	FTR->NoSlices = NoSlices;

	Ret = FindNewPoint(FTR, RS);
	return Ret;
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

	MapRatesToTree(Tree, Trees->NoOfSites, FTR);
	
	SetNodeAns(Trees, Tree->Root, 0.0);
	
	MapTreeToRates(Tree, Trees->NoOfSites, FTR);
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

	MapRatesToTree(Tree, Trees->NoOfSites, FTR);

//	PrintAnsVect(Tree);
	
	TSite = RandUSInt(Rates->RS) % Trees->NoOfSites;
	TNode = GetSliceSampleNode(Tree, Rates->RS);

	
	SetStableDist(FTR->SDList[TSite], FTR->Alpha[TSite], FTR->Scale[TSite]);
	
	CAns = TNode->FatTailNode->Ans[TSite];
	PCAns = AnsStateLh(CAns, TSite, TNode, FTR->SDList[TSite]);
	
	GetSiteMinMax(FTR, TSite, &Min, &Max);

//	TNode->FatTailNode->Ans[TSite] = Min + (RandDouble(Rates->RS) * (Max - Min));
//	MapTreeToRates(Tree, Trees->NoOfSites, FTR);

	SetXPosVect(Min, Max, NO_SLICE_STEPS, FTR->SliceX); 
	
	Valid = SetYPosVect(Opt, TNode, TSite, FTR->SDList[TSite], NO_SLICE_STEPS, FTR->SliceX, FTR->SliceY);

	if(Valid == FALSE)
		NAns = Min + (RandDouble(Rates->RS) * (Max - Min));
	else
		NAns = GetSlices(FTR, Rates->RS, PCAns);

	TNode->FatTailNode->Ans[TSite] = NAns;
	
	MapTreeToRates(Tree, Trees->NoOfSites, FTR);
	
	return;
}

void	NodeSliceSampleFatTail(NODE N, int SiteNo, OPTIONS *Opt, TREES *Trees, RATES *Rates) 
{
	int Changed, Valid, Index;
	double CLh, CAns, NAns, NLh, Min, Max;
	FATTAILRATES *FTR;

	FTR = Rates->FatTailRates;
	
	CAns = N->FatTailNode->Ans[SiteNo];
	CLh = AnsStateLh(CAns, SiteNo, N, FTR->SDList[SiteNo]);

	FTR = Rates->FatTailRates;
	
	GetSiteMinMax(FTR, SiteNo, &Min, &Max);

	SetXPosVect(Min, Max, NO_SLICE_STEPS, FTR->SliceX); 
	
	Valid = SetYPosVect(Opt, N, SiteNo, FTR->SDList[SiteNo], NO_SLICE_STEPS, FTR->SliceX, FTR->SliceY);
	
	Changed = FALSE;
	do
	{
		if(Valid == FALSE)
			NAns = Min + (RandDouble(Rates->RS) * (Max - Min));
		else
			NAns = GetSlices(FTR, Rates->RS, CLh);

		NLh = AnsStateLh(NAns, SiteNo, N, FTR->SDList[SiteNo]);

		if(log(RandDouble(Rates->RS)) < (NLh - CLh))
		{
		//	printf("New:\t%f\tOld\t%f\tDiffer\t%f\t%f\t%f\n", NLh, CLh, NLh - CLh, NAns, CAns);fflush(stdout);
			Changed = TRUE;
		}
	} while(Changed == FALSE);

	N->FatTailNode->Ans[SiteNo] = NAns;
}


void	TestSample(OPTIONS *Opt, TREES *Trees, RATES *Rates)
{
	double Lh;
	TREE *Tree;
	NODE N;
	FATTAILRATES *FTR;

	Rates->Rates[0] = 1.657206;
	Rates->Rates[1] = 0.036151;
	
	Tree = Trees->Tree[Rates->TreeNo];

	FTR = Rates->FatTailRates;
	
	Lh = Likelihood(Rates, Trees, Opt);

	printf("Lh:\t%f\n", Lh);

	MapRatesToTree(Tree, Trees->NoOfSites, FTR);


	SliceSampleFatTail(Opt, Trees, Rates);


	MapTreeToRates(Tree, Trees->NoOfSites, FTR);

	Lh = Likelihood(Rates, Trees, Opt);
	
	printf("Lh:\t%f\n", Lh);

	exit(0);
}

void	AllSliceSampleFatTail(OPTIONS *Opt, TREES *Trees, RATES *Rates)
{
	int SIndex, NIndex;
	double Lh;
	TREE *Tree;
	NODE N;
	FATTAILRATES *FTR;

//	TestSample(Opt, Trees, Rates);

	Tree = Trees->Tree[Rates->TreeNo];

	FTR = Rates->FatTailRates;

//	Lh = Likelihood(Rates, Trees, Opt);
//	printf("Lh:\t%f\n", Lh);
//	exit(0);

	MapRatesToTree(Tree, Trees->NoOfSites, FTR);

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

	MapTreeToRates(Tree, Trees->NoOfSites, FTR);

//	Lh = Likelihood(Rates, Trees, Opt);
//	printf("Lh:\t%f\n", Lh);
//	exit(0);
}

void MutateFatTailRates(OPTIONS *Opt, TREES* Trees, RATES* Rates, SCHEDULE*	Shed)
{
	int Pos;
	double NewR, OldR, Dev;

	Shed->PNo = RandUSInt(Rates->RS) % Shed->NoParm;

	Pos = Shed->PNo;
	Dev = Opt->RateDevList[Shed->PNo];

	OldR = Rates->Rates[Pos];

	do
	{
		NewR = OldR + (RandDouble(Rates->RS) * Dev) - (Dev / 2.0);
	} while(NewR < 0.0);
	
	Rates->Rates[Pos] = NewR;
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

	for(Index=0;Index<Trees->NoOfSites;Index++)
		fprintf(Opt->LogFatTail, "Alpha %d\tScale %d\t", Index+1, Index+1);

	for(Index=0;Index<NID;Index++)
	{
		for(SIndex=0;SIndex<Trees->NoOfSites;SIndex++)
			fprintf(Opt->LogFatTail, "Node-%05d - %d\t",Index, SIndex+1);
	}
	
	fprintf(Opt->LogFatTail, "\n");

	fflush(Opt->LogFatTail);
}

void	OutputFatTail(long long Itter, OPTIONS *Opt, TREES *Trees, RATES *Rates)
{
	int Index, SIndex;
	NODE N;
	TREE *Tree;

	Tree = Trees->Tree[Rates->TreeNo];

	MapRatesToTree(Tree, Trees->NoOfSites, Rates->FatTailRates);

	fprintf(Opt->LogFatTail, "%lld\t%f\t", Itter, Rates->Lh);
	for(Index=0;Index<Rates->NoOfRates;Index++)
		fprintf(Opt->LogFatTail, "%f\t", Rates->Rates[Index]);

	for(Index=0;Index<Tree->NoNodes;Index++)
	{
		N = Tree->NodeList[Index];
		if(N->Tip == FALSE)
		{
			for(SIndex=0;SIndex<Trees->NoOfSites;SIndex++)
				fprintf(Opt->LogFatTail, "%f\t", N->FatTailNode->Ans[SIndex]);
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

void	SetRandFatTail(RATES *Rates, int SiteNo)
{
	int Pos;
	PRIORS *P;

	Pos = SiteNo * 2;
		
	P = Rates->Prios[Pos];
	Rates->Rates[Pos] = RandUniDouble(Rates->RS, P->DistVals[0], P->DistVals[1]);
	Pos++;

	P = Rates->Prios[Pos];
	Rates->Rates[Pos] = RandUniDouble(Rates->RS, P->DistVals[0], P->DistVals[1]);
	Pos++;

}

void	InitFatTailRates(OPTIONS *Opt, TREES *Trees, RATES *Rates)
{
	double Lh;
	int Index;
	
	do
	{
		for(Index=0;Index<Trees->NoOfSites;Index++)
			SetRandFatTail(Rates, Index);

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