/*
*  BayesTriats 4.0
*
*  copyright 2022
*
*  Andrew Meade
*  School of Biological Sciences
*  University of Reading
*  Reading
*  Berkshire
*  RG6 6BX
*
* BayesTriats is free software: you can redistribute it and/or modify
* it under the terms of the GNU General Public License as published by
* the Free Software Foundation, either version 3 of the License, or
* (at your option) any later version.
* 
* This program is distributed in the hope that it will be useful,
* but WITHOUT ANY WARRANTY; without even the implied warranty of
* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
* GNU General Public License for more details.
* 
* You should have received a copy of the GNU General Public License
* along with this program.  If not, see <http://www.gnu.org/licenses/>
*
*/



#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>

#include "FatTail.h"
#include "TypeDef.h"
#include "GenLib.h"
#include "StableDist.h"
#include "Likelihood.h"
#include "Part.h"
#include "Trees.h"
#include "Praxis.h"
#include "Geo.h"
#include "SliceSampler.h"
#include "MCMC.h"
#include "DistData.h"
#include "Threaded.h"
#include "IntraNode.h"



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

double*	SetPartialLikelihoodMem(TREES *Trees)
{
	assert(Trees->NoTrees == 1);
	return (double*)SMalloc(sizeof(double) * Trees->Tree[0]->NoInternalNodes);
}

SLICESAMPLER**	AllocSliceSamplers(int NoT, int NoSteps)
{
	SLICESAMPLER** Ret;
	int Index;

	Ret = (SLICESAMPLER**)SMalloc(sizeof(SLICESAMPLER*) * NoT);

	for(Index=0;Index<NoT;Index++)
		Ret[Index] = CrateSliceSampler(NoSteps);

	return Ret;
}

FATTAILRATES*	AllocFatTailRates(OPTIONS *Opt, TREES *Trees)
{
	FATTAILRATES* Ret;
	int NoSites;

	NoSites = Trees->NoSites;

	Ret = (FATTAILRATES*)SMalloc(sizeof(FATTAILRATES));

	Ret->SiteMin = (double*)SMalloc(sizeof(double) * NoSites);
	Ret->SiteMax = (double*)SMalloc(sizeof(double) * NoSites);
	Ret->SiteSD  = (double*)SMalloc(sizeof(double) * NoSites);


	Ret->SliceSamplers = AllocSliceSamplers(GetMaxThreads(), Opt->NoSliceSampleSteps);
	
	Ret->AnsVect = (double*)SMalloc(sizeof(double) * NoSites * Trees->MaxNodes);

	Ret->PartialLh = SetPartialLikelihoodMem(Trees);
	
	return Ret;
}

void			GetSiteInfo(int SiteNo, TREES *Trees, FATTAILRATES* FTR)
{
	int Index;
	double Data, Mean, SD;


	FTR->SiteMax[SiteNo] = Trees->Taxa[0]->ConData[SiteNo];
	FTR->SiteMin[SiteNo] = Trees->Taxa[0]->ConData[SiteNo];

	Mean = 0;

	for(Index=0;Index<Trees->NoTaxa;Index++)
	{
		Data = Trees->Taxa[Index]->ConData[SiteNo];
				
		Mean += Data;

		if(Data > FTR->SiteMax[SiteNo])
			FTR->SiteMax[SiteNo] = Data;

		if(Data < FTR->SiteMin[SiteNo])
			FTR->SiteMin[SiteNo] = Data;
	}

	Mean = Mean / Trees->NoTaxa;

	SD = 0;
	for(Index=0;Index<Trees->NoTaxa;Index++)
	{
		Data = Trees->Taxa[Index]->ConData[SiteNo];
		SD += (Data - Mean) * (Data - Mean);
	}
	SD = SD / Trees->NoTaxa;
	SD = sqrt(SD);

	FTR->SiteSD[SiteNo] = SD;
}

FATTAILRATES*	CreateFatTailRates(RATES *Rates, OPTIONS *Opt, TREES *Trees)
{
	FATTAILRATES* Ret;
	int	Index;

	Ret = AllocFatTailRates(Opt, Trees);

	for(Index=0;Index<Rates->NoOfRates;Index++)
	{
		Rates->Rates[Index] = 0.5;
	}
	
	for(Index=0;Index<Trees->NoSites;Index++)
		GetSiteInfo(Index, Trees, Ret);

	SetInitAnsStates(Opt, Trees, Trees->Tree[0]);
	FatTailGetAnsSates(Trees->Tree[0], Trees->NoSites, Ret);
	
	CheckRestictedMaps(Trees, Ret);

	return Ret;
}

void	FreeFatTailRates(FATTAILRATES* FTR, int NoSites, int NoCores)
{
	int Index;

	free(FTR->SiteMin);
	free(FTR->SiteMax);
	free(FTR->SiteSD);

	for(Index=0;Index<NoCores;Index++)
		FreeSliceSampler(FTR->SliceSamplers[Index]);
	free(FTR->SliceSamplers);
		
	free(FTR->AnsVect);

	free(FTR->PartialLh);
	
	free(FTR);
}

void			CopyFatTailRates(TREES *Trees, FATTAILRATES *A, FATTAILRATES *B)
{
	// Alpha and Scale should be mapped from globabal rates
	
	memcpy(A->AnsVect, B->AnsVect, sizeof(double) * Trees->MaxNodes * Trees->NoSites);
}

FATTAILNODE*	InitFatTailNode(int NoSites, NODE N, double *AnsVect)
{
	FATTAILNODE*	Ret;
	int				Index;
	size_t			Pos;

	Ret = (FATTAILNODE*)SMalloc(sizeof(FATTAILNODE));
	
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

void	AllocFatTailTreeGroups(TREE *Tree)
{
	FATTAILTREE *FTT;
	size_t	Index, MaxNoNodes, MaxGroups;
	size_t i,j;

	MaxNoNodes = Tree->NoInternalNodes;
	MaxGroups = Tree->NoParallelGroups;

	FTT = Tree->FatTailTree;

	FTT->NoParallelGroups = 0;
	
	FTT->ParallelNodeListLength = (int*)SMalloc(sizeof(int) * MaxGroups);
	for(Index=0;Index<MaxGroups;Index++)
		FTT->ParallelNodeListLength[Index] = 0;

	FTT->ParallelNodeList = (NODE**)SMalloc(sizeof(NODE*) * MaxGroups);
	FTT->ParallelNodeList[0] = (NODE*)SMalloc(sizeof(NODE) * MaxGroups * MaxNoNodes);

	for(Index=1;Index<MaxGroups;Index++)
		FTT->ParallelNodeList[Index] = FTT->ParallelNodeList[0] + Index * MaxNoNodes;


	for(i=0;i<MaxGroups;i++)
		for(j=0;j<MaxNoNodes;j++)
			FTT->ParallelNodeList[i][j] = NULL;
}

int		ValidInGroup(NODE Node, NODE *NodeList, int GroupSize)
{
	int NIndex, DecIndex;
	NODE Dec;

	for(NIndex=0;NIndex<GroupSize;NIndex++)
	{
		Dec = NodeList[NIndex];

		for(DecIndex=0;DecIndex<Node->NoNodes;DecIndex++)
			if(Dec == Node->NodeList[DecIndex])
				return FALSE;
	}

	return TRUE;
}

int		FindNodeGroup(NODE N, FATTAILTREE *FTT)
{
	int GIndex;

	for(GIndex=0;GIndex<FTT->NoParallelGroups;GIndex++)
	{
		if(ValidInGroup(N, FTT->ParallelNodeList[GIndex], FTT->ParallelNodeListLength[GIndex]) == TRUE)
			return GIndex;
	}

	return -1;
}

void	AddNodeToGroup(NODE N, FATTAILTREE *FTT, int Group)
{
	int Pos;

	if(Group == -1)
	{
		Group = FTT->NoParallelGroups;
		FTT->NoParallelGroups++;
	}

	Pos = FTT->ParallelNodeListLength[Group];

	FTT->ParallelNodeList[Group][Pos] = N;
	FTT->ParallelNodeListLength[Group]++;
}

void	SetFatTailTreeGroups(TREE *Tree)
{
	FATTAILTREE *FTT;
	int GIndex, NIndex, Group;
	NODE Node;

	FTT = Tree->FatTailTree;

	AllocFatTailTreeGroups(Tree);

	memcpy(FTT->ParallelNodeList[0], Tree->ParallelNodes[0], sizeof(NODE) * Tree->ParallelGroupSize[0]);
	memcpy(FTT->ParallelNodeList[1], Tree->ParallelNodes[1], sizeof(NODE) * Tree->ParallelGroupSize[1]);

	FTT->NoParallelGroups = 2;
	FTT->ParallelNodeListLength[0] = Tree->ParallelGroupSize[0];
	FTT->ParallelNodeListLength[1] = Tree->ParallelGroupSize[1];

	for(GIndex=2;GIndex<Tree->NoParallelGroups;GIndex++)
	{
		for(NIndex=0;NIndex<Tree->ParallelGroupSize[GIndex];NIndex++)
		{
			Node = Tree->ParallelNodes[GIndex][NIndex];
			Group = FindNodeGroup(Node, FTT);
			AddNodeToGroup(Node, FTT, Group);
		}

	}
}

void	PrintFatTailGroups(TREE *Tree)
{
	FATTAILTREE *FTT;
	int Index, NIndex;
	NODE Node;

	FTT = Tree->FatTailTree;

	printf("FTT Group\n");
	for(Index=0;Index<FTT->NoParallelGroups;Index++)
		printf("%d\t%d\n", Index, FTT->ParallelNodeListLength[Index]);
	
	printf("Parallel Group\n");
	for(Index=0;Index<Tree->NoParallelGroups;Index++)
		printf("%d\t%d\n", Index, Tree->ParallelGroupSize[Index]);

	printf("IDs and Taxa\n");
	for(Index=0;Index<Tree->NoInternalNodes;Index++)
	{
		Node = Tree->NodeList[Index];
		printf("Node\t%d\t", Node->ID);
		RecPRintNodeTaxa(Node, ' ');
		printf("\n");
	}


	printf("FTT groups\n");
	for(Index=0;Index<FTT->NoParallelGroups;Index++)
	{
		printf("%d\t", Index);
		for(NIndex=0;NIndex<FTT->ParallelNodeListLength[Index];NIndex++)
		{
			Node = FTT->ParallelNodeList[Index][NIndex];

			printf("%d\t", Node->ID);
		}

		printf("\n");
	}

	printf("Parallel Groups\n");
	for(Index=0;Index<Tree->NoParallelGroups;Index++)
	{
		printf("%d\t", Index);
		for(NIndex=0;NIndex<Tree->ParallelGroupSize[Index];NIndex++)
		{
			Node = Tree->ParallelNodes[Index][NIndex];

			printf("%d\t", Node->ID);
		}

		printf("\n");
	}
	exit(0);
}

void	InitFatTailTree(OPTIONS *Opt, TREE *Tree, TREES *Trees)
{
	int Index;
	NODE N;

	for(Index=0;Index<Tree->NoNodes;Index++)
	{
		N = Tree->NodeList[Index];
		N->FatTailNode = InitFatTailNode(Trees->NoSites, N, Tree->FatTailTree->AnsVect);
	}

	SetFatTailTreeGroups(Tree);
//	PrintFatTailGroups(Tree);
}

FATTAILTREE*	AllocFatTailTree(TREE *Tree, int NoSites)
{
	FATTAILTREE* Ret;

	Ret = (FATTAILTREE*)SMalloc(sizeof(FATTAILTREE));

	Ret->AnsVect = (double*)SMalloc(sizeof(double) * Tree->NoNodes * NoSites);

	Ret->ParallelNodeList = NULL;
	Ret->ParallelNodeListLength = NULL;
	Ret->NoParallelGroups = -1;

	return Ret;
}

void			FreeFatTailTree(FATTAILTREE *FatTailTree)
{
	free(FatTailTree->AnsVect);
		
	free(FatTailTree->ParallelNodeList[0]);
	free(FatTailTree->ParallelNodeList);
	free(FatTailTree->ParallelNodeListLength);

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
	// Values is only set from the first two nodes, 
	// its only an inishal aproximuation. 

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

	for(SIndex=0;SIndex<Trees->NoSites;SIndex++)
	{
		SetInitAnsStateNodes(SIndex, Tree);
		SetContrastAnsStates(Tree->Root);
		GetInitAnsStateNodes(SIndex, Tree);
	}

	if(Opt->Model == M_GEO)
		CorrectIntGeoNodes(Tree);
}


void	InitFatTailTrees(OPTIONS *Opt, TREES *Trees)
{
	int Index;
	TREE *Tree;

	for(Index=0;Index<Trees->NoTrees;Index++)
	{
		Tree = Trees->Tree[Index];
		Tree->FatTailTree = AllocFatTailTree(Tree, Trees->NoSites);
		InitFatTailTree(Opt, Tree, Trees);

	}

	SetTreesRNG(Trees, Opt->Seed);
}

double	CalcNodeStableLh(NODE N, int NoSites, double *ScaleList, int UseGeoModel)
{
	int Index, SIndex;
	double Ret;
	double L, x;
	double Scale;

	Ret = 0;
	Scale = ScaleList[0];
	for(SIndex=0;SIndex<NoSites;SIndex++)
	{
		if(UseGeoModel == FALSE)
			Scale = ScaleList[SIndex];


		for(Index=0;Index<N->NoNodes;Index++)
		{
			x = N->FatTailNode->Ans[SIndex]- N->NodeList[Index]->FatTailNode->Ans[SIndex];
			
			L = StableDistTPDF(Scale, x , N->NodeList[Index]->Length);
							
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
	int	UseGeoModel;
	
	Tree = Trees->Tree[Rates->TreeNo];
	NoSites = Trees->NoSites;
	FTR = Rates->FatTailRates;

#ifdef FAT_TAIL_ML_PARAM
	FTR->AnsVect[0] = FAT_TAIL_ML_ROOT;
	FTR->SDList[0]->Scale = FAT_TAIL_ML_SIG2;
	Rates->Rates[0] = FAT_TAIL_ML_SIG2;
#endif

	UseGeoModel = FALSE;
	
	if(Opt->Model == M_GEO)
		UseGeoModel = TRUE;
	
	FatTailSetAnsSates(Tree, NoSites, FTR);

	if(Opt->UseDistData == TRUE)
		SetTreeDistData(Rates, Opt, Trees);

	if(Opt->UseIntraNode == FALSE)
	{
		#ifdef OPENMP_THR
			#pragma omp parallel for num_threads(Opt->Cores)
		#endif
		for(Index=0;Index<Tree->NoInternalNodes;Index++)
			FTR->PartialLh[Index] = CalcNodeStableLh(Tree->InternalNodesList[Index], NoSites, Rates->Rates, UseGeoModel);

		Ret = 0;
		for(Index=0;Index<Tree->NoInternalNodes;Index++)
			Ret += FTR->PartialLh[Index];
	}
	else
		Ret = CalcIntraNodeLh(Trees, Rates);

	if(ValidLh(Ret, Opt->ModelType) == FALSE)
		return ERRLH;


	return Ret;
}

NODE	GetSliceSampleNode(TREE *Tree, gsl_rng *RNG)
{
	NODE Ret;
	int Pos;

	do
	{
//		Pos = RandUSInt(RS) % Tree->NoNodes;
		Pos = (int)gsl_rng_uniform_int(RNG, Tree->NoNodes);
		Ret = Tree->NodeList[Pos];
	}while(Ret->Tip == TRUE);

	return Ret;
}

double	AnsStateLh(double X, int SiteNo, NODE N, double Scale)
{
	double Ret, Val;
	int Index;

	Ret = 0;

	for(Index=0;Index<N->NoNodes;Index++)
	{
		Val = X - N->NodeList[Index]->FatTailNode->Ans[SiteNo];
		Ret += StableDistTPDF(Scale, Val, N->NodeList[Index]->Length);
	}

	if(N->Ans != NULL)
	{
		Val = X - N->Ans->FatTailNode->Ans[SiteNo];
		Ret += StableDistTPDF(Scale, Val, N->Length);
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

int		FatTailSetYPosVect(SLICESAMPLER *SS, OPTIONS *Opt, NODE N, int SiteNo, double Scale)
{
	int Index; 

	for(Index=0;Index<SS->NoSteps;Index++)
		SS->SliceY[Index] = AnsStateLh(SS->SliceX[Index], SiteNo, N, Scale);
	
	return TRUE;
//	Used to test if the lh differnce between slices is to big, posible lh err
/*	for(Index=0;Index<SS->NoSteps;Index++)
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

	return TRUE;*/
}

void	GetSiteMinMax(FATTAILRATES *FTR, int SiteNo, double *Min, double *Max)
{
	if(FTR->SiteMin[SiteNo] > 0)
		*Min = FTR->SiteMin[SiteNo] - (FTR->SiteMin[SiteNo] * 0.1);
	else
		*Min = FTR->SiteMin[SiteNo] + (FTR->SiteMin[SiteNo] * 0.1);

	if(FTR->SiteMax[SiteNo] > 0)
		*Max = FTR->SiteMax[SiteNo] + (FTR->SiteMax[SiteNo] * 0.1);
	else
		*Max = FTR->SiteMax[SiteNo] - (FTR->SiteMax[SiteNo] * 0.1);
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
	for(Index=0;Index<Trees->NoSites;Index++)
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

	FatTailSetAnsSates(Tree, Trees->NoSites, FTR);
	
	SetNodeAns(Trees, Tree->Root, 0.0);
	
	FatTailGetAnsSates(Tree, Trees->NoSites, FTR);
	Lh = Likelihood(Rates, Trees, Opt);
	printf("Lh:\t%f\n", Lh);
		PrintAnsStates(Trees, Tree->Root);

	exit(0);
}

void	SSNodeFatTail(NODE N, int SiteNo, OPTIONS *Opt, TREES *Trees, RATES *Rates, SLICESAMPLER *SS) 
{
	int Changed, Valid;
	double CLh, CAns, NAns, NLh, Min, Max;
	FATTAILRATES *FTR;
	gsl_rng *RNG;

	RNG = N->RNG;

	FTR = Rates->FatTailRates;
	
	CAns = N->FatTailNode->Ans[SiteNo];
	CLh = AnsStateLh(CAns, SiteNo, N, Rates->Rates[SiteNo]);

	FTR = Rates->FatTailRates;
	
	GetSiteMinMax(FTR, SiteNo, &Min, &Max);

	SSSetXPosVect(SS, Min, Max); 
	
	Valid = FatTailSetYPosVect(SS, Opt, N, SiteNo, Rates->Rates[SiteNo]);

	Changed = FALSE;
	do
	{
		if(Valid == FALSE)
			NAns = Min + (gsl_rng_uniform(RNG) * (Max - Min));
		else
			NAns = SSGetNewPoint(SS, RNG, CLh);

		NLh = AnsStateLh(NAns, SiteNo, N, Rates->Rates[SiteNo]);

		if(log(gsl_rng_uniform(RNG)) < (NLh - CLh))
			Changed = TRUE;

	} while(Changed == FALSE);

	N->FatTailNode->Ans[SiteNo] = NAns;
}

void	PrintAllAnsStates(TREES *Trees)
{
	int NIndex, SIndex;
	TREE *Tree;
	NODE N;

	Tree = Trees->Tree[0];

	for(NIndex=0;NIndex<Tree->NoNodes;NIndex++)
	{
		N = Tree->NodeList[NIndex];
		for(SIndex=0;SIndex<Trees->NoSites;SIndex++)
			printf("%f\t", N->FatTailNode->Ans[SIndex]);
	}
	printf("\n");
	fflush(stdout);
}


void	SSAllAnsStatesFatTailSite(OPTIONS *Opt, TREES *Trees, RATES *Rates, FATTAILRATES *FTR, TREE *Tree, int SiteNo)
{
	int FIndex, NIndex;
	NODE N;
	int TNo;

	for(FIndex=0;FIndex<Tree->NoParallelGroups;FIndex++)
	{
#ifdef OPENMP_THR
	#pragma omp parallel for num_threads(Opt->Cores) private(TNo, N) schedule(dynamic, 1)
#endif
		for(NIndex=0;NIndex<Tree->ParallelGroupSize[FIndex];NIndex++)
		{			
			N = Tree->ParallelNodes[FIndex][NIndex];
			TNo = GetThreadNo();
			SSNodeFatTail(N, SiteNo, Opt, Trees, Rates, FTR->SliceSamplers[TNo]);
		}
	}
}

void	SSAllAnsStatesFatTail(OPTIONS *Opt, TREES *Trees, RATES *Rates)
{
	int SIndex;
	TREE *Tree;
	
	FATTAILRATES *FTR;

	LhTransformTree(Rates, Trees, Opt);

	Tree = Trees->Tree[Rates->TreeNo];

	FTR = Rates->FatTailRates;

	FatTailSetAnsSates(Tree, Trees->NoSites, FTR);
	
	for(SIndex=0;SIndex<Trees->NoSites;SIndex++)
		SSAllAnsStatesFatTailSite(Opt, Trees, Rates, FTR, Tree, SIndex);

	FatTailGetAnsSates(Tree, Trees->NoSites, FTR);

	Rates->Lh = Likelihood(Rates, Trees, Opt);

	Rates->AutoAccept = TRUE;
}

void	SSAnsStatesFatTail(OPTIONS *Opt, TREES *Trees, RATES *Rates)
{
	int SIndex;
	NODE N;
	TREE *Tree;
	FATTAILRATES *FTR;

	LhTransformTree(Rates, Trees, Opt);

	Tree = Trees->Tree[Rates->TreeNo];

	FTR = Rates->FatTailRates;

	FatTailSetAnsSates(Tree, Trees->NoSites, FTR);
	
	SIndex = (int)gsl_rng_uniform_int(Rates->RNG, Trees->NoSites);

	N = GetSliceSampleNode(Tree, Rates->RNG);

	SSNodeFatTail(N, SIndex, Opt, Trees, Rates, FTR->SliceSamplers[0]);

	FatTailGetAnsSates(Tree, Trees->NoSites, FTR);

	Rates->Lh = Likelihood(Rates, Trees, Opt);

	Rates->AutoAccept = TRUE;
}

int	GetMutateFatTailRatesPos(OPTIONS *Opt, TREES* Trees, RATES* Rates, SCHEDULE* Shed)
{
	int Pos;

	if(Opt->Model == M_GEO)
		return 0;
	
	if(Opt->FatTailNormal == FALSE)
		return (int)gsl_rng_uniform_int(Rates->RNG, Shed->NoParm);


	do
	{
		Pos = (int)gsl_rng_uniform_int(Rates->RNG, Shed->NoParm);
	}while(Rates->Rates[Pos] == FAT_TAIL_NORMAL_VAL);

	return Pos;
}

void MutateFatTailRates(OPTIONS *Opt, TREES* Trees, RATES* Rates, SCHEDULE*	Shed)
{
	int Pos;
	double NewR, OldR, Dev;
		
//	Shed->PNo = RandUSInt(Rates->RS) % Shed->NoParm;
	Shed->PNo = GetMutateFatTailRatesPos(Opt, Trees, Rates, Shed);

	Pos = Shed->PNo;

	Shed->CurrentAT = Shed->RateDevATList[Shed->PNo];

	Dev = Shed->CurrentAT->CDev;
	OldR = Rates->Rates[Pos];

	do
	{
		NewR = OldR + (gsl_rng_uniform_pos(Rates->RNG) * Dev) - (Dev / 2.0);
	} while(NewR < 0.0);
	
	Rates->Rates[Pos] = NewR;
}

void	InitFattailFile(OPTIONS *Opt, TREES *Trees)
{
	
	TREE	*Tree;
	NODE	N;
	int		Index, SIndex, TIndex, NID;
	PART	*Part;
	TAXA	*Taxa;

	Opt->LogFatTail = OpenWithExt(Opt->CheckPointAppendFiles, Opt->BaseOutputFN, OUTPUT_EXT_ANC);

	if(Opt->CheckPointAppendFiles == TRUE)
		return;


	fprintf(Opt->LogFatTail, "Node Name\tBrach length\tHeight\tRestriction Map\tTaxa\n");

	Tree = Trees->Tree[0];
	NID = 0;
	for(Index=0;Index<Tree->NoNodes;Index++)
	{
		N = Tree->NodeList[Index];
		//		if(N->Tip == FALSE)
		{
			Part = N->Part;

			fprintf(Opt->LogFatTail, "Node-%05d\t", NID++);

			fprintf(Opt->LogFatTail, "%f\t%f\t", N->UserLength, N->Height);

			if(N->NodeResMap == NULL)
				fprintf(Opt->LogFatTail, "None\t");
			else
				fprintf(Opt->LogFatTail, "%s\t", N->NodeResMap->ResMap->FileName);

			
			for(TIndex=0;TIndex<Part->NoTaxa;TIndex++)
			{
				Taxa = Trees->Taxa[Part->Taxa[TIndex]];
				fprintf(Opt->LogFatTail, "%s\t", Taxa->Name);
			}
			fprintf(Opt->LogFatTail, "\n");
		}
	}

	fprintf(Opt->LogFatTail, "Itter\tLh\t");

	if(Opt->Model == M_GEO)
		fprintf(Opt->LogFatTail, "Scale\t");
	else
	{
		for(Index=0;Index<Trees->NoSites;Index++)
			fprintf(Opt->LogFatTail, "Sig2 %d\t", Index+1);
	}

	for(Index=0;Index<NID;Index++)
	{
		fprintf(Opt->LogFatTail, "Node-%05d - Branch Length\t", Index);
		if(Opt->Model == M_GEO)
			fprintf(Opt->LogFatTail, "Node-%05d - Long\tNode-%05d - Lat\t", Index, Index);
		else
		{
			for(SIndex=0;SIndex<Trees->NoSites;SIndex++)
				fprintf(Opt->LogFatTail, "Node-%05d - %d\t",Index, SIndex+1);
		}
	}
	
	fprintf(Opt->LogFatTail, "\n");

	fflush(Opt->LogFatTail);
}

void	OutputFatTail(size_t Itter, OPTIONS *Opt, TREES *Trees, RATES *Rates)
{
	int Index, SIndex;
	NODE N;
	TREE *Tree;

	double Long, Lat;
	
	Tree = Trees->Tree[Rates->TreeNo];

	FatTailSetAnsSates(Tree, Trees->NoSites, Rates->FatTailRates);

	fprintf(Opt->LogFatTail, "%zu\t%f\t", Itter, Rates->Lh);

	for(Index=0;Index<Rates->NoOfRates;Index++)
		fprintf(Opt->LogFatTail, "%f\t", Rates->Rates[Index]);

	for(Index=0;Index<Tree->NoNodes;Index++)
	{
		N = Tree->NodeList[Index];

		fprintf(Opt->LogFatTail, "%f\t", N->Length);

		if(Opt->Model == M_GEO)
		{
			NodeToLongLat(N, &Long, &Lat);
			fprintf(Opt->LogFatTail, "%f\t%f\t", Long, Lat);
		}
		else
		{				
			for(SIndex=0;SIndex<Trees->NoSites;SIndex++)
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



void	SetRandFatTail(OPTIONS *Opt, RATES *Rates, int SiteNo)
{
	int Pos;
	PRIOR *P;

	Pos = SiteNo * 2;
		
	P = Rates->Priors[Pos];
	Rates->Rates[Pos] = FAT_TAIL_NORMAL_VAL;
	Pos++;

	P = Rates->Priors[Pos];
	Rates->Rates[Pos] = gsl_rng_uniform_pos(Rates->RNG) * 100;
}



/*
void	InitFatTailRates(OPTIONS *Opt, TREES *Trees, RATES *Rates)
{
	int Index;
	
	do
	{
		for(Index=0;Index<Rates->FatTailRates->NoSD;Index++)
			SetRandFatTail(Opt, Rates, Index);

		MapRatesToFatTailRate(Rates, Rates->FatTailRates);

	} while(ValidMCMCParameters(Opt, Trees, Rates) == ERRLH);
	
	return;
}
*/

void CheckFatTailBL(TREES *Trees)
{
	int TIndex, NIndex;
	TREE *Tree;
	NODE N;

	for(TIndex=0;TIndex<Trees->NoTrees;TIndex++)
	{
		Tree = Trees->Tree[TIndex];

		for(NIndex=0;NIndex<Tree->NoNodes;NIndex++)
		{
			N = Tree->NodeList[NIndex];
			if(N->Ans != NULL)
			{
				if(N->Length < 1.00E-06)
				{
					printf("Branch length (%f) too short for fat tail model, in tree %d. Please resolve as a hard polytomy", N->Length, TIndex+1);
					exit(1);
				}
			}
		}
	}
}
// "5000000	288.328853	0.009159	2.488875	2.714179	3.828256	2.287889	2.436114	1.973251	2.610829	2.462855	2.366684	1.973528	1.915911	1.513107	1.949874	1.981611	1.845254	2.421675	2.438663	2.637304	2.716137	2.750677	2.692969	2.716249	2.718206	1.473751	1.444184	1.546933	1.345488	1.506435	1.191772	1.394676	1.453141	1.585412	1.311838	1.480633	2.239537	1.863841	1.582355	2.521889	2.460154	2.592014	2.575806	2.816637	2.858764	2.494994	3.068337	2.861372	3.34047	2.175996	1.696778	1.860562	1.654537	1.285238	1.210652	1.015188	0.91011	1.11409	1.286758	1.047825	0.800071	1.307353	1.138519	1.350845	1.455858	1.441877	1.311041	1.401308	1.314259	1.336229	1.41349	1.399126	1.454459	1.489199	1.652721	1.350327	1.808023	1.742999	1.832743	1.622685	1.78563	1.693727	1.79398	1.647075	1.837716	1.79414	1.878611	1.94758	1.834736	1.972279	1.970635	2.264778	1.881545	2.178631	2.602169	2.500089	2.837382	2.825085	2.806964	2.923083	2.894581	2.968518	2.064428	1.891093	2.165527	2.83135	3.769385	3.981791	4.156853	2.683093	2.564882	2.822815	2.318354	1.914291	1.42983	3.518993	3.577619	3.091506	3.271363	3.53834	3.276928	3.606499	3.251272	3.655934	3.390693	3.570037	3.391614	3.273107	3.223141	2.174581	2.105966	1.705942	1.871201	2.141587	2.349143	2.355046	2.485614	2.405754	2.652292	2.709332	2.750151	2.754016	3.07484	3.280935	2.975086	3.013273	3.07358	2.989121	3.230679	3.315549	2.995441	3.405316	3.280307	3.696253	3.830474	3.465136	3.868947	3.861395	3.935801	4.05411	4.006014	4.039399	3.982164	3.791263	3.83816	3.968226	3.875255	3.903599	3.847653	3.488136	3.791121	3.377113	3.229324	3.777927	3.830893	3.82772	3.90478	3.27615	3.885445	3.873099	3.773435	3.904537	3.988921	4.110446	4.488499	4.245838	4.355263	4.255057	3.970921	3.90133	3.888987	4.005717	2.950052	3.208981	3.899755	6.328724	3.824202	4.513047	5.455294	5.224986	3.640851	3.613381	3.642334	3.12964	2.874457	2.384268	2.422745	2.03584	1.881973	1.890797	1.559959	1.547085	1.384456	1.421835	1.3689	1.572464	2.564552	2.034348	1.938359	1.754484	2.049914	1.997939	2.004397	1.755421	2.033655	1.092865	2.430967	2.233602	2.566693	2.408008	2.947896	3.01911	3.504402	4.020505	3.413684	3.695372	3.608074	4.180622	4.092932	1.94134	3.116863	3.425187	3.186558	2.725957	2.717289	2.743878	3.33938	3.269981	3.027702	3.255405	2.910217	2.71203	2.410046	2.513022	2.43887	1.903971	1.520301	1.883268	1.520064	1.347288	1.324283	1.098534	1.32977	1.142124	1.381428	1.866591	1.476635	1.621764	1.448901	1.984072	2.31183	2.245032	2.202679	2.178377	2.29683	2.290049	2.528582	2.405509	2.581221	2.50146	2.521912	2.704181	1.541403	1.581988	1.32717	0.944778	1.319603	1.079948	1.883243	1.290012	1.510336	1.179025	1.329177	1.399776	1.146988	0.727066	1.028973	1.568123	1.299919	0.552584	0.88313	1.122093	1.145754	0.952377	0.841528	1.117876	1.192532	0.994264	1.092281	0.802539	0.966701	0.945302	0.839669	0.82117	0.8118	0.70935	1.416082	0.825442	1.127654	0.870107	1.03057	0.661831	1.19419	1.074857	2.965292	2.994679	3.483155	3.977105	3.527319	3.584723	3.779551	3.875942	3.600567	3.855506	4.534721	4.543124	4.197179	4.188201	4.354411	4.148646	3.67944	3.592273	3.583951	4.143906	3.768068	3.811982	3.534825	3.46987	3.622445	3.519138	3.546955	4.152086	4.233698	4.098076	4.226609	4.10048	4.016219	4.231931	3.951036	4.017774	4.015374	4.239259	4.360552	4.527007	4.537584	4.571744	4.599536	3.439116	3.449352	3.406891	3.350054	3.177855	3.336719	3.501829	3.510148	3.267031	3.14	3.255208	3.396066	4.010964	3.391588	3.446204	3.430027	3.410819	3.560547	3.827218	3.961372	2.940951	2.956639	3.011285	2.909924	3.327	3.394304	3.325797	3.419205	3.3964	3.467836	3.381985	3.08594	3.396044	2.802666	2.878521	2.958037	3.064197	3.073763	2.943718	3.111696	3.165858	3.158301	3.067134	3.040442	3.089068	2.958902	3.195707	3.401997	2.888825	3.070416	2.985022	2.905622	3.903234	3.663285	3.618768	3.543705	3.4646	3.534107	3.69838	3.439259	3.083101	3.433154	3.3138	3.498844	3.440372	3.608005	3.698932	4.018339	4.131731	4.180791	4.148828	4.155728	3.691199	3.756576	3.860385	4.221993	4.489002	4.708629	4.67741	5.184019	4.75228	3.56355	4.530603	5.220056	5.064104	5.098528	4.942626	4.936727	4.952022	5.250718	5.340725	5.268505	5.43427	5.251498	5.516991	5.683141	5.762604	6.281444	5.481649	5.443539	5.517257	4.900764	5.084162	4.996756	4.98568	5.07354	5.075486	4.97012	4.998746	3.631919	3.313503	3.220856	3.364994	3.061171	3.172841	2.971009	3.217174	3.683464	3.499517	3.555442	3.42756	3.349525	3.165912	2.552132	2.534502	2.46022	2.372138	2.391688	2.279393	2.893875	2.913217	2.81115	2.858261	3.182086	3.089216	3.344557	3.422959	3.371829	3.351497	3.351032	3.399583	3.361952	3.071977	2.916587	2.745957	3.567932	4.057832	3.058394	3.355602	3.294619	3.29417	3.331011	3.337424	3.37387	2.625779	3.97401	3.817033	3.944383	3.69299	3.89864	3.901337	3.699108	3.979422	3.114123	4.940323	5.648905	5.663221	5.642607	5.979863	5.225889	5.31036	5.310935	5.329918	6.191	6.107279	5.81668	4.174387	5.627681	5.595859	5.100356	4.228232	3.959525	4.234557	4.250259	4.906137	4.624911	4.962006	4.640088	4.827511	5.02429	4.748401	5.337503	5.697029	5.708959	6.742147	7.025272	7.365861	7.42683	7.414737	7.407938	7.396991	5.719552	6.127467	5.590704	5.991929	5.024087	5.906706	6.325128	5.228112	5.258957	5.183383	5.28181	5.796586	5.237657	4.74738	4.625186	5.339725	4.954325	5.116866	4.959014	4.926243	4.832664	4.777801	4.986087	4.865049	5.192405	5.69623	5.773444	5.776611	5.818697	5.169925	4.47837	3.504204	3.69271	4.805485	5.292095	5.492819	4.748109	4.838423	4.281321	4.136986	4.753604	4.787786	4.919221	4.548926	4.558351	4.484581	4.513301	4.51593	4.558095	4.261139	4.505914	4.203131	4.81465	5.064637	4.399539	5.064236	4.729761	5.073453	5.15251	4.716583	5.175491	5.142026	4.870511	4.998788	5.683374	5.580037	5.084833	5.264724	5.177861	5.327499	5.445266	5.745471	5.896691	5.818511	5.737836	5.808963	5.831892	4.289287	4.302827	4.685351	5.101524	5.368913	5.034078	5.312283	4.810819	5.299002	5.301533	4.834968	4.764399	4.429332	4.386843	4.569099	4.751589	4.568147	4.680156	4.641494	5.033071	4.801645	4.686751	4.917881	4.993577	4.665361	4.664663	4.898144	4.593645	4.44995	4.527467	4.653373	4.534402	4.360135	4.427435	4.349719	4.364284	4.270164	4.287885	4.384796	4.714857	4.717974	4.063407	3.903238	4.083903	4.296171	4.281324	4.29868	4.392782	4.554072	4.085237	4.106291	3.923147	3.941795	4.064781	4.681461	4.749017	4.764411	4.521055	4.889703	4.804172	4.819617	4.752538	4.333462	2.001899	1.57877	1.470016	1.487342	1.737041	1.657038	1.700079	1.574166	1.700285	1.494494	1.389988	1.16325	1.563042	1.093147	1.285999	1.455413	1.451257	1.397538	1.46887	1.615289	1.803784	1.523043	1.609957	2.243954	1.76831	1.900234	1.92181	1.63602	1.498977	1.543104	1.814916	2.081153	1.989369	1.842064	1.868173	1.861609	1.504473	1.537714	2.168964	1.523438	1.093632	1.540229	1.481677	1.567672	1.798434	1.558982	1.964692	1.565388	1.938698	2.358935	2.276719	1.779453	1.810623	2.580789	2.49404	2.506407	2.56216	2.852338	2.789829	2.558808	2.512607	2.892889	2.51855	2.488443	2.559369	2.604994	2.453356	2.543013	2.394414	1.929648	1.896476	1.453427	1.432928	1.320796	1.3122	1.189344	1.115706	1.083463	1.669057	1.626044	1.166139	0.666385	0.706729	0.770933	1.206218	0.750048	1.016491	0.857923	0.920604	0.95844	1.124106	1.328427	0.956711	0.50114	0.81074	0.494205	1.122775	1.398888	0.888205	1.231949	1.340765	1.113406	1.654064	1.346832	1.089813	1.003166	1.049264	0.897672	0.789974	1.152601	0.869533	0.895426	0.91173	0.722766	0.963131	1.054248	1.020831	1.000112	1.022657	0.49443	1.060554	0.995433	1.081958	0.793358	0.957781	0.60583	0.658096	1.03418	1.120607	1.023411	1.089842	1.046433	0.89623	1.288804	0.976767	0.811891	0.722923	0.789269	1.204334	0.80793	1.164508	0.816647	1.177643	1.075499	0.764959	0.959785	1.056175	1.586049	1.084024	1.02754	1.022406	1.056642	1.047421	0.918487	0.9361	0.903725	1.068459	0.640138	0.797024	1.541001	1.675339	1.579323	1.439218	1.46436	1.354629	1.519582	1.692756	1.931315	1.096204	0.907506	1.154111	1.325412	1.180731	0.89295	0.992595	1.092533	0.877584	0.985866	1.142519	1.019727	1.169729	1.007333	0.952033	1.04608	0.965102	1.13795	1.341343	1.256206	0.941407	0.96224	1.350961	1.342724	1.239198	1.528144	1.135076	1.067619	1.242077	1.37191	1.147791	0.964038	1.211392	0.952493	1.29469	0.95408	1.000118	1.40606	1.305875	1.192694	1.468268	1.510399	1.56984	1.594298	1.231853	1.343726	1.30935	1.477012	1.279015	1.284323	1.137793	0.992166	1.378889	1.199859	1.081696	1.306382	1.334317	1.033191	1.477541	1.4033	1.503211	1.179529	1.376544	1.452638	0.980626	0.93014	1.305832	1.325908	1.434419	1.571439	1.2883	1.428538	1.301554	1.390685	1.395323	1.412595	1.378376	1.304097	1.305792	1.273205	1.718516	1.822697	1.591531	1.58188	1.351864	1.461052	1.555077	1.485254	1.450432	1.473261	1.383119	1.588889	0.914861	0.805829	1.062221	0.521035	1.01823	1.193314	1.011273	1.12176	0.56462	0.792729	1.033569	0.885001	0.846469	1.231627	1.022295	1.066897	1.093819	1.231856	0.649191	0.550323	0.646316	1.362901	0.955253	1.021388	1.469559	1.319453	1.619879	1.264651	0.674415	0.635855	0.849385	0.914083	1.083027	1.262708	1.071772	1.228542	1.061983	1.051157	1.137339	1.275581	1.24506	1.266296	1.271727	1.231516	1.220919	0.602261	0.410016	0.61797	0.658997	0.645386	0.742233	0.80434	0.902178	0.841594	1.062508	1.060002	1.093648	1.006119	0.855254	0.99606	1.002667	0.76643	0.95201	1.261533	0.454808	0.953101	0.976167	0.869797	0.677704	1.003734	0.620628	0.575415	0.935557	1.051374	0.968535	1.338	0.877809	0.935406	1.095641	3.044181	3.00726	1.641647	1.539915	2.042289	1.585062	1.830448	1.625608	3.107925	3.12666	3.223167	3.193486	2.557767	2.427467	2.584942	2.353067	2.420533	2.82967	2.100901	2.100189	2.608456	2.611631	2.302992	2.618632	3.24161	2.922991	2.119153	2.052788	2.50678	2.004037	1.899648	1.53904	3.159504	3.423602	3.552795	3.123471	3.58185	2.69558	3.346004	3.347626	3.321325	3.313583	3.362911	3.549689	3.55765	2.038482	3.69057	3.21492	3.263426	3.937601	3.98244	4.017853	4.048299	4.425612	3.896843	3.816529	3.791627	3.762379	3.995593	3.826526	3.335443	3.12408	3.366642	3.407069	3.53001	3.000332	3.2943	3.041053	2.966684	2.96513	2.940803	3.029888	3.096614	3.010787	2.951308	2.924965	3.057771	3.051717	3.341321	3.498719	3.521909	2.832775	2.878584	2.884819	3.021115	3.024264	2.822809	2.586378	2.877696	2.560375	2.659782	2.631145	2.739778	2.718856	2.610422	2.753738	2.935463	2.946579	2.940594	2.641076	2.3463	2.344883	2.692081	2.563931	2.565752	2.522621	2.557485	2.480369	3.600757	3.991126	4.408179	4.126422	4.418662	4.547348	3.786165	3.695181	3.700932	3.743701	3.79342	3.760141	3.90296	3.962509	3.932102	4.036678	4.029628	3.997979	3.734698	3.997257	4.149939	3.998958	3.837906	3.980096	3.994007	4.147309	4.140488	3.922917	3.88538	3.739087	3.678093	3.843902	3.942858	3.8569	3.69058	3.726701	3.782284	3.741762	3.691323	3.73779	3.760424	3.643098	3.654814	3.598694	3.743131	3.478845	3.539205	3.35236	3.685176	3.83251	3.830079	4.098456	4.032622	4.014288	4.013793	4.129258	4.052956	4.03432	3.825639	3.851393	3.907384	4.009233	4.051068	4.10391	4.009329	3.960562	3.67539	3.785662	3.852302	3.080547	2.478201	2.450115	2.130141	1.765028	1.694884	1.852551	2.315806	1.818285	2.158008	2.176325	2.160224	3.057682	3.47028	3.409992	3.031156	3.010566	3.145917	3.00383	2.972189	2.991963	3.021807	2.653776	2.987536	2.975598	2.900719	2.918501	3.228798	3.161463	3.227897	3.122151	3.23143	3.352306	3.280025	3.240888	3.016151	3.239878	3.378292	3.27722	3.191641	3.393824	3.336473	3.399181	3.342429	3.389784	3.431328	3.395857	3.45507	3.386907	3.42297	3.391699	3.39347	3.392601	3.575548	3.533943	3.567446	2.641432	2.44262	2.114553	1.777585	1.678287	1.838757	1.64987	1.69321	1.810777	1.806837	1.552074	2.916905	2.749161	2.698883	2.398962	2.411459	2.544742	2.745947	2.83678	2.775417	2.687545	2.524639	2.52196	2.76374	3.03309	2.816177	3.143458	2.858815	3.103029	2.384924	2.137008	2.009938	2.106576	2.096122	2.023238	2.007612	2.016415	2.085326	2.422149	2.259172	2.071547	2.59429	2.106004	2.096224	2.124499	2.585988	2.594737	2.577883	2.584221	2.570674	2.553359	2.544423	2.535455	2.551197	2.547351	2.517564	2.548633	2.568474	2.564636	2.618986	2.645132	2.55106	2.571009	2.526495	2.553031	2.619064	2.672814	2.748129	2.60241	2.525513	2.390789	2.610879	2.648698	2.665768	2.599887	2.626445	2.616826	1.786561	1.752758	1.650547	2.054615	1.970385	2.164949	1.883635	1.451598	2.094596	2.662185	3.04403	2.972689	2.490982	2.034446	2.126305	2.494663	2.175247	2.898361	2.610798	2.283942	2.100763	2.172736	2.575074	2.072492	1.698915	1.679644	1.758548	1.540766	1.777746	1.964229	1.811303	1.88227	1.751036	1.985624	2.188818	1.954941	1.759871	1.720797	1.801659	1.663939	1.939096	1.893995	1.801717	1.832966	1.848196	1.804633	2.538405	2.037189	2.821399	3.378882	3.402447	3.390795	3.616603	3.618228	2.849498	2.898342	2.938736	2.967116	2.622414	2.463651	2.480098	2.530987	2.487429	2.431122	2.512463	2.436239	2.459522	2.713746	2.694735	2.714152	2.636504	2.45573	2.265047	2.389633	2.469485	2.43473	2.403421	2.400821	2.359274	2.441696	2.481539	2.589617	2.912863	2.767105	2.217683	2.55904	2.407748	2.429751	2.705071	1.816898	2.767092	2.463932	2.726691	2.271639	2.527388	2.182041	2.528708	3.168074	3.258669	3.122671	3.196541	3.498005	3.323163	2.675779	2.584372	2.113143	2.738973	2.371364	2.374592	2.730027	2.445243	2.80879	2.533891	2.641076	3.069448	3.945128	3.48608	3.508487	3.927548	3.382665	2.805639	2.705021	2.181479	2.138632	2.721134	2.292844	3.279215	3.210883	3.139264	3.132057	3.554256	3.188822	2.699931	2.295414	2.800597	2.81466	4.433173	1.397184	1.719093	1.913127	2.138264	1.993216	1.903645	2.14794	2.287846	2.247758	2.379382	2.524266	2.308939	2.480742	2.499216	2.799759	2.381136	2.356775	2.11956	1.776742	1.931844	1.982697	1.953078	1.929628	1.771919	2.015756	1.846199	1.124471	1.071664	0.882485	1.45261	1.316747	1.43067	0.943144	1.10183	1.107278	1.082378	1.006319	1.312942	1.174358	1.2695	1.030932	1.909992	1.914772	1.959812	1.855748	1.895497	2.023818	1.907408	1.89978	1.585737	1.610695	1.951128	1.738108	1.84509	1.842109	1.86458	2.880286	2.737796	2.036271	1.762114	2.024987	2.483432	2.348363	1.662598	1.92248	1.633804	1.64897	1.274123	2.170293	2.175575	1.426258	1.318585	1.805702	1.823975	2.08184	1.699847	1.702106	1.530561	1.598902	1.723258	2.303382	1.761552	3.063341	1.64465	2.085041	1.086545	1.868912	1.268485	3.206292	1.976815	2.758743	1.540398	1.143778	1.702466	1.429001	1.779498	1.669594	1.953022	1.992236	2.246294	1.955615	2.596272	1.754435	2.106335	1.959561	1.491381	1.900123	2.086418	2.556344	2.418479	2.102166	1.978019	2.663501	1.932341	1.531493	1.727647	1.594064	2.02783	2.024775	1.489725	1.601669	1.681863	1.962146	2.316045	1.862737	1.304942	1.829519	2.220902	1.897085	1.967048	1.898283	1.663221	1.625396	1.861048	0.979001	1.247642	0.90893	2.023488	1.941846	1.872489	1.92293	1.794402	2.195959	1.995927	2.09712	1.952561	3.131706	1.863816	1.240044	2.342719	2.302001	2.349208	2.089456	1.696148	2.236627	1.767458	1.894641	1.738399	2.454647	2.515546	2.684496	2.546375	2.129842	2.364465	1.640314	1.409969	1.322446	1.102633	0.973488	1.472323	1.113738	1.375691	1.481579	1.210797	1.018274	1.170368	1.00272	1.652912	1.607633	1.391675	1.668125	1.815265	1.889737	1.632391	1.59918	1.60146	2.081546	2.033069	2.018792	2.17997	2.178135	1.741803	1.903483	1.998012	1.564835	1.815908	1.459024	1.637448	1.758174	1.74428	1.798407	1.772692	1.754942	1.570446	1.424266	1.360086	1.517918	1.528473	1.594205	1.360996	1.406422	1.746113	1.71965	1.637867	1.44489	1.672059	1.720909	1.794521	1.891561	1.992588	1.940881	1.770101	1.556848	1.696256	1.772899	1.169367	2.002531	1.924383	2.198287	1.873368	1.385447	2.073173	1.800007	1.466019	1.822308	1.823853	1.740455	1.553563	1.898027	1.726659	1.928933	1.441906	1.609996	1.576838	1.428858	1.565301	1.332564	1.430038	1.240117	1.468241	1.27026	1.369215	1.495837	2.223442	1.470386	1.648409	1.628353	1.508848	1.418949	1.164689	1.532683	1.758491	1.612448	1.610978	1.693305	1.900156	1.501576	1.604071	1.797371	1.599293	1.591632	1.645477	1.300615	1.307882	1.707649	1.355027	1.417227	1.52338	1.411784	1.555471	1.637536	1.897719	1.749276	1.377248	1.320646	1.495031	1.817703	2.106496	1.784606	2.524566	2.802686	2.065375	1.946157	1.973803	1.873459	1.753765	2.319649	2.681222	2.642835	2.685164	2.581282	3.187069	3.193465	0.83629	0.913617	1.026765	1.612573	1.580452	1.923644	2.234628	2.127584

void	LoadRatesFromStr(char *Str, RATES *Rates, OPTIONS *Opt, TREES *Trees)
{
	FATTAILRATES *FTR;
	TREE *Tree;
	int Index, Tokes, Pos;
	char *Buffer;
	char **Passed;
	NODE N;


	Tree = Trees->Tree[Rates->TreeNo];

	if(Trees->NoSites != 1)
	{
		printf("LoadRatesFromStr is only avalable for 1 site.\n");
		exit(0);
	}
	
	Buffer = StrMake(Str);
	Passed = (char**)SMalloc(sizeof(char*) * strlen(Str));

	Tokes = MakeArgv(Buffer, Passed, (int)strlen(Str));

	Rates->Rates[0] = atof(Passed[2]);

	FTR = Rates->FatTailRates;
	Pos = 3;
	for(Index=0;Index<Tree->NoNodes;Index++)
	{
		N = Tree->NodeList[Index];
		if(N->Tip == FALSE)
			N->FatTailNode->Ans[0] = atof(Passed[Pos++]);
		
	}

	FatTailGetAnsSates(Tree, 1, FTR);
	
//	for(Index=0;Index<)

	free(Buffer);
	free(Passed);

	
	Rates->Lh = Likelihood(Rates, Trees, Opt);
}

