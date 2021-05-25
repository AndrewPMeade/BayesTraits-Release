#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

#include "TypeDef.h"
#include "GenLib.h"
#include "LandscapeRJGroups.h"

void FreeRateGroup(RATE_GROUP *RateGroup)
{
	free(RateGroup->PramList);
	free(RateGroup);
}

RATE_GROUP* AllocRateGroup(int NoRates)
{
	RATE_GROUP*	Ret;
	
	Ret = (RATE_GROUP*)SMalloc(sizeof(RATE_GROUP));
	Ret->NoRates = 0;
	Ret->Rate = -1;
	Ret->PramList = (int*)SMalloc(sizeof(int) * NoRates);
	
	return Ret;
}

LAND_RATE_GROUPS*	AllocLandRateGroup(int NoRates)
{
	LAND_RATE_GROUPS* Ret;

	Ret = (LAND_RATE_GROUPS*)SMalloc(sizeof(LAND_RATE_GROUPS));

	Ret->NoGroups = 0;
	Ret->RateGroupList = (RATE_GROUP**)SMalloc(sizeof(RATE_GROUP*) * NoRates);
	
	return Ret;
}

void	SetOneGroup(LAND_RATE_GROUPS* LandRateGroups, double RateValue)
{
	RATE_GROUP* RateG;
	int Index;

	LandRateGroups->NoGroups = 1;

	RateG = LandRateGroups->RateGroupList[0];

	RateG->NoRates = LandRateGroups->NoRates;
	RateG->Rate = RateValue;

	for(Index=0;Index<RateG->NoRates;Index++)
		RateG->PramList[Index] = Index;
	
}

void	FreeLandRateGroups(LAND_RATE_GROUPS* LandRateG)
{
	int Index;
	for(Index=0;Index<LandRateG->NoRates;Index++)
		FreeRateGroup(LandRateG->RateGroupList[Index]);
	free(LandRateG->RateGroupList);
	free(LandRateG);
}

LAND_RATE_GROUPS*	CreateLandRateGroups(int NoRates)
{
	LAND_RATE_GROUPS* Ret;
	int Index;

	Ret = AllocLandRateGroup(NoRates);

	Ret->NoRates = NoRates;
	
	for(Index=0;Index<Ret->NoRates;Index++)
		Ret->RateGroupList[Index] = AllocRateGroup(NoRates);

	Ret->NoGroups = 1;

	SetOneGroup(Ret, 1.0);

	return Ret;
}

void	CopyRateGroup(RATE_GROUP *A, RATE_GROUP *B)
{
	A->Rate = B->Rate;
	A->NoRates = B->NoRates;

	memcpy(A->PramList, B->PramList, sizeof(int) * A->NoRates);
}

void				CopyLandRateGroups(LAND_RATE_GROUPS* A, LAND_RATE_GROUPS* B)
{
	int Index;

	A->NoGroups = B->NoGroups;
	A->NoRates = B->NoRates;

	for(Index=0;Index<B->NoGroups;Index++)
		CopyRateGroup(A->RateGroupList[Index], B->RateGroupList[Index]);
}

void	MapLandRateGroup(TREE *Tree, RATE_GROUP *RateGroup)
{
	NODE N;
	int Index, Pos;

	for(Index=0;Index<RateGroup->NoRates;Index++)
	{
		Pos = RateGroup->PramList[Index];
		Pos += 1;

		N = Tree->NodeList[Pos];
		N->LandscapeBeta = RateGroup->Rate;
	}
}

void	MapLandRateGroups(OPTIONS *Opt, TREES *Trees, RATES *Rates)
{
	TREE *Tree;
	int		GIndex;
	LAND_RATE_GROUPS *LandRateGroups;

	LandRateGroups = Rates->RJLandscapeRateGroups;
	Tree = Trees->Tree[0];

	Tree->Root->LandscapeBeta = 0.0;

	for(GIndex=0;GIndex<LandRateGroups->NoGroups;GIndex++)
		MapLandRateGroup(Tree, LandRateGroups->RateGroupList[GIndex]);
}

void	LandRateGroupsSplit(RATES *Rates, LAND_RATE_GROUPS *LandRateGroups)
{
	int GroupNo, Pos;
	RATE_GROUP	*RateGroup;

	GroupNo = RandInt(Rates->RS) % LandRateGroups->NoGroups;
	RateGroup = LandRateGroups->RateGroupList[GroupNo];

	ShuffleIntList(Rates->RS, RateGroup->PramList, RateGroup->NoRates);

	Pos = RandUSInt(Rates->RS) % (RateGroup->NoRates - 1);
	Pos += 1;

	
}

void	LandRateGroupsMurge(RATES *Rates, LAND_RATE_GROUPS *LandRateGroups)
{

}


void	LandRateGroupsSplitMurge(OPTIONS *Opt, TREES *Trees, RATES *Rates)
{
	LAND_RATE_GROUPS *LandRateGroups;

	LandRateGroups = Rates->RJLandscapeRateGroups;

	if(LandRateGroups->NoGroups == 1)
	{
		LandRateGroupsSplit(LandRateGroups);
		return;
	}

	if(LandRateGroups->NoGroups == LandRateGroups->NoRates)
	{
		LandRateGroupsMurge(LandRateGroups);
		return;
	}
	
	if(RandDouble(Rates->RS) < 0.5)
		LandRateGroupsSplit(Rates, LandRateGroups);
	else
		LandRateGroupsMurge(Rates, LandRateGroups);

	
}

