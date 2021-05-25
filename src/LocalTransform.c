#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "typedef.h"
#include "VarRates.h"
#include "genlib.h"
#include "trees.h"


LOCAL_TRANSFORM*	CreateLocalTransforms(char *Name, TAG **TagList, int NoTags, TRANSFORM_TYPE Type, int Est, double Scale)
{
	LOCAL_TRANSFORM* Ret;

	Ret = (LOCAL_TRANSFORM*)SMalloc(sizeof(LOCAL_TRANSFORM));


	Ret->Name = StrMake(Name);
	Ret->TagList = CloneMem(sizeof(TAG**) * NoTags, TagList);
	Ret->NoTags = NoTags;	
	Ret->Type = Type;
	Ret->Est = Est;
	Ret->Scale = Scale;

	return Ret;
}

void		FreeLocalTransforms(LOCAL_TRANSFORM* LTrans)
{
	free(LTrans->Name);
	free(LTrans->TagList);
	free(LTrans);
}

void		CopyLocalTransforms(LOCAL_TRANSFORM* ATrans, LOCAL_TRANSFORM* BTrans)
{
	ATrans->Scale	= BTrans->Scale;
}

LOCAL_TRANSFORM*	CloneLocalTransform(LOCAL_TRANSFORM* LTrans)
{
	return CreateLocalTransforms(LTrans->Name, LTrans->TagList, LTrans->NoTags, LTrans->Type, LTrans->Est, LTrans->Scale);
}

int			NoEstLocalTransform(LOCAL_TRANSFORM** List, int NoRates)
{
	int Index, Ret;

	Ret = 0;
	for(Index=0;Index<NoRates;Index++)
		if(List[Index]->Est == TRUE)
			Ret++;

	return Ret;
}

int			EstLocalTransforms(LOCAL_TRANSFORM** List, int NoRates)
{

	if(NoEstLocalTransform(List, NoRates) > 0)
		return TRUE;

	return FALSE;
}

void		PrintLocalTransform(FILE *Str, LOCAL_TRANSFORM* Trans)
{
	fprintf(Str, "    %s %s ", Trans->Name, VarRatesTypeToStr(Trans->Type));	
	if(Trans->Est == TRUE)
		fprintf(Str, "Estimate");
	else
		fprintf(Str, "%f", Trans->Scale);

	fprintf(Str, "\n");
}

void		PrintLocalTransforms(FILE *Str, LOCAL_TRANSFORM** List, int NoRates)
{
	int Index;

	if(NoRates == 0)
		return;

	fprintf(Str, "Local Rates:\n");
	for(Index=0;Index<NoRates;Index++)
		PrintLocalTransform(Str, List[Index]);
}

void	ApplyLocalTransforms(RATES *Rates, TREES *Trees, OPTIONS *Opt, int Norm)
{
	LOCAL_TRANSFORM *LRate;
	int TIndex, Index;
	TREE *Tree;
	NODE N;

	Tree = Trees->Tree[Rates->TreeNo];
	
	for(Index=0;Index<Rates->NoLocalTransforms;Index++)
	{
		LRate = Rates->LocalTransforms[Index];

		for(TIndex=0;TIndex<LRate->NoTags;TIndex++)
		{
			N = LRate->TagList[TIndex]->NodeList[Rates->TreeNo];
			VarRatesNode(Trees, Tree, N, LRate->Scale, LRate->Type);
		}
	}

//	SaveTrees("sout.trees", Trees);exit(0);
}

double	ChangeLocalScale(RANDSTATES	*RS, double Scale, double Dev)
{
	double		Ret;
	
	do
	{
		Ret = ((RandDouble(RS) * Dev) - (Dev / 2.0)) + Scale; 
	} while(Ret <= 0);
	
	return Ret;
}

LOCAL_TRANSFORM*	GetEstRate(RATES *Rates)
{
	int Ret;
	do
	{
		Ret = RandUSInt(Rates->RS) % Rates->NoLocalTransforms;
	}while(Rates->LocalTransforms[Ret]->Est == FALSE);

	return Rates->LocalTransforms[Ret];
}

void		ChangeLocalTransform(OPTIONS *Opt, TREES *Trees, RATES *Rates, SCHEDULE *Shed)
{
	double NRate, Dev;
	LOCAL_TRANSFORM *LRate;

	Shed->CurrentAT = Shed->LocalRatesAT;
	Dev = Shed->CurrentAT->CDev;
	
	LRate = GetEstRate(Rates); 

	NRate = ChangeLocalScale(Rates->RS, LRate->Scale, Dev);

	if(NRate > MAX_LOCAL_RATE || NRate < MIN_LOCAL_RATE)
		NRate = LRate->Scale;

	LRate->Scale = NRate;
}

int	GetNoTransformType(TRANSFORM_TYPE TType, RATES *Rates)
{
	int Ret, Index;

	Ret = 0;

	for(Index=0;Index<Rates->VarRates->NoNodes;Index++)
		if(Rates->VarRates->NodeList[Index]->Type == TType)
			Ret++;

	return Ret;
}