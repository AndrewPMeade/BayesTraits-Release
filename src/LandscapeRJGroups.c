#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_cdf.h>

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

RATE_GROUP*	GetRateGroup(LAND_RATE_GROUPS* LandRGroup, int PNo, int *Pos)
{
	int GRIndex, Index;
	RATE_GROUP*	RateG;

	for(GRIndex=0;GRIndex<LandRGroup->NoGroups;GRIndex++)
	{
		RateG = LandRGroup->RateGroupList[GRIndex];

		for(Index=0;Index<RateG->NoRates;Index++)
		{
			if(RateG->PramList[Index] == PNo)
			{
				*Pos = Index;
				return RateG;
			}
		}
	}

	Pos = NULL;
	return NULL;
}


void	FixedBetaDiscretised(LAND_RATE_GROUPS* LandRGroup, RANDSTATES *RS)
{
	int NoGroup, Index, GNo;
	double CDF_P, X, P, Sig;
	RATE_GROUP*	RateG;


	Sig = sqrt(2.0);
	NoGroup = 11;
	LandRGroup->NoGroups = NoGroup;

	for(Index=0;Index<NoGroup;Index++)
	{
		CDF_P = ((Index+1) - 0.5) / (double)NoGroup;
	
		X = gsl_cdf_gaussian_Pinv(CDF_P, Sig );
		P = gsl_ran_gaussian_pdf(X, Sig);

		LandRGroup->RateGroupList[Index]->Rate = X;
		LandRGroup->RateGroupList[Index]->NoRates = 0;
	}

	for(Index=0;Index<LandRGroup->NoRates;Index++)
	{
		GNo = RandUSInt(RS) % LandRGroup->NoGroups;

		GNo = 5;
		
		RateG = LandRGroup->RateGroupList[GNo];

		RateG->PramList[RateG->NoRates++] = Index;
	}
}

LAND_RATE_GROUPS*	CreateLandRateGroups(RATES *Rates, int NoRates)
{
	LAND_RATE_GROUPS* Ret;
	int Index;

	Ret = AllocLandRateGroup(NoRates);

	Ret->NoRates = NoRates;
	
	for(Index=0;Index<Ret->NoRates;Index++)
		Ret->RateGroupList[Index] = AllocRateGroup(NoRates);

	Ret->NoGroups = 1;

	SetOneGroup(Ret, 1.0);

	FixedBetaDiscretised(Ret, Rates->RS);

//	PrintLandRateGroups(stdout, Ret);
//	exit(0);

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

	LandRateGroups = Rates->LandscapeRateGroups;
	Tree = Trees->Tree[0];

	Tree->Root->LandscapeBeta = 0.0;

	for(GIndex=0;GIndex<LandRateGroups->NoGroups;GIndex++)
		MapLandRateGroup(Tree, LandRateGroups->RateGroupList[GIndex]);
}

void	RemovePramFromGroup(RATE_GROUP *RGroup, int Pos)
{
	RGroup->PramList[Pos] = RGroup->PramList[RGroup->NoRates-1];
	RGroup->NoRates--;
}

void	AddPramToGroup(RATE_GROUP *RGroup, int Pram)
{
	RGroup->PramList[RGroup->NoRates] = Pram;
	RGroup->NoRates += 1;
}

void	SwapRateGroupParam(RATES *Rates)
{
	int Rate, Pos;
	LAND_RATE_GROUPS *LandRateGroups;
	RATE_GROUP	*RGOld, *RGNew; 
	
	LandRateGroups = Rates->LandscapeRateGroups;


	Rate = RandUSInt(Rates->RS) % LandRateGroups->NoRates;
	RGOld = GetRateGroup(LandRateGroups, Rate, &Pos);

	RemovePramFromGroup(RGOld, Pos);
	
	do
	{
		Pos = RandUSInt(Rates->RS) % LandRateGroups->NoGroups;
		RGNew = LandRateGroups->RateGroupList[Pos];
	}while(RGOld == RGNew);

	AddPramToGroup(RGNew, Rate);
}

void	PrintLandRateGroups(FILE *File, LAND_RATE_GROUPS* LandRateGroups)
{
	int GIndex, Index;
	RATE_GROUP	*RateG;

	fprintf(File, "NoRateGroups:\t%d\n", LandRateGroups->NoGroups);
	fprintf(File, "NoRates:\t%d\n", LandRateGroups->NoRates);

	for(GIndex=0;GIndex<LandRateGroups->NoGroups;GIndex++)
	{
		RateG = LandRateGroups->RateGroupList[GIndex];

		fprintf(File, "%d\t%d\t%f\t", GIndex, RateG->NoRates, RateG->Rate);
		for(Index=0;Index<RateG->NoRates;Index++)
			fprintf(File, "%d\t", RateG->PramList[Index]);
		fprintf(File, "\n");
	}
}

void	RecPrintNodeTaxa(FILE *Str, NODE N)
{
	int Index;

	if(N->Tip == TRUE)
	{
		fprintf(Str, "%s\t", N->Taxa->Name);
		return;
	}

	for(Index=0;Index<N->NoNodes;Index++)
		RecPrintNodeTaxa(Str, N->NodeList[Index]);

}

void	InishLandRateGroupsOutput(OPTIONS *Opt, TREES *Trees, RATES *Rates)
{
	int Index;
	TREE *Tree;
	LAND_RATE_GROUPS *LandRateG;
	NODE N;

	LandRateG = Rates->LandscapeRateGroups;

		
	Opt->LogLandscapeGroups = OpenWriteWithExt(Opt->BaseOutputFN, OUTPUT_EXT_LAND_GROUP);

	fprintf(Opt->LogLandscapeGroups, "NoRates:\t%d\n", LandRateG->NoRates);

	Tree = Trees->Tree[0];
	for(Index=1;Index<Tree->NoNodes;Index++)
	{
		N = Tree->NodeList[Index];
		fprintf(Opt->LogLandscapeGroups, "%d\t%f\t", Index-1, N->UserLength);
		RecPrintNodeTaxa(Opt->LogLandscapeGroups, N);		
		fprintf(Opt->LogLandscapeGroups, "\n");
	}

	fprintf(Opt->LogLandscapeGroups, "It\tLh\tNoGroups\tGroupID\tBeta\tNoNodes\tNodes\t");
	fprintf(Opt->LogLandscapeGroups, "\n");

	fflush(Opt->LogLandscapeGroups);
}

void	PrintRateGroupPramList(FILE *Out, RATE_GROUP *RateG)
{
	int Index;
	
	if(RateG->NoRates == 0)
	{
		fprintf(Out, "|\t");
		return;
	}

	for(Index=0;Index<RateG->NoRates-1;Index++)
		fprintf(Out, "%d|", RateG->PramList[Index]);
	fprintf(Out, "%d\t", RateG->PramList[Index]);
}

void	OutputLandRateGroupsOutput(long long Itter, OPTIONS *Opt, TREES *Trees, RATES *Rates)
{
	LAND_RATE_GROUPS *LandRateG;
	FILE *Out;
	RATE_GROUP	*RateG;
	int GIndex, Index;

	Out = Opt->LogLandscapeGroups;
	
	LandRateG = Rates->LandscapeRateGroups;
	fprintf(Out, "%lld\t%f\t%d\t", Itter, Rates->Lh, LandRateG->NoGroups);
	

	for(GIndex=0;GIndex<LandRateG->NoGroups;GIndex++)
	{
		RateG = LandRateG->RateGroupList[GIndex];
		fprintf(Out, "%d\t%f\t%d\t", GIndex, RateG->Rate, RateG->NoRates);
		PrintRateGroupPramList(Out, RateG);

	}

	fprintf(Out, "\n");	
	fflush(Out);
}