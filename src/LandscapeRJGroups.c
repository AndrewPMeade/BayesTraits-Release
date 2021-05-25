#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_cdf.h>

#include "TypeDef.h"
#include "GenLib.h"
#include "LandscapeRJGroups.h"
#include "Priors.h"
#include "Rates.h"

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
	Ret->Sig = -1;

	return Ret;
}

void	EmptyLandRateGroup(LAND_RATE_GROUPS* LandRGroup)
{
	int Index;

	for(Index=0;Index<LandRGroup->NoRates;Index++)
		LandRGroup->RateGroupList[Index]->NoRates = 0;
}

void	SetOneGroup(LAND_RATE_GROUPS* LandRateGroups, double RateValue)
{
	RATE_GROUP* RateG;
	int Index;

	EmptyLandRateGroup(LandRateGroups);

	LandRateGroups->NoGroups = 1;

	RateG = LandRateGroups->RateGroupList[0];

	RateG->NoRates = LandRateGroups->NoRates;
	RateG->Rate = RateValue;

	for(Index=0;Index<RateG->NoRates;Index++)
		RateG->PramList[Index] = LandRateGroups->RateList[Index];
	
}

void	FreeLandRateGroups(LAND_RATE_GROUPS* LandRateG)
{
	int Index;
	for(Index=0;Index<LandRateG->NoRates;Index++)
		FreeRateGroup(LandRateG->RateGroupList[Index]);
	free(LandRateG->RateGroupList);
	free(LandRateG->RateList);
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

void	SetDiscretisedDist(LAND_RATE_GROUPS* LandRGroup)
{
	int Index;
	double CDF_P, X;

	for(Index=0;Index<LandRGroup->NoGroups;Index++)
	{
		CDF_P = ((Index+1) - 0.5) / (double)LandRGroup->NoGroups;

		X = gsl_cdf_gaussian_Pinv(CDF_P, LandRGroup->Sig);

		LandRGroup->RateGroupList[Index]->Rate = X;
	}
}

void	FixedBetaDiscretised(RATES *Rates, LAND_RATE_GROUPS* LandRGroup, RANDSTATES *RS, int NoGroup)
{
	int Index, GNo;
	double P, Sig;
	RATE_GROUP*	RateG;
	PRIOR *Prior;

	EmptyLandRateGroup(LandRGroup);
	
	Prior = GetPriorFromName("RJ_Landscape_Rate_Group", Rates->Priors, Rates->NoPriors);
	Sig = Prior->DistVals[1];

	LandRGroup->NoGroups = NoGroup;

	LandRGroup->Sig = Sig;
	SetDiscretisedDist(LandRGroup);


	for(Index=0;Index<LandRGroup->NoRates;Index++)
	{
		GNo = RandUSInt(RS) % LandRGroup->NoGroups;

//		GNo = NoGroup / 2;
		
		RateG = LandRGroup->RateGroupList[GNo];

		RateG->PramList[RateG->NoRates++] = LandRGroup->RateList[Index];
	}
}


void SetNoFixedGroups(RATES *Rates, LAND_RATE_GROUPS* LandRGroup, int NoGroup)
{
	int Index, GNo;
	RATE_GROUP*	RateG;
	PRIOR *Prior;
	
	EmptyLandRateGroup(LandRGroup);

	LandRGroup->NoGroups = NoGroup;
	

	for(Index=0;Index<LandRGroup->NoRates;Index++)
	{
		GNo = RandUSInt(Rates->RS) % LandRGroup->NoGroups;
		RateG = LandRGroup->RateGroupList[GNo];
		RateG->PramList[RateG->NoRates++] = LandRGroup->RateList[Index];
	}
	
	RateG = LandRGroup->RateGroupList[0];
	RateG->Rate = 0.0;

	Prior = GetPriorFromName("RJ_Landscape_Rate_Group", Rates->Priors, Rates->NoPriors);
		
	for(Index=1;Index<LandRGroup->NoGroups;Index++)
	{ 
		RateG = LandRGroup->RateGroupList[Index];
		RateG->Rate = RandFromPrior(Rates->RNG, Prior);
	}
}

int*	CreateRateList(TREE *Tree, int *NoRates)
{
	int *Ret, Index;
	NODE N;

	Ret = (int*)SMalloc(sizeof(int) * Tree->NoNodes);
	*NoRates = 0;

	for(Index=0;Index<Tree->NoNodes;Index++)
	{
		N = Tree->NodeList[Index];
		if(N->Tip == FALSE && N != Tree->Root)
		{
			Ret[*NoRates] = Index;
			(*NoRates)++;
		}
	}

	return Ret;
}

LAND_RATE_GROUPS*	CreateLandRateGroups(RATES *Rates, TREES *Trees)
{
	LAND_RATE_GROUPS* Ret;
	int Index, NoRates, *RateList;
		
	RateList = CreateRateList(Trees->Tree[0], &NoRates);

	Ret = AllocLandRateGroup(NoRates);

	Ret->NoRates = NoRates;
	Ret->RateList = RateList;
	
	for(Index=0;Index<Ret->NoRates;Index++)
		Ret->RateGroupList[Index] = AllocRateGroup(NoRates);

	Ret->NoGroups = 1;
	SetOneGroup(Ret, 0.0);

	FixedBetaDiscretised(Rates, Ret, Rates->RS, 11);
	
//	SetNoFixedGroups(Rates, Ret, 2);
	
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

	A->Sig = B->Sig;

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
	Rate = LandRateGroups->RateList[Rate];
	RGOld = GetRateGroup(LandRateGroups, Rate, &Pos);

	RemovePramFromGroup(RGOld, Pos);
	
//	do
//	{
		Pos = RandUSInt(Rates->RS) % LandRateGroups->NoGroups;
		RGNew = LandRateGroups->RateGroupList[Pos];
//	}while(RGOld == RGNew);

	AddPramToGroup(RGNew, Rate);
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
	for(Index=0;Index<Tree->NoNodes;Index++)
	{
		N = Tree->NodeList[Index];
		fprintf(Opt->LogLandscapeGroups, "%d\t%f\t", Index, N->UserLength);
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

double CalcPriorLandRateGoup(RATES *Rates)
{
	double Ret, BetaPrior;
	PRIOR *Prior;
	RATE_GROUP	*RateG;
	int Index;
	LAND_RATE_GROUPS* LandRateGroups;
	
	return 0;

	LandRateGroups = Rates->LandscapeRateGroups;

	Prior = GetPriorFromName("RJ_Landscape_Rate_Group", Rates->Priors, Rates->NoPriors);
	
//	LandRateGroups->Sig = 3.4761097485661408e-17;
	Prior->DistVals[1] = LandRateGroups->Sig;

	
	Ret = 0;
	for(Index=0;Index<LandRateGroups->NoGroups;Index++)
	{
		RateG = LandRateGroups->RateGroupList[Index];
		BetaPrior = CalcLhPriorP(RateG->Rate, Prior);

//		printf("%d\t%f\t%f\n", Index, RateG->Rate, BetaPrior);

		if(BetaPrior == ERRLH)
			return ERRLH;

		if(RateG->Rate != 0)
	//		Ret += BetaPrior * RateG->NoRates;
			Ret += BetaPrior;
	}
//	exit(0);
	return Ret;
}

void ChangeLandRateGoup(RATES *Rates,  SCHEDULE* Shed)
{
	LAND_RATE_GROUPS* LandRateGroups;
	RATE_GROUP	*RateG;
	int RateGNo;
	double Dev;

	Shed->CurrentAT = Shed->LandscapeRateChangeAT; 

	Dev = Shed->CurrentAT->CDev;

	LandRateGroups = Rates->LandscapeRateGroups;

	if(LandRateGroups->NoGroups < 2)
		return;

	RateGNo = RandUSInt(Rates->RS) % (LandRateGroups->NoGroups - 1);
	RateGNo += 1;

	RateG = LandRateGroups->RateGroupList[RateGNo];

	RateG->Rate += RandNormal(Rates->RS, 0, Dev);
}

void	ChangeLandRateGroupSigDist(RATES *Rates,  SCHEDULE* Shed)
{
	LAND_RATE_GROUPS* LandRateGroups;
	double Dev;

	Shed->CurrentAT = Shed->LandscapeRateChangeAT; 

	Dev = Shed->CurrentAT->CDev;

	LandRateGroups = Rates->LandscapeRateGroups;

	LandRateGroups->Sig = ChangeRateExp(LandRateGroups->Sig, Dev, Rates->RS, &Rates->LnHastings);

	SetDiscretisedDist(LandRateGroups);
}