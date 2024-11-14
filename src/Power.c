



#include <math.h>
#include <stdlib.h>
#include <string.h>

#include "Power.h"
#include "TypeDef.h"
#include "GenLib.h"
#include "Priors.h"

int		GetNoPowerSites(OPTIONS *Opt)
{
	int Index, No;

	No = 0;
	for(Index=0;Index<Opt->NoOfSites;Index++)
		if(Opt->PowerSites[Index] == TRUE)
			No += 1;

	return No;
}

void	PrintPowerOpt(FILE* Str, OPTIONS *Opt)
{
	int Index, No;

	No = GetNoPowerSites(Opt);

	if(No == 0)
		return;

	fprintf(Str, "Power Estimates For Sites:       ");
	for(Index=0;Index<Opt->NoOfSites;Index++)
		if(Opt->PowerSites[Index] == TRUE)
			fprintf(Str, "%d ", Index);
	fprintf(Str, "\n");
}

SITE_POWER*	AllocSitePowers(OPTIONS *Opt)
{
	int NoSites;
	SITE_POWER*	Ret;

	NoSites = GetNoPowerSites(Opt);
	
	Ret = (SITE_POWER*)SMalloc(sizeof(SITE_POWER));

	Ret->NoSites = GetNoPowerSites(Opt);
	Ret->SiteIndex = (int*)SMalloc(sizeof(int) * Ret->NoSites);
	Ret->Powers = (double*)SMalloc(sizeof(double) * Ret->NoSites);

	return Ret;
}

SITE_POWER*	CrateSitePowers(OPTIONS *Opt)
{
	int NoSites, Index;
	SITE_POWER*	SitePower;

	NoSites = GetNoPowerSites(Opt);
	if(NoSites == 0)
		return NULL;

	SitePower = AllocSitePowers(Opt);

	NoSites = 0;
	for(Index=0;Index<Opt->NoOfSites;Index++)
	{
		if(Opt->PowerSites[Index] == TRUE)
		{
			SitePower->SiteIndex[NoSites] = Index;
			SitePower->Powers[NoSites] = 1.0;

			NoSites++;
		}
	}

	return SitePower;
}

void	FreeSitePowers(SITE_POWER* SitePower)
{
	free(SitePower->Powers);
	free(SitePower->SiteIndex);
	free(SitePower);
}



void	PrintDataVals(RATES* Rates, TREES *Trees, OPTIONS *Opt)
{
	int Index, SIndex, NoSites;
	NODE Node;
	CONDATA* Con;
	TREE *Tree;
	SITE_POWER* SitePowers;
//	return;

	SitePowers = Rates->SitePowers;
	NoSites = Trees->NoSites;
	Tree = Trees->Tree[Rates->TreeNo];

	for(Index=0;Index<Tree->NoNodes;Index++)
	{
		Node = Tree->NodeList[Index];
		if(Node->Tip == TRUE)
		{
			printf("Taxa:\t%s\t", Node->Taxa->Name);
			Con = Node->ConData;
			for(SIndex=0;SIndex<Trees->NoSites;SIndex++)
				printf("%f\t",Con->Contrast[0]->Data[SIndex]);
			printf("\n");
		}
	}
}


void SetSitePower(RATES* Rates, TREES *Trees, OPTIONS *Opt)
{
	int Index, SIndex, NoSites, SiteIndex;
	NODE Node;
	CONDATA* Con;
	TREE *Tree;
	SITE_POWER* SitePowers;

	SitePowers = Rates->SitePowers;
	NoSites = Trees->NoSites;
	Tree = Trees->Tree[Rates->TreeNo];

	for(Index=0;Index<Tree->NoNodes;Index++)
	{
		Node = Tree->NodeList[Index];
		if(Node->Tip == TRUE)
		{
			Con = Node->ConData;
			for(SIndex=0;SIndex<SitePowers->NoSites;SIndex++)
			{
				SiteIndex = SitePowers->SiteIndex[SIndex];
				Con->Contrast[0]->Data[SiteIndex] = pow(Con->Contrast[0]->Data[SiteIndex], SitePowers->Powers[SIndex]);
			}
		}
	}
}

void CopySitePowers(SITE_POWER* A, SITE_POWER* B)
{
	memcpy(A->Powers, B->Powers, sizeof(double) * A->NoSites);
}


void	ChangeSitePowers(RATES *Rates, SCHEDULE* Shed)
{
	SITE_POWER* SitePowers;
	int Pos;

	SitePowers = Rates->SitePowers;

	Shed->CurrentAT = Shed->SitePower;
		
	Pos = (int)gsl_rng_uniform_int(Rates->RNG, SitePowers->NoSites);

	SitePowers->Powers[Pos] += gsl_ran_gaussian(Rates->RNG, Shed->CurrentAT->CDev);

//	printf("%f\n", SitePowers->Powers[Pos]);

	Rates->LnHastings = 0.0;
}

double	CalcSitePowersPrior(RATES *Rates)
{
	int Index, Pos;
	double PLh, Lh;
	char *PriorStr;
	SITE_POWER* SitePowers;
	PRIOR *Prior;

	SitePowers = Rates->SitePowers;

	PriorStr = (char*)SMalloc(sizeof(char) * 128);
	Lh = 0;
	for(Index=0;Index<SitePowers->NoSites;Index++)
	{
		Pos = SitePowers->SiteIndex[Index];
		sprintf(PriorStr, "Power-%d", Pos);
		Prior = GetPriorFromName(PriorStr, Rates->Priors, Rates->NoPriors);
		PLh = CalcLhPriorP(SitePowers->Powers[Index], Prior);

		if(PLh == ERRLH)
		{
			free(PriorStr);
			return ERRLH;
		}
		Lh += PLh;
	}

	free(PriorStr);
	return Lh;



	
	

}