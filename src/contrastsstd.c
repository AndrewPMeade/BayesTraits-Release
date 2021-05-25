#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

#include "typedef.h"
#include "genlib.h"

double	CaclStdContrastLhMLSite(OPTIONS *Opt, TREES *Trees, RATES *Rates, int SiteNo)
{
	int			Index, CIndex;
	TREE		*T;
	NODE		N;
	CONTRAST	*Con;
	CONTRASTR	*ConRates;
	double		GlobalVar;
	double		SumLogVar;
	double		T1, Ret;
	int			NoCon;

	ConRates = Rates->Contrast;
	T = Trees->Tree[Rates->TreeNo];

	NoCon = 0;
	GlobalVar = 0;
	SumLogVar = 0;

	for(Index=0;Index<T->NoNodes;Index++)
	{
		N = T->NodeList[Index];
		if(N->Tip == FALSE)
		{
			for(CIndex=0;CIndex<N->ConData->NoContrast;CIndex++)
			{
				Con = N->ConData->Contrast[CIndex];

				GlobalVar += (Con->Cont[SiteNo] * Con->Cont[SiteNo]) / Con->Var[SiteNo];

				SumLogVar += log(Con->Var[SiteNo]);
				NoCon++;
			}
		}
	}
	NoCon = T->NoContrast;
//	GlobalVar = T->Root->ConData->GVar[SiteNo];
//	SumLogVar = T->Root->ConData->SumLogVar[SiteNo];

//	exit(0);
	SumLogVar += log(T->Root->ConData->Contrast[0]->Err[SiteNo]);

	T1 = GlobalVar;
	GlobalVar = GlobalVar / (NoCon+1);
	Ret = (NoCon+1) * log(6.28318530717958647692528676655900576839 * GlobalVar);
	Ret += SumLogVar + (T1 / GlobalVar);
	Ret *= -0.5;

	Rates->Contrast->Alpha[SiteNo] = T->Root->ConData->Contrast[0]->Data[SiteNo];
	Rates->Contrast->Sigma[SiteNo] = GlobalVar;

	return Ret;
}

double	CaclStdContrastLhML(OPTIONS *Opt, TREES *Trees, RATES *Rates)
{
	double Ret;
	int Index;

	Ret = 0;

	for(Index=0;Index<Trees->NoOfSites;Index++)
		Ret += CaclStdContrastLhMLSite(Opt, Trees, Rates, Index);

	return Ret;
}

double	CaclAlphaErr(NODE N, double EstAlpha, int SiteNo)
{
	double		Ret;
	CONTRAST	*Con;

	Con = N->ConData->Contrast[0];

	Ret = (EstAlpha - Con->Data[SiteNo]) * (EstAlpha - Con->Data[SiteNo]);
	Ret = Ret / Con->Err[SiteNo];

	return Ret;
}

double CaclStdContrastLhMCMCSite(OPTIONS *Opt, TREES* Trees, RATES* Rates, int SiteNo)
{
	int			Index, CIndex;
	TREE		*T;
	NODE		N;
	CONTRAST	*Con;
	CONTRASTR	*ConRates;
	double		GlobalVar;
	double		SumLogVar;
	double		T1, T2;
	int			NoCon;
	double		Ret;

	ConRates = Rates->Contrast;
	T = Trees->Tree[Rates->TreeNo];

	NoCon = 0;
	GlobalVar = 0;
	SumLogVar = 0;

#ifdef RES_SIGMA
	ConRates->Sigma[SiteNo] = RES_SIGMA;
#endif

#ifdef RES_ALPHA
	ConRates->Alpha[SiteNo] = RES_ALPHA;
#endif

	for(Index=0;Index<T->NoNodes;Index++)
	{
		N = T->NodeList[Index];
		if(N->Tip == FALSE)
		{
			for(CIndex=0;CIndex<N->ConData->NoContrast;CIndex++)
			{
				Con = N->ConData->Contrast[CIndex];
				GlobalVar += (Con->Cont[SiteNo] * Con->Cont[SiteNo]) / Con->Var[SiteNo];

				SumLogVar += log(Con->Var[SiteNo]);
				NoCon++;
			}
		}
	}
	
	NoCon = T->NoContrast;
//	GlobalVar = T->Root->ConData->GVar[SiteNo];
//	SumLogVar = T->Root->ConData->SumLogVar[SiteNo];

	SumLogVar += log(T->Root->ConData->Contrast[0]->Err[SiteNo]);

	T1 = GlobalVar;
	GlobalVar = GlobalVar / NoCon;

	Ret = (NoCon+1) * log(6.28318530717958647692528676655900576839 * ConRates->Sigma[SiteNo]);

	T2 = ((NoCon+1) * GlobalVar) +  CaclAlphaErr(T->Root, ConRates->Alpha[SiteNo], SiteNo);
	Ret += SumLogVar + (T2 / ConRates->Sigma[SiteNo]);
	Ret *= -0.5;

	return Ret;
}

double	CaclStdContrastLhMCMC(OPTIONS *Opt, TREES *Trees, RATES *Rates)
{
	double Ret;
	int Index;

	Ret = 0;

	for(Index=0;Index<Trees->NoOfSites;Index++)
		Ret += CaclStdContrastLhMCMCSite(Opt, Trees, Rates, Index);

	return Ret;
}

double	CaclStdContrastLh(OPTIONS *Opt, TREES *Trees, RATES *Rates)
{
	if(Opt->Analsis == ANALML)
		return CaclStdContrastLhML(Opt, Trees, Rates);
	
	if(Opt->Analsis == ANALMCMC)
		return CaclStdContrastLhMCMC(Opt, Trees, Rates);

	return ERRLH;
}
