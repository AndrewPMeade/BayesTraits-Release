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
#include <stdlib.h>
#include <math.h>
#include <string.h>

#include "TypeDef.h"
#include "GenLib.h"
#include "Likelihood.h"
void		PrintContrast(RATES *Rates, TREES *Trees);


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

//	PrintContrast(Rates, Trees);exit(0);

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

				GlobalVar += (Con->Cont[SiteNo] * Con->Cont[SiteNo]) / Con->Var;

				SumLogVar += log(Con->Var);
				NoCon++;
			}
		}
	}
	NoCon = T->NoContrast;
//	GlobalVar = T->Root->ConData->GVar[SiteNo];
//	SumLogVar = T->Root->ConData->SumLogVar[SiteNo];

//	exit(0);
	SumLogVar += log(T->Root->ConData->Contrast[0]->Err);

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

	for(Index=0;Index<Trees->NoSites;Index++)
		Ret += CaclStdContrastLhMLSite(Opt, Trees, Rates, Index);

	return Ret;
}

double	CaclAlphaErr(NODE N, double EstAlpha, int SiteNo)
{
	double		Ret;
	CONTRAST	*Con;

	Con = N->ConData->Contrast[0];

	Ret = (EstAlpha - Con->Data[SiteNo]) * (EstAlpha - Con->Data[SiteNo]);
	Ret = Ret / Con->Err;

	return Ret;
}


double CaclStdContrastLhMCMC(OPTIONS *Opt, TREES* Trees, RATES* Rates)
{
	int			Index, CIndex, SIndex;
	TREE		*T;
	NODE		N;
	CONTRAST	*Con;
	CONTRASTR	*ConRates;
	double		*GlobalVar;
	double		SumLogVar;
	double		T1, T2;
	int			NoCon, NoSites;
	double		Ret, GRet;


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

	NoSites = Trees->NoSites;

	GlobalVar = (double*)SMalloc(sizeof(double) * NoSites);
	for(SIndex=0;SIndex<NoSites;SIndex++)
		GlobalVar[SIndex] = 0;

	for(Index=0;Index<T->NoNodes;Index++)
	{
		N = T->NodeList[Index];
		if(N->Tip == FALSE)
		{
			for(CIndex=0;CIndex<N->ConData->NoContrast;CIndex++)
			{
				Con = N->ConData->Contrast[CIndex];
				SumLogVar += log(Con->Var);

				for(SIndex=0;SIndex<NoSites;SIndex++)
					GlobalVar[SIndex] += (Con->Cont[SIndex] * Con->Cont[SIndex]) / Con->Var;
				
				NoCon++;
			}
		}
	}
	
	NoCon = T->NoContrast;

	SumLogVar += log(T->Root->ConData->Contrast[0]->Err);

	GRet = 0;
	for(SIndex=0;SIndex<NoSites;SIndex++)
	{

		T1 = GlobalVar[SIndex];
		GlobalVar[SIndex] = GlobalVar[SIndex] / NoCon;

		Ret = (NoCon+1) * log(6.28318530717958647692528676655900576839 * ConRates->Sigma[SIndex]);

		T2 = ((NoCon+1) * GlobalVar[SIndex]) +  CaclAlphaErr(T->Root, ConRates->Alpha[SIndex], SIndex);
		Ret += SumLogVar + (T2 / ConRates->Sigma[SIndex]);
		Ret *= -0.5;

		GRet += Ret;


	}

	free(GlobalVar);
	return GRet;
}


double	CaclStdContrastLh(OPTIONS *Opt, TREES *Trees, RATES *Rates)
{
	if(Opt->Analsis == ANALML)
		return CaclStdContrastLhML(Opt, Trees, Rates);
	
	if(Opt->Analsis == ANALMCMC)
		return CaclStdContrastLhMCMC(Opt, Trees, Rates);

	return ERRLH;
}
