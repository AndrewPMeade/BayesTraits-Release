#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "typedef.h"
#include "Rates.h"
#include "genlib.h"
#include "RandLib.h"
#include "trees.h"
#include "continuous.h"
#include "revjump.h"
#include "priors.h"
#include "likelihood.h"
#include "data.h"
#include "matrix.h"
#include "randdists.h"
#include "contrasts.h"
#include "VarRates.h"
#include "BigLh.h"
#include "ml.h"
#include "schedule.h"
#include "modelfile.h"
#include "SimData.h"
#include "RJDummy.h"
#include "ContrastRegOutput.h"
#include "FatTail.h"
#include "Geo.h"
#include "Threaded.h"
#include "LocalTransform.h"
#include "DistData.h"

//double**	LoadModelFile(RATES* Rates, OPTIONS *Opt);
//void		SetFixedModel(RATES *Rates, OPTIONS *Opt);


void	SetRegBetaZero(int NoSites, RATES *Rates)
{
	int Index;

	for(Index=0;Index<NoSites;Index++)
		Rates->Beta[Index] = 0;
}

int		FindNoConRates(OPTIONS *Opt)
{
	switch(Opt->Model)
	{
		case M_CONTINUOUS_RR:
			return Opt->Trees->NoUserSites;
		break;
			
		case M_CONTINUOUS_DIR:
			return Opt->Trees->NoUserSites  * 2;
		break;

		case M_CONTINUOUS_REG:
			return Opt->Trees->NoUserSites; 
		break;

		case M_CONTRAST_CORREL:
			return Opt->Trees->NoUserSites;
		break;

		case M_CONTRAST_REG:
			return Opt->Trees->NoUserSites - 1; 
		break;

		case M_CONTRAST:
			return Opt->Trees->NoUserSites * 2;
		break;

		case M_FATTAIL:
			return Opt->Trees->NoUserSites * 2;
		break;

		case M_GEO:
			return 2;
		break;
	}

	printf("Unkonwn model %s::%d\n", __FILE__, __LINE__);
	exit(0);
	return 0;
}

int		FindNoOfRates(OPTIONS *Opt)
{
	int	Index;
	int	Ret;

	Ret = 0;

	if(Opt->UseRModel == TRUE)
	{
		if(Opt->RModelP == -1)
			Ret++;
	}
	else
	{
		for(Index=0;Index<Opt->NoOfRates;Index++)
			if(Opt->ResTypes[Index] == RESNONE)
				Ret++;	
	}

	return Ret;
}

char**	GetRateNames(OPTIONS *Opt)
{
	char **Ret;
	int No, Index, Pos;

	No = FindNoOfRates(Opt);

	Ret = (char**)SMalloc(sizeof(char*) * No);

	Pos = 0;
	for(Index=0;Index<Opt->NoOfRates;Index++)
	{
		if(Opt->ResTypes[Index] == RESNONE)
			Ret[Pos++] = StrMake(Opt->RateName[Index]);
	}

	return Ret;
}


double	FindRateVal(int Pos, RATES *Rates, OPTIONS *Opt)
{
	int	RateIndex=0;
	int	OptIndex=0;

	if(Opt->ResTypes[Pos] == RESCONST)
		return Opt->ResConst[Pos];

	for(;;)
	{
		if(OptIndex==Pos)
			return Rates->Rates[RateIndex];

		
		if(Opt->ResTypes[OptIndex] == RESNONE)
			RateIndex++;
		OptIndex++;
	}
	return -1;
}

void	MapMCMCConRates(RATES* Rates, OPTIONS *Opt)
{
	int	Index, NoSites;

	NoSites = Opt->Trees->NoOfSites;

	if(Opt->ModelType == MT_FATTAIL)
	{
		MapRatesToFatTailRate(Rates, Rates->FatTailRates);
		return;
	}

	if(Opt->ModelType == MT_CONTRAST)
	{
		MapRatesToConVals(Opt, Rates, Rates->Contrast);
		return;
	}

	if(Opt->Model == M_CONTINUOUS_REG)
	{
		
		Rates->Means[0] = Rates->Rates[0];

//		memcpy(Rates->Beta, Rates->Rates, sizeof(double) * Rates->NoOfRates);
		for(Index=1;Index<Rates->NoOfRates;Index++)
			Rates->Beta[Index - 1] = Rates->Rates[Index];

		return;
	}
	
	for(Index=0;Index<Opt->Trees->NoOfSites;Index++)
		Rates->Means[Index] = Rates->Rates[Index];
	

	if(Opt->Model == M_CONTINUOUS_RR)
		return;

	if(Opt->Model == M_CONTINUOUS_DIR)
	{
		for(;Index<Rates->NoOfRates;Index++)
			Rates->Beta[Index - Opt->Trees->NoOfSites] = Rates->Rates[Index];
		return;
	}


	
}

int		FindRatePos(int Rate, OPTIONS *Opt)
{
	int	Pos;

	if((Opt->ResTypes[Rate] == RESNONE) ||(Opt->ResTypes[Rate] == RESCONST))
		return Rate;

	Pos = Rate;
	do
	{
		Pos = Opt->ResNo[Pos];
		if((Opt->ResTypes[Pos] == RESNONE) || (Opt->ResTypes[Pos] == RESCONST))
			return Pos;

	}while(1==1);
	return Pos;
}


void	MapRates(RATES* Rates, OPTIONS *Opt)
{
	int	Index;
	int	Pos;

	if(Opt->DataType == CONTINUOUS)
	{
		if(Opt->Analsis == ANALMCMC)
			MapMCMCConRates(Rates, Opt);
	
		return;
	}

	if(Opt->UseRJMCMC == TRUE)
	{
//		MapRJRates(Rates->Rates, Rates->MappingVect, Rates->NoOfFullRates, Rates->FullRates);
		MapRJRates(Opt, Rates);
		return;
	}

	for(Index=0;Index<Rates->NoOfRates;Index++)
	{
		if((Rates->Rates[Index] < MINRATE) || (IsNum(Rates->Rates[Index]) == FALSE))
			Rates->Rates[Index] = MINRATE;

		if(Rates->Rates[Index] > MAXRATE)
			Rates->Rates[Index] = MAXRATE;
	}

	if(Opt->UseRModel == TRUE)
	{
		if(Opt->RModelP != -1)
			Rates->FullRates[0] = Opt->RModelP;
		else
			Rates->FullRates[0] = Rates->Rates[0];
	}
	else
	{
		for(Index=0;Index<Rates->NoOfFullRates;Index++)
		{
			Pos = FindRatePos(Index, Opt);
			Rates->FullRates[Index] = FindRateVal(Pos, Rates, Opt);

			if(Rates->FullRates[Index] < MINRATE)
				Rates->FullRates[Index] = MINRATE;

			if(Rates->FullRates[Index] > MAXRATE)
				Rates->FullRates[Index] = MAXRATE;
		}
	}
	

	Pos = Rates->NoOfRates;
	if((Opt->UseCovarion == TRUE) && (Opt->Analsis == ANALML))
		Pos = Pos - 1;

	if((Opt->EstKappa == TRUE) && (Opt->Analsis == ANALML))
		Pos = Pos - 1;

	if((Opt->EstGamma == TRUE) && (Opt->Analsis == ANALML))
		Pos = Pos - 1;

	if((Opt->UseCovarion == TRUE) && (Opt->Analsis == ANALML))
	{
		Rates->OnToOff = Rates->Rates[Pos++];
	//	Rates->OffToOn = Rates->Rates[Pos++];
		Rates->OffToOn = Rates->OnToOff; 
	}

	if((Opt->EstKappa == TRUE) && (Opt->Analsis == ANALML))
	{
		Rates->Kappa = Rates->Rates[Pos++];
	
		if(Rates->Kappa < 0)
			Rates->Kappa = 0.0000001;

		if(Rates->Kappa > 5)
			Rates->Kappa = 5;
	}

	if((Opt->EstGamma == TRUE) && (Opt->Analsis == ANALML))
	{
		Rates->Gamma = Rates->Rates[Pos++];
	
		if(Rates->Gamma < MIN_GAMMA)
			Rates->Gamma = MIN_GAMMA;

		if(Rates->Gamma > MAX_GAMMA)
			Rates->Gamma = MAX_GAMMA;
	}
}

void	FindEmpPis(RATES *Rates, OPTIONS *Opt)
{
	TREES	*Trees;
	double	*TempPis;
	int		State;
	double	Weight;
	double	Total;
	int		SIndex,TIndex,SymbolIndex;
	TAXA	*Taxa;

	Trees = Opt->Trees;

	TempPis = (double*)malloc(sizeof(double)*Trees->NoOfStates);
	if(TempPis == NULL)
		MallocErr();
	for(SIndex=0;SIndex<Trees->NoOfStates;SIndex++)
		TempPis[SIndex] = 0;

	for(SIndex=0;SIndex<Trees->NoOfSites;SIndex++)
	{
		for(TIndex=0;TIndex<Trees->NoOfTaxa;TIndex++)
		{
			Taxa = Trees->Taxa[TIndex];

			if(SiteHadUnKnownState(Taxa->DesDataChar[SIndex]) == FALSE)
			{
				Weight = (double)1/(double)strlen(Taxa->DesDataChar[SIndex]);

				for(SymbolIndex=0;SymbolIndex<(int)strlen(Taxa->DesDataChar[SIndex]);SymbolIndex++)
				{
					State = SymbolToPos(Taxa->DesDataChar[SIndex][SymbolIndex], Trees->SymbolList);
					TempPis[State] += Weight;
				}
			}
			else
			{
				Weight = (double)1/(double)Trees->NoOfStates;
				for(SymbolIndex=0;SymbolIndex<Trees->NoOfStates;SymbolIndex++)
					TempPis[SymbolIndex] += Weight;
			}
		}
	}

	Total = (double)Trees->NoOfSites * (double)Trees->NoOfTaxa;

	for(SIndex=0;SIndex<Trees->NoOfStates;SIndex++)
		Rates->Pis[SIndex] = TempPis[SIndex] / Total;

	free(TempPis);
}

void	SetPiValues(RATES *Rates, OPTIONS *Opt)
{
	int		Index;
	TREES	*Trees;

	Trees = Opt->Trees;

	if(Opt->PiTypes == PIUNI)
	{
		for(Index=0;Index<Trees->NoOfStates;Index++)
			Rates->Pis[Index] = (double)1/Trees->NoOfStates;

		return;
	}

	if(Opt->PiTypes == PINONE)
	{
		for(Index=0;Index<Trees->NoOfStates;Index++)
			Rates->Pis[Index] = 1;

		return;
	}

	FindEmpPis(Rates, Opt);
}

double	GetHMean(OPTIONS *Opt, RATES *Rates)
{
#ifndef BIG_LH
	return log(Rates->HMeanCount / Rates->HMeanSum);
#else
	mpfr_t	t1, t2;
	double Ret;

	mpfr_init2(t1, Opt->Precision);
	mpfr_init2(t2, Opt->Precision);
	
	mpfr_si_div(t1, Rates->HMeanCount, Rates->HMeanSum, DEF_ROUND);
	mpfr_log(t2, t1, DEF_ROUND);

	Ret  = mpfr_get_d(t2, DEF_ROUND);
	
	mpfr_clears(t1, t2, NULL);

	return Ret;
#endif
}

int		FindNoEstDataPoint(OPTIONS *Opt, TREES *Trees)
{
	int Index, SIndex, Ret;
	TAXA *Taxa;

	Ret = 0;

	for(Index=0;Index<Trees->NoOfTaxa;Index++)
	{
		Taxa = Trees->Taxa[Index];
		for(SIndex=0;SIndex<Trees->NoOfSites;SIndex++)
			if(Taxa->EstDataP[SIndex] == TRUE)
				Ret++;
		
		if(Taxa->EstDepData == TRUE)
			Ret++;
	}

	return Ret;
}


void	CreatCRates(OPTIONS *Opt, RATES *Rates)
{
	int		Index;
	TREES	*Trees;

	Rates->NoOfRates = 0;

	Rates->Means = NULL;
	Rates->Rates = NULL;
	Rates->Beta = NULL;

	Rates->Delta = 1;
	Rates->Kappa = 1;
	Rates->Lambda= 1;

	if(Opt->EstOU == TRUE)
		Rates->OU = MIN_OU;
	else
		Rates->OU = 0;

	if(Opt->Analsis == ANALMCMC)
	{
		Rates->NoOfRates = FindNoConRates(Opt);
	
		Rates->NoOfFullRates = Rates->NoOfRates;

		Rates->Rates = (double*)malloc(sizeof(double) * Rates->NoOfRates);
		if(Rates->Rates == NULL)
			MallocErr();

		for(Index=0;Index<Rates->NoOfRates;Index++)
			Rates->Rates[Index] = 0;
		
		if(Opt->ModelType == MT_CONTRAST)
		{
			Rates->Means = NULL;
			Rates->Beta	 = NULL;
		}
		else
		{
			if(Opt->Model == M_CONTINUOUS_REG)
				Rates->Means = (double*)malloc(sizeof(double));
			else
				Rates->Means = (double*)malloc(sizeof(double) * Opt->Trees->NoOfSites);
			
			if(Rates->Means == NULL)
				MallocErr();

			if((Opt->Model == M_CONTINUOUS_DIR) || (Opt->Model == M_CONTINUOUS_REG))
			{
				Rates->Beta = (double*)malloc(sizeof(double) * Opt->Trees->NoOfSites);
				if(Rates->Beta== NULL)
					MallocErr();
			}
			else
				Rates->Beta = NULL;
		}
	}

	Trees = Opt->Trees;

	Rates->UseEstData	=	FALSE;
	Rates->EstData		=	NULL;

	Rates->NoEstData	=	FindNoEstDataPoint(Opt, Trees);

	if(Rates->NoEstData > 0)
	{
		Rates->UseEstData = TRUE;
		Rates->EstData = (double*)malloc(sizeof(double) * Rates->NoEstData);
		if(Rates->EstData == NULL)
			MallocErr();
		for(Index=0;Index<Rates->NoEstData;Index++)
			Rates->EstData[Index] = 0;
	}
	
	if(Opt->ModelType == MT_CONTRAST)
		Rates->Contrast = CreatContrastRates(Opt, Rates);

	if(Opt->LoadModels == TRUE)
	{
		Rates->ModelFile = LoadModelFile(Opt->LoadModelsFN, Opt, Opt->Trees, Rates);
		ChangeModelFile(Rates, Rates->RS);
	}

	if(Opt->ModelType == MT_FATTAIL)
	{
		Rates->FatTailRates = CreateFatTailRates(Opt, Trees);
		MapFatTailRateToRates(Rates, Rates->FatTailRates);
	}
}

int		FindNoEstData(TREES *Trees, OPTIONS *Opt)
{
	int		TIndex;
	int		SIndex;
	TAXA	*Taxa;
	int		Ret;
	int		NOS;
	
	Ret = 0;

	for(TIndex=0;TIndex<Trees->NoOfTaxa;TIndex++)
	{
		Taxa = Trees->Taxa[TIndex];
		
		if(Opt->Model == M_MULTISTATE)
			NOS = Trees->NoOfSites;
		else
			NOS = 2;

		for(SIndex=0;SIndex<NOS;SIndex++)
		{
			if(Taxa->EstDataP[SIndex] == TRUE)
				Ret++;
		}		
	}

	return Ret;
}

void	InitHMean(RATES* Rates, OPTIONS *Opt)
{
	Rates->HMeanCount = 0;

#ifndef BIG_LH
	Rates->HMeanSum		= 0.0;
#else
	mpfr_init2(Rates->HMeanSum, Opt->Precision);
	mpfr_set_d(Rates->HMeanSum, 0.0, DEF_ROUND);
#endif
}

INVINFO**	CreatInvInfo(int NOS,  int NoM)
{
	int	 Index;
	INVINFO** Ret;

	Ret = (INVINFO**)malloc(sizeof(INVINFO**) * NoM);
	if(Ret == NULL)
		MallocErr();

	for(Index=0;Index<NoM;Index++)
		Ret[Index] = AllocInvInfo(NOS);

	return Ret;
}

void	FreeHetero(HETERO* Hetero)
{
	int Index;

	for(Index=0;Index<Hetero->NoModels;Index++)
		FreeInvInfo(Hetero->ModelInv[Index]);

	free(Hetero->ModelInv);
	free(Hetero->MList);
	free(Hetero);
}

HETERO*	CreatHetero(OPTIONS *Opt, RATES* Rates)
{
	HETERO *Ret;
	TREES *Trees;
	int		Index;

	Trees = Opt->Trees;

	Ret = (HETERO*)malloc(sizeof(HETERO));
	if(Ret == NULL)
		MallocErr();

	Ret->NoModels = 2;
	Ret->ModelInv = CreatInvInfo(Trees->NoOfStates, Ret->NoModels);
	
	Ret->MListSize = Trees->MaxNodes;

	Ret->MList = (int*)malloc(sizeof(int) * Ret->MListSize);
	if(Ret->MList == NULL)
		MallocErr();
	
	for(Index=0;Index<Ret->MListSize;Index++)
		Ret->MList[Index] = RandUSInt(Rates->RS) % Ret->NoModels;

	return Ret;
}

void	CopyHetero(HETERO *A, HETERO *B)
{
	memcpy(A->MList, B->MList, sizeof(int) * A->MListSize);
}

void	 MutateHetero(RATES *Rates)
{
	int No, New;
	HETERO *Hetero;

	Hetero = Rates->Hetero;
	No = RandUSInt(Rates->RS) % Hetero->MListSize;

	do
	{
		New = RandUSInt(Rates->RS) % Hetero->NoModels;
	}while(New == Hetero->MList[No]);

	Hetero->MList[No] = New;
}

void	SetRatesLocalRates(RATES *Rates, OPTIONS *Opt)
{	
	int Index;

	Rates->UseLocalTransforms	= FALSE;
	Rates->EstLocalTransforms	= FALSE;
	Rates->NoLocalTransforms		= Opt->NoLocalTransforms;
	Rates->LocalTransforms	= NULL;

	if(Opt->NoLocalTransforms == 0)
		return;

	Rates->LocalTransforms = (LOCAL_TRANSFORM**)malloc(sizeof(LOCAL_TRANSFORM*) * Opt->NoLocalTransforms);
	if(Rates->LocalTransforms == NULL)
		MallocErr();

	for(Index=0;Index<Rates->NoLocalTransforms;Index++)
		Rates->LocalTransforms[Index] = CloneLocalTransform(Opt->LocalTransforms[Index]);
	
	Rates->EstLocalTransforms = EstLocalTransforms(Opt->LocalTransforms, Opt->NoLocalTransforms);
	Rates->UseLocalTransforms = TRUE;

	for(Index=0;Index<Rates->NoLocalTransforms;Index++)
		if(Rates->LocalTransforms[Index]->Est == TRUE)
		{
			if(Rates->LocalTransforms[Index]->Type == VR_OU)
				Rates->LocalTransforms[Index]->Scale = MIN_OU;
			else
				Rates->LocalTransforms[Index]->Scale = 1.0;
		}
}

void	SetDiscretisedPrior(PRIOR *Prior, OPTIONS *Opt)
{
	if(Prior == NULL)
		return;

	Prior->Discretised = TRUE;
	Prior->Width = 1.0 / Opt->PriorCats;
}

void	SetDiscretisedPriors(RATES *Rates, OPTIONS *Opt)
{
	PRIOR *Prior;
	int Index;
	
	Prior = GetPriorFromName("RJRates", Rates->Priors, Rates->NoPriors);
	SetDiscretisedPrior(Prior, Opt);

	for(Index=0;Index<Rates->NoOfRates;Index++)
	{
		Prior = GetPriorFromName(Rates->RateNames[Index], Rates->Priors, Rates->NoPriors);
		SetDiscretisedPrior(Prior, Opt);
	}
}

RATES*	CreatRates(OPTIONS *Opt)
{
	RATES*	Ret;
	int		Index;

	Ret = (RATES*)SMalloc(sizeof(RATES));

	Ret->NoOfFullRates	= Opt->NoOfRates;

	if(Opt->UseRModel == TRUE)
		Ret->NoOfFullRates = 1;

	Ret->NoOfRates		= FindNoOfRates(Opt);
	Ret->RateNames		= GetRateNames(Opt);
	Ret->NoOfRJRates	= -1;
	Ret->TreeNo			= 0;
	Ret->Rates			= NULL;

	Ret->Pis			= NULL;
	Ret->FullRates		= NULL;
	Ret->Means			= NULL;
	Ret->Beta			= NULL;
	Ret->MappingVect	= NULL;
	Ret->LhPrior		= 0;
	Ret->LnHastings		= 0;
	Ret->LnJacobion		= 0;

	Ret->Gamma			= -1;
	Ret->GammaCats		= 1;
	Ret->GammaMults		= NULL;

	Ret->Lh				= 0;

	InitHMean(Ret, Opt);

	Ret->NoEstData		=	0;
	Ret->EstData		=	NULL;
	Ret->ModelFile		=	NULL;
	Ret->ModelNo		=	-1;

	Ret->EstData		=	NULL;
	Ret->EstDescData	=	NULL;
	Ret->UseEstData		=	FALSE;
	Ret->NoEstData		=	0;

	Ret->Kappa			=	-1;
	Ret->Lambda			=	-1;
	Ret->Delta			=	-1;
	Ret->OU				=	-1;

	Ret->VarRates		=	NULL;
	Ret->Hetero			=	NULL;
	Ret->ModelFile		=	NULL;

	Ret->Contrast		=	NULL;
	Ret->RJDummy		=	NULL;
	Ret->FatTailRates	=	NULL;
	
	Ret->DistDataRates	=	NULL;

	Ret->NoPriors		=	0;
	Ret->Priors			=	NULL;
	Ret->AutoAccept		=	FALSE;
	Ret->CalcLh			=	TRUE;
	
	Ret->RS				=	CreateSeededRandStates(Opt->Seed);
	Ret->RSList			=	CreateRandStatesList(Ret->RS, GetMaxThreads());
	Ret->RNG			=	gsl_rng_alloc(gsl_rng_mt19937);
	gsl_rng_set(Ret->RNG, Opt->Seed);

	SetRatesLocalRates(Ret, Opt);

	if(Opt->UseDistData == TRUE)
		Ret->DistDataRates = CreateDistDataRates(Opt->DistData, Ret->RS);

	if(Opt->Analsis == ANALMCMC)
	{
		CrateRatePriors(Opt, Ret);
		if(Opt->ModelType == MT_DISCRETE)
			SetDiscretisedPriors(Ret, Opt);
	}

	
	
	if(Opt->UseGamma == TRUE)
	{
		Ret->GammaMults= (double*)SMalloc(sizeof(double) * Opt->GammaCats);

		Ret->GammaCats = Opt->GammaCats;
		
		if(Opt->EstGamma == FALSE)
			Ret->Gamma = Opt->FixGamma;
		else
			Ret->Gamma = 1;
	}

	if(Opt->UseKappa == FALSE)
		Ret->Kappa = -1;
	else
	{
		if(Opt->EstKappa == FALSE)
			Ret->Kappa = Opt->FixKappa;
		else
			Ret->Kappa = 1;
	}

	if(UseNonParametricMethods(Opt) == TRUE)
		Ret->VarRates = CreatVarRates(Ret, Opt->Trees, Opt);
	
	if(Opt->DataType == CONTINUOUS)
	{
		CreatCRates(Opt, Ret);
		return Ret;
	}

	if(Ret->NoOfRates > 0)
		Ret->Rates = (double*)SMalloc(sizeof(double)*Ret->NoOfRates);

	Ret->FullRates = (double*)SMalloc(sizeof(double)*Ret->NoOfFullRates);
	
	for(Index=0;Index<Ret->NoOfRates;Index++)
		Ret->Rates[Index] = 1;

	Ret->Pis = (double*)SMalloc(sizeof(double)*Opt->Trees->NoOfStates);

	SetPiValues(Ret, Opt);

	if(Opt->UseCovarion == TRUE)
	{
		Ret->OffToOn = 1;
		Ret->OnToOff = 1;
		Ret->CoVarPis[0] = 1;
		Ret->CoVarPis[1] = 1;
	}

	if(Opt->UseRJMCMC == TRUE)
	{	
		Ret->NoOfRJRates	= Ret->NoOfRates;

		Ret->MappingVect = (int*)SMalloc(sizeof(int) * Ret->NoOfRates);

		/* Inishal all rates to be in unique rate classes */
		if(Opt->CapRJRatesNo == -1)
		{
			for(Index=0;Index<Ret->NoOfRates;Index++)
				Ret->MappingVect[Index] = Index;
			Ret->NoOfRJRates = Ret->NoOfRates;
		}
		else
		{
		/* Inishal all rates to be in the same classes */
			for(Index=0;Index<Ret->NoOfRates;Index++)
				Ret->MappingVect[Index] = 0;
			Ret->NoOfRJRates = 1;
		}
	}
	else
		Ret->MappingVect = NULL;

	Ret->NoEstData	= FindNoEstData(Opt->Trees, Opt);
	if(Ret->NoEstData > 0)
	{
		Ret->UseEstData = TRUE;
		Ret->EstDescData = (int*)SMalloc(sizeof(int) * Ret->NoEstData);

		for(Index=0;Index<Ret->NoEstData;Index++)
		{
			if(Opt->Model == M_MULTISTATE)
				Ret->EstDescData[Index] = RandUSLong(Ret->RS) % Opt->Trees->NoOfStates;
			else
				Ret->EstDescData[Index] = RandUSLong(Ret->RS) % 2;
		}
	}

	if(Opt->Model == M_DESCHET)
		Ret->Hetero = CreatHetero(Opt, Ret);

	MapRates(Ret, Opt);

	if(Opt->LoadModels == TRUE)
	{
		Ret->ModelFile = LoadModelFile(Opt->LoadModelsFN, Opt, Opt->Trees, Ret);
		ChangeModelFile(Ret, Ret->RS);
	}
	
	SimData(Opt, Opt->Trees, Ret);

	return Ret;
}

void	PrintConRegVarCoVarHeadder(FILE* Str, int NoOfSites)
{
	int	x;

	NoOfSites++;

	fprintf(Str, "s.e. Alpha\t");

	for(x=1;x<NoOfSites;x++)
		fprintf(Str, "s.e. Beta-%d\t", x);
	
/*
	for(x=0;x<NoOfSites;x++)
		for(y=x+1;y<NoOfSites;y++)
			if((y != DepSiteNo) && (x != DepSiteNo))
				fprintf(Str, "Trait %d,%d CoVar\t", x+1,y+1);
*/	
}

void	PrintEstDataHeader(FILE *Str, OPTIONS *Opt)
{
	TREES	*Trees;
	int		Index;
	int		NOS;
	TAXA	*Taxa;
	int		x;
	
	Trees = Opt->Trees;

	NOS = Trees->NoOfSites;

	if((Opt->Model == M_DESCINDEP) || (Opt->Model == M_DESCDEP))
		NOS = 2;

	for(Index=0;Index<Trees->NoOfTaxa;Index++)
	{
		Taxa = Trees->Taxa[Index];
		if(Taxa->EstData == TRUE)
		{
			for(x=0;x<NOS;x++)
			{
				if(Taxa->EstDataP[x] == TRUE)
					fprintf(Str, "Est %s - %d\t", Taxa->Name, x+1);
				
			}
			if(Taxa->EstDepData == TRUE)
				fprintf(Str, "Est %s - Dep\t", Taxa->Name);
		}
	}
}

void	PrintConRecNodesHeadder(FILE *Str, OPTIONS *Opt)
{
	int		Index, SiteIndex;
	RECNODE	*RNode;
		
	if(Opt->ModelType != MT_CONTRAST)
		return;
	
	for(Index=0;Index<Opt->NoOfRecNodes;Index++)
	{
		RNode = Opt->RecNodeList[Index];

		for(SiteIndex=0;SiteIndex<Opt->Trees->NoOfSites;SiteIndex++)
		{
			if(Opt->Trees->NoOfSites == 1)
				fprintf(Str, "%s Alpha\t", RNode->Name);
			else
				fprintf(Str, "%s %d Alpha\t", RNode->Name, SiteIndex + 1);
		}
	}
}


void	FreeParamNames(int No, char **PName)
{
	int Index;

	for(Index=0;Index<No;Index++)
		free(PName[Index]);

	free(PName);
}

void	PrintLocalRateHeader(FILE *Str, OPTIONS *Opt)
{
	int Index;
	LOCAL_TRANSFORM *LRate;

	for(Index=0;Index<Opt->NoLocalTransforms;Index++)
	{
		LRate = Opt->LocalTransforms[Index];
		if(LRate->Est == TRUE)
			fprintf(Str, "%s - %s\t", LRate->Name, VarRatesTypeToStr(LRate->Type));
	}
}

void	PrintLocalTransformHeadder(FILE *Str, OPTIONS *Opt)
{
	if(Opt->UseVarRates == TRUE)
	{
		fprintf(Str, "No RJ Local Branch\t");
		fprintf(Str, "No RJ Local Node\t");
	}

	if(Opt->UseRJLocalScalar[VR_KAPPA] == TRUE)
		fprintf(Str, "No RJ Local Kappa\t");

	if(Opt->UseRJLocalScalar[VR_LAMBDA] == TRUE)
		fprintf(Str, "No RJ Local Lambda\t");

	if(Opt->UseRJLocalScalar[VR_DELTA] == TRUE)
		fprintf(Str, "No RJ Local Delta\t");

	if(Opt->UseRJLocalScalar[VR_OU] == TRUE)
		fprintf(Str, "No RJ Local OU\t");
}

void	PrintRatesHeadderCon(FILE *Str, OPTIONS *Opt)
{
	int		Index, NOS;
	int		x,y;

	NOS = Opt->Trees->NoOfSites;

	if(Opt->LoadModels == TRUE)
		fprintf(Str, "Model No\t");

	if(Opt->Model == M_CONTRAST_CORREL)
	{
		for(Index=0;Index<Opt->Trees->NoOfSites;Index++)
			fprintf(Str, "Alpha %d\t", Index+1);

		for(Index=0;Index<Opt->Trees->NoOfSites;Index++)
			fprintf(Str, "Sigma^2 %d\t", Index+1);

		for(x=0;x<NOS;x++)
		{
			for(y=0;y<x;y++)
				fprintf(Str, "R Trait %d %d\t", y+1, x+1);
		}
	}

	if(Opt->Model == M_CONTRAST)
	{
		for(Index=0;Index<Opt->Trees->NoOfSites;Index++)
			fprintf(Str, "Alpha %d\t", Index+1);

		for(Index=0;Index<Opt->Trees->NoOfSites;Index++)
			fprintf(Str, "Sigma^2 %d\t", Index+1);
	}

	if(Opt->Model == M_CONTRAST_REG)
	{

		if(Opt->RJDummy == TRUE)
			fprintf(Str, "No Dummy Codes\t");

		fprintf(Str, "Alpha\t");
		for(Index=1;Index<Opt->Trees->NoOfSites;Index++)
			fprintf(Str, "Beta %d\t", Index);

		fprintf(Str, "Var\t");
		fprintf(Str, "R^2\t");

		fprintf(Str, "s.e. Alpha\t");
		for(Index=1;Index<Opt->Trees->NoOfSites;Index++)
			fprintf(Str, "s.e. Beta-%d\t", Index);
	}

	if((Opt->Model == M_CONTINUOUS_DIR) || (Opt->Model == M_CONTINUOUS_RR))
	{
		for(Index=0;Index<Opt->Trees->NoOfSites;Index++)
		{
			fprintf(Str, "Alpha %d\t", Index+1);
		}
	}
	
	if(Opt->Model == M_CONTINUOUS_DIR)
	{
		for(Index=0;Index<Opt->Trees->NoOfSites;Index++)
			fprintf(Str, "Beta %d\t", Index+1);
	}

	if(Opt->Model == M_CONTINUOUS_REG)
	{
		fprintf(Str, "Alpha\t");

		for(Index=1;Index<Opt->NoOfRates;Index++)
		{
			fprintf(Str, "Beta %d\t", Index);
		}

		fprintf(Str, "Var\t");
		fprintf(Str, "R^2\tSSE\tSST\t");

		if(Opt->Analsis == ANALML)
			fprintf(Str, "Error Ratio\t");

		PrintConRegVarCoVarHeadder(Str, Opt->Trees->NoOfSites);
	}
	
	if((Opt->Model == M_CONTINUOUS_DIR) || (Opt->Model == M_CONTINUOUS_RR))
	{
		for(Index=0;Index<Opt->Trees->NoOfSites;Index++)
			fprintf(Str, "Sigma^2 %d \t", Index+1);

		for(x=0;x<Opt->Trees->NoOfSites;x++)
			for(y=x+1;y<Opt->Trees->NoOfSites;y++)
				fprintf(Str, "R Trait %d %d\t", x+1, y+1);
	}		

	if(Opt->ModelType == MT_FATTAIL)
	{
		if(Opt->Model == M_GEO)
			fprintf(Str, "Alpha\tScale\t");
		else
		{
			for(Index=0;Index<Opt->Trees->NoOfSites;Index++)
				fprintf(Str, "Alpha %d\tScale %d\t", Index+1, Index+1);
		}
	}


	if((Opt->EstKappa == TRUE) || (Opt->FixKappa != -1))
		fprintf(Str, "Kappa\t");

	if((Opt->EstDelta == TRUE) || (Opt->FixDelta != -1))
		fprintf(Str, "Delta\t");

	if((Opt->EstLambda == TRUE) || (Opt->FixLambda != -1))
		fprintf(Str, "Lambda\t");

	if((Opt->EstOU == TRUE) || (Opt->FixOU != -1))
		fprintf(Str, "OU\t");

	PrintLocalRateHeader(Str, Opt);

	if(Opt->NodeBLData == TRUE)
	{
		fprintf(Str, "Slope Nodes\tSlope  Root to Tip\t");
		fprintf(Str, "Min Nodes\tMax Nodes\t");
	}

	PrintEstDataHeader(Str, Opt);

	PrintConRecNodesHeadder(Str, Opt);

	PrintLocalTransformHeadder(Str, Opt);


	if(Opt->UseDistData == TRUE)
		OutputDataDistHeadder(Str, Opt);

	if(Opt->Analsis == ANALML)
		fprintf(Str, "\n");
}

void	PrintRecNodeHeadder(FILE* Str, OPTIONS *Opt, char* Name, int SiteNo)
{
	int		Index;
	int		NOS;
	TREES	*Trees;

	if((Opt->Model == M_DESCDEP) || (Opt->Model == M_DESCINDEP))
	{
		fprintf(Str, "%s - P(0,0)\t", Name);
		fprintf(Str, "%s - P(0,1)\t", Name);
		fprintf(Str, "%s - P(1,0)\t", Name);
		fprintf(Str, "%s - P(1,1)\t", Name);
		return;
	}

	if(Opt->Model == M_DESCCV)
	{
		fprintf(Str, "%s - I P(0,0)\t", Name);
		fprintf(Str, "%s - I P(0,1)\t", Name);
		fprintf(Str, "%s - I P(1,0)\t", Name);
		fprintf(Str, "%s - I P(1,1)\t", Name);
		fprintf(Str, "%s - D P(0,0)\t", Name);
		fprintf(Str, "%s - D P(0,1)\t", Name);
		fprintf(Str, "%s - D P(1,0)\t", Name);
		fprintf(Str, "%s - D P(1,1)\t", Name);

		return;
	}

	if(Opt->UseCovarion == TRUE)
		NOS = (Opt->Trees->NoOfStates / 2);
	else
		NOS = Opt->Trees->NoOfStates;

	if(Opt->NOSPerSite == FALSE)
	{
		for(Index=0;Index<NOS;Index++)
		{
			if(SiteNo != -1)
				fprintf(Str, "%s - S(%d) - P(%c)\t", Name, SiteNo,Opt->Trees->SymbolList[Index]);
			else
				fprintf(Str, "%s P(%c)\t", Name, Opt->Trees->SymbolList[Index]);
		}
	}
	else
	{
		Trees	= Opt->Trees;
		NOS		= Trees->NOSList[SiteNo];
		
		for(Index=0;Index<NOS;Index++)
			fprintf(Str, "%s - S(%d) - P(%c)\t", Name, SiteNo,Trees->SiteSymbols[SiteNo][Index]);
	}
}



void	PrintRatesHeadder(FILE* Str, OPTIONS *Opt)
{
	int			Index;
	int			SiteIndex;
	RECNODE		*RNode;

	if(Opt->Analsis == ANALMCMC)
		fprintf(Str, "Iteration\tLh\tHarmonic Mean\tTree No\t");
	else
		fprintf(Str, "Tree No\tLh\t");


	if(Opt->DataType == CONTINUOUS)
	{
		PrintRatesHeadderCon(Str, Opt);
		return;
	}	

	if(Opt->UseRJMCMC == TRUE)
	{
		fprintf(Str, "No Off Parmeters\t");
		fprintf(Str, "No Off Zero\t");
		fprintf(Str, "Model string\t");
		if(Opt->Model == M_DESCDEP)
			fprintf(Str, "Dep / InDep\t");
	}

	if(Opt->LoadModels == TRUE)
		fprintf(Str, "Model No\t");
	
	if(Opt->UseRModel == FALSE)
	{
		if(Opt->NOSPerSite == FALSE)
		{
			for(Index=0;Index<Opt->NoOfRates;Index++)
				fprintf(Str, "%s\t", Opt->RateName[Index]);
		}
		else
			fprintf(Str, "Mue\t");
	}
	else
		fprintf(Str, "R Model\t");

	if(Opt->UseCovarion == TRUE)
		fprintf(Str, "Covar On to Off\t Covar Off to On\t");

	if(Opt->Model == M_DESCHET)
		fprintf(Str, "No Indep\tNo Dep\tMap\t");
	

	if(Opt->UseKappa == TRUE)
		fprintf(Str, "Kappa\t");

	if(Opt->UseGamma == TRUE)
		fprintf(Str, "Gamma\t");

	PrintLocalRateHeader(Str, Opt);

	PrintEstDataHeader(Str, Opt);

	PrintLocalTransformHeadder(Str, Opt);

	for(SiteIndex=0;SiteIndex<Opt->Trees->NoOfSites;SiteIndex++)
	{
		if((Opt->Trees->NoOfSites == 1) && (Opt->NOSPerSite == FALSE))
			PrintRecNodeHeadder(Str, Opt, "Root", -1);
		else
			PrintRecNodeHeadder(Str, Opt, "Root", SiteIndex);
	}

	for(Index=0;Index<Opt->NoOfRecNodes;Index++)
	{
		RNode = Opt->RecNodeList[Index];

		for(SiteIndex=0;SiteIndex<Opt->Trees->NoOfSites;SiteIndex++)
		{
			if(Opt->Trees->NoOfSites == 1)
				PrintRecNodeHeadder(Str, Opt, RNode->Name, -1);
			else
				PrintRecNodeHeadder(Str, Opt, RNode->Name, SiteIndex);
		}
	}	



	if(Opt->Analsis == ANALML)
		fprintf(Str, "\n");
}

double	TransVarCoVar(int N, double x)
{
	return (x * (double)N) / (double)(N-1);
}

double	CalcR(double CV, double VT1, double VT2)
{
	return CV / (sqrt(VT1) * sqrt(VT2));
}

double	FindMean(double *List, int WSize)
{
	double	Ret=0;
	int		Index;

	for(Index=0;Index<WSize;Index++)
		Ret += List[Index];

	return Ret / WSize;
}

double	FindCorrelation(double *ListX, double *ListY, int WSize)
{
	double	Top;
	double	BotX;
	double	BotY;
	double	XMean;
	double	YMean;
	int		Index;

	XMean = FindMean(ListX, WSize);
	YMean = FindMean(ListY, WSize);

	Top = 0;
	for(Index=0;Index<WSize;Index++)
		Top += (ListX[Index] - XMean) * (ListY[Index] - YMean);

	BotX = 0;
	BotY = 0;
	for(Index=0;Index<WSize;Index++)
	{
		BotX += (ListX[Index] - XMean) * (ListX[Index] - XMean);
		BotY += (ListY[Index] - YMean) * (ListY[Index] - YMean);
	}

	BotX = sqrt(BotX);
	BotY = sqrt(BotY);

	return Top / (BotX * BotY);
}

double	FindYEst(double Alpha, double *Beta, double *Sites, int NoSites)
{
	double	Ret;
	int		Index;

	Ret = Alpha;
	for(Index=0;Index<NoSites;Index++)
		Ret += (Beta[Index] * Sites[Index]);

	return Ret;
}

double	FindERatio(RATES* Rates, OPTIONS *Opt)
{
	TREES	*Trees;
	TREE	*Tree;
	TAXA	*Taxa;
	CONVAR	*CV;
	double	*Y;
	double	*YP;
	double	*TempV;
	int		Index;
	double	Ret;
	double	SSy;
	double	SSe;
	double	Mean;

	Trees	= Opt->Trees;
	Tree	= Trees->Tree[Rates->TreeNo];
	CV		= Tree->ConVars;
	
	Y		= (double*)malloc(sizeof(double) * Trees->NoOfTaxa);
	YP		= (double*)malloc(sizeof(double) * Trees->NoOfTaxa);
	TempV	= (double*)malloc(sizeof(double) * Trees->NoOfTaxa);

	if((Y == NULL) || (YP == NULL) || (TempV == NULL))
		MallocErr();

	for(Index=0;Index<Trees->NoOfTaxa;Index++)
	{
		Taxa = Trees->Taxa[Index];

		Y[Index] = Taxa->Dependant;
		YP[Index] = Taxa->Dependant - FindYEst(CV->Alpha[0], CV->Beta, Taxa->ConData, Trees->NoOfSites);
	}

	Mean = MLFindAlphaReg(Trees, Tree, Y);
	for(Index=0;Index<Trees->NoOfTaxa;Index++)
		Y[Index] -= Mean;

	VectByMatrixMult(Y, Tree->ConVars->InvV, TempV);
	SSy = VectByVectMult(Y, TempV, Trees->NoOfTaxa);

	VectByMatrixMult(YP, Tree->ConVars->InvV, TempV);
	SSe = VectByVectMult(YP, TempV, Trees->NoOfTaxa);

	Ret = (SSy - SSe) / SSy;

 	free(Y);
	free(YP);
	free(TempV);

	return Ret;
}


void FindRSquared(RATES* Rates, OPTIONS *Opt, double *R2, double *SSE, double *SST)
{
	TREES	*Trees;
	TREE	*Tree;
	TAXA	*Taxa;
	CONVAR	*CV;
	double	*Y;
	double	*YP;
	double	*TempV;
	int		Index;
	double	MeanY;
	double	MeanYP;
	double	T,B1,B2;

	Trees	= Opt->Trees;
	Tree	= Trees->Tree[Rates->TreeNo];
	CV		= Tree->ConVars;
	
	Y		= (double*)malloc(sizeof(double) * Trees->NoOfTaxa);
	YP		= (double*)malloc(sizeof(double) * Trees->NoOfTaxa);
	TempV	= (double*)malloc(sizeof(double) * Trees->NoOfTaxa);

	if((Y == NULL) || (YP == NULL) || (TempV == NULL))
		MallocErr();

	for(Index=0;Index<Trees->NoOfTaxa;Index++)
	{
		Taxa = Trees->Taxa[Index];

		Y[Index] = Taxa->Dependant;
		YP[Index] = FindYEst(CV->Alpha[0], CV->Beta, Taxa->ConData, Trees->NoOfSites);
	}

	MeanY = MLFindAlphaReg(Trees, Tree, Y);
	MeanYP= MLFindAlphaReg(Trees, Tree, YP);
	for(Index=0;Index<Trees->NoOfTaxa;Index++)
	{
		Y[Index] -= MeanY;
		YP[Index] -= MeanYP;
	}

	VectByMatrixMult(Y, Tree->ConVars->InvV, TempV);
	T = VectByVectMult(YP, TempV, Trees->NoOfTaxa);
	T = T * T;

	VectByMatrixMult(Y, Tree->ConVars->InvV, TempV);
	B1 = VectByVectMult(Y, TempV, Trees->NoOfTaxa);

	VectByMatrixMult(YP, Tree->ConVars->InvV, TempV);
	B2 = VectByVectMult(YP, TempV, Trees->NoOfTaxa);

	(*R2) = T / (B1 * B2);
	
	(*SSE) = (1-(*R2)) * B1;

	(*SST) = B1;

 	free(Y);
	free(YP);
	free(TempV);
}


void	PrintRegVarCoVar(FILE* Str, RATES *Rates, OPTIONS *Opt)
{
	MATRIX	*Var;
	int		Index;
	TREES	*Trees;

	Trees = Opt->Trees;

	Var = FindRegVar(Opt->Trees, Rates, Opt->AlphaZero);

	if(Opt->AlphaZero == FALSE)
	{
		for(Index=0;Index<Trees->NoOfSites+1;Index++)
			fprintf(Str, "%f\t", sqrt(Var->me[Index][Index]));
	}
	else
	{
		fprintf(Str, "0\t");
		for(Index=0;Index<Trees->NoOfSites;Index++)
			fprintf(Str, "%f\t", sqrt(Var->me[Index][Index]));
	}

	FreeMatrix(Var);
}

void	PrintConRecNodes(FILE *Str, RATES* Rates, OPTIONS *Opt)
{
	int			Index, SiteIndex;
	RECNODE		*RNode;
	NODE		N;
	CONTRAST	*Con;
	double		Alpha; //, Sigma, Lh;
		
	for(Index=0;Index<Opt->NoOfRecNodes;Index++)
	{
		RNode = Opt->RecNodeList[Index];
		N = RNode->TreeNodes[Rates->TreeNo];

		if(N == NULL)
		{
			for(SiteIndex=0;SiteIndex<Opt->Trees->NoOfSites;SiteIndex++)
				fprintf(Str, "--\t--\t--\t");
		}
		else
		{
			Con = N->ConData->Contrast[0];

			for(SiteIndex=0;SiteIndex<Opt->Trees->NoOfSites;SiteIndex++)
			{
				RecIntNode(N, SiteIndex, &Alpha);
				fprintf(Str, "%f\t", Alpha);
			}
		}
	}
}

void	PrintLocalTransformNo(FILE* Str, RATES* Rates, OPTIONS *Opt)
{
	if(Opt->UseVarRates == TRUE)
	{
		fprintf(Str, "%d\t", GetNoTransformType(VR_BL, Rates));
		fprintf(Str, "%d\t", GetNoTransformType(VR_NODE, Rates));
	}

	if(Opt->UseRJLocalScalar[VR_KAPPA] == TRUE)
		fprintf(Str, "%d\t", GetNoTransformType(VR_KAPPA, Rates));

	if(Opt->UseRJLocalScalar[VR_LAMBDA] == TRUE)
		fprintf(Str, "%d\t", GetNoTransformType(VR_LAMBDA, Rates));

	if(Opt->UseRJLocalScalar[VR_DELTA] == TRUE)
		fprintf(Str, "%d\t", GetNoTransformType(VR_DELTA, Rates));

	if(Opt->UseRJLocalScalar[VR_OU] == TRUE)
		fprintf(Str, "%d\t", GetNoTransformType(VR_OU, Rates));
}

void	PrintRateLocalTransform(FILE *Str, RATES *Rates)
{
	int Index;

	for(Index=0;Index<Rates->NoLocalTransforms;Index++)
		if(Rates->LocalTransforms[Index]->Est == TRUE)
			fprintf(Str, "%0.12f\t", Rates->LocalTransforms[Index]->Scale);
}

void	PrintRatesCon(FILE* Str, RATES* Rates, OPTIONS *Opt)
{
	int		Index;
	int		x,y, NOS;
	CONVAR	*ConVar;
	int		MinNodes, MaxNodes;
	TAXA	*Taxa;
	double	R2, SSE, SST;
	TREES	*Trees;

	Trees = Opt->Trees;
	NOS = Trees->NoOfSites;
	ConVar = Opt->Trees->Tree[Rates->TreeNo]->ConVars;

	if(Opt->LoadModels == TRUE)
		fprintf(Str, "%d\t", Rates->ModelNo);

	if((Opt->Model == M_CONTINUOUS_RR) || (Opt->Model == M_CONTINUOUS_DIR))
	{
		for(Index=0;Index<Trees->NoOfSites;Index++)
		{
			fprintf(Str, "%0.12f\t", ConVar->Alpha[Index]);
			if(Opt->Model == M_CONTINUOUS_DIR)
				fprintf(Str, "%0.12f\t", ConVar->Beta[Index]);
		}	
			for(Index=0;Index<Trees->NoOfSites;Index++)
				fprintf(Str, "%0.12f\t", TransVarCoVar(Opt->Trees->NoOfTaxa, ConVar->Sigma->me[Index][Index]));
		
			for(x=0;x<Trees->NoOfSites;x++)
				for(y=x+1;y<Trees->NoOfSites;y++)
					fprintf(Str, "%0.12f\t", CalcR(ConVar->Sigma->me[x][y], ConVar->Sigma->me[x][x], ConVar->Sigma->me[y][y]));
	}

	if(Opt->Model == M_CONTRAST_CORREL)
	{
//		PrintMatrix(Rates->Contrast->SigmaMat, "Sig = ", stdout);exit(0);

		for(Index=0;Index<NOS;Index++)
			fprintf(Str, "%0.12f\t", Rates->Contrast->Alpha[Index]);

		for(Index=0;Index<NOS;Index++)
			fprintf(Str, "%0.12f\t", Rates->Contrast->SigmaMat->me[Index][Index]);

		for(x=0;x<NOS;x++)
		{
			for(y=0;y<x;y++)
				fprintf(Str, "%0.12f\t", CalcR(Rates->Contrast->SigmaMat->me[x][y], Rates->Contrast->SigmaMat->me[x][x], Rates->Contrast->SigmaMat->me[y][y]));
			//	fprintf(Str, "%0.12f\t", Rates->Contrast->SigmaMat->me[x][y]);
		}

	}
	
	if(Opt->Model == M_CONTRAST)
	{
		for(Index=0;Index<NOS;Index++)
			fprintf(Str, "%0.12f\t", Rates->Contrast->Alpha[Index]);

		for(Index=0;Index<NOS;Index++)
			fprintf(Str, "%0.12f\t", Rates->Contrast->Sigma[Index]);
	}

	if(Opt->Model == M_CONTRAST_REG)
		OutputConReg(Str, Opt, Trees, Rates);

	if(Opt->Model == M_CONTINUOUS_REG)
	{
		fprintf(Str, "%0.12f\t", ConVar->Alpha[0]);
		for(Index=0;Index<Trees->NoOfSites;Index++)
			fprintf(Str, "%0.12f\t", ConVar->Beta[Index]);

		fprintf(Str, "%0.12f\t", ConVar->Sigma->me[0][0]);

		FindRSquared(Rates, Opt, &R2, &SSE, &SST);
		fprintf(Str, "%0.12f\t%0.12f\t%0.12f\t", R2, SSE, SST);
	
		if(Opt->Analsis == ANALML)
			fprintf(Str, "%0.12f\t", FindERatio(Rates, Opt));

 		PrintRegVarCoVar(Str, Rates, Opt);
	}

	if(Opt->ModelType == MT_FATTAIL)
	{
		if(Opt->Model == M_GEO)
			fprintf(Str, "%0.12f\t%0.12f\t", Rates->FatTailRates->Alpha[0], Rates->FatTailRates->Scale[0]);
		else
		{
			for(Index=0;Index<Trees->NoOfSites;Index++)
				fprintf(Str, "%0.12f\t%0.12f\t", Rates->FatTailRates->Alpha[Index], Rates->FatTailRates->Scale[Index]);
		}
	}


	if(Opt->EstKappa == TRUE)
		fprintf(Str, "%0.12f\t", Rates->Kappa);

	if(Opt->FixKappa != -1)
		fprintf(Str, "%0.12f\t", Opt->FixKappa);

	if(Opt->EstDelta == TRUE)
		fprintf(Str, "%0.12f\t", Rates->Delta);

	if(Opt->FixDelta != -1)
		fprintf(Str, "%0.12f\t", Opt->FixDelta);

	if(Opt->EstLambda == TRUE)
		fprintf(Str, "%0.12f\t", Rates->Lambda);

	if(Opt->FixLambda != -1)
		fprintf(Str, "%0.12f\t", Opt->FixLambda);

	if(Opt->EstOU == TRUE)
		fprintf(Str, "%0.12f\t", Rates->OU);

	if(Opt->FixOU != -1)
		fprintf(Str, "%0.12f\t", Opt->FixOU);

	PrintRateLocalTransform(Str, Rates);

	if(Opt->NodeBLData == TRUE)
	{
		fprintf(Str, "%f\t", ConVar->Sigma->me[0][1] / ConVar->Sigma->me[0][0]);
		fprintf(Str, "%f\t", ConVar->Sigma->me[0][1] / ConVar->Sigma->me[1][1]);
		
		MinNodes = MaxNodes = (int)Opt->Trees->Taxa[0]->ConData[0];
		for(Index=1;Index<Opt->Trees->NoOfTaxa;Index++)
		{
			Taxa = Opt->Trees->Taxa[Index];
			if((int)Taxa->ConData[0] > MaxNodes)
				MaxNodes = (int)Taxa->ConData[0];

			if((int)Taxa->ConData[0] < MinNodes)
				MinNodes = (int)Taxa->ConData[0];
		}

		fprintf(Str, "%d\t%d\t", MinNodes, MaxNodes);
	}

	for(Index=0;Index<Rates->NoEstData;Index++)
		fprintf(Str, "%0.12f\t", Rates->EstData[Index]);

	PrintConRecNodes(Str, Rates, Opt);


	PrintLocalTransformNo(Str, Rates, Opt);

	if(Opt->UseDistData == TRUE)
		OutputDataDist(Str, Rates, Opt);
}

double	GetPartailPi(RATES *Rates, NODE N, int StateNo, int SiteNo)
{
	return N->Partial[SiteNo][StateNo] * Rates->Pis[StateNo];
}


void	PrintNodeRecDep(RATES *Rates, OPTIONS *Opt, FILE *Str, double Total, NODE Node)
{
		if(	(Node->Partial[0][0]/Total > 1) ||
			(Node->Partial[0][1]/Total > 1) ||
			(Node->Partial[0][2]/Total > 1) ||
			(Node->Partial[0][3]/Total > 1))
		{
			Total = Total;
//			printf("Err\n");
		}
			
		if(Opt->UseCovarion == FALSE)
		{
			fprintf(Str, "%f\t", (Node->Partial[0][0])/Total);
			fprintf(Str, "%f\t", (Node->Partial[0][1])/Total);
			fprintf(Str, "%f\t", (Node->Partial[0][2])/Total);
			fprintf(Str, "%f\t", (Node->Partial[0][3])/Total);
		}
		else
		{
			fprintf(Str, "%f\t", (Node->Partial[0][0] + Node->Partial[0][4])/Total);
			fprintf(Str, "%f\t", (Node->Partial[0][1] + Node->Partial[0][5])/Total);
			fprintf(Str, "%f\t", (Node->Partial[0][2] + Node->Partial[0][6])/Total);
			fprintf(Str, "%f\t", (Node->Partial[0][3] + Node->Partial[0][7])/Total);
		}
}

void	PrintNodeRec(FILE *Str, NODE Node, int NOS, int NoOfSites, RATES* Rates, OPTIONS *Opt)
{
	double	Tot=0;
	int		Index;
	int		SiteIndex;
	int		TrueNOS;
	TREES	*Trees;

	Trees = Opt->Trees;



	if(Node == NULL)
	{
		for(SiteIndex=0;SiteIndex<NoOfSites;SiteIndex++)
			for(Index=0;Index<NOS;Index++)
				fprintf(Str, "--\t");
		return;
	}

#ifdef BIG_LH
	SetBigLhNodeRec(Node, NOS, NoOfSites, Rates, Opt);
#endif

#ifdef QUAD_DOUBLE
	SetQuadDoubleNodeRec(Node, NOS, NoOfSites, Rates, Opt);
#endif

	for(SiteIndex=0;SiteIndex<NoOfSites;SiteIndex++)
	{
		if(Opt->NOSPerSite == FALSE)
			NOS = Trees->NoOfStates;
		else
		{
			NOS = Trees->NOSList[SiteIndex];
			for(Index=0;Index<NOS;Index++)
				Rates->Pis[Index] = (double)1/NOS;
		}

		Tot=0;
		for(Index=0;Index<NOS;Index++)
			Tot += (Node->Partial[SiteIndex][Index] * Rates->Pis[Index]);
//			Tot += (Node->Partial[SiteIndex][Index]);

	
	/*
		if(Opt->Model == DESCINDEP)
		{
			if(Opt->UseCovarion == FALSE)
			{
				fprintf(Str, "%f\t", (Node->Partial[0][0] + Node->Partial[0][1])/Tot);
				fprintf(Str, "%f\t", (Node->Partial[0][2] + Node->Partial[0][3])/Tot);
				fprintf(Str, "%f\t", (Node->Partial[0][0] + Node->Partial[0][2])/Tot);
				fprintf(Str, "%f\t", (Node->Partial[0][1] + Node->Partial[0][3])/Tot);
			}
			else
			{
				fprintf(Str, "%f\t", (Node->Partial[0][0] + Node->Partial[0][1] + Node->Partial[0][4] + Node->Partial[0][5])/Tot);
				fprintf(Str, "%f\t", (Node->Partial[0][2] + Node->Partial[0][3] + Node->Partial[0][6] + Node->Partial[0][7])/Tot);
				fprintf(Str, "%f\t", (Node->Partial[0][0] + Node->Partial[0][2] + Node->Partial[0][4] + Node->Partial[0][6])/Tot);
				fprintf(Str, "%f\t", (Node->Partial[0][1] + Node->Partial[0][3] + Node->Partial[0][5] + Node->Partial[0][7])/Tot);
			}
		}

		if(Opt->Model == DESCDEP)
			PrintNodeRecDep(Rates, Opt, Str, Tot, Node);
		

		if(Opt->Model == DESCCV)
		{
			for(Index=0;Index<NOS;Index++)
				fprintf(Str, "%f\t", Node->Partial[0][Index] / Tot);
		}

		if(Opt->Model == DESCHET)
		{
			for(Index=0;Index<NOS;Index++)
				fprintf(Str, "%f\t", Node->Partial[0][Index] / Tot);
		}
		*/
//		if(Opt->Model == MULTISTATE)
		{
			if(Opt->UseCovarion == FALSE)
			{
				for(Index=0;Index<NOS;Index++)
					fprintf(Str, "%f\t", (Node->Partial[SiteIndex][Index] *  Rates->Pis[Index]) / Tot);
//					fprintf(Str, "%f\t", (Node->Partial[SiteIndex][Index]) / Tot);
			}
			else
			{
				TrueNOS = NOS / 2;
				for(Index=0;Index<TrueNOS;Index++)
					fprintf(Str, "%f\t", ((Node->Partial[SiteIndex][Index] *  Rates->Pis[Index]) +
										  (Node->Partial[SiteIndex][Index+TrueNOS] *  Rates->Pis[Index+TrueNOS]))/ Tot);
			}

		}
	}
}

char	RJModelType(int *ModelStr)
{

	if(	(ModelStr[0] == ModelStr[5]) &&
		(ModelStr[1] == ModelStr[3]) &&
		(ModelStr[2] == ModelStr[7]) &&
		(ModelStr[4] == ModelStr[6])
		)
		return 'I';
	else
		return 'D';
}

int		NoZeroRate(RATES *Rates)
{
	int Ret, Index;
	Ret = 0;

	for(Index=0;Index<Rates->NoOfRates;Index++)
		if(Rates->MappingVect[Index] == ZERORATENO)
			Ret++;


	return Ret;
}

void	PrintHetro(FILE* Str, RATES *Rates)
{
	int NoI, NoD, Index;


	NoI = NoD = 0;
	for(Index=1;Index<Rates->Hetero->MListSize;Index++)
	{
		if(Rates->Hetero->MList[Index] == 0)
			NoI++;

		if(Rates->Hetero->MList[Index] == 1)
			NoD++;
	}

	fprintf(Str, "%d\t%d\t", NoI, NoD);
	for(Index=1;Index<Rates->Hetero->MListSize-1;Index++)
		fprintf(Str, "%d,", Rates->Hetero->MList[Index]);
	fprintf(Str, "%d\t", Rates->Hetero->MList[Index]);
}

void	PrintAutoTuneDevAcc(FILE *Str,AUTOTUNE *AT)
{
	fprintf(Str,"%f\t", AT->CDev);
	fprintf(Str,"%d\t",AT->NoTried);
	fprintf(Str,"%f\t", AutoTuneCalcAcc(AT));
}

void	PrintAutoTune(FILE* Str, OPTIONS *Opt, SCHEDULE* Shed)
{
	int Index;

	for(Index=0;Index<Shed->NoFullATList;Index++)
		PrintAutoTuneDevAcc(Str, Shed->FullATList[Index]);
}

void	PrintRates(FILE* Str, RATES* Rates, OPTIONS *Opt, SCHEDULE* Shed)
{
	int		Index;
	

	if(Opt->Analsis == ANALML)
		fprintf(Str, "%d\t%f\t", Rates->TreeNo+1,Rates->Lh);

	if(Opt->DataType == CONTINUOUS)
	{
		PrintRatesCon(Str, Rates, Opt);
		return;
	}
	
	if(Opt->UseRJMCMC == TRUE)
	{
		fprintf(Str, "%d\t", NoOfPramGroups(Rates, NULL, NULL));
		fprintf(Str, "%d\t", NoZeroRate(Rates));
		fprintf(Str, "'");
//		for(Index=0;Index<Rates->NoOfFullRates;Index++)
		for(Index=0;Index<Rates->NoOfRates;Index++)
		{
			if(Rates->MappingVect[Index] == ZERORATENO)
				fprintf(Str, "Z ");
//				fprintf(Str, "Z");
			else
			{
				// TODO Phoneim remove. 
				fprintf(Str, "%d ",  Rates->MappingVect[Index]);
//				if(Rates->MappingVect[Index] <= 9)
//					fprintf(Str, "%d", Rates->MappingVect[Index]);
//				else
//					fprintf(Str, "%c", Rates->MappingVect[Index] + 'A');
			}
		}
		fprintf(Str, "\t");

		if(Opt->Model == M_DESCDEP)
			fprintf(Str, "%c\t", RJModelType(Rates->MappingVect));
	}

	if(Opt->LoadModels == TRUE)
		fprintf(Str, "%d\t", Rates->ModelNo);

	if(Opt->UseRModel == TRUE)
		fprintf(Str, "%f\t", Rates->FullRates[0]);
	else
	{
		if(Opt->NOSPerSite == FALSE)
		{
			for(Index=0;Index<Opt->NoOfRates;Index++)
				fprintf(Str, "%f\t", Rates->FullRates[Index]);
		}
		else
			fprintf(Str, "%f\t", Rates->FullRates[0]);
	}

	if(Opt->UseCovarion == TRUE)
		fprintf(Str, "%f\t%f\t", Rates->OffToOn, Rates->OnToOff);

	if(Opt->Model == M_DESCHET)
		PrintHetro(Str, Rates);

	if(Opt->UseKappa == TRUE)
		fprintf(Str, "%f\t", Rates->Kappa);

	if(Opt->UseGamma == TRUE)
		fprintf(Str, "%f\t", Rates->Gamma);

	PrintRateLocalTransform(Str, Rates);
	
	for(Index=0;Index<Rates->NoEstData;Index++)
		fprintf(Str, "%c\t", Opt->Trees->SymbolList[Rates->EstDescData[Index]]);

	PrintLocalTransformNo(Str, Rates, Opt);

	PrintNodeRec(Str, Opt->Trees->Tree[Rates->TreeNo]->Root, Opt->Trees->NoOfStates, Opt->Trees->NoOfSites, Rates, Opt);

	for(Index=0;Index<Opt->NoOfRecNodes;Index++)
		PrintNodeRec(Str, Opt->RecNodeList[Index]->TreeNodes[Rates->TreeNo], Opt->Trees->NoOfStates, Opt->Trees->NoOfSites, Rates, Opt);

}

void	CopyRJRtaes(RATES *A, RATES *B, OPTIONS *Opt)
{
	A->NoOfRJRates = B->NoOfRJRates ;

	memcpy(A->Rates, B->Rates, sizeof(double)*B->NoOfRJRates);
	memcpy(A->MappingVect, B->MappingVect, sizeof(int)*B->NoOfRates);
}

void	CopyRates(RATES *A, RATES *B, OPTIONS *Opt)
{
	int	Index;

	A->Delta = B->Delta;
	A->Kappa = B->Kappa;
	A->Lambda= B->Lambda;
	A->OU	 = B->OU;

	A->HMeanCount	= B->HMeanCount;
	
#ifndef BIG_LH
	A->HMeanSum		= B->HMeanSum;
#else
	mpfr_set(A->HMeanSum, B->HMeanSum, DEF_ROUND);
#endif

	if(Opt->UseDistData == TRUE)
		CopyDistDataRates(A->DistDataRates, B->DistDataRates);

	A->ModelNo		= B->ModelNo;

	if(Opt->UseRJMCMC == FALSE)
	{
		if(A->Rates != NULL)
			memcpy(A->Rates, B->Rates, sizeof(double) * A->NoOfRates);
	}
	else
		CopyRJRtaes(A, B, Opt);

	A->TreeNo		= B->TreeNo;
	A->Lh			= B->Lh;
	A->LhPrior		= B->LhPrior;
	A->LnJacobion	= 0;
	A->LnHastings	= 0;

	if(Opt->UseCovarion==TRUE)
	{
		A->OffToOn	= B->OffToOn;
		A->OnToOff	= B->OnToOff;
	}

	if(Opt->Analsis == ANALMCMC)
	{
		A->NoPriors = A->NoPriors;
		for(Index=0;Index<A->NoPriors;Index++)
			CopyPrior(A->Priors[Index], B->Priors[Index]);
	}

	if(Opt->UseGamma == TRUE)
	{
		A->Gamma	= B->Gamma;
		A->GammaCats= B->GammaCats;

		memcpy(A->GammaMults, B->GammaMults, sizeof(double) * A->GammaCats);
	}

	if(B->UseEstData == TRUE)
	{
		A->NoEstData = B->NoEstData;

		if(Opt->DataType == CONTINUOUS)
			memcpy(A->EstData, B->EstData, sizeof(double)*A->NoEstData);
		else
			memcpy(A->EstDescData, B->EstDescData, sizeof(int)*A->NoEstData);
	}

	if(Opt->ModelType == MT_CONTRAST)
		CopyContrastRates(Opt, A, B, Opt->Trees->NoOfSites);
	
	if(UseNonParametricMethods(Opt) == TRUE)
		VarRatesCopy(A, B);

	if(A->Hetero != NULL)
		CopyHetero(A->Hetero, B->Hetero);

	if(Opt->RJDummy == TRUE)
		RJDummyCopy(A, B);

	if(Opt->ModelType == MT_FATTAIL)
		CopyFatTailRates(Opt->Trees, A->FatTailRates, B->FatTailRates);

	for(Index=0;Index<B->NoLocalTransforms;Index++)
		CopyLocalTransforms(A->LocalTransforms[Index], B->LocalTransforms[Index]);
} 

double ChangeRatesTest(RATES *Rates, double RateV, double dev)
{
	int Index;
	double Scale;

	for(Index=0;Index<10000;Index++)
	{
		Scale = exp(dev * (RandDouble(Rates->RS) - 0.5));
		
		printf("%d\t%f\t%f\n", Index, dev, Scale);
	}
	exit(0);
}
/*
double ChangeRate(RATES *Rates, double RateV, double dev)
{
	int		Exit;
	double	Ret, Scale;

//	ChangeRatesTest(Rates, 1, dev);

	if(RateV >= MAXRATE)
		return  MAXRATE;

	do
	{
		Exit = TRUE;

		Scale = exp(dev * (RandDouble(Rates->RS) - 0.5));

		Ret = RateV * Scale;

		if(Ret > MAXRATE)
			Exit = FALSE;

		if(Ret < MINRATE)
			Exit = FALSE;

	} while(Exit == FALSE);


	// Working ish with 1. 
	Rates->LnHastings += log(Ret / RateV);

//	Rates->LnHastings += log(Ret/ RateV) / dev;

	return Ret;
}
*/

double ChangeRate(RATES *Rates, double RateV, double dev)
{
	int		Exit;
	double	Ret;

	if(RateV >= MAXRATE)
		return  MAXRATE;

	do
	{
		Exit = TRUE;
#ifdef RATE_CHANGE_UNI
		Ret = (RandDouble(Rates->RS) * dev) - (dev / 2.0); 
		Ret += RateV;
#endif

#ifdef RATE_CHANGE_NORM
		Ret = RandNormal(Rates->RS, RateV, dev);
#endif
		if(Ret > MAXRATE)
			Exit = FALSE;

		if(Ret < MINRATE)
			Exit = FALSE;

	} while(Exit == FALSE);

	return Ret;
}

void	ChangeGammaRates(RATES *Rates, SCHEDULE* Shed)
{
	double Dev;

	Shed->CurrentAT = Shed->GammaAT;
	Dev = Shed->CurrentAT->CDev;

	Rates->Gamma =  ChangeRate(Rates, Rates->Gamma, Dev);
}

double	MultePram(RATES *Rates, double Val, double Min, double Max, double Dev)
{
	double	Ret;
	int		Exit;

	do
	{
		Exit = TRUE;

		Ret = (RandDouble(Rates->RS) * Dev) - (Dev / 2.0); 
		Ret += Val;
//		Ret = RandNormal(Rates->RS, Val, Dev);

		if(Ret > Max)
			Exit = FALSE;

		if(Ret < Min)
			Exit = FALSE;

	} while(Exit == FALSE);

	return Ret;		
}


void	TestMult(RATES *Rates, double Val, double Min, double Max, double Dev)
{
	int Index;
	double Ret;

	Val = 10;
	Dev = 5;

	for(Index=0;Index<10000;Index++)
	{
		Ret = MultePram(Rates, Val, Min, Max, Dev);
		printf("%d\t%f\n", Index, Ret);
	}

	exit(0);
}


int		ValidMove(RATES *Rates, int No)
{
	VARRATES *PP;

	if(No == SPPMOVE || No == SPPCHANGESCALE)
	{
		PP = Rates->VarRates;
		if(PP->NoNodes == 0)
			return FALSE;
	}

	if(No == SRJDUMMYMOVE || No == SRJDUMMYCHANGEBETA)
	{
		if(Rates->RJDummy->NoDummyCode == 0)
			return FALSE;
	}


	return TRUE;
}

int		PickFromVect(RANDSTATES *RS, double *Vect, int Size)
{
	double	Val, Sum;
	int		Index;

	Val = RandDouble(RS);
	Sum = 0;
	for(Index=0;Index<Size;Index++)
	{
		if((Val > Sum) && (Val <= (Sum + Vect[Index])))
			return Index;

		Sum += Vect[Index];
	}

	printf("Error in PickACat %s::%d\n", __FILE__, __LINE__);
	exit(0);
	return -1;
}

int	PickACat(RATES *Rates, double *Vect, int Size)
{
	int		Cat;

	do
	{
		Cat = PickFromVect(Rates->RS, Vect, Size);
	}while(ValidMove(Rates, Cat) == FALSE);

	return Cat;
}

int		NumInList(int *List, int No, int Size)
{
	int	i;

	for(i=0;i<Size;i++)
		if(List[i] == No)
			return TRUE;
	
	return FALSE;
}

int*	PickEstChangeSites(int No, int Max, RANDSTATES *RS)
{
	int *Ret;
	int	Index;
	int	Pick;

	Ret = (int*)malloc(sizeof(int) * No);
	if(Ret == NULL)
		MallocErr();
	for(Index=0;Index<No;Index++)
		Ret[Index] = -1;

	if(No >= Max)
	{
		for(Index=0;Index<Max;Index++)
			Ret[Index] = Index;

		return Ret;
	}

	for(Index=0;Index<No;Index++)
	{
		do
		{
			Pick = RandUSLong(RS) % Max;
		} while (NumInList(Ret, Pick, Index) == TRUE);

		Ret[Index] = Pick;
	}

	return Ret;
}

double*	GetMultVarChanges(RATES *Rates, OPTIONS *Opt, SCHEDULE* Shed)
{
	double	*Ret;
	TREE	*Tree;
	TREES	*Trees;
	int		Index;

	Shed->CurrentAT = Shed->DataDevAT;

	Trees = Opt->Trees;
	Tree = Trees->Tree[Rates->TreeNo];
		
	Ret = (double*)malloc(sizeof(double) * Trees->NoOfTaxa);
	if(Ret == NULL)
		MallocErr();

	genmn(Tree->ConVars->MultiVarNormState, Ret, Tree->ConVars->MultiVarNormTemp);

	for(Index=0;Index<Trees->NoOfTaxa;Index++)
		Ret[Index] = Shed->CurrentAT->CDev * Ret[Index];
	
	return Ret;
}
/*
double*	GetMultVarChanges(RATES *Rates, OPTIONS *Opt)
{
	double	*Ret;
	TREE	*Tree;
	TREES	*Trees;
	MATRIX	*VarCo;
	MATRIX	*Changes;
	int		Index;

	Trees = Opt->Trees;
	Tree = Trees->Tree[Rates->TreeNo];

	VarCo = AllocMatrix(Trees->NoOfTaxa, Trees->NoOfTaxa);
	
	CopyMatrix(VarCo, Tree->ConVars->V);

	Changes = MultivariateNormal(1, VarCo);
	
	Ret = (double*)malloc(sizeof(double) * Trees->NoOfTaxa);
	if(Ret == NULL)
		MallocErr();

	for(Index=0;Index<Trees->NoOfTaxa;Index++)
		Ret[Index] = Opt->EstDataDev * Changes->me[0][Index];

	FreeMatrix(Changes);
	FreeMatrix(VarCo);

	return Ret;
}
*/

void	MutateEstRatesDiscrete(OPTIONS* Opt, RATES* Rates, SCHEDULE* Shed)
{
	int		Old, New;
	int		Site;

	Site = RandUSLong(Rates->RS) % Rates->NoEstData;
	Old = Rates->EstDescData[Site];
//	do
//	{
		if(Opt->Model == M_MULTISTATE)
			New = RandUSLong(Rates->RS) % Opt->Trees->NoOfStates;
		else
			New = RandUSLong(Rates->RS) % 2;
//	} while(New == Old);
	Rates->EstDescData[Site] = New;
}

void	MutateEstRates(OPTIONS* Opt, RATES* Rates, SCHEDULE* Shed)
{
	double	*Changes;
	
	int		RIndex, SIndex, TIndex;
	TAXA	*Taxa;
	TREES	*Trees;
	int		Site;
	
	Trees = Opt->Trees;

	if(Opt->DataType == DISCRETE)
	{
		MutateEstRatesDiscrete(Opt, Rates, Shed);
		return;
	}

//	printf("Hello\n");

//	Change1EstData(Opt, Rates); return;

	Changes	=	GetMultVarChanges(Rates, Opt, Shed);

//	Changes =	GetPhyChanges(Trees, Trees->Tree[Rates->TreeNo], Opt->EstDataDev, Rates->RS);
			
	Site	=	Opt->EstDataSites[RandUSLong(Rates->RS) % Opt->NoEstDataSite];
	RIndex	=	0;
	for(TIndex=0;TIndex<Trees->NoOfTaxa;TIndex++)
	{
		Taxa = Trees->Taxa[TIndex];

		for(SIndex=0;SIndex<Trees->NoOfSites;SIndex++)
		{
			if(Taxa->EstDataP[SIndex] == TRUE)
			{
				if(SIndex != Site)
					RIndex++;
				else
					Rates->EstData[RIndex++] += Changes[TIndex];
			}
		}

		if(Taxa->EstDepData == TRUE)
		{
			if(Site != -1)
				RIndex++;
			else
				Rates->EstData[RIndex++] += Changes[TIndex];
		}
	}

	free(Changes);
}

int		ForceMerge(OPTIONS *Opt, RATES *Rates, int NoOfGroups)
{
	if(NoOfGroups == Rates->NoOfRates)
		return TRUE;

	if(Opt->CapRJRatesNo != -1)
	{
		if(NoOfGroups == Opt->CapRJRatesNo)
			return TRUE;
	}

	return FALSE;
}

int		TryRJMove(OPTIONS* Opt, RATES* Rates, SCHEDULE* Shed)
{
	int NoOfGroups;

	NoOfGroups = NoOfPramGroups(Rates, NULL, NULL);

	if(RandDouble(Rates->RS) < 0.25)
	{
		if(RandDouble(Rates->RS) < 0.5)
			return RJAugment(Rates, Opt);
		else
			return RJReduce(Rates, Opt);
	}

	if(RandDouble(Rates->RS) < 0.5)
		return RJSplit(Rates, Opt);
	
	return RJMerge(Rates, Opt);
}


void	RJMove(OPTIONS* Opt, RATES* Rates, SCHEDULE* Shed)
{
	int Success;

	do
	{
		Success = TryRJMove(Opt, Rates, Shed);
	} while(Success == FALSE);
}

void	ChangeConRates(OPTIONS* Opt, RATES* Rates, SCHEDULE* Shed)
{
	int Index, Pos;
	double Dev;

	Pos = RandUSLong(Rates->RS) % Rates->NoOfRates;
	Shed->PNo = Pos;
	Shed->CurrentAT = Shed->RateDevATList[Pos];

	Dev = Shed->CurrentAT->CDev;

//	Uniform
	Rates->Rates[Shed->PNo] += (RandDouble(Rates->RS) * Dev) - (Dev / 2.0);

//		Normal Does not seem to work well. 
//		Rates->Rates[Shed->PNo] = RandNormal(Rates->RS, Rates->Rates[Shed->PNo], Dev); 

	if(Opt->AlphaZero == TRUE)
	{
		if(Opt->Model == M_CONTINUOUS_REG)
			Rates->Rates[0] = 0;
		else
		{
			for(Index=0;Index<Opt->Trees->NoOfSites;Index++)
				Rates->Rates[Index] = 0;
		}
	}
}

void	ChangeRates(OPTIONS* Opt, RATES* Rates, SCHEDULE* Shed, long long It)
{
	int Index, NoOfRates;

	if(Opt->LoadModels == TRUE)
	{
		ChangeModelFile(Rates, Rates->RS);
		return;
	}

	if(Opt->DataType == DISCRETE)
	{
		NoOfRates = Rates->NoOfRates;
		if(Opt->UseRJMCMC == TRUE)
			NoOfRates = Rates->NoOfRJRates;
		
		Shed->CurrentAT = Shed->RateDevATList[0];

#ifdef RATE_CHANGE_ONE
		Index = RandUSLong(Rates->RS) % NoOfRates;
		Rates->Rates[Index] = ChangeRates(Rates, Rates->Rates[Index], Shed->CurrentAT->CDev);
#else
		for(Index=0;Index<NoOfRates;Index++)
			Rates->Rates[Index] = ChangeRate(Rates, Rates->Rates[Index], Shed->CurrentAT->CDev);
#endif
		Shed->PNo = 0;
	}
	else
	{
		if(Opt->ModelType == MT_FATTAIL)
		{
			MutateFatTailRates(Opt, Opt->Trees, Rates, Shed);
			return;
		}

		if(Opt->ModelType == MT_CONTRAST)
			MutateContrastRates(Opt, Opt->Trees, Rates, Shed);
		else
			ChangeConRates(Opt, Rates, Shed);
	}
}

void	ChangeCVRates(OPTIONS* Opt, RATES* Rates, SCHEDULE* Shed)
{
	double Dev;

	Shed->CurrentAT = Shed->RateDevATList[0];
	Dev = Shed->CurrentAT->CDev;
	
	Rates->OffToOn = ChangeRate(Rates, Rates->OffToOn, Dev);
	Rates->OnToOff = ChangeRate(Rates, Rates->OnToOff, Dev);

	// Set the off / on rate to the same. 
	Rates->OffToOn = Rates->OnToOff;
}


void	ChangeKappa(OPTIONS* Opt, RATES* Rates, SCHEDULE* Shed)
{
	double Dev;

	Shed->CurrentAT = Shed->KappaAT;
	Dev = Shed->CurrentAT->CDev;

	if(Dev > MAX_KAPPA)
		Rates->Kappa = RandUniDouble(Rates->RS, MIN_KAPPA, MAX_KAPPA);
	else
		Rates->Kappa = MultePram(Rates, Rates->Kappa, MIN_KAPPA, MAX_KAPPA, Dev);
}

void	ChangeDelta(OPTIONS* Opt, RATES* Rates, SCHEDULE* Shed)
{
	double Dev;

	Shed->CurrentAT = Shed->DeltaAT;
	Dev = Shed->CurrentAT->CDev;

	if(Dev > MAX_DELTA)
		Rates->Delta = RandUniDouble(Rates->RS, MIN_DELTA, MAX_DELTA);
	else
		Rates->Delta = MultePram(Rates, Rates->Delta, MIN_DELTA, MAX_DELTA, Dev);
}

void	ChangeLambda(OPTIONS* Opt, RATES* Rates, SCHEDULE* Shed)
{
	double Dev;

	Shed->CurrentAT = Shed->LambdaAT;
	Dev = Shed->CurrentAT->CDev;

	if(Dev > MAX_LAMBDA)
		Rates->Lambda = RandUniDouble(Rates->RS, MIN_LAMBDA, MAX_LAMBDA);
	else
		Rates->Lambda = MultePram(Rates, Rates->Lambda, MIN_LAMBDA, MAX_LAMBDA, Dev);
}

void	ChangeOU(OPTIONS* Opt, RATES* Rates, SCHEDULE* Shed)
{
	double Dev;

	Shed->CurrentAT = Shed->OUAT;
	Dev = Shed->CurrentAT->CDev;

	if(Dev > MAX_OU)
		Rates->OU = RandUniDouble(Rates->RS, MIN_OU, MAX_OU);
	else
		Rates->OU = MultePram(Rates, Rates->OU, MIN_OU, MAX_OU, Dev);
}

void	MutateRates(OPTIONS* Opt, RATES* Rates, SCHEDULE* Shed, long long It)
{
	Shed->Op = PickACat(Rates, Shed->OptFreq, Shed->NoOfOpts);

	Shed->CurrentAT = NULL;
	Rates->AutoAccept = FALSE;
	Rates->CalcLh = TRUE;
	
	switch(Shed->Op)
	{
		case SRATES:
			ChangeRates(Opt, Rates, Shed, It);
		break;

		case SCV:
			ChangeCVRates(Opt, Rates, Shed);
		break;

		case SKAPPA:
			ChangeKappa(Opt, Rates, Shed);
		break;

		case SDELTA:
			ChangeDelta(Opt, Rates, Shed);
		break;

		case SLABDA:
			ChangeLambda(Opt, Rates, Shed);
		break;

		case SOU:
			ChangeOU(Opt, Rates, Shed);
		break;

		case SJUMP:
			RJMove(Opt, Rates ,Shed);
		break;

		case SPPROR:
			MutatePriorsNormal(Rates, Rates->Priors, Rates->NoPriors, Opt->HPDev);
		break;

		case SESTDATA:
			MutateEstRates(Opt, Rates, Shed);
		break;

		case SSOLOTREEMOVE:
			Rates->TreeNo = RandUSLong(Rates->RS) % Opt->Trees->NoOfTrees;
		break;

		case SPPADDREMOVE:
			VarRatesAddRemove(Rates, Opt->Trees, Opt, Shed, It);
		break;

		case SPPMOVE:
			VarRatesMoveNode(Rates, Opt->Trees, Opt);
		break;
		
		case SPPCHANGESCALE:
			ChangeVarRatesScale(Rates, Opt->Trees, Opt, Shed);
		break;

		case SPPHYPERPRIOR:
			ChangeVarRatesHyperPrior(Rates, Opt);
		break;

		case SHETERO:
			MutateHetero(Rates);
		break;

		case STREEMOVE:
			Rates->TreeNo = RandUSLong(Rates->RS) % Opt->Trees->NoOfTrees;
		break;
		
		case SGAMMAMOVE:
			ChangeGammaRates(Rates, Shed);
		break;

		case SRJDUMMY:
			RJDummyMove(It, Opt, Opt->Trees, Rates);
		break;

		case SRJDUMMYMOVE:
			RJDummyMoveNode(Opt, Opt->Trees, Rates);
		break;

		case SRJDUMMYCHANGEBETA:
			RJDummyChange(Opt, Opt->Trees, Rates);
		break;

		case SFATTAILANS:
			if(Opt->Model == M_FATTAIL)
			{
			//	SliceSampleFatTail(Opt, Opt->Trees, Rates);
				AllSliceSampleFatTail(Opt, Opt->Trees, Rates);
			}
			else
			//	GeoUpDateAllAnsStates(Opt, Opt->Trees, Rates);
				GeoUpDateAnsStates(Opt, Opt->Trees, Rates);
		break;

		case SLOCALRATES:
			ChangeLocalTransform(Opt, Opt->Trees, Rates, Shed);
		break;

		case SDATADIST:
			ChangeTreeDistData(Opt, Rates);
		break;
	}
}

void	FreeRates(RATES *Rates, TREES *Trees)
{
	int MaxT, Index;

	if(Rates->Priors != NULL)
		FreePriors(Rates);
	
	if(Rates->DistDataRates != NULL)
		FreeDistDataRates(Rates->DistDataRates);

	if(Rates->FatTailRates != NULL)
		FreeFatTailRates(Rates->FatTailRates, Trees->NoOfSites);

	if(Rates->Contrast != NULL)
		FreeContrastRates(Rates);

	if(Rates->Rates != NULL)
		free(Rates->Rates);

	for(Index=0;Index<Rates->NoOfRates;Index++)
		free(Rates->RateNames[Index]);
	free(Rates->RateNames);

	if(Rates->Pis != NULL)
		free(Rates->Pis);

	if(Rates->FullRates != NULL)
		free(Rates->FullRates);

	if(Rates->Means != NULL)
		free(Rates->Means);

	if(Rates->Beta	!= NULL)
		free(Rates->Beta);

	if(Rates->MappingVect != NULL)
		free(Rates->MappingVect);

	if(Rates->GammaMults != NULL)
		free(Rates->GammaMults);

	if(Rates->EstData != NULL)
		free(Rates->EstData);

	if(Rates->EstDescData != NULL)
		free(Rates->EstDescData);

	if(Rates->VarRates != NULL)
		FreeVarRates(Rates->VarRates);
		
#ifdef BIG_LH
	mpfr_clear(Rates->HMeanSum);
#endif

	FreeRandStates(Rates->RS);

	MaxT = GetMaxThreads();
	for(Index=0;Index<MaxT;Index++)
		FreeRandStates(Rates->RSList[Index]);
	free(Rates->RSList);
	
	gsl_rng_free(Rates->RNG);
	
	if(Rates->Hetero != NULL)
		FreeHetero(Rates->Hetero);

	if(Rates->ModelFile != NULL)
		FreeModelFile(Rates->ModelFile);

	if(Rates->RJDummy != NULL)
		FreeRJDummyCode(Rates->RJDummy);

	for(Index=0;Index<Rates->NoLocalTransforms;Index++)
		FreeLocalTransforms(Rates->LocalTransforms[Index]);

	if(Rates->LocalTransforms != NULL)
		free(Rates->LocalTransforms);

	free(Rates);
}

void	UpDataSummaryNo(SUMMARYNO *SumNo, double No)
{
	SumNo->Sum		+= No;
	SumNo->SumSqrs	+= No * No;
	SumNo->N		+= 1;
}

double	GetSummaryAve(SUMMARYNO *SumNo)
{
	return (SumNo->Sum / SumNo->N);
}

double	GetSummaryVar(SUMMARYNO *SumNo)
{
	double	Ret;

	Ret = (SumNo->Sum * SumNo->Sum) / (double)SumNo->N;
	Ret = SumNo->SumSqrs - Ret;

	return Ret / (SumNo->N - 1);	
}

void	UpDataSummary(SUMMARY *Summary, RATES* Rates, OPTIONS *Opt)
{
	int		Index;
	double	Pct;
	double	*RootP;

	RootP = Opt->Trees->Tree[Rates->TreeNo]->Root->Partial[0];

	UpDataSummaryNo(&Summary->Lh, Rates->Lh);

	for(Index=0;Index<Opt->NoOfRates;Index++)
		UpDataSummaryNo(&Summary->Rates[Index], Rates->Rates[Index]);

	for(Index=0;Index<Opt->Trees->NoOfStates;Index++)
	{
		Pct = GetStateProbPct(Index, Opt->Trees->NoOfStates, RootP);
		UpDataSummaryNo(&Summary->Root[Index], Pct);
	}
}

void	InitSummaryNo(SUMMARYNO *SumNo)
{
	SumNo->N		=	0;
	SumNo->Sum		=	0;
	SumNo->SumSqrs	=	0;
}

void	FreeSummary(SUMMARY*	Summary)
{
	free(Summary->Rates);
	free(Summary->Root);
}

SUMMARY*	CreatSummary(OPTIONS *Opt)
{
	SUMMARY* Ret=NULL;
	int		 Index;

	Ret = (SUMMARY*)malloc(sizeof(SUMMARY));
	if(Ret==NULL)
		MallocErr();

	InitSummaryNo(&Ret->Lh);
	
	Ret->Rates = (SUMMARYNO*)malloc(sizeof(SUMMARYNO)*Opt->NoOfRates);
	if(Ret->Rates == NULL)
		MallocErr();
	for(Index=0;Index<Opt->NoOfRates;Index++)
		InitSummaryNo(&Ret->Rates[Index]);


	Ret->Root = (SUMMARYNO*)malloc(sizeof(SUMMARYNO)*Opt->Trees->NoOfStates);
	if(Ret->Root== NULL)
		MallocErr();
	for(Index=0;Index<Opt->Trees->NoOfStates;Index++)
		InitSummaryNo(&Ret->Root[Index]);

	return Ret;
}

void	PrintSummaryHeadder(FILE* Str, SUMMARY	*Summary, OPTIONS *Opt)
{
	int	Index;

	fprintf(Str, "Tree File\tData File\tNOS\t");

	fprintf(Str, "Lh - Ave\tLh - Var\t");


	fprintf(Str, "%s - Ave\t%s - Var\t", Opt->RateName[0], Opt->RateName[0]);
	fprintf(Str, "%s - Ave\t%s - Var\t", Opt->RateName[1], Opt->RateName[1]);

	
	for(Index=0;Index<Opt->Trees->NoOfStates;Index++)
		fprintf(Str, "%d - Ave\t%d - Var\t", Index, Index);

	fprintf(Str, "\n");
}

void	PrintSummary(FILE* Str, SUMMARY	*Summary, OPTIONS *Opt)
{
	int	Index;

	fprintf(Str, "%s\t%s\t%d\t", Opt->TreeFN, Opt->DataFN, Opt->Trees->NoOfStates);

	fprintf(Str, "%f\t%f\t", GetSummaryAve(&Summary->Lh), GetSummaryVar(&Summary->Lh));

	fprintf(Str, "%f\t%f\t", GetSummaryAve(&Summary->Rates[0]), GetSummaryVar(&Summary->Rates[0]));
	fprintf(Str, "%f\t%f\t", GetSummaryAve(&Summary->Rates[1]), GetSummaryVar(&Summary->Rates[1]));

	for(Index=0;Index<Opt->Trees->NoOfStates;Index++)
		fprintf(Str, "%f\t%f\t", GetSummaryAve(&Summary->Root[Index]), GetSummaryVar(&Summary->Root[Index]));

	fprintf(Str, "\n");
}
/*
int		GetNoExpRatesModelFile(RATES *Rates, OPTIONS *Opt)
{
	int		Ret;

	Ret = Rates->NoOfRates;

	if(Opt->EstKappa == TRUE)
		Ret++;

	if(Opt->EstDelta == TRUE)
		Ret++;

	if(Opt->EstLambda == TRUE)
		Ret++;

	return Ret;	
}

*/
