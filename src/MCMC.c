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
#include <string.h>
#include <math.h>


#include "TypeDef.h"
#include "Trees.h"
#include "Rates.h"
#include "Priors.h"
#include "Likelihood.h"
#include "GenLib.h"
#include "RandLib.h"
#include "Options.h"
#include "RevJump.h"
#include "Data.h"
#include "Gamma.h"
#include "ML.h"
#include "VarRates.h"
#include "Threaded.h"
#include "Schedule.h"
#include "ModelFile.h"
#include "Stones.h"
#include "RJDummy.h"
#include "Contrasts.h"
#include "FatTail.h"
#include "Geo.h"
#include "TransformTree.h"
#include "Prob.h"
#include "MCMC.h"
#include "VarRatesMLLocalTransfom.h"
#include "IntraNode.h"
#include "Output.h"

#include <gsl/gsl_rng.h>




void	UpDateHMean(OPTIONS *Opt, RATES *Rates)
{
#ifndef BIG_LH
	Rates->HMeanCount++;
	Rates->HMeanSum += 1/exp(Rates->Lh);
#else
	mpfr_t	t1, t2;

	Rates->HMeanCount++;

	mpfr_init2(t1, Opt->Precision);
	mpfr_init2(t2, Opt->Precision);

	mpfr_set_d(t1, Rates->Lh, DEF_ROUND);

	mpfr_exp(t2, t1, DEF_ROUND);

	mpfr_d_div(t1, 1.0, t2, DEF_ROUND);

	mpfr_add(t2, t1, Rates->HMeanSum, DEF_ROUND);

	mpfr_set(Rates->HMeanSum, t2, DEF_ROUND);

	mpfr_clears(t1, t2, NULL);
#endif
}

void	PrintPrior(FILE* Str, PRIOR *Prior)
{
	int	Index;

	for(Index=0;Index<DISTPRAMS[Prior->Dist];Index++)
		fprintf(Str, "%f\t", Prior->DistVals[Index]);
}

void	PrintMCMCSample(size_t Itters, SCHEDULE* Shed, OPTIONS *Opt, RATES *Rates, FILE* Str, TREES* Trees)
{
	int		PIndex;
	double	HMean;
	PRIOR	*Prior;
	

	HMean = GetHMean(Opt, Rates);
	fprintf(Str, "%zu\t%f\t%d\t", Itters, Rates->Lh, Rates->TreeNo+1);

	PrintRates(Str, Rates, Opt, Shed, Trees);

	for(PIndex=0;PIndex<Rates->NoPriors;PIndex++)
	{
		Prior = Rates->Priors[PIndex];
		if(Prior->UseHP == TRUE)
			PrintPrior(Str, Rates->Priors[PIndex]);
	}

	fprintf(Str, "\n");

	fflush(stdout);
}

void	PrintTest(int Itters, RATES* Rates)
{
	char	MType;
	int		Index;

	MType = RJModelType(Rates->MappingVect);

	printf("%d\t%d\t'", Itters, NoOfPramGroups(Rates, NULL, NULL));

	for(Index=0;Index<Rates->NoOfFullRates;Index++)
	{
		if(Rates->MappingVect[Index] == ZERO_RATE_NO)
			printf("Z");
		else
		{
			if(Rates->MappingVect[Index] <= 9)
				printf("%d", Rates->MappingVect[Index]);
			else
				printf("%c", Rates->MappingVect[Index] + 'A');
		}
	}
	printf("\n");
}

void	TestInitPho(OPTIONS *Opt, TREES *Trees, RATES *Rates)
{
	int Index, RIndex;

	for(Index=0;Index<1000;Index++)
	{
		for(RIndex=0;RIndex<Rates->NoOfRates;RIndex++)
		{
			Rates->Rates[RIndex] = gsl_rng_uniform_pos(Rates->RNG) * 0.01;
		}

		Rates->Lh = Likelihood(Rates, Trees, Opt);

		printf("Lh\t%d\t%f\n", Index, Rates->Lh);

		fflush(stdout);
	}

	exit(0);
}

double	GetRandValFromType(TRANSFORM_TYPE Type, RATES *Rates, gsl_rng *RNG)
{
	PRIOR *Prior;

	Prior = NULL;

	if(Type == VR_KAPPA)
		Prior = GetPriorFromName("Kappa", Rates->Priors, Rates->NoPriors);

	if(Type == VR_LAMBDA)
		Prior = GetPriorFromName("Lambda", Rates->Priors, Rates->NoPriors);

	if(Type == VR_DELTA)
		Prior = GetPriorFromName("Delta", Rates->Priors, Rates->NoPriors);

	if(Type == VR_OU)
		Prior = GetPriorFromName("OU", Rates->Priors, Rates->NoPriors);

	if(Type == VR_NODE)
		Prior = GetPriorFromName("VRNode", Rates->Priors, Rates->NoPriors);

	if(Type == VR_BL)
		Prior = GetPriorFromName("VRBL", Rates->Priors, Rates->NoPriors);

	if(Type == VR_FABRIC_BETA)
		Prior = GetPriorFromName("FabricBeta", Rates->Priors, Rates->NoPriors);

	assert(Prior != NULL);

	return RandFromPrior(RNG, Prior);
}


void	SetDefMCMCParameters(OPTIONS *Opt, TREES *Trees, RATES *Rates, 	gsl_rng *RNG)
{
	int Index;
	LOCAL_TRANSFORM *LR;
	PRIOR *Prior;
	double PriorVal;


	if(Opt->EstKappa == TRUE)
		Rates->Kappa = GetRandValFromType(VR_KAPPA, Rates, RNG);

	if(Opt->EstLambda == TRUE)
		Rates->Lambda = GetRandValFromType(VR_LAMBDA, Rates, RNG);

	if(Opt->EstDelta == TRUE)
		Rates->Delta = GetRandValFromType(VR_DELTA, Rates, RNG);

	if(Opt->EstOU == TRUE)
		Rates->OU = GetRandValFromType(VR_OU, Rates, RNG);

	if(Opt->UseGlobalTrend == TRUE)
	{
		Prior = GetPriorFromName("GlobalTrend", Rates->Priors, Rates->NoPriors);
		Rates->GlobalTrend = RandFromPrior(RNG, Prior);
//		Rates->GlobalTrend = 0.0;
	}

	if(Opt->EstGamma == TRUE)
	{
		Prior = GetPriorFromName("Gamma", Rates->Priors, Rates->NoPriors);
		Rates->Gamma = RandFromPrior(RNG, Prior);
	}

	for(Index=0;Index<Rates->NoLocalTransforms;Index++)
	{
		LR = Rates->LocalTransforms[Index];
		if(LR->Est == TRUE)
			LR->Scale = GetRandValFromType(LR->Type, Rates, RNG);
	}

	if(Rates->NoEstData > 0 && Opt->DataType == CONTINUOUS)
		SetEstDataFromPrior(Rates);

	if(Opt->UseCovarion == TRUE)
	{
		Prior = GetPriorFromName("CVSwichRate",Rates->Priors,Rates->NoPriors);
		PriorVal = RandFromPrior(RNG, Prior);
		Rates->OffToOn = Rates->OnToOff = PriorVal;
	}
}

double		ValidMCMCParameters(OPTIONS *Opt, TREES *Trees, RATES *Rates)
{

	Rates->Lh = Likelihood(Rates, Trees, Opt);

	if(Rates->Lh == ERRLH)
		return ERRLH;

	CalcPriors(Rates, Opt);

	if(Rates->LhPrior == ERRLH)
		return ERRLH;

	return Rates->Lh + Rates->LhPrior;
}

void	SetAllMCMCRates(double Val, RATES *Rates)
{
	int Index;

	for(Index=0;Index<Rates->NoOfRates;Index++)
		Rates->Rates[Index] = Val;
}

int	FindValidStartRateAllSame(OPTIONS *Opt, TREES *Trees, RATES *Rates)
{
	double CRate, BRate;
	double CLh, BLh;

	SetAllMCMCRates(1.0, Rates);
	CLh = ValidMCMCParameters(Opt, Trees, Rates);


	CRate = RATE_MAX;
	BRate = -1.0;
	BLh = ERRLH;
	do
	{
		SetAllMCMCRates(CRate, Rates);
		CLh = ValidMCMCParameters(Opt, Trees, Rates);
		if(CLh != ERRLH && CLh > BLh)
		{
			BRate = CRate;
			BLh = CLh;
		}

		CRate = CRate * 0.1;

	} while(CRate > RATE_MIN);

	if(BLh == ERRLH)
		return FALSE;

	SetAllMCMCRates(BRate, Rates);
	Rates->Lh = ValidMCMCParameters(Opt, Trees, Rates);

	return TRUE;
}


double	RandFromPriorPosition(int Pos, OPTIONS *Opt, TREES *Trees, RATES *Rates, gsl_rng *RNG)
{
	PRIOR *P;

	if(Opt->UseRJMCMC == TRUE)
		P = GetPriorFromName("RJRates", Rates->Priors, Rates->NoPriors);
	else
		P =  GetPriorFromName(Rates->RateNames[Pos], Rates->Priors, Rates->NoPriors);;

	return RandFromPrior(RNG, P);
}

void	RandRatesFromPrior(OPTIONS *Opt, TREES *Trees, RATES *Rates, gsl_rng *RNG)
{
	int TNo, RIndex;

	for(TNo=0;TNo<NO_RAND_START_TRIES;TNo++)
	{
		if(gsl_rng_uniform(RNG) < 0.1)
			FindValidStartRateAllSame(Opt, Trees, Rates);
		{
			for(RIndex=0;RIndex<Rates->NoOfRates;RIndex++)
				Rates->Rates[RIndex] = RandFromPriorPosition(RIndex, Opt, Trees, Rates, RNG);
		}

		SetDefMCMCParameters(Opt, Trees, Rates, RNG);

		if(ValidMCMCParameters(Opt, Trees, Rates) != ERRLH)
			return;
	}

	printf("Cannot find initial starting set of parameters.");
	exit(0);
}

void FindValidStartLh(OPTIONS *Opt, TREES *Trees, RATES *Rates, gsl_rng *RNG)
{
	RandRatesFromPrior(Opt, Trees, Rates, RNG);
}

void	InitMCMC(OPTIONS *Opt, TREES *Trees, RATES *Rates)
{
	Rates->TreeNo = 0;

	SetDefMCMCParameters(Opt, Trees, Rates, Rates->RNG);

	if(Opt->MCMCMLStart == TRUE)
	{
		printf("Starting MCMC form ML is not avalable, please e-mail support.\n");
		exit(0);
		MLTree(Opt, Trees, Rates);
	}

	FindValidStartLh(Opt, Trees, Rates, Rates->RNG);

	if(Opt->VarRatesCheckPoint != NULL)
		SetVarRatesFromStr(Rates, Opt, Opt->VarRatesCheckPoint, Trees);
}

void	ShowTimeSec(double StartT, double EndT)
{
	printf("Sec:\t%f\n", EndT - StartT);
}



int		ExitMCMC(OPTIONS *Opt, STONES *Stones, size_t Itters)
{
	if(Stones != NULL)
		return StoneExit(Stones, Itters);
	

	if((Opt->Itters == Itters) && (Opt->Itters != 0))
		return TRUE;

	return FALSE;
}


void	TestArea(OPTIONS *Opt, TREES *Trees, RATES *Rates)
{
//	int Index;
	double Lh, X;

	TestDummyCodeSig(Opt, Trees, Rates);
	exit(0);

	Lh = Likelihood(Rates, Trees, Opt);
	printf("Lh:\t%f\n", Lh);

	for(X=-1;X<2;X += 0.0001)
	{
		Rates->Rates[0] = X;
		Lh = Likelihood(Rates, Trees, Opt);
		printf("%f\t%f\n", X, Lh);
	}

	exit(0);
}

void	ReSetAccFlags(RATES *Rates)
{
	Rates->AutoAccept = FALSE;
	Rates->CalcLh = TRUE;
}

int		MCMCAccept(size_t Itters, OPTIONS *Opt, TREES *Trees, SCHEDULE* Shed, RATES *CRates, RATES *NRates, STONES *Stones, double Temperature)
{
	double Heat;
	
	if(NRates->AutoAccept == TRUE)
	{
		ReSetAccFlags(NRates);
		return TRUE;
	}

	ReSetAccFlags(NRates);

	if((Shed->Op == S_FAT_TAIL_ANS_ALL || Shed->Op == S_FAT_TAIL_ANS) && Opt->Model == M_FATTAIL)
		return TRUE;		

	Heat = (NRates->Lh - CRates->Lh) * Temperature;

	if(Stones != NULL)
		Heat = GetStoneHeat(Stones, Itters, Heat);

	Heat += NRates->LhPrior - CRates->LhPrior;
	Heat += NRates->LnHastings;

	if(log(gsl_rng_uniform_pos(CRates->RNG)) <= Heat)
		return TRUE;

	return FALSE;
}

void PrintStdOutHeader(OPTIONS *Opt, TREES *Trees)
{
	PrintOptions(stdout, Opt, Trees);
	printf("Iteration\tLh\tLh Prior\tElapsed Seconds\tState");
	printf("\n");
	fflush(stdout);
}

void PrintChainState(CHAIN_STATE ChainState)
{
	switch(ChainState)
	{
		case STATE_BURN_IN:
			printf("Burn In"); break;

		case STATE_SAMPLING:
			printf("Sampling"); break;

		case STATE_STEPPING_STONES:
			printf("Stepping Stone Sampler"); break;
	}
}

void PrintStdOut(OPTIONS *Opt, RATES *Rates, size_t Itters, double Seconds, CHAIN_STATE ChainState)
{
	printf("%zu\t%f\t%f\t%f\t", Itters, Rates->Lh, Rates->LhPrior, Seconds);

	PrintChainState(ChainState);

	printf("\n");
	fflush(stdout);
}

int GetNodeDepth(NODE N)
{
	int Ret;

	Ret = 0;
	while(N->Ans != NULL)
	{
		N = N->Ans;
		Ret++;
	}

	return Ret;
}

void 	MCMCTest(OPTIONS *Opt, TREES *Trees, RATES*	Rates, SCHEDULE* Shed)
{
	double		StartT, EndT, Lh;
	TREE *Tree;
	int Index;

	Tree = Trees->Tree[0];

	Lh = Likelihood(Rates, Trees, Opt);

	printf("%f\n", Lh);

	StartT = GetSeconds();

	for(Index=0;Index<2500;Index++)
	{
		GeoUpDateAllAnsStates(Opt, Trees, Rates);
		Lh = Likelihood(Rates, Trees, Opt);
	}

	EndT = GetSeconds();

	Lh = Likelihood(Rates, Trees, Opt);
	
	printf("%d\t%f\n", Opt->Cores, Lh);

	printf("Total Time:\t%f\n", EndT - StartT);

	exit(0);


	for(Index=0;Index<Tree->NoParallelGroups;Index++)
	{
		printf("%d\t%d\n", Index, Tree->ParallelGroupSize[Index]);
	}
	exit(0);



	for(Index=0;Index<1000;Index++)
	{
	//	Likelihood(Rates, Trees, Opt);
		GeoUpDateAllAnsStates(Opt, Trees, Rates);
	}

	printf("Total Time:\t%f\n", EndT - StartT);

	exit(0);
}

CHAIN_STATE	SetChainState(OPTIONS *Opt, size_t Itters, int EqualTreeBurntIn)
{

	if(Opt->UseEqualTrees == TRUE)
	{
		if(EqualTreeBurntIn == FALSE)
			return STATE_BURN_IN;

		if(Itters < Opt->EqualTreesBI)
			return STATE_BURN_IN;

		if(Itters > Opt->Itters)
			return STATE_STEPPING_STONES;

		return STATE_SAMPLING;

	}

	if(Itters < Opt->BurnIn)
		return STATE_BURN_IN;

	if(Itters > Opt->Itters && Opt->Itters != -1)
		return STATE_STEPPING_STONES;

	return STATE_SAMPLING;
}

void BugTest(OPTIONS *Opt, TREES *Trees, RATES *Rates)
{
	double X;

	Rates->Rates[0] = 0.604366105;
	Rates->Rates[1] = 0.273410018;

	printf("Lh:\t%f\n", Likelihood(Rates, Trees, Opt));

	for(X=-5;X<10;X+=0.001)
	{
		Rates->SitePowers->Powers[0] = X;
		printf("%f\t%f\n", X, Likelihood(Rates, Trees, Opt));
	}

	exit(0);
}

void	MCMC(OPTIONS *Opt, TREES *Trees)	
{
	RATES		*CRates;
	RATES		*NRates;
	size_t		Itters;
	double		StartT;
	SCHEDULE	*Shed;
	int			EqualTreeBurntIn;
	CHAIN_STATE	ChainState;
	STONES		*Stones;

	CRates	=	NULL;
	NRates	=	NULL;
	Shed = NULL;
	Stones = NULL;

	EqualTreeBurntIn = FALSE;

	Itters = 0;

	if(Opt->LoadCheckPointFile == TRUE)
	{
		Itters++;
	}
	else
	{
		CRates	=	CreatRates(Opt, Trees, Opt->Seed);
		NRates	=	CreatRates(Opt, Trees, Opt->Seed);
		Shed = CreatSchedule(Opt, CRates->RNG, Trees);
		InitMCMC(Opt, Trees, CRates);
	}

#ifdef SIM_GEO_DATA
	SimGeoData(Opt, Trees, CRates);
#endif

	if(Opt->StoneOptions != NULL)
	{
		Stones = CratesStones(Opt->StoneOptions->NoStones, Opt->StoneOptions->ItPerStone, Opt->StoneOptions->Alpha, Opt->StoneOptions->Beta);
		Stones->ItStart = Opt->Itters + 1;
	}

	SetOutputFile(Opt, Trees, CRates, Shed, Stones);

	CRates->Lh	=	Likelihood(CRates, Trees, Opt);
	CalcPriors(CRates, Opt);

	

	StartT = GetSeconds();

	for(;;Itters++)
	{
		ChainState = SetChainState(Opt, Itters, EqualTreeBurntIn);

		SetCustomSchedule(Opt, Opt->ShedFile, Itters, Shed);

 		CopyRates(NRates, CRates, Opt, Trees);

		MutateRates(Opt, NRates, Trees, Shed, Itters);

		if(Opt->NodeData == TRUE)
			SetTreeAsData(Opt, Trees, NRates->TreeNo);

		NRates->Lh = Likelihood(NRates, Trees, Opt);

		if(NRates->Lh == ERRLH)
			Itters--;
		else
		{
			CalcPriors(NRates, Opt);

			if(MCMCAccept(Itters, Opt, Trees, Shed, CRates, NRates, Stones, 1.0) == TRUE)
			{
				Swap((void**)&NRates, (void**)&CRates);
				UpDateShedAcc(TRUE, Shed);
			}
			else
				UpDateShedAcc(FALSE, Shed);

			if(Itters % Opt->Sample == 0)
			{
				// Needed to make sure the tree has correct values set
				// Also a sandity check.
				if(CRates->Lh != Likelihood(CRates, Trees, Opt))
				{
					printf("Itteration %zu: Likelihood Error: %d::%s\n", Itters, __LINE__, __FILE__);
					exit(1);
				}
				PrintStdOut(Opt, CRates, Itters, GetSeconds() - StartT, ChainState);
			}

			if(Itters == Opt->BurnIn)
			{
				if(	Opt->UseEqualTrees == TRUE && 
					EqualTreeBurntIn == FALSE)
				{
					EqualTreeBurntIn = TRUE;
					Itters = -1;
				}
			}


			if(Itters % Opt->Sample == 0 && ChainState == STATE_SAMPLING)
			{
				UpDateHMean(Opt, CRates);

				PrintMCMCSample(Itters, Shed, Opt, CRates, Opt->LogFile, Trees);
				fflush(Opt->LogFile);

				if(Opt->UseSchedule == TRUE)
					PrintShed(Opt, Shed, Opt->ShedFile);

				if(UseNonParametricMethods(Opt) == TRUE)
					PrintVarRatesOutput(Opt, Trees, CRates, Itters);

				if(Opt->ModelType == MT_FATTAIL)
				{
					OutputFatTail(Itters, Opt, Trees, CRates);
					if(Opt->UseIntraNode == TRUE)
						OutputIntraNode(Itters, Opt, Trees, CRates, Opt->LogIntraNode);
				}

				if(Opt->RJDummy == TRUE)
					PrintRJDummy(Itters, Opt, Trees, CRates);

				if(Opt->SaveModels == TRUE)
					SaveModelFile(Opt->SaveModelFile, Opt, Trees, CRates);

				if(Opt->SaveTrees == TRUE)
					OutputTree(Opt, Trees, CRates, Itters, Opt->OutTrees);
			}

			if(Itters % MCMC_SCHEDULE_UPDATE == 0)
			{
				// The schedule should be updated even when stones are running
				UpDateSchedule(Opt, Shed, CRates->RNG);
				BlankSchedule(Shed);
			}

			if(ExitMCMC(Opt, Stones, Itters) == TRUE)
			{
				if( (Opt->UseEqualTrees == FALSE) ||
					(CRates->TreeNo == Trees->NoTrees - 1))
				{

					FreeRates(CRates, Trees);
					FreeRates(NRates, Trees);

					FreeeSchedule(Shed);

					if(Stones != NULL)
						FreeStones(Stones);

					if(Opt->StoneFile != NULL)
						fclose(Opt->StoneFile);

					if(Opt->ShedFile != NULL)
						fclose(Opt->ShedFile);

					if(Opt->VarRatesLog != NULL)
						fclose(Opt->VarRatesLog);

					if(Opt->SaveModelFile != NULL)
						fclose(Opt->SaveModelFile);

					return;
				}

				if(EqualTreeBurntIn == TRUE)
				{
					CRates->TreeNo++;
					CRates->Lh = Likelihood(CRates, Trees, Opt);
					Itters = 0;
					BlankSchedule(Shed);
					Shed->GNoAcc = Shed->GNoTried = 0;
				}
			}

			if(Stones != NULL)
				StoneItter(Stones, Itters, CRates->Lh, Opt->StoneFile);
			
			if(Opt->Itters == Itters && Opt->Itters != 0 && Opt->RJLockModel == TRUE)
				SetRJLockedModel(Shed);
		}
	}
}


