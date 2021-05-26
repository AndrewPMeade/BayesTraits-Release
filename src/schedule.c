/*
*  BayesTriats 3.0
*
*  copyright 2017
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

#include "TypeDef.h"
#include "GenLib.h"
#include "Schedule.h"
#include "Data.h"
#include "AutoTune.h"
#include "Rates.h"
#include "VarRates.h"
#include "LocalTransform.h"
#include "TimeSlices.h"

void	AddToFullATList(SCHEDULE* Shed, AUTOTUNE *AT)
{
	AUTOTUNE **NList;

	NList = (AUTOTUNE**)SMalloc(sizeof(AUTOTUNE*) * (Shed->NoFullATList + 1));

	if(Shed->FullATList != NULL)
	{
		memcpy(NList, Shed->FullATList, sizeof(AUTOTUNE*) * Shed->NoFullATList);
		free(Shed->FullATList);
	}

	NList[Shed->NoFullATList] = AT;

	Shed->FullATList = NList;

	Shed->NoFullATList++;
}

void	UpDateShedAcc(int Acc, SCHEDULE* Shed)
{
	Shed->Tryed[Shed->Op]++;

	if(Shed->CurrentAT != NULL)
		Shed->CurrentAT->NoTried++;

	Shed->GNoTried++;
	Shed->SNoTried++;

	if(Acc == FALSE)
		return;



	Shed->GNoAcc++;
	Shed->SNoAcc++;
	Shed->Accepted[Shed->Op]++;

	if(Shed->CurrentAT != NULL)
		Shed->CurrentAT->NoAcc++;
}

int		UsingHP(OPTIONS *Opt)
{
	int Index;

	for(Index=0;Index<Opt->NoAllPriors;Index++)
	{
		if(Opt->AllPriors[Index]->UseHP == TRUE)
			return TRUE;
	}

	return FALSE;
}

void	BlankSchedule(SCHEDULE*	Shed)
{
	int	Index;

	Shed->SNoAcc = Shed->SNoTried = 0;

	for(Index=0;Index<Shed->NoOfOpts;Index++)
	{
		Shed->Accepted[Index] = 0;
		Shed->Tryed[Index] = 0;
	}
}


int		MultiTree(OPTIONS *Opt)
{
	if(Opt->UseEqualTrees == TRUE)
		return FALSE;

	if(Opt->Trees->NoTrees == 1)
		return FALSE;

	return TRUE;
}

void	SetVarRatesShed(OPTIONS *Opt, SCHEDULE *Shed)
{
	int Index, Max, No;

	Max = NO_RJ_LOCAL_SCALAR + 2;

	Shed->FreqVarRatesOp	= (double*)SMalloc(sizeof(double) * Max);
	Shed->VarRatesOp		= (TRANSFORM_TYPE*)SMalloc(sizeof(TRANSFORM_TYPE) * Max);

	No = 0;

	for(Index=0;Index<NO_RJ_LOCAL_SCALAR;Index++)
	{
		if(Opt->UseRJLocalScalar[Index] == TRUE)
		{
			Shed->VarRatesOp[No] = (TRANSFORM_TYPE)Index;
			Shed->FreqVarRatesOp[No] = 0.1;
			No++;
		}
	}

	Shed->NoVarRatesOp = No;

	NormaliseVector(Shed->FreqVarRatesOp, Shed->NoVarRatesOp);
}

void	SetCustomSchdule(SCHEDULE *Shed, OPTIONS *Opt)
{
	return;
}

void	SetSchedule(SCHEDULE *Shed, OPTIONS *Opt)
{
	int		Rates, Index;
	
	for(Index=0;Index<Shed->NoOfOpts;Index++)
		Shed->OptFreq[Index] = 0.0;

	if(Opt->UseCovarion == TRUE && Opt->LoadModels == FALSE)
		Shed->OptFreq[S_CV] = 0.2;

	if(Opt->EstKappa == TRUE && Opt->LoadModels == FALSE)
		Shed->OptFreq[S_KAPPA] = 0.1;

	if(Opt->EstDelta == TRUE && Opt->LoadModels == FALSE)
		Shed->OptFreq[S_DELTA] = 0.1;

	if(Opt->EstLambda == TRUE  && Opt->LoadModels == FALSE)
		Shed->OptFreq[S_LABDA] = 0.1;

	if(Opt->UseRJMCMC == TRUE && Opt->LoadModels == FALSE)
		Shed->OptFreq[S_JUMP] = 0.1;

	if(UsingHP(Opt) == TRUE)
		Shed->OptFreq[S_PPROR] = 0.1;

	if(EstData(Opt->Trees) == TRUE)
		Shed->OptFreq[S_EST_DATA] = 0.5;

	if(Opt->EstOU == TRUE && Opt->LoadModels == FALSE)
		Shed->OptFreq[S_OU] = 0.1;

	if(Opt->EstGamma == TRUE)
		Shed->OptFreq[S_GAMMA_MOVE] = 0.1;

	if(Opt->Model == M_GEO)
	{
		Shed->OptFreq[S_GEO_MOVE_ALL] = 0.9;
	}

	if(Opt->Model == M_FATTAIL)
		Shed->OptFreq[S_FAT_TAIL_ANS_ALL] = 0.8;

	if(MultiTree(Opt) == TRUE)
		Shed->OptFreq[S_TREE_MOVE] = 0.1;
	
	if(EstLocalTransforms(Opt->LocalTransforms, Opt->NoLocalTransforms) == TRUE && Opt->LoadModels == FALSE)
		Shed->OptFreq[S_LOCAL_RATES] = 0.1;

	if(UseRJLocalScalar(Opt) == TRUE)
	{
		Shed->OptFreq[S_VARRATES_ADD_REMOVE] = 0.5;
		Shed->OptFreq[S_VARRATES_MOVE] = 0.05;
		Shed->OptFreq[S_VARRATES_CHANGE_SCALE] = 0.4;
		
		SetVarRatesShed(Opt, Shed);
	}

	if(Opt->Model == M_DISC_HET)
		Shed->OptFreq[S_HETERO] = 0.4;

	Rates = 0;
	if(Opt->DataType == CONTINUOUS)
		Rates = Opt->Trees->NoSites;
	else
		for(Index=0;Index<Opt->NoOfRates;Index++)
			if(Opt->ResTypes[Index] == RESNONE)
				Rates++;

	if(Rates == 0)
		Shed->OptFreq[S_RATES] = 0;
	else
	{
		Shed->OptFreq[S_RATES] = 0.5;
		if(Opt->ModelType == MT_FATTAIL)
			Shed->OptFreq[S_RATES] = 0.1;
	}
	
#ifdef CONTRAST_ML_PARAM
	if(Opt->ModelType == MT_CONTRAST)
		Shed->OptFreq[0] = 0;
#endif

	if(Opt->LoadModels == TRUE)
		Shed->OptFreq[S_RATES] = 0.4;
	
	if(Opt->RJDummy == TRUE)
	{
		Shed->OptFreq[S_RATES] = 0.5;
		Shed->OptFreq[S_RJ_DUMMY] = 0.2;
		Shed->OptFreq[S_RJ_DUMMY_MOVE] = 0.2;
		Shed->OptFreq[S_RJ_DUMMY_CHANG_EBETA] = 0.2;
	}

	if(Opt->UseDistData == TRUE)
		Shed->OptFreq[S_DATA_DIST] = 0.2;

	if(TimeSliceEstTime(Opt->TimeSlices) == TRUE)
		Shed->OptFreq[S_TIME_SLICE_TIME] = 0.1;
	
	if(TimeSliceEstScale(Opt->TimeSlices) == TRUE)
		Shed->OptFreq[S_TIME_SLICE_SCALE] = 0.1;

	if(Opt->NormQMat == TRUE)
		Shed->OptFreq[S_GLOBAL_RATE] = 0.1;

	if(Opt->UseGlobalTrend == TRUE)
		Shed->OptFreq[S_GLOBAL_TREND] = 0.1;
	   
	NormaliseVector(Shed->OptFreq, Shed->NoOfOpts);

	memcpy(Shed->DefShed, Shed->OptFreq, sizeof(double) * Shed->NoOfOpts);
}

void	PrintATHeader(FILE *Str,AUTOTUNE *AT)
{
	fprintf(Str,"%s - Dev\t",AT->Name);
	fprintf(Str,"%s - Tried\t",AT->Name);
	fprintf(Str,"%s - Accepted\t",AT->Name);
}

void	PrintAutoTuneHeader(FILE* Str,SCHEDULE* Shed)
{
	int Index;

	for(Index=0;Index<Shed->NoFullATList;Index++)
		PrintATHeader(Str,Shed->FullATList[Index]);
}

void	PrintCustomShedHeadder(OPTIONS* Opt, SCHEDULE* Shed, FILE* Str)
{
	int ShedNo, Index;
	CUSTOM_SCHEDULE *CSched;

	for(ShedNo=0;ShedNo<Shed->NoCShed;ShedNo++)
	{
		CSched = Shed->CShedList[ShedNo];

		fprintf(Str, "\nCustom schedule %d\tStarting\t%lld\n", ShedNo, CSched->Iteration);
	//	fprintf(Str, "It:\t%lld\t", CSched->Iteration);

		for(Index=0;Index<Shed->NoOfOpts;Index++)
		{
			if(CSched->Frequencies[Index] != 0)
				fprintf(Str, "%s\t%2.2f\n", SHEDOP[Index], CSched->Frequencies[Index]*100);
		}

		fprintf(Str, "\n\n");
	}
}

void	PrintShedHeadder(OPTIONS* Opt, SCHEDULE* Shed, FILE* Str)
{
	int	Index;
	
	fprintf(Str, "Default schedule\n");

	for(Index=0;Index<Shed->NoOfOpts;Index++)
	{
		if(Shed->OptFreq[Index] != 0)
			fprintf(Str, "%s\t%2.2f\n", SHEDOP[Index], Shed->OptFreq[Index]*100);
	}
		

	PrintCustomShedHeadder(Opt, Shed, Str);

	for(Index=0;Index<Shed->NoOfOpts;Index++)
	{
		if(Shed->OptFreq[Index] != 0)
			fprintf(Str, "%s Tried\t%% Accepted\t", SHEDOP[Index]);
	}

	PrintAutoTuneHeader(Str, Shed);

	fprintf(Str, "Sample Ave Acceptance\tTotal Ave Acceptance\t");

	fprintf(Str, "\n");	
	fflush(Str);
}

void	PrintShed(OPTIONS* Opt, SCHEDULE* Shed, FILE* Str)
{
	int	Index;
	double Last;
	
	Last = 0;
	for(Index=0;Index<Shed->NoOfOpts;Index++)
	{
		if(Shed->OptFreq[Index] != 0)
		{
			if(Shed->Tryed[Index] != 0)
				fprintf(Str, "%d\t%f\t",Shed->Tryed[Index], (double)Shed->Accepted[Index] / (double)Shed->Tryed[Index] );
			else
				fprintf(Str, "0.0\t0\t");
		}
		Last = Shed->OptFreq[Index];
	}
	
	PrintAutoTune(Str, Opt, Shed);

	fprintf(Str, "%f\t", (double)Shed->SNoAcc / Shed->SNoTried);
	fprintf(Str, "%f\t", (double)Shed->GNoAcc / Shed->GNoTried);
	fprintf(Str, "\n");

	fflush(Str);
}

SCHEDULE*	AllocSchedule()
{
	SCHEDULE*	Ret;
	
	Ret = (SCHEDULE*)SMalloc(sizeof(SCHEDULE));
	
	Ret->NoOfOpts = NO_SCHEDULE_OPT;

	Ret->OptFreq	=	(double*)SMalloc(sizeof(double) * Ret->NoOfOpts);
	Ret->Tryed		=	(int*)SMalloc(sizeof(int) * Ret->NoOfOpts);
	Ret->Accepted	=	(int*)SMalloc(sizeof(int) * Ret->NoOfOpts);
	Ret->DefShed	=	(double*)SMalloc(sizeof(double) * Ret->NoOfOpts);
	
	Ret->CShedList	=	NULL;
	Ret->NoCShed	=	0;

	Ret->RateDevATList	= NULL;
	Ret->DataDevAT		= NULL;
	Ret->VarRateAT		= NULL;

	Ret->NoParm			= -1;
	Ret->PTried			= NULL;
	Ret->PAcc			= NULL;

	Ret->KappaAT		= NULL;
	Ret->DeltaAT		= NULL;
	Ret->LambdaAT		= NULL;
	Ret->OUAT			= NULL;

	Ret->GammaAT		= NULL;

	Ret->RJDummyBetaAT	= NULL;

	Ret->LocalRatesAT	= NULL;
	
	Ret->NoVarRatesOp	= 0;
	Ret->FreqVarRatesOp = NULL;
	Ret->VarRatesOp		= NULL;

	Ret->NoFullATList	= 0;
	Ret->FullATList		= 0;

	Ret->TimeSliceTimeAT = NULL;
	Ret->TimeSliceScaleAT= NULL;

	Ret->GlobalRateAT	= NULL;
	Ret->GlobalTrendAT	= NULL;

	Ret->LandscapeRateChangeAT	 = NULL;

	Ret->StochasticBeta	=	NULL;
	Ret->StochasticBetaPrior =	NULL;
	
	return Ret;
}

int		FindNoOfAutoCalibRates(OPTIONS *Opt)
{
	if(Opt->ModelType == MT_DISCRETE)
		return 1;

	return FindNoConRates(Opt);
}

char**	GetAutoParamNames(OPTIONS *Opt)
{
	char	**Ret,*Buffer;
	int		NoP,PIndex,Index,NoS;
	
	PIndex = 0;

	NoP = FindNoOfAutoCalibRates(Opt);
	Buffer = (char*)SMalloc(sizeof(char) * BUFFERSIZE);
	Ret = (char**)SMalloc(sizeof(char*) * NoP);

	if(Opt->DataType == DISCRETE)
	{
		sprintf(Buffer,"Rates");
		Ret[PIndex++] = StrMake(Buffer);
		free(Buffer);
		return Ret;
	}

	if(Opt->ModelType == MT_FATTAIL)
	{
		NoS = Opt->Trees->NoSites;
		if(Opt->Model == M_GEO)
		{
			sprintf(Buffer,"Alpha");
			Ret[0] = StrMake(Buffer);

			sprintf(Buffer,"Scale");
			Ret[1] = StrMake(Buffer);
			free(Buffer);
			return Ret;
		}
			

		for(Index=0;Index<NoS;Index++)
		{
/*			sprintf(Buffer,"Alpha %d",Index+1);
			Ret[PIndex++] = StrMake(Buffer);

			sprintf(Buffer,"Scale %d",Index+1);
			Ret[PIndex++] = StrMake(Buffer);*/

			sprintf(Buffer,"Sig2 %d",Index+1);
			Ret[PIndex++] = StrMake(Buffer);
		}

		free(Buffer);
		return Ret;
	}

	if(Opt->Model == M_CONTRAST_CORREL)
	{
		for(Index=0;Index<Opt->Trees->NoSites;Index++)
		{
			sprintf(Buffer,"Alpha %d",Index+1);
			Ret[PIndex++] = StrMake(Buffer);
		}

		free(Buffer);
		return Ret;
	}

	if(Opt->Model == M_CONTRAST_REG)
	{
		for(Index=1;Index<Opt->Trees->NoSites;Index++)
		{
			sprintf(Buffer,"Beta %d", Index);
			Ret[PIndex++] = StrMake(Buffer);
		}

		free(Buffer);
		return Ret;
	}

	if(Opt->Model == M_CONTRAST)
	{
		for(Index=0;Index<Opt->Trees->NoSites;Index++)
		{
			sprintf(Buffer,"Alpha %d",Index+1);
			Ret[PIndex++] = StrMake(Buffer);
		}

		for(Index=0;Index<Opt->Trees->NoSites;Index++)
		{
			sprintf(Buffer,"Sigma^2 %d",Index+1);
			Ret[PIndex++] = StrMake(Buffer);
		}

		free(Buffer);
		return Ret;
	}

	if(Opt->Model == M_CONTINUOUS_REG)
	{
		sprintf(Buffer,"Alpha");
		Ret[PIndex++] = StrMake(Buffer);

		for(Index=1;Index<Opt->Trees->NoSites+1;Index++)
		{
			sprintf(Buffer,"Beta Trait %d",Index);
			Ret[PIndex++] = StrMake(Buffer);
		}
		free(Buffer);
		return Ret;
	}

	for(Index=0;Index<Opt->Trees->NoSites;Index++)
	{
		sprintf(Buffer,"Alpha Trait %d",Index+1);
		Ret[PIndex++] = StrMake(Buffer);
	}

	if(Opt->Model == M_CONTINUOUS_DIR)
	{
		for(Index=0;Index<Opt->Trees->NoSites;Index++)
		{
			sprintf(Buffer,"Beta Trait %d",Index+1);
			Ret[PIndex++] = StrMake(Buffer);
		}
	}


	free(Buffer);

	return Ret;
}

void	SetRateDevPerParm(SCHEDULE* Shed, OPTIONS *Opt, RANDSTATES *RS)
{
	int Index;
	char **PNames;

	Shed->NoParm			= FindNoOfAutoCalibRates(Opt);
	
	Shed->RateDevATList		= (AUTOTUNE**)SMalloc(sizeof(AUTOTUNE*) * Shed->NoParm);

	PNames = GetAutoParamNames(Opt);
	
	for(Index=0;Index<Shed->NoParm;Index++)
	{
		Shed->RateDevATList[Index] = CreatAutoTune(PNames[Index], RandDouble(RS) * 10, MIN_VALID_ACC, MAX_VALID_ACC);
		AddToFullATList(Shed, Shed->RateDevATList[Index]);

		if(Opt->ModelType == MT_DISCRETE)
			SetMaxDev(Shed->RateDevATList[Index], 5.0);
	}

	for(Index=0;Index<Shed->NoParm;Index++)
		free(PNames[Index]);
	free(PNames);
}

CUSTOM_SCHEDULE* CloneCustomSchedule(CUSTOM_SCHEDULE* CShed)
{
	CUSTOM_SCHEDULE* Ret;

	Ret = AllocCustomSchedule();

	Ret->Default = CShed->Default;
	Ret->Iteration = CShed->Iteration;

	memcpy(Ret->Frequencies, CShed->Frequencies, sizeof(double) * NO_SCHEDULE_OPT);
	return Ret;

}

CUSTOM_SCHEDULE**		CloneCustomScheduleList(int NoCShed, CUSTOM_SCHEDULE **CShedList)
{
	int Index;
	CUSTOM_SCHEDULE**	Ret;

	if(NoCShed == 0)
		return NULL;

	Ret = (CUSTOM_SCHEDULE**)SMalloc(sizeof(CUSTOM_SCHEDULE*) * NoCShed);

	for(Index=0;Index<NoCShed;Index++)
		Ret[Index] = CloneCustomSchedule(CShedList[Index]);

	return Ret;
}

SCHEDULE*	CreatSchedule(OPTIONS *Opt, RANDSTATES *RS)
{
	SCHEDULE*	Ret;

	Ret = AllocSchedule();
	
	BlankSchedule(Ret);

	SetSchedule(Ret, Opt);

	Ret->GNoAcc = Ret->GNoTried = 0;
	Ret->SNoAcc = Ret->SNoTried = 0;

	SetRateDevPerParm(Ret, Opt, RS);

	Ret->NoCShed = Opt->NoCShed;
	Ret->CShedList = CloneCustomScheduleList(Opt->NoCShed, Opt->CShedList);

	// Set Auto tune Data Dev
	if(Opt->EstData == TRUE)
	{
		Ret->DataDevAT = CreatAutoTune("Data", RandDouble(RS) * 10, MIN_VALID_ACC, MAX_VALID_ACC);
		AddToFullATList(Ret,Ret->DataDevAT);
	}

	// Set VarRates Auto Tune
	if(UseNonParametricMethods(Opt) == TRUE)
	{
		Ret->VarRateAT = CreatAutoTune("VarRates", RandDouble(RS), MIN_VALID_ACC, MAX_VALID_ACC);
		SetMaxDev(Ret->VarRateAT, 200.0);
		AddToFullATList(Ret,Ret->VarRateAT);
	}

	if(Opt->EstKappa == TRUE)
	{
		Ret->KappaAT = CreatAutoTune("Kappa", RandDouble(RS) * 10, MIN_VALID_ACC, MAX_VALID_ACC);
		SetMaxDev(Ret->KappaAT, 10.0);
		AddToFullATList(Ret, Ret->KappaAT);
	}

	if(Opt->EstLambda == TRUE)
	{
		Ret->LambdaAT = CreatAutoTune("Lambda", RandDouble(RS) * 10, MIN_VALID_ACC, MAX_VALID_ACC);
		SetMaxDev(Ret->LambdaAT, 10.0);
		AddToFullATList(Ret, Ret->LambdaAT);
	}

	if(Opt->EstDelta == TRUE)
	{
		Ret->DeltaAT = CreatAutoTune("Delta", RandDouble(RS) * 10, MIN_VALID_ACC, MAX_VALID_ACC);
		SetMaxDev(Ret->DeltaAT, 10.0);
		AddToFullATList(Ret, Ret->DeltaAT);
	}

	if(Opt->EstOU == TRUE)
	{
		Ret->OUAT = CreatAutoTune("OU", RandDouble(RS) * 10, MIN_VALID_ACC, MAX_VALID_ACC);
		SetMaxDev(Ret->OUAT, 10.0);
		AddToFullATList(Ret, Ret->OUAT);
	}

	if(Opt->RJDummy == TRUE)
	{
		Ret->RJDummyBetaAT = CreatAutoTune("Dummy", RandDouble(RS) * 10, MIN_VALID_ACC, MAX_VALID_ACC);
		AddToFullATList(Ret, Ret->RJDummyBetaAT);
	}
	
	if(Opt->EstGamma == TRUE)
	{
		Ret->GammaAT = CreatAutoTune("Gamma", RandDouble(RS) * 10, MIN_VALID_ACC, MAX_VALID_ACC);
		SetMaxDev(Ret->GammaAT, 10.0);
		AddToFullATList(Ret, Ret->GammaAT);
	}

	if(EstLocalTransforms(Opt->LocalTransforms, Opt->NoLocalTransforms) == TRUE)
	{
		Ret->LocalRatesAT = CreatAutoTune("LocalTransform", RandDouble(RS), MIN_VALID_ACC, MAX_VALID_ACC);
		SetMaxDev(Ret->LocalRatesAT, 10.0);
		AddToFullATList(Ret, Ret->LocalRatesAT);
	}

	if(TimeSliceEstTime(Opt->TimeSlices) == TRUE)
	{
		Ret->TimeSliceTimeAT = CreatAutoTune("Time Slice Time", RandDouble(RS), MIN_VALID_ACC, MAX_VALID_ACC);
		SetMaxDev(Ret->TimeSliceTimeAT, 1.0);
		AddToFullATList(Ret, Ret->TimeSliceTimeAT);
	}

	if(TimeSliceEstScale(Opt->TimeSlices) == TRUE)
	{
		Ret->TimeSliceScaleAT = CreatAutoTune("Time Slice Scale", RandDouble(RS), MIN_VALID_ACC, MAX_VALID_ACC);
		SetMaxDev(Ret->TimeSliceScaleAT, 10.0);
		AddToFullATList(Ret, Ret->TimeSliceScaleAT);
	}

	if(Opt->UseGlobalTrend == TRUE)
	{
		Ret->GlobalTrendAT = CreatAutoTune("Global Trend", RandDouble(RS), MIN_VALID_ACC, MAX_VALID_ACC);
		SetMaxDev(Ret->GlobalTrendAT, 100.0);
		AddToFullATList(Ret, Ret->GlobalTrendAT);
	}

	if(Opt->NormQMat == TRUE)
	{
		Ret->GlobalRateAT = CreatAutoTune("Global Rate", RandDouble(RS), MIN_VALID_ACC, MAX_VALID_ACC);
		SetMaxDev(Ret->GlobalRateAT, 10000.0);
		AddToFullATList(Ret, Ret->GlobalRateAT);
	}

	return Ret;
}


void		FreeeSchedule(SCHEDULE* Sched)
{
	int Index;

	free(Sched->OptFreq);
	free(Sched->Tryed);
	free(Sched->Accepted);

	if(Sched->RateDevATList != NULL)
	{
		for(Index=0;Index<Sched->NoParm;Index++)
			FreeAutoTune(Sched->RateDevATList[Index]);
		free(Sched->RateDevATList);
		free(Sched->PAcc);
		free(Sched->PTried);
	}

	if(Sched->DataDevAT != NULL)
		FreeAutoTune(Sched->DataDevAT);

	if(Sched->VarRateAT != NULL)
		FreeAutoTune(Sched->VarRateAT);

	if(Sched->KappaAT != NULL)
		FreeAutoTune(Sched->KappaAT);

	if(Sched->LambdaAT != NULL)
		FreeAutoTune(Sched->LambdaAT);

	if(Sched->DeltaAT != NULL)
		FreeAutoTune(Sched->DeltaAT);

	if(Sched->OUAT != NULL)
		FreeAutoTune(Sched->OUAT);

	if(Sched->LocalRatesAT != NULL)
		FreeAutoTune(Sched->LocalRatesAT);

	if(Sched->GammaAT != NULL)
		FreeAutoTune(Sched->GammaAT);

	if(Sched->RJDummyBetaAT != NULL)
		FreeAutoTune(Sched->RJDummyBetaAT);
	
	if(Sched->VarRatesOp != NULL)
		free(Sched->VarRatesOp);

	if(Sched->FreqVarRatesOp != NULL)
		free(Sched->FreqVarRatesOp);

	if(Sched->FullATList != NULL)
		free(Sched->FullATList);

	if(Sched->TimeSliceScaleAT != NULL)
		FreeAutoTune(Sched->TimeSliceScaleAT);

	if(Sched->TimeSliceTimeAT != NULL)
		FreeAutoTune(Sched->TimeSliceTimeAT);

	if(Sched->GlobalRateAT != NULL)
		FreeAutoTune(Sched->GlobalRateAT);

	if(Sched->LandscapeRateChangeAT != NULL)
		FreeAutoTune(Sched->LandscapeRateChangeAT);

	if(Sched->GlobalTrendAT != NULL)
		FreeAutoTune(Sched->GlobalTrendAT);

	if(Sched->StochasticBeta != NULL)
		FreeAutoTune(Sched->StochasticBeta);

	if(Sched->StochasticBetaPrior != NULL)
		FreeAutoTune(Sched->StochasticBetaPrior);

	if(Sched->NoCShed > 0)
	{
		for(Index=0;Index<Sched->NoCShed;Index++)
			FreeCustomSchedule(Sched->CShedList[Index]);

		free(Sched->CShedList);
	}



	free(Sched->DefShed);

	free(Sched);
}

double	GetAccRate(int Op, SCHEDULE* Shed)
{
	int	Acc, Tried;

	Tried = Shed->Tryed[Op];
	Acc = Shed->Accepted[Op]; 

	return Acc / (double)Tried;
}

double	GetRDDecAccRate(OPTIONS *Opt, SCHEDULE* Shed)
{
	double Ret;
	int	Acc, Tried;

	if(Opt->UseCovarion == FALSE)
		Ret = (double)Shed->PAcc[0] / Shed->PTried[0];
	else
	{
		Acc = Shed->PAcc[0] + Shed->Accepted[S_CV];
		Tried = Shed->PTried[0] + Shed->Tryed[S_CV];
		Ret = (double)Acc/Tried;
	}

	return Ret;
}

void	UpDateSchedule(OPTIONS *Opt, SCHEDULE* Shed, RANDSTATES *RS)
{
	int		Index;
	
	for(Index=0;Index<Shed->NoFullATList;Index++)
		AutoTuneUpDate(Shed->FullATList[Index], RS);
}

void		SetCustomShed(SCHEDULE* Shed)
{
	int Index;

	for(Index=0;Index<Shed->NoOfOpts;Index++)
		Shed->OptFreq[Index] = 0.0;

	Shed->OptFreq[S_VARRATES_MOVE] = 1.0;

	NormaliseVector(Shed->OptFreq, Shed->NoOfOpts);
}


void	SetShedOpFreq(SCHEDULE*	Shed, int No, double Val)
{
	Shed->OptFreq[No] = Val;
	NormaliseVector(Shed->OptFreq, Shed->NoOfOpts);
}

CUSTOM_SCHEDULE*	AllocCustomSchedule(void)
{
	CUSTOM_SCHEDULE* Ret;
	int Index;
	
	Ret = (CUSTOM_SCHEDULE*)SMalloc(sizeof(CUSTOM_SCHEDULE));

	Ret->Default = FALSE;
	Ret->Iteration = -1;

	Ret->Frequencies = (double*)SMalloc(sizeof(double) * NO_SCHEDULE_OPT);

	for(Index=0;Index<NO_SCHEDULE_OPT;Index++)
		Ret->Frequencies[Index] = 0.0;

	return Ret;
}


void	FreeCustomSchedule(CUSTOM_SCHEDULE*	CShed)
{
	free(CShed->Frequencies);
	free(CShed);
}

void PrintCustomSchedule(FILE *Str, int NoCShed, CUSTOM_SCHEDULE **ShedList)
{
	int Index, FIndex;
	CUSTOM_SCHEDULE *Shed;

	for(Index=0;Index<NoCShed;Index++)
	{
		Shed = ShedList[Index];

		fprintf(Str, "\t%lld\t", Shed->Iteration);

		if(Shed->Default == TRUE)
			fprintf(Str, "Default\n");
		else
		{
			for(FIndex=0;FIndex<NO_SCHEDULE_OPT;FIndex++)
				fprintf(Str, "%f\t", Shed->Frequencies[FIndex]);
			fprintf(Str, "\n");
		}
	}
}

void SetCustomSchedule(OPTIONS* Opt, FILE* ShedFile, long long Itters, SCHEDULE* Shed)
{
	int Index;
	CUSTOM_SCHEDULE *NShed;

	if(Shed->NoCShed == 0)
		return;

	for(Index=0;Index<Shed->NoCShed;Index++)
	{
		NShed = Shed->CShedList[Index];
		if(Itters == NShed->Iteration)
		{
			if(NShed->Default == TRUE)
				memcpy(Shed->OptFreq, Shed->DefShed, sizeof(double) * NO_SCHEDULE_OPT);
			else
				memcpy(Shed->OptFreq, NShed->Frequencies, sizeof(double) * NO_SCHEDULE_OPT);

			PrintShedHeadder(Opt, Shed, ShedFile);

			return;
		}
	}
}