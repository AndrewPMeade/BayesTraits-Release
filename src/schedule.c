#include <stdio.h>
#include <math.h>
#include <stdlib.h>

#include "typedef.h"
#include "genlib.h"
#include "schedule.h"
#include "data.h"
#include "AutoTune.h"
#include "rates.h"

void	UpDateShedAcc(int Acc, SCHEDULE* Shed)
{
	Shed->Tryed[Shed->Op]++;
	if((Shed->Op == SRATES) && (Shed->RateDevPerParm == TRUE))
		Shed->PTried[Shed->PNo]++;
	
	Shed->GNoTried++;
	Shed->SNoTried++;

	if(Acc == FALSE)
		return;

	Shed->GNoAcc++;
	Shed->SNoAcc++;
	Shed->Accepted[Shed->Op]++;
	if((Shed->Op == SRATES) && (Shed->RateDevPerParm == TRUE))
		Shed->PAcc[Shed->PNo]++;
}

int		UsingHP(OPTIONS *Opt)
{
	int Index;

	for(Index=0;Index<Opt->NoOfRates;Index++)
	{
		if(Opt->Priors[Index]->UseHP == TRUE)
			return TRUE;
	}

	if(Opt->RJPrior->UseHP == TRUE)
		return TRUE;

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

	if(Shed->RateDevPerParm == TRUE)
	{
		for(Index=0;Index<Shed->NoParm;Index++)
			Shed->PAcc[Index] = Shed->PTried[Index] = 0;
	}
}

void	ScaleSchedVect(SCHEDULE * Sched)
{
	double	SF=0;
	int		Index=0;

	for(Index=0;Index<Sched->NoOfOpts;Index++)
		SF += Sched->OptFreq[Index];

	SF = 1 / SF;

	for(Index=0;Index<Sched->NoOfOpts;Index++)
		Sched->OptFreq[Index] *= SF;
}
/*
void		ScaleVect(double *Vect, int VectSize)
{
	double	SF=0;
	int		Index=0;

	for(Index=0;Index<VectSize;Index++)
		SF += Vect[Index];


	SF = 1 / SF;

	for(Index=0;Index<VectSize;Index++)
		Vect[Index] *= SF;

	SF = 0;
	for(Index=0;Index<VectSize;Index++)
	{
		Vect[Index] = Vect[Index] + SF;
		SF = Vect[Index];
	}
}
*/
int		MultiTree(OPTIONS *Opt)
{
	if(Opt->UseEqualTrees == TRUE)
		return FALSE;

	if(Opt->Trees->NoOfTrees == 1)
		return FALSE;

	return TRUE;
}


void	SetSchedule(SCHEDULE*	Shed, OPTIONS *Opt)
{
	double	Left;
	int		Rates, Index;
	
	Left = 1.0;

	for(Index=0;Index<Shed->NoOfOpts;Index++)
		Shed->OptFreq[Index] = 0.0;

	if((Opt->UseCovarion == TRUE) && (Opt->LoadModels == FALSE))
	{
		Shed->OptFreq[SCV] = 0.2;
		Left = Left - Shed->OptFreq[SCV];
	}

	if((Opt->EstKappa == TRUE) && (Opt->LoadModels == FALSE))
	{
		Shed->OptFreq[SKAPPA] = 0.1;
		Left = Left - Shed->OptFreq[SKAPPA];
	}

	if((Opt->EstDelta == TRUE) && (Opt->LoadModels == FALSE))
	{
		Shed->OptFreq[SDELTA] = 0.1;
		Left = Left - Shed->OptFreq[SDELTA];
	}

	if((Opt->EstLambda == TRUE)  && (Opt->LoadModels == FALSE))
	{
		Shed->OptFreq[SLABDA] = 0.1;
		Left = Left - Shed->OptFreq[SLABDA];
	}


	if((Opt->UseRJMCMC == TRUE) && (Opt->LoadModels == FALSE))
	{
		Shed->OptFreq[SJUMP] = 0.1;
		Left = Left - Shed->OptFreq[SJUMP];
	}

	if(UsingHP(Opt) == TRUE)
	{
		Shed->OptFreq[SPPROR] = 0.1;
		Left = Left - Shed->OptFreq[SPPROR];
	}

	if(EstData(Opt->Trees) == TRUE)
	{
		Shed->OptFreq[SESTDATA] = 0.5;
		Left = Left - Shed->OptFreq[SESTDATA];
	}

	if(Opt->UseVarData == TRUE)
	{
		Shed->OptFreq[SVARDATA] = 0.1;
		Left = Left - Shed->OptFreq[SVARDATA];
	}

	if((Opt->EstOU == TRUE)  && (Opt->LoadModels == FALSE))
	{
		Shed->OptFreq[SOU] = 0.1;
		Left = Left - Shed->OptFreq[SOU];
	}

	if(Opt->EstGamma == TRUE)
	{
		Shed->OptFreq[SGAMMA] = 0.1;
		Left = Left - Shed->OptFreq[SGAMMA];
	}

	if(Opt->Model == M_FATTAIL)
	{
		Shed->OptFreq[SFATTAILANS] = Left * 0.9;
		Left = Left - Shed->OptFreq[SFATTAILANS];
	}

/*
	if(Opt->SoloTreeMove == TRUE)
	{
		Shed->OptFreq[9] = 0.2;
		Left = Left - Shed->OptFreq[9];
	}
*/

	if(MultiTree(Opt) == TRUE)
		Shed->OptFreq[STREEMOVE] = 0.1;
	else
		Shed->OptFreq[STREEMOVE] = 0.0;
	

	if(Opt->UseVarRates == TRUE)
	{
/*
		SPPADDREMOVE=10,
		SPPMOVE=11,
		SPPCHANGESCALE=12,
		SPPHYPERPRIOR=13,
*/
		Shed->OptFreq[10] = 0;
		Shed->OptFreq[11] = 0;
		Shed->OptFreq[12] = 0;
		Shed->OptFreq[13] = 0;

		Shed->OptFreq[10] = 0.5;
		Shed->OptFreq[11] = 0.05;
		Shed->OptFreq[12] = 0.4;
//		Shed->OptFreq[13] = 0.05;

		Left = Left - (Shed->OptFreq[10] + Shed->OptFreq[11] + Shed->OptFreq[12] + Shed->OptFreq[13]);
	}

	if(Opt->Model == M_DESCHET)
	{
		Shed->OptFreq[14] = 0.4;
		Left = Left - Shed->OptFreq[14];
	}

	Rates = 0;
	if(Opt->DataType == CONTINUOUS)
	{
		Rates = Opt->Trees->NoOfSites;
	}
	else
		for(Index=0;Index<Opt->NoOfRates;Index++)
			if(Opt->ResTypes[Index] == RESNONE)
				Rates++;

	if(Rates == 0)
		Shed->OptFreq[SRATES] = 0;
	else
	{
		if(Left < 0.05)
			Left = 0.05;

		Shed->OptFreq[SRATES] = Left;
	}

#ifdef CONTRAST_ML_PARAM
	if(Opt->ModelType == MT_CONTRAST)
		Shed->OptFreq[0] = 0;
#endif

	if(Opt->LoadModels == TRUE)
		Shed->OptFreq[SRATES] = 0.3;
	
	if(Opt->RJDummy == TRUE)
	{
		Shed->OptFreq[SRATES]		= 0.1;
		Shed->OptFreq[SRJDUMMY]		= 0.4;
		Shed->OptFreq[SRJDUMMYMOVE] = 0.1;
		Shed->OptFreq[SRJDUMMYCHANGEBETA] = 0.5;
	}

	ScaleSchedVect(Shed);
}

void	PrintShedHeadder(OPTIONS* Opt, SCHEDULE* Shed, FILE* Str)
{
	int	Index;
	double	Last;

	Last = 0;
	
	for(Index=0;Index<Shed->NoOfOpts;Index++)
	{
		if(Shed->OptFreq[Index] != 0)
			fprintf(Str, "%s\t%2.2f\n", SHEDOP[Index], Shed->OptFreq[Index]*100);
		Last = Shed->OptFreq[Index];
	}
	
	Last = 0;
	for(Index=0;Index<Shed->NoOfOpts;Index++)
	{
		if(Shed->OptFreq[Index] != 0)
			fprintf(Str, "%s Tried\t%% Accepted\t", SHEDOP[Index]);
		Last = Shed->OptFreq[Index];
	}

	if(UseAutoTune(Opt) == TRUE)
		PrintAutoTuneHeader(Str, Opt);

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

	if(UseAutoTune(Opt) == TRUE)
		PrintAutoTune(Str, Opt, Shed);

	fprintf(Str, "%f\t", (double)Shed->SNoAcc / Shed->SNoTried);
	fprintf(Str, "%f\t", (double)Shed->GNoAcc / Shed->GNoTried);
	fprintf(Str, "\n");

	fflush(Str);
}

SCHEDULE*	AllocSchedule()
{
	SCHEDULE*	Ret;
	
	Ret = (SCHEDULE*)malloc(sizeof(SCHEDULE));
	if(Ret==NULL)
		MallocErr();
	
	Ret->NoOfOpts = NOOFOPERATORS;

	Ret->OptFreq	=	(double*)malloc(sizeof(double) * Ret->NoOfOpts);
	Ret->Tryed		=	(int*)malloc(sizeof(int) * Ret->NoOfOpts);
	Ret->Accepted	=	(int*)malloc(sizeof(int) * Ret->NoOfOpts);
	if((Ret->OptFreq == NULL) || (Ret->Tryed == NULL) || (Ret->Accepted == NULL))
		MallocErr();

	Ret->RateDevATList	= NULL;
	Ret->DataDevAT		= NULL;
	Ret->VarRateAT		= NULL;

	Ret->RateDevPerParm	= FALSE;
	Ret->NoParm			= -1;
	Ret->PTried			= NULL;
	Ret->PAcc			= NULL;

	Ret->KappaAT		= NULL;
	Ret->DeltaAT		= NULL;
	Ret->LambdaAT		= NULL;
	Ret->OUAT			= NULL;

	Ret->RJDummyBetaAT	=	NULL;
	
	return Ret;
}

int		FindNoOfAutoCalibRates(OPTIONS *Opt)
{
	if(Opt->ModelType == MT_DISCRETE)
		return 1;

	return FindNoConRates(Opt);
}

void	SetRateDevPerParm(SCHEDULE* Shed, OPTIONS *Opt, RANDSTATES *RS)
{
	int Index;

	Shed->RateDevPerParm	= TRUE;
	Shed->NoParm			= FindNoOfAutoCalibRates(Opt);
	Shed->PAcc				= (int*)malloc(sizeof(int) * Shed->NoParm);
	Shed->PTried			= (int*)malloc(sizeof(int) * Shed->NoParm);
	Shed->RateDevATList		= (AUTOTUNE**)malloc(sizeof(AUTOTUNE*) * Shed->NoParm);

	if((Shed->PAcc == NULL) || (Shed->PTried == NULL) || (Shed->RateDevATList == NULL))
		MallocErr();

	for(Index=0;Index<Shed->NoParm;Index++)
	{
		Shed->PTried[Index] = Shed->PAcc[Index] = 0;
		Shed->RateDevATList[Index] = CreatAutoTune(0.2, 0.4);
		Opt->RateDevList[Index] = RandDouble(RS) * 10;
	}
}

SCHEDULE*	CreatSchedule(OPTIONS *Opt, RANDSTATES *RS)
{
	SCHEDULE*	Ret;

	Ret = AllocSchedule();
	
	BlankSchedule(Ret);

	SetSchedule(Ret, Opt);

	Ret->GNoAcc = Ret->GNoTried = 0;
	Ret->SNoAcc = Ret->SNoTried = 0;

	if(Opt->AutoTuneRD == TRUE)
		SetRateDevPerParm(Ret, Opt, RS);
		
	if(Opt->AutoTuneDD == TRUE)
	{
		Ret->DataDevAT = CreatAutoTune(0.2, 0.4);
		Opt->EstDataDev= RandDouble(RS) * 10;
	}

	if(Opt->AutoTuneVarRates == TRUE)
	{
		Ret->VarRateAT = CreatAutoTune(0.2, 0.4);
		Opt->VarRatesScaleDev = RandDouble(RS) * 100;
	}

	if(Opt->EstKappa == TRUE)
	{
		Ret->KappaAT = CreatAutoTune(0.2, 0.4);
		Opt->RateDevKappa = RandDouble(RS) * 10;
	}

	if(Opt->EstLambda == TRUE)
	{
		Ret->LambdaAT = CreatAutoTune(0.2, 0.4);
		Opt->RateDevLambda = RandDouble(RS) * 10;
	}

	if(Opt->EstDelta == TRUE)
	{
		Ret->DeltaAT = CreatAutoTune(0.2, 0.4);
		Opt->RateDevDelta = RandDouble(RS) * 10;

	}

	if(Opt->EstOU == TRUE)
	{
		Ret->OUAT = CreatAutoTune(0.2, 0.4);
		Opt->RateDevOU = RandDouble(RS) * 10;
	}

	if(Opt->RJDummy == TRUE)
	{
		Ret->RJDummyBetaAT = CreatAutoTune(0.2, 0.4);
		Opt->RJDummyBetaDev= RandDouble(RS);
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

	if(Sched->RJDummyBetaAT != NULL)
		FreeAutoTune(Sched->RJDummyBetaAT);

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
		Acc = Shed->PAcc[0] + Shed->Accepted[SCV];
		Tried = Shed->PTried[0] + Shed->Tryed[SCV];
		Ret = (double)Acc/Tried;
	}

	return Ret;
}

void	UpDateRateDevPerP(OPTIONS *Opt, SCHEDULE* Shed, RANDSTATES *RS)
{
	int Index;
	double Acc;

	if(Opt->DataType == DISCRETE)
	{
		if(Shed->PTried[0] > 4)
		{
			Acc = GetRDDecAccRate(Opt, Shed);
			Opt->RateDevList[0] = AutoTuneNextRD(Shed->RateDevATList[0], RS, Opt->RateDevList[0], Acc);
		}
	}
	else
	{
		for(Index=0;Index<Shed->NoParm;Index++)
		{
			if(Shed->PTried[Index] > 4)
			{
				Acc = (double)Shed->PAcc[Index] / Shed->PTried[Index];
				Opt->RateDevList[Index] = AutoTuneNextRD(Shed->RateDevATList[Index], RS, Opt->RateDevList[Index], Acc);
			}
		}
	}

	Opt->RateDev = Opt->RateDevList[0];
}

void	UpDateSchedule(OPTIONS *Opt, SCHEDULE* Shed, RANDSTATES *RS)
{
	int		Tried;
	double	Acc;

	if((Opt->AutoTuneRD == TRUE) && (Shed->RateDevPerParm == TRUE))
		UpDateRateDevPerP(Opt, Shed, RS);
	
	if(Shed->DataDevAT != NULL)
	{
		Tried = Shed->Tryed[SESTDATA];
		Acc = GetAccRate(SESTDATA, Shed);

		if(Tried > 2)
			Opt->EstDataDev = AutoTuneNextRD(Shed->DataDevAT, RS, Opt->EstDataDev, Acc);
	}

	if(Shed->VarRateAT != NULL)
	{
		Tried = Shed->Tryed[SPPCHANGESCALE];
		Acc = GetAccRate(SPPCHANGESCALE, Shed);

		if(Tried > 2)
			Opt->VarRatesScaleDev = AutoTuneNextRD(Shed->VarRateAT, RS, Opt->VarRatesScaleDev, Acc);
	}

	if(Shed->KappaAT != NULL)
	{
		Tried = Shed->Tryed[SKAPPA];
		Acc = GetAccRate(SKAPPA, Shed);

		if(Tried > 2)
		{
			Opt->RateDevKappa = AutoTuneNextRD(Shed->KappaAT, RS, Opt->RateDevKappa, Acc);
			if(Opt->RateDevKappa > MAX_KAPPA)
				Opt->RateDevKappa = MAX_KAPPA;
		}
	}

	if(Shed->DeltaAT != NULL)
	{
		Tried = Shed->Tryed[SDELTA];
		Acc = GetAccRate(SDELTA, Shed);

		if(Tried > 2)
		{
			Opt->RateDevDelta = AutoTuneNextRD(Shed->DeltaAT, RS, Opt->RateDevDelta, Acc);
			if(Opt->RateDevDelta > MAX_DELTA)
				Opt->RateDevDelta = MAX_DELTA;
		}
	}

	if(Shed->LambdaAT != NULL)
	{
		Tried = Shed->Tryed[SLABDA];
		Acc = GetAccRate(SLABDA, Shed);

		if(Tried > 2)
		{
			Opt->RateDevLambda = AutoTuneNextRD(Shed->LambdaAT, RS, Opt->RateDevLambda, Acc);
			if(Opt->RateDevLambda > MAX_LAMBDA)
				Opt->RateDevLambda = MAX_LAMBDA;
		}
	}

	if(Shed->OUAT != NULL)
	{
		Tried = Shed->Tryed[SOU];
		Acc = GetAccRate(SOU, Shed);

		if(Tried > 2)
		{
			Opt->RateDevOU = AutoTuneNextRD(Shed->OUAT, RS, Opt->RateDevOU, Acc);
			if(Opt->RateDevOU  > MAX_OU)
				Opt->RateDevOU = MAX_OU;
		}
	}

	if(Shed->RJDummyBetaAT != NULL)
	{
		Tried = Shed->Tryed[SRJDUMMYCHANGEBETA];
		Acc = GetAccRate(SRJDUMMYCHANGEBETA, Shed);

		if(Tried > 4)
		{
			Opt->RJDummyBetaDev = AutoTuneNextRD(Shed->RJDummyBetaAT, RS, Opt->RJDummyBetaDev, Acc);
			if(Opt->RJDummyBetaDev > 100)
				Opt->RJDummyBetaDev = 100;
		}

	}
}

int		UseAutoTune(OPTIONS *Opt)
{
	if(Opt->AutoTuneRD == TRUE)
		return TRUE;

	if(Opt->AutoTuneDD == TRUE)
		return TRUE;

	if(Opt->AutoTuneVarRates == TRUE)
		return TRUE;

	if(Opt->EstKappa == TRUE)
		return TRUE;

	if(Opt->EstDelta == TRUE)
		return TRUE;

	if(Opt->EstLambda == TRUE)
		return TRUE;

	if(Opt->EstOU == TRUE)
		return TRUE;

	return FALSE;
}

