#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>


#include "typedef.h"
#include "trees.h"
#include "rates.h"
#include "priors.h"
#include "likelihood.h"
#include "genlib.h"
#include "rand.h"
#include "options.h"
#include "revjump.h"
#include "data.h"
#include "gamma.h"
#include "ml.h"
#include "phyloplasty.h"

#ifdef	 JNIRUN
//	extern void	SetProgress(JNIEnv *Env, jobject Obj, int Progress);
	#include "BayesTraitsJNI.h"
#endif

void	PrintPriorHeadder(FILE* Str, OPTIONS *Opt, RATES* Rates)
{
	int		PIndex;
	PRIORS	*Prios;

	if(Opt->UseRJMCMC == FALSE)
	{
		for(PIndex=0;PIndex<Rates->NoOfRates;PIndex++)
		{
			Prios = Rates->Prios[PIndex];

			if(Prios->UseHP == TRUE)
			{
				switch(Prios->Dist)
				{
					case UNIFORM:
						fprintf(Str, "%s-Min\t%s-Max\t", Opt->RateName[Prios->RateNo], Opt->RateName[Prios->RateNo]);
						break;
				
					case EXP:
						fprintf(Str, "%s-Mean\t", Opt->RateName[Prios->RateNo]);
						break;

					default:
						fprintf(Str, "%s-Mean\t%s-Var\t", Opt->RateName[Prios->RateNo], Opt->RateName[Prios->RateNo]);
				}
			}
		}
	}
	else
	{
		Prios = Opt->RJPrior;

		if(Prios->UseHP == TRUE)
		{
			switch(Prios->Dist)
			{
				case UNIFORM:
					fprintf(Str, "RJ Prior Min\tRJ Prior Max\t");
				break;

				case EXP:
					fprintf(Str, "RJ Prior Mean\t");
				break;

				case GAMMA:
					fprintf(Str, "RJ Alpha\tRJ Beta\t");
				break;

				default:
					fprintf(Str, "RJ Mean\tRJ-Var\t");
				break;
			}
		}
	}

	fprintf(Str, "Acceptance\n");
}

void	PrintPrior(FILE* Str, PRIORS *Prior)
{
	int	Index;

	if(Prior->UseHP == FALSE)
		return;

	for(Index=0;Index<DISTPRAMS[Prior->Dist];Index++)
		fprintf(Str, "%f\t", Prior->DistVals[Index]);
}

void	PrintMCMCSample(int Itters, int Acc, OPTIONS *Opt, RATES *Rates, FILE* Str)
{
	TREES*	Trees;
	int		PIndex;

	Trees = Opt->Trees;
	Rates->Lh = Likelihood(Rates, Trees, Opt);

	fprintf(Str, "%d\t", Itters);

	Rates->HMeanCount++;
	Rates->HMeanSum += 1/exp(Rates->Lh);

	PrintRates(Str, Rates, Opt);

	if(Opt->UseRJMCMC == FALSE)
	{
		for(PIndex=0;PIndex<Rates->NoOfRates;PIndex++)
			PrintPrior(Str, Rates->Prios[PIndex]);
	}
	else
		PrintPrior(Str, Rates->Prios[0]);

	fprintf(Str, "%f\n", (double)Acc/(double)Opt->Sample);
}


void	PrintTest(int Itters, RATES* Rates)
{
	char	MType;
	int		Index;

	MType = RJModelType(Rates->MappingVect);

	printf("%d\t%d\t'", Itters, NoOfPramGroups(Rates, NULL, NULL));
		
	for(Index=0;Index<Rates->NoOfFullRates;Index++)
	{
		if(Rates->MappingVect[Index] == ZERORATENO)
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

void	InitMCMC(OPTIONS *Opt, TREES *Trees, RATES *Rates)
{
	double	*Vec;

	if(Likelihood(Rates, Trees, Opt) != ERRLH)
		return;

	if(Opt->DataType == CONTINUOUS)
	{
		printf("Err (%s,%d): Could not get inish Lh:\n", __FILE__, __LINE__);
		exit(0);
	}

	Vec = (double*)malloc(sizeof(double) * Rates->NoOfRates);
	if(Vec == NULL)
		MallocErr();

	FindValidStartSet(Vec, Rates, Trees, Opt);
	memcpy(Rates->Rates, Vec, sizeof(double) * Rates->NoOfRates);

	free(Vec);
}


#ifdef	 JNIRUN
	void	MCMC(OPTIONS *Opt, TREES *Trees, JNIEnv *Env, jobject Obj)
#else
	void	MCMC(OPTIONS *Opt, TREES *Trees)
#endif
{
	RATES*		CRates=NULL;
	RATES*		NRates=NULL;
	int			Itters;
	double		Heat;
	int			Acc=0;
	FILE		*SumOut;
	SUMMARY		*Summary;
	char		Buffer[1024];
	SCHEDULE*	Shed=NULL;
	FILE*		ShedFile=NULL;
#ifdef	JNIRUN
	long		FP;
#endif

	if(Opt->UseSchedule == TRUE)
		ShedFile = OpenWrite(Opt->ScheduleFile);
	
	Shed = CreatSchedule(Opt);

	if(Opt->UseSchedule == TRUE)
		PrintShedHeadder(Opt, Shed, ShedFile);
	
	if(Opt->Summary == TRUE)
		Summary = CreatSummary(Opt);

	CRates	=	CreatRates(Opt);
	NRates	=	CreatRates(Opt);
	
	CreatPriors(Opt, CRates);
	CreatPriors(Opt, NRates);

	SetRatesToPriors(Opt, CRates);
	SetRatesToPriors(Opt, NRates);

	if(Opt->UsePhyloPlasty == TRUE)
		InitPPFiles(Opt, Trees, CRates);


	#ifndef JNIRUN
		PrintOptions(stdout, Opt);
		PrintRatesHeadder(stdout, Opt);
		PrintPriorHeadder(stdout, Opt, CRates);
		fflush(stdout);
	#endif

	PrintOptions(Opt->LogFile, Opt);

	#ifdef JNIRUN
		fflush(Opt->LogFile);
		FP = ftell(Opt->LogFile);	
/*		GotoFileEnd(Opt->LogFileRead, Opt->LogFileBuffer, LOGFILEBUFFERSIZE); */
	#endif

	PrintRatesHeadder(Opt->LogFile, Opt);
	PrintPriorHeadder(Opt->LogFile, Opt, CRates);

	#ifdef JNIRUN
		fflush(Opt->LogFile);
		fseek(Opt->LogFileRead, FP, SEEK_SET);
		fgets(Opt->LogFileBuffer, LOGFILEBUFFERSIZE, Opt->LogFileRead);
		ProcessHeaders(Env, Obj, Opt);
	#endif

	InitMCMC(Opt, Trees, CRates);
	

	CRates->Lh	=	Likelihood(CRates, Trees, Opt);
	CalcPriors(CRates, Opt);

/*	FindNoBL(Opt, Trees, CRates);  */


	for(Itters=0;;Itters++)
	{ 
		CopyRates(NRates, CRates, Opt);
		MutateRates(Opt, NRates, Shed);

		if(Opt->NodeData == TRUE)
			SetTreeAsData(Opt, Trees, NRates->TreeNo);

		NRates->Lh = Likelihood(NRates, Trees, Opt);

		Heat = NRates->Lh - CRates->Lh;
		CalcPriors(NRates, Opt);
		
		Heat += NRates->LhPrior - CRates->LhPrior;
		Heat += NRates->LogHRatio;
		
		if((log(GenRandState(CRates->RandStates)) <= Heat) && (NRates->Lh != ERRLH))
		{
			Swap(&NRates, &CRates);
			Acc++;
			Shed->Accepted[Shed->Op]++;	
		}

		if(NRates->Lh == ERRLH)
			Itters--;
		else
		{
			if(Itters%Opt->Sample==0)
			{
				if(Itters >= Opt->BurnIn)
				{
					#ifndef JNIRUN
						PrintMCMCSample(Itters, Acc, Opt, CRates, stdout);
						fflush(stdout);
					#endif

					PrintMCMCSample(Itters, Acc, Opt, CRates, Opt->LogFile);
					fflush(Opt->LogFile);

					if(Opt->Summary == TRUE)
						UpDataSummary(Summary, CRates, Opt);

					if(Opt->UseSchedule == TRUE)
						PrintShed(Opt, Shed, ShedFile);

					if(Opt->UsePhyloPlasty == TRUE)
						PrintPPOutput(Opt, Trees, CRates, Itters);
					
					BlankSchedule(Shed);
					#ifdef JNIRUN
						fgets(Opt->LogFileBuffer, LOGFILEBUFFERSIZE, Opt->LogFileRead);
						AddResults(Env, Obj, Opt);
					#endif
				}

				Acc=0;
				#ifdef JNIRUN
					SetProgress(Env, Obj, Itters);
				#endif
			}
		}

		CRates->LhPrior = NRates->LhPrior = 0;
		CRates->LogHRatio = NRates->LogHRatio = 0;
		CRates->LogJacobion = NRates->LogJacobion = 0;

		if((Opt->Itters == Itters) && (Opt->Itters != -1))
		{

			FreePriors(CRates);
			FreePriors(NRates);

			FreeRates(CRates);
			FreeRates(NRates);

			free(Shed);

			if(Opt->Summary == TRUE)
			{
				sprintf(&Buffer[0], "%s.%s", Opt->DataFN, SUMMARYFILEEXT);
				SumOut = OpenWrite(Buffer);

				PrintSummaryHeadder(SumOut, Summary, Opt);
				PrintSummary(SumOut, Summary, Opt);

				fclose(SumOut);

				FreeSummary(Summary);
			}

			if(Opt->UseSchedule == TRUE)
				fclose(ShedFile);

			return;
		}

		#ifdef JNIRUN
			if(Itters%100==0)
			{
				CheckStop(Env, Obj, Trees);
				if(Trees->JStop == TRUE)
				{
					FreePriors(CRates);
					FreePriors(NRates);

					FreeRates(CRates);
					FreeRates(NRates);

					free(Shed);

					return;
				}
			}
		#endif
	}
}

void	LhOverAllModels(OPTIONS *Opt, TREES *Trees)
{
	RATES	*Rates;
	int		Index;

	Rates = CreatRates(Opt);

	CreatPriors(Opt, Rates);
	SetRatesToPriors(Opt, Rates);

	Rates->HMeanCount = 0;
	Rates->HMeanSum = 0;
	printf("\n");
	
	for(Index=0;Index<Rates->NoOfModels;Index++)
	{
		Rates->ModelNo = Index;
		memcpy(Rates->Rates, Rates->FixedModels[Index], sizeof(double) * Rates->NoOfRates);
		Rates->Lh = Likelihood(Rates, Trees, Opt);
	
		if(Rates->Lh != ERRLH)
			PrintMCMCSample(Index, 1, Opt, Rates, stdout);
	}

	exit(0);
}