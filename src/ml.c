#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "typedef.h"
#include "rand.h"
#include "genlib.h"
#include "rates.h"
#include "praxis.h"
#include "likelihood.h"
#include "options.h"
#include "1dopt.h"
#include "data.h"
#include "continuous.h"
#include "contrasts.h"

double	Min1D(RATES* Rates, TREES *Trees, OPTIONS *Opt, double From, double To, int Steps);

#ifdef	 JNIRUN
	#include "BayesTraitsJNI.h"
//	extern void	SetProgress(JNIEnv *Env, jobject Obj, int Progress);
#endif


void	FindValidStartSet(double *Vec, RATES *Rates, TREES *Trees, OPTIONS *Opt)
{
	double		Lh;
	int			Index;
	int			i;
	PRAXSTATE	*PState;

	if(Rates->NoOfRates == 0)
		return;

	PState			= IntiPraxis(LhPraxis, NULL, Rates->NoOfRates, 0, 1, 4, 5000);
	PState->Opt		= Opt;
	PState->Trees	= Trees;
	PState->Rates	= Rates;

	i = 0;
	do
	{ 
		for(Index=0;Index<Rates->NoOfRates;Index++)
		{
			if(Opt->DataType == CONTINUOUS)
			{
				if(i == 0)		
					Vec[Index] = 1;
				else
					Vec[Index] = GenRandState(Rates->RandStates);
			}
			else
			{
				if(GenRandState(Rates->RandStates)<0.1)
					Vec[Index] = GenRandState(Rates->RandStates) * 10;
				else
					Vec[Index] = GenRandState(Rates->RandStates) * 0.1;
			}
		}	
	
		Lh = LhPraxis(PState, Vec);
		
		i++;
	} while(Lh == -ERRLH);

	FreePracxStates(PState);
}

double	PraxisGo(OPTIONS *Opt, RATES *Rates, TREES *Trees)
{
	double		*TempVec;
	int			Index;
	double		Ret;
	PRAXSTATE*	PState;

	TempVec = (double*)malloc(sizeof(double)*Rates->NoOfRates);
	if(TempVec == NULL)
		MallocErr();

	FindValidStartSet(TempVec, Rates, Trees, Opt);

	if(Rates->NoOfRates > 0)
	{
		if(Rates->NoOfRates > 1)
			PState = IntiPraxis(LhPraxis, TempVec, Rates->NoOfRates, 0, 1, 4, 5000);	
		else
			PState = IntiPraxis(LhPraxis, TempVec, Rates->NoOfRates, 0, 1, 4, 250);

		PState->Opt		= Opt;
		PState->Trees	= Trees;
		PState->Rates	= Rates;

		Ret = praxis(PState);
		FreePracxStates(PState);		
	}

	for(Index=0;Index<Rates->NoOfRates;Index++)
		if(TempVec[Index] < MINRATE)
			Rates->Rates[Index] = MINRATE;
		else
			Rates->Rates[Index] = TempVec[Index];

	Rates->Lh = Ret;
	free(TempVec);

	return Ret;
}

void	TestCL(OPTIONS *Opt, TREES* Trees, RATES *Rates)
{
	double D;

	for(D=0.001;D<12;D+=0.001)
	{
		Rates->Delta = D;
		printf("%f\t%f\n", D, Likelihood(Rates, Trees, Opt));
	}
}

void	Test(OPTIONS *Opt, TREES* Trees, RATES* Rates)
{
	double d;
	int		i;


	Rates->Rates[0] = 0.833651;

	printf("\n");
	d = 0.0001;
	for(i=0;i<1000;i++)
	{		
		Rates->Rates[1] = d;
		printf("%d\t%f\t%f\t", i, Rates->Rates[0], Rates->Rates[1]);
		Rates->Lh = Likelihood(Rates, Trees, Opt);
		printf("%f\n", Rates->Lh);
		fflush(stdout);
		d += 0.1;
	}

	return;
}
#ifdef	 JNIRUN
	void	FindML(OPTIONS *Opt, TREES *Trees, JNIEnv *Env, jobject Obj)
#else
	void	FindML(OPTIONS *Opt, TREES *Trees)
#endif
{
	int		TIndex;
	int		OIndex;
	RATES	*Rates;
	double	*BRates=NULL;
	double	CLh;
	double	BLh;
	FILE	*SumOut;
	SUMMARY	*Summary;
	char	Buffer[1024];
	int		ti;
	long	FP;

	if(Opt->Summary == TRUE)
		Summary = CreatSummary(Opt);

	Rates = CreatRates(Opt);

/*	Some Test code */
/*	
	Test(Opt, Trees, Rates);
	FreeRates(Rates);
	return;
*/
/*	End of some test code */

	BRates = (double*)malloc(sizeof(double)*Rates->NoOfRates);
	if(BRates==NULL)
		MallocErr();

	#ifndef JNIRUN
		PrintOptions(stdout, Opt);
		PrintRatesHeadder(stdout, Opt);
	#endif

	if(Opt->Headers == TRUE)
	{
		PrintOptions(Opt->LogFile, Opt);
		
		#ifdef JNIRUN
			fflush(Opt->LogFile);
			FP = ftell(Opt->LogFile);
		/*	GotoFileEnd(Opt->LogFileRead, Opt->LogFileBuffer, LOGFILEBUFFERSIZE); */
		#endif	

		PrintRatesHeadder(Opt->LogFile, Opt);

		#ifdef JNIRUN
			fflush(Opt->LogFile);
			fseek(Opt->LogFileRead, FP, SEEK_SET);
			fgets(Opt->LogFileBuffer, LOGFILEBUFFERSIZE, Opt->LogFileRead);
			ProcessHeaders(Env, Obj, Opt);
		#endif
	}

	for(TIndex=0;TIndex<Trees->NoOfTrees;TIndex++)
	{
		if(Opt->DataType == CONTINUOUS)
			InitContinusTree(Opt, Trees, TIndex);

		Rates->TreeNo = TIndex;

		if((Opt->NodeData == TRUE) || (Opt->NodeBLData == TRUE))
			SetTreeAsData(Opt, Trees, Rates->TreeNo);
		
		BLh = ERRLH;

		if(Opt->Model == CONTRASTM)
			CalcContrastLh(Opt, Trees, Rates);
		else
		{
			if(Rates->NoOfRates > 1)
			{		
				for(OIndex=0;OIndex<Opt->MLTries;OIndex++)
				{
					CLh = PraxisGo(Opt, Rates, Trees);
				
					if((CLh < BLh) || (BLh == ERRLH))
					{
						BLh = CLh;
						for(ti=0;ti<Rates->NoOfRates;ti++)
							BRates[ti] = Rates->Rates[ti];

						if(Opt->DataType == CONTINUOUS)
							LHRandWalk(Opt, Trees, Rates);
					}

					#ifdef JNIRUN
						CheckStop(Env, Obj, Trees);
						if(Trees->JStop == TRUE)
						{
							FreeRates(Rates);
							free(BRates);
							return;
						}
					#endif
				}

				for(ti=0;ti<Rates->NoOfRates;ti++)
					Rates->Rates[ti] = BRates[ti];

				if(Opt->DataType == CONTINUOUS)
				//	LhPraxisCon(BRates);
					LHRandWalk(Opt, Trees, Rates);
			}
			else
			{
				if(Rates->NoOfRates == 1)
					Opt1D(Opt, Rates, Trees);
			}
		}

		Rates->Lh = Likelihood(Rates, Trees, Opt);
  		
		if(Opt->Summary == TRUE)
			UpDataSummary(Summary, Rates, Opt);

		PrintRates(Opt->LogFile, Rates, Opt);
		fprintf(Opt->LogFile, "\n");
		fflush(Opt->LogFile);

		#ifndef JNIRUN
			PrintRates(stdout, Rates, Opt);
			printf("\n");
			fflush(stdout);
		#else
			fgets(Opt->LogFileBuffer, LOGFILEBUFFERSIZE, Opt->LogFileRead);
			AddResults(Env, Obj, Opt);

			SetProgress(Env, Obj, TIndex);
		#endif

		if(Opt->DataType == CONTINUOUS)
		{
			if(Opt->Model == CONTRASTM)
			{
				
			}
			else
			{
				FreeConVar(Trees->Tree[TIndex].ConVars, Trees->NoOfTaxa);
				Trees->Tree[TIndex].ConVars = NULL;
			}
		}
	}

	FreeRates(Rates);
	free(BRates);

	if(Opt->Summary == TRUE)
	{
		sprintf(&Buffer[0], "%s.%s", Opt->DataFN, SUMMARYFILEEXT);
		SumOut = OpenWrite(Buffer);

		PrintSummaryHeadder(SumOut, Summary, Opt);
		PrintSummary(SumOut, Summary, Opt);

		fclose(SumOut);

		FreeSummary(Summary);
	}
}



double	Min1D(RATES* Rates, TREES *Trees, OPTIONS *Opt, double From, double To, int Steps)
{
	double	StepSize;
	int		Itters;
	double	Val;
	double	Ret=To;
	double	Best=-999999;
	double	New;

	StepSize = (To - From) / (double)Steps;
	
	for(Itters=0,Val=From;Itters<Steps;Itters++,Val+=StepSize)
	{
		Rates->Rates[0] = Val;

		New = Likelihood(Rates, Trees, Opt);

		if((New > Best) && (New < 0))
		{
			Ret = Val;
			Best = New;
		}

	}

	Rates->Rates[0] = Ret;
	Rates->Lh		= Best;

	return Best;
}

