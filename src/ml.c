#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "typedef.h"
#include "RandLib.h"
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

void	SetRatesVecFixed(double *Vec, double Val, int Size)
{
	int Index;

	for(Index=0;Index<Size;Index++)
		Vec[Index] = Val;
}

void	SetRatesVec(double *Vec, double Mag, RANDSTATES *RS, int Size)
{
	int Index;

	for(Index=0;Index<Size;Index++)
	{
		Vec[Index] = RandDouble(RS);
	//	Vec[Index] = RandDouble(RS) * Mag;

	}
}

double	RandMag(double Min, double Max, RANDSTATES *RS)
{
	double Ret;
	int		No, Index;

	Ret = Max / Min;
	No = (int)log10(Ret);
	No++;

	No = RandUSInt(RS) % No;

	Ret = Min;
	for(Index=0;Index<No;Index++)
		Ret = Ret * 10;

	return Ret;
}

int		SetRandStart(OPTIONS *Opt, TREES *Trees, RATES *Rates)
{
	int Index;
	double Mag;


	for(Index=0;Index<10000;Index++)
	{
		Mag = RandMag(0.001, 1000, Rates->RS);
		
		SetRatesVec(Rates->Rates, Mag, Rates->RS, Rates->NoOfRates);

		Rates->Lh = Likelihood(Rates, Trees, Opt);

		if(Rates->Lh != ERRLH)
			return TRUE;

	}

	return FALSE;	
}

int		Set1RandRate(OPTIONS *Opt, TREES *Trees, RATES *Rates)
{
	int Index;
	double Mag;

	for(Index=0;Index<10000;Index++)
	{
		Mag = RandMag(0.001, 1000, Rates->RS);

		SetRatesVecFixed(Rates->Rates, RandDouble(Rates->RS) * Mag, Rates->NoOfRates);

		Rates->Lh = Likelihood(Rates, Trees, Opt);
		
		if(Rates->Lh != ERRLH)
			return TRUE;

	}
	return FALSE;	
}


void	FindValidStartSet(OPTIONS *Opt, TREES *Trees, RATES *Rates)
{
	if(Rates->NoOfRates == 0)
		return;

	/* Set a fully random set of rates */
	if(SetRandStart(Opt, Trees, Rates) == TRUE)
		return;
	
	if(Set1RandRate(Opt, Trees, Rates) == TRUE)
		return;

	printf("%s::%d Cannot find a valid set of starting parameters, likelihood calculation may be singular.\n", __FILE__, __LINE__);
	exit(0);
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

	FindValidStartSet(Opt, Trees, Rates);

	memcpy(TempVec, Rates->Rates, sizeof(double)*Rates->NoOfRates);
	/*
	printf("Start Vect:\t%d\t%f\t", Rates->TreeNo, Rates->Lh);
	for(Index=0;Index<Rates->NoOfRates;Index++)
		printf("%f\t", Rates->Rates[Index]);
	printf("\n");
	*/
	if(Rates->NoOfRates > 0)
	{
		if(Rates->NoOfRates > 1)
//			PState = IntiPraxis(LhPraxis, TempVec, Rates->NoOfRates, 0, 1, 4, 5000);	
			PState = IntiPraxis(LhPraxis, TempVec, Rates->NoOfRates, 0, 0, 4, 5000);	
		else
			PState = IntiPraxis(LhPraxis, TempVec, Rates->NoOfRates, 0, 0, 4, 250);

		PState->Opt		= Opt;
		PState->Trees	= Trees;
		PState->Rates	= Rates;

		Ret = praxis(PState);

		FreePracxStates(PState);		
	}


	CheckRatesVec(TempVec, Rates->NoOfRates);
	memcpy(Rates->Rates, TempVec, sizeof(double)*Rates->NoOfRates);

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

	InitContinusTree(Opt, Trees, 0);


	Rates->Lambda = 0;
	printf("\n");
	d = 0.000001;
	for(i=0;i<1001;i++)
	{		
		Rates->Lambda = d;
		printf("%d\t%f\t", i, Rates->Lambda);
		Rates->Lh = Likelihood(Rates, Trees, Opt);
		printf("%f\n", Rates->Lh);
		fflush(stdout);
		d += 0.001;
	}

	exit(0);
	return;
}

void	MLTree(OPTIONS *Opt, TREES *Trees, RATES *Rates)
{
	int OIndex, Index;
	double CLh, BLh;
	double *BRates;
	
	if(Rates->NoOfRates == 0)
	{
		Rates->Lh = Likelihood(Rates, Trees, Opt);
		return;
	}

	if(Rates->NoOfRates == 1)
	{
		Opt1D(Opt, Rates, Trees);
		Rates->Lh = Likelihood(Rates, Trees, Opt);
		return;
	}

	FindValidStartSet(Opt, Trees, Rates);

	BRates = (double*)malloc(sizeof(double)*Rates->NoOfRates);
	if(BRates == NULL)
		MallocErr();
	memcpy(BRates, Rates->Rates, sizeof(double)*Rates->NoOfRates);
	BLh = -Likelihood(Rates, Trees, Opt);;
	
	for(OIndex=0;OIndex<Opt->MLTries;OIndex++)
	{
		CLh = PraxisGo(Opt, Rates, Trees);

	//	printf("BestLh\t%f\n", CLh);
	//	exit(0);		

		if(CLh < BLh)
		{
			BLh = CLh;
			memcpy(BRates, Rates->Rates, sizeof(double)*Rates->NoOfRates);
		}

		#ifdef JNIRUN
			CheckStop(Env, Obj, Trees);
			if(Trees->JStop == TRUE)
			{
				free(BRates);
				return;
			}
		#endif
	}

	memcpy(Rates->Rates, BRates, sizeof(double)*Rates->NoOfRates);
	free(BRates);

	Rates->Lh = Likelihood(Rates, Trees, Opt);
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
	FILE	*SumOut;
	SUMMARY	*Summary;
	char	Buffer[1024];
	int		ti;
#ifdef JNIRUN
	long	FP;
#endif

	if(Opt->Summary == TRUE)
		Summary = CreatSummary(Opt);

	Rates = CreatRates(Opt);

		
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

		MLTree(Opt, Trees, Rates);
		
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
			CheckStop(Env, Obj, Trees);
			if(Trees->JStop == TRUE)
			{
				FreeRates(Rates);
				return;
			}

			fgets(Opt->LogFileBuffer, LOGFILEBUFFERSIZE, Opt->LogFileRead);
			AddResults(Env, Obj, Opt);

			SetProgress(Env, Obj, TIndex);
		#endif

		if(Opt->DataType == CONTINUOUS)
		{
			if(Opt->Model != CONTRASTM)
			{
				FreeConVar(Trees->Tree[TIndex]->ConVars, Trees->NoOfTaxa);
				Trees->Tree[TIndex]->ConVars = NULL;
			}
		}
	}

	FreeRates(Rates);

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

