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
#include "threaded.h"

double	Min1D(RATES* Rates, TREES *Trees, OPTIONS *Opt, double From, double To, int Steps);

#ifdef	 JNIRUN
	#include "BayesTraitsJNI.h"
//	extern void	SetProgress(JNIEnv *Env, jobject Obj, int Progress);
#endif

#ifdef	NLOPT_BT
//	#include "nlopt.h"
	#include "NLOptBT.h"
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

void SetOtherRandRates(OPTIONS *Opt, TREES *Trees, RATES *Rates, double Mag)
{
	if(Opt->UseCovarion == TRUE)
	{
		Rates->OffToOn = RandDouble(Rates->RS) * Mag;
		Rates->OnToOff = Rates->OffToOn;
	}

	if(Opt->EstKappa == TRUE)
	{
		Rates->Kappa = RandDouble(Rates->RS);
	}
}


int		SetRandStart(OPTIONS *Opt, TREES *Trees, RATES *Rates)
{
	int Index;
	double Mag;


	for(Index=0;Index<10000;Index++)
	{
		Mag = RandMag(0.001, 1000, Rates->RS);
		
		SetRatesVec(Rates->Rates, Mag, Rates->RS, Rates->NoOfRates);

		SetOtherRandRates(Opt, Trees, Rates, Mag);

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

void	SetRandStaes(OPTIONS *Opt, TREES *Trees, RATES *Rates)
{
	int NoRates;
	double Mag;

	NoRates = Rates->NoOfRates;
	if(Opt->UseRJMCMC == TRUE)
		NoRates = Rates->NoOfRJRates;

	Mag = RandMag(0.001, 1000, Rates->RS);

	SetRatesVecFixed(Rates->Rates, RandDouble(Rates->RS) * Mag, NoRates);
}

void	FindRandConSet(OPTIONS *Opt, TREES *Trees, RATES *Rates)
{
	int Pos;
	
	do
	{
		Pos = 0;
		if(Opt->EstKappa == TRUE)
			Rates->Rates[Pos++] = RandDouble(Rates->RS) * (MAX_KAPPA - MIN_KAPPA);

		if(Opt->EstDelta == TRUE)
			Rates->Rates[Pos++] = RandDouble(Rates->RS) * (MAX_DELTA - MIN_DELTA);

		if(Opt->EstLambda == TRUE)
			Rates->Rates[Pos++] = RandDouble(Rates->RS) * (MAX_LAMBDA - MIN_LAMBDA);

		if(Opt->EstOU == TRUE)
//			Rates->Rates[Pos++] = RandDouble(Rates->RS) * (MAX_OU - MIN_OU);
			Rates->Rates[Pos++] = RandDouble(Rates->RS) * 1.0;
//			Rates->Rates[Pos++] = 0;

		MapConParams(Opt, Rates, Rates->Rates);
		Rates->Lh = Likelihood(Rates, Trees, Opt);

	}while(Rates->Lh == ERRLH);
}

void	FindValidStartSet(OPTIONS *Opt, TREES *Trees, RATES *Rates)
{
	double Lh;
	if(Rates->NoOfRates == 0)
		return;

	if(Opt->ModelType == MT_CONTRAST || Opt->ModelType == MT_FATTAIL)
	{
		Lh = Likelihood(Rates, Trees, Opt);
		if(Lh == ERRLH)
		{
			printf("Starting likelihood for contrast is invalid.\n");
			exit(0);
		}
		return;
	}

	if(Opt->DataType == CONTINUOUS)
	{
		FindRandConSet(Opt, Trees, Rates);
		return;
	}

	/* Set a fully random set of rates */
	if(SetRandStart(Opt, Trees, Rates) == TRUE)
		return;
	
	if(Set1RandRate(Opt, Trees, Rates) == TRUE)
		return;

	printf("Cannot find a valid set of starting parameters. Please see manual for solutions.\n", __FILE__, __LINE__);
	exit(0);
}

double	PraxisGo(OPTIONS *Opt, RATES *Rates, TREES *Trees)
{
	double		*TempVec;
	double		Ret;
	PRAXSTATE*	PState;

	TempVec = (double*)malloc(sizeof(double)*Rates->NoOfRates);
	if(TempVec == NULL)
		MallocErr();
	
	FindValidStartSet(Opt, Trees, Rates);

	memcpy(TempVec, Rates->Rates, sizeof(double)*Rates->NoOfRates);

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

//		Ret = praxis(PState);

#ifndef NLOPT_BT
		Ret = praxis(PState);
		Ret = LhPraxis((void*)PState, TempVec);
#else
		Ret = NLOptBT(TempVec, PState);
#endif
				
		FreePracxStates(PState);		
	}

#ifndef NLOPT_BT
	CheckRatesVec(TempVec, Rates->NoOfRates);
	memcpy(Rates->Rates, TempVec, sizeof(double)*Rates->NoOfRates);
	
	Rates->Lh = Likelihood(Rates, Trees, Opt);
#endif

	free(TempVec);
	return Rates->Lh;
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


void	MLTree(OPTIONS *Opt, TREES *Trees, RATES *Rates)
{
	int Index;
	double CLh, BLh, TLh;
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
	BLh = Likelihood(Rates, Trees, Opt);
	
	for(Index=0;Index<Opt->MLTries;Index++)
	{
		CLh = PraxisGo(Opt, Rates, Trees);
		TLh = Likelihood(Rates, Trees, Opt);
		
		if(CLh > BLh)
		{
			BLh = CLh;
			memcpy(BRates, Rates->Rates, sizeof(double)*Rates->NoOfRates);

		}

//		printf("%d\tCLh:\t%f\tBLh:\t%f\n", Index, CLh, BLh);fflush(stdout);

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
	if(Opt->DataType == CONTINUOUS)
		MapConParams(Opt, Rates, BRates);
	free(BRates);

	Rates->Lh = Likelihood(Rates, Trees, Opt);
/*	printf("My Lh:\t%f\n", Rates->Lh);
	printf("%f\t%f\n", Rates->Rates[0], Rates->Rates[1]);
	fflush(stdout); exit(0);*/
}

void	Test(OPTIONS *Opt, TREES* Trees, RATES* Rates)
{
	double Rate;

	for(Rate=0.00000000001;Rate<20;Rate+=0.001)
	{
		Rates->Rates[0] = Rate;
		Rates->Lh = Likelihood(Rates, Trees, Opt);
		printf("%f\t%f\n", Rates->Lh, Rate);
	}

	exit(0);
	return;
}

void	InitML(OPTIONS *Opt, TREES *Trees, RATES *Rates)
{
	if(Opt->Model == M_CONTRAST_REG)
	{
		NormaliseReg(Opt, Trees, Rates);
	}
}

#ifdef	 JNIRUN
	void	FindML(OPTIONS *Opt, TREES *Trees, JNIEnv *Env, jobject Obj)
#else
	void	FindML(OPTIONS *Opt, TREES *Trees)
#endif
{
	int		TIndex;
	RATES	*Rates;
	FILE	*SumOut;
	SUMMARY	*Summary;
	char	Buffer[1024];
	double	StartT, EndT;
#ifdef JNIRUN
	long	FP;
#endif

	if(Opt->Summary == TRUE)
		Summary = CreatSummary(Opt);

	StartT = GetSeconds();	

	Rates = CreatRates(Opt);
	
	InitML(Opt, Trees, Rates);

//	Test(Opt, Trees, Rates);
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

	fflush(stdout);
	
//	Test(Opt, Trees, Rates);
	for(TIndex=0;TIndex<Trees->NoOfTrees;TIndex++)
	{
		if(Opt->ModelType == MT_CONTINUOUS)
			InitContinusTree(Opt, Trees, TIndex);

		Rates->TreeNo = TIndex;

		if((Opt->NodeData == TRUE) || (Opt->NodeBLData == TRUE))
			SetTreeAsData(Opt, Trees, Rates->TreeNo);

		MLTree(Opt, Trees, Rates);
		
		Rates->Lh = Likelihood(Rates, Trees, Opt);
  		
		if(Opt->Summary == TRUE)
			UpDataSummary(Summary, Rates, Opt);

		PrintRates(Opt->LogFile, Rates, Opt, NULL);
		fprintf(Opt->LogFile, "\n");
		fflush(Opt->LogFile);
	
		#ifndef JNIRUN
			PrintRates(stdout, Rates, Opt, NULL);
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

		if(Opt->ModelType == MT_CONTINUOUS)
		{
			FreeConVar(Trees->Tree[TIndex]->ConVars, Trees->NoOfTaxa);
			Trees->Tree[TIndex]->ConVars = NULL;
		}
	}

	EndT = GetSeconds();

	printf("Sec:\t%f\n", EndT - StartT);

	FreeRates(Rates, Trees);

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

