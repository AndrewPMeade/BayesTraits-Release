#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>


#include "typedef.h"
#include "trees.h"
#include "Rates.h"
#include "priors.h"
#include "likelihood.h"
#include "genlib.h"
#include "RandLib.h"
#include "options.h"
#include "revjump.h"
#include "data.h"
#include "gamma.h"
#include "ml.h"
#include "VarRates.h"
#include "Threaded.h"
#include "schedule.h"
#include "modelfile.h"
#include "Stones.h"
#include "RJDummy.h"
#include "contrasts.h"
#include "FatTail.h"
#include "Geo.h"
#include "TransformTree.h"
#include "Prob.h"

#include <gsl/gsl_rng.h>


#ifdef	 JNIRUN
//	extern void	SetProgress(JNIEnv *Env, jobject Obj, int Progress);
	#include "BayesTraitsJNI.h"
#endif

void	PrintPriorHeadder(FILE* Str, OPTIONS *Opt, RATES* Rates)
{
	int		PIndex;
	PRIOR	*Prior;

	for(PIndex=0;PIndex<Rates->NoPriors;PIndex++)
	{
		Prior = Rates->Priors[PIndex];

		if(Prior->UseHP == TRUE)
		{
			switch(Prior->Dist)
			{
				case EXP:
					fprintf(Str, "%s - Mean\t", Prior->Name);
				break;

				case GAMMA:
					fprintf(Str, "%s - Shape\t%s - Scale\t", Prior->Name, Prior->Name);
				break;
	
				default:
					printf("%s::%d Hyper Prior not supported.", __FILE__, __LINE__);
					exit(0);
			}
		}
	}

	fprintf(Str, "\n");
}


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

void	PrintMCMCSample(long long Itters, SCHEDULE* Shed, OPTIONS *Opt, RATES *Rates, FILE* Str)
{
	TREES*	Trees;
	int		PIndex;
	double	HMean;
	PRIOR	*Prior;

	Trees = Opt->Trees;

	HMean = GetHMean(Opt, Rates);
	fprintf(Str, "%lld\t%f\t%f\t%d\t", Itters, Rates->Lh, HMean, Rates->TreeNo+1);
		
	PrintRates(Str, Rates, Opt, Shed);

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

void	TestInitPho(OPTIONS *Opt, TREES *Trees, RATES *Rates)
{
	int Index, RIndex;

	for(Index=0;Index<1000;Index++)
	{
		for(RIndex=0;RIndex<Rates->NoOfRates;RIndex++)
		{
			Rates->Rates[RIndex] = RandDouble(Rates->RS) * 0.01;
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

	return RandFromPrior(RNG, Prior);
}


void	SetDefMCMCParameters(OPTIONS *Opt, TREES *Trees, RATES *Rates, 	gsl_rng *RNG)
{
	int Index;
	LOCAL_TRANSFORM *LR;
	PRIOR *Prior;

	if(Opt->EstKappa == TRUE)
		Rates->Kappa = GetRandValFromType(VR_KAPPA, Rates, RNG);

	if(Opt->EstLambda == TRUE)
		Rates->Lambda = GetRandValFromType(VR_LAMBDA, Rates, RNG);

	if(Opt->EstDelta == TRUE)
		Rates->Delta = GetRandValFromType(VR_DELTA, Rates, RNG);
	
	if(Opt->EstOU == TRUE)
		Rates->OU = GetRandValFromType(VR_OU, Rates, RNG);
	
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


	CRate = 100000;
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

	} while(CRate > 1E-20);

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


	for(TNo=0;TNo<10000;TNo++)
	{
		for(RIndex=0;RIndex<Rates->NoOfRates;RIndex++)
			Rates->Rates[RIndex] = RandFromPriorPosition(RIndex, Opt, Trees, Rates, RNG);

		if(ValidMCMCParameters(Opt, Trees, Rates) != ERRLH)
			return;
	}

	printf("Cannot find initial starting set of parameters.");
	exit(0);
}

void FindValidStartLh(OPTIONS *Opt, TREES *Trees, RATES *Rates, gsl_rng *RNG)
{

	if(Opt->ModelType == MT_DISCRETE)
	{
		if(FindValidStartRateAllSame(Opt, Trees, Rates) == TRUE)
			return;
	}
	
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
}

void	ShowTimeSec(double StartT, double EndT)
{
	printf("Sec:\t%f\n", EndT - StartT);
}


FILE*	SetScheduleFile(OPTIONS *Opt, SCHEDULE*	Shed)
{
	FILE *Ret;
	char	*Buffer;

	Buffer = (char*)malloc(sizeof(char) * BUFFERSIZE);
	if(Buffer == NULL)
		MallocErr();

	sprintf(Buffer, "%s.Schedule.txt", Opt->LogFN);

	Ret = OpenWrite(Buffer);

	PrintShedHeadder(Opt, Shed, Ret);

	free(Buffer);

	return Ret;
}

FILE*	CreatStoneOuput(OPTIONS *Opt)
{
	FILE*	Ret;
	char*	Buffer;

	Buffer = (char*)malloc(sizeof(char) * BUFFERSIZE);
	if(Buffer == NULL)
		MallocErr();

	sprintf(Buffer, "%s.Stones.txt", Opt->LogFN);
	Ret = OpenWrite(Buffer);

	OutputStoneHeadder(Ret, Opt->Stones);

	free(Buffer);

	return Ret;
}

int		ExitMCMC(OPTIONS *Opt, long long Itters)
{
	if(Opt->Stones != NULL)
		return StoneExit(Opt->Stones, Itters);
	
	if((Opt->Itters == Itters) && (Opt->Itters != -1))
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

int		MCMCAccept(long long Itters, OPTIONS *Opt, TREES *Trees, SCHEDULE* Shed, RATES *CRates, RATES *NRates)
{
	double Heat;

	if(Shed->Op == SFATTAILANS && Opt->Model == M_FATTAIL)
		return TRUE;
	
	Heat = NRates->Lh - CRates->Lh;
						
	if(Opt->Stones != NULL)
		Heat = GetStoneHeat(Opt->Stones, Itters, Heat);
			
	Heat += NRates->LhPrior - CRates->LhPrior;
	Heat += NRates->LnHastings;

	if(log(RandDouble(CRates->RS)) <= Heat)
		return TRUE;

	return FALSE;
}

void	MCMCTest(OPTIONS *Opt, TREES *Trees, RATES *Rates)
{
	double Lh, X;

	Rates->Rates[0] = FAT_TAIL_NORMAL_VAL;
	Rates->Rates[1] = 10.0;

	Lh = Likelihood(Rates, Trees, Opt);

	X = 0.000001;
	while(X<1000)
	{
		Rates->Rates[1] = X;
		Lh = Likelihood(Rates, Trees, Opt);
		printf("%f\t%f\n", X, Lh);
		X += 0.01;
	}

	exit(0);
}

#ifdef	 JNIRUN
	void	MCMC(OPTIONS *Opt, TREES *Trees, JNIEnv *Env, jobject Obj)
#else
	void	MCMC(OPTIONS *Opt, TREES *Trees)
#endif
{
	RATES*		CRates;
	RATES*		NRates;
	long long	Itters;
	double		StartT, EndT;
	SCHEDULE*	Shed;
	FILE*		ShedFile;
	FILE*		SaveModelF;
	FILE*		StoneF;
	int			BurntIn, GBurntIn;

#ifdef	JNIRUN
	long		FP;
#endif
	
	

	ShedFile	= NULL;
	SaveModelF	= NULL;
			
	CRates	=	CreatRates(Opt);
	NRates	=	CreatRates(Opt);

	Shed = CreatSchedule(Opt, CRates->RS);
	
	if(Opt->ModelType == MT_FATTAIL)
	{
		InitFatTailRates(Opt, Trees, CRates);
		InitFattailFile(Opt, Trees);
	}

	if(UseNonParametricMethods(Opt) == TRUE)
		InitVarRatesFiles(Opt, Trees, CRates);
		
	if(Opt->RJDummy == TRUE)
		InitRJDummyFile(Opt);
	
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

	//	STest(Opt, Trees, CRates);

//	if(Opt->LoadModels == TRUE)
//		TestModelFile(Opt, Trees, CRates);

	InitMCMC(Opt, Trees, CRates);
	
//	SetVarRatesFromStr("64680000	456.345426	338.643858	16	2.631978	0.003743	1.100000	11	3.101924	23747469	Kappa	20	4.834260	23802931	Branch	23	3.525218	23803910	Branch	15	633.971682	23803921	Kappa	164	9.358599	64526669	Kappa	84	3.034389	64677737	Node	431	6.568727	64678104	Node	181	5.733636	64678397	Node	244	5.479540	64679268	Branch	26	1.634304	64679343	Node	100	1.204039	64679581	Node	233	7.079552	64679638	Branch	59	6.473202	64679780	Node	195	5.609137	64679832	Branch	325	9.689859	64679865	Branch	45	0.975597	64679910	Kappa	\0", CRates, Opt);
//	SetVarRatesFromStr("1286940000	462.782556	341.393881	32	3.408467	0.003828	1.100000	164	19.825595	1285507302	Kappa	10	3.641360	1285543613	Kappa	16	14.043385	1285676718	Node	20	6.311118	1285677152	Branch	23	61.092926	1285677605	Branch	15	24.752145	1285677829	Kappa	11	16.287350	1285834783	Kappa	12	28.561212	1285919910	Node	14	33.714321	1286035856	Branch	19	0.127157	1286248466	Branch	13	44.946621	1286777888	Branch	395	11.139661	1286922071	Branch	161	0.396306	1286931542	Kappa	87	5.562612	1286934136	Branch	433	10.474905	1286935635	Node	120	6.248020	1286938212	Branch	284	3.690947	1286938845	Branch	374	3.322525	1286939083	Kappa	116	3.098829	1286939439	Kappa	221	2.047578	1286939463	Node	256	20.416036	1286939504	Node	107	0.358568	1286939650	Kappa	136	0.897210	1286939678	Branch	377	1.832391	1286939727	Kappa	398	0.119313	1286939755	Branch	57	16.939028	1286939782	Branch	287	2.412661	1286939786	Node	194	7.603891	1286939858	Branch	125	0.126994	1286939920	Kappa	220	8.858728	1286939933	Branch	320	7.863254	1286939974	Node	208	0.625917	1286939977	Kappa", CRates, Opt);
//	LoadGeoData("11807000	-60.858809	1.663948	0.317928	-1.808946	-1.155580	-1.755331	-1.377126	1.149021	1.370203	-1.076007	1.017187	2.139111	0.541804	0.251747	0.508075	0.182317	", Opt, Trees, CRates);
	
	CRates->Lh	=	Likelihood(CRates, Trees, Opt);
	CalcPriors(CRates, Opt);


	if(Opt->UseSchedule == TRUE)
		ShedFile = SetScheduleFile(Opt, Shed);
	
	if(Opt->SaveModels == TRUE)
		SaveModelF = InitSaveModelFile(Opt->SaveModelsFN, Opt, Trees, CRates);

	StoneF = NULL;
	if(Opt->Stones != NULL)
		StoneF = CreatStoneOuput(Opt);
	
	GBurntIn = BurntIn = FALSE;
	if(Opt->BurnIn == 0)
		BurntIn = TRUE;

//	MCMCTest(Opt, Trees, CRates);

	fflush(stdout);
	StartT = GetSeconds();	
	for(Itters=1;;Itters++)
	{ 
 		CopyRates(NRates, CRates, Opt);

		MutateRates(Opt, NRates, Shed, Itters);

		if(Opt->NodeData == TRUE)
			SetTreeAsData(Opt, Trees, NRates->TreeNo);

		if(!(Shed->Op == SFATTAILANS && Opt->Model == M_GEO))
			NRates->Lh = Likelihood(NRates, Trees, Opt);
	
		if(NRates->Lh == ERRLH)
			Itters--;
		else
		{
			CalcPriors(NRates, Opt);
			
			if(MCMCAccept(Itters, Opt, Trees, Shed, CRates, NRates) == TRUE)
			{
				Swap((void**)&NRates, (void**)&CRates);
				UpDateShedAcc(TRUE, Shed);
			}
			else
				UpDateShedAcc(FALSE, Shed);

			if( (Itters % Opt->Sample) == 0 && 
				BurntIn == TRUE &&
				StonesStarted(Opt->Stones, Itters) == FALSE)
			{
				UpDateHMean(Opt, CRates);
				CRates->Lh = Likelihood(CRates, Trees, Opt);

				#ifndef JNIRUN
					PrintMCMCSample(Itters, Shed, Opt, CRates, stdout);
					fflush(stdout);
				#endif

				PrintMCMCSample(Itters, Shed, Opt, CRates, Opt->LogFile);
				fflush(Opt->LogFile);

				if(Opt->UseSchedule == TRUE)
					PrintShed(Opt, Shed, ShedFile);

				if(UseNonParametricMethods(Opt) == TRUE)
					PrintVarRatesOutput(Opt, Trees, CRates, Itters);

				if(Opt->ModelType == MT_FATTAIL)
					OutputFatTail(Itters, Opt, Trees, CRates);

				if(Opt->RJDummy == TRUE)
					PrintRJDummy(Itters, Opt, Trees, CRates);

				#ifdef JNIRUN
					fgets(Opt->LogFileBuffer, LOGFILEBUFFERSIZE, Opt->LogFileRead);
					AddResults(Env, Obj, Opt);
				#endif
				
			#ifdef JNIRUN
				SetProgress(Env, Obj, Itters);
			#endif

				if(Opt->SaveModels == TRUE)
					SaveModelFile(SaveModelF, Opt, Trees, CRates);
			}

			if(Itters % Opt->Sample == 0)
			{
				if(StonesStarted(Opt->Stones, Itters) == FALSE)
				{
					UpDateSchedule(Opt, Shed, CRates->RS);
					BlankSchedule(Shed);
				}
			}

			if(ExitMCMC(Opt, Itters) == TRUE)
			{

				if( (Opt->UseEqualTrees == FALSE) || 
					(CRates->TreeNo == Trees->NoOfTrees - 1))
				{	
					EndT = GetSeconds();
					printf("Sec:\t%f\n", EndT - StartT);

					FreeRates(CRates, Trees);
					FreeRates(NRates, Trees);

					FreeeSchedule(Shed);

					if(StoneF != NULL)
						fclose(StoneF);
									
					if(Opt->UseSchedule == TRUE)
						fclose(ShedFile);

					if(UseNonParametricMethods(Opt) == TRUE)
						FinishVarRatesFiles(Opt);

					if(SaveModelF != NULL)
						fclose(SaveModelF);
					return;
				}

				if(GBurntIn == TRUE)
				{
					CRates->TreeNo++;
					CRates->Lh = Likelihood(CRates, Trees, Opt);
					Itters = 0;
					BurntIn = FALSE;
					BlankSchedule(Shed);
					Shed->GNoAcc = Shed->GNoTried = 0;
				}
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
			
			if(Opt->UseEqualTrees == TRUE)
			{
				if((Itters == Opt->ETreeBI) && (GBurntIn == TRUE))
					BurntIn = TRUE;
			}

			if(Itters == Opt->BurnIn)
			{
				if((Opt->UseEqualTrees == TRUE) && (GBurntIn == FALSE))
				{
					GBurntIn = TRUE;
					Itters = 1;
					BurntIn = FALSE;
				}
				else
					BurntIn = TRUE;
			}

			if(Opt->Stones != NULL)
				StoneItter(Opt->Stones, Itters, CRates->Lh, StoneF);
		}
	}
}

