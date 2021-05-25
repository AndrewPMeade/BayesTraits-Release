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
				case PDIST_EXP:
					fprintf(Str, "%s - Mean\t", Prior->Name);
				break;

				case PDIST_GAMMA:
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

	if(Type == VR_NODE)
		Prior = GetPriorFromName("VRNode", Rates->Priors, Rates->NoPriors);

	if(Type == VR_BL)
		Prior = GetPriorFromName("VRBL", Rates->Priors, Rates->NoPriors);

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

	Buffer = (char*)SMalloc(sizeof(char) * BUFFERSIZE);

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

void	ReSetAccFlags(RATES *Rates)
{
	Rates->AutoAccept = FALSE;
	Rates->CalcLh = TRUE;
}

int		MCMCAccept(long long Itters, OPTIONS *Opt, TREES *Trees, SCHEDULE* Shed, RATES *CRates, RATES *NRates)
{
	double Heat;

	if(NRates->AutoAccept == TRUE)
	{
		ReSetAccFlags(NRates);
		return TRUE;
	}

	ReSetAccFlags(NRates);

	if(Shed->Op == S_FAT_TAILANS && Opt->Model == M_FATTAIL)
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

	X = Rates->LocalTransforms[0]->Scale ;
	Rates->Lh = Likelihood(Rates, Trees, Opt);

	printf("%f\t%f\n", X, Rates->Lh);

	Rates->LocalTransforms[0]->Scale *= 10;
	X = Rates->LocalTransforms[0]->Scale ;
	Rates->Lh = Likelihood(Rates, Trees, Opt);

	printf("%f\t%f\n", X, Rates->Lh);

//	exit(0);
	printf("====\n");

	for(X=0;X<10;X+=0.01)
	{
		Rates->LocalTransforms[0]->Scale = X;
		Rates->Lh = Likelihood(Rates, Trees, Opt);

		printf("%f\t%f\n", X, Rates->Lh);
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
//	SetShedOpFreq(Shed, S_VARRATES_ADD_REMOVE, 0);

	if(Opt->ModelType == MT_FATTAIL)
	{
//		InitFatTailRates(Opt, Trees, CRates);
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



	InitMCMC(Opt, Trees, CRates);
	
//	SetVarRatesFromStr(CRates, Opt, "75700000	-25329.349319	-26266.275528	222	NA	NA	NA	767	28.018929	75678795	Branch	141	0.201699	26323	Node	882	0.101275	6282	Node	893	0.051791	18063	Node	258	47.255215	54420207	Branch	538	21.539686	66939157	Branch	542	0.485291	449747	Node	144	0.117957	21507362	Node	328	36.977394	75023473	Branch	334	62.452925	269255	Branch	39	7.786450	75589854	Branch	1009	24.744504	12975034	Node	1011	0.010496	532029	Node	175	0.479850	74147515	Node	181	21.107805	49981506	Node	183	28.571179	52412343	Branch	362	26.195650	301487	Branch	628	32.385695	94814	Branch	1066	25.764678	522667	Node	99	0.040256	67831	Node	973	4.904898	39677542	Node	634	0.016185	48093	Node	406	10.151407	58892704	Node	900	47.658934	633206	Branch	67	0.028830	3845	Node	646	24.382860	74682733	Node	744	0.029419	261935	Node	71	1.886235	75697006	Branch	1042	41.771942	49576028	Node	560	0.047869	66819297	Node	1083	22.145681	68003090	Branch	1044	0.064360	74920461	Node	212	0.003491	510	Node	854	0.010728	39958	Node	1085	1.798938	75693571	Node	930	42.069486	6814178	Node	993	13.488178	70706241	Node	832	0.178693	45996959	Node	1089	20.311537	37520019	Node	932	0.001889	31698	Node	857	22.294719	63487011	Branch	424	0.067895	57455011	Node	147	2.219527	74893329	Node	903	4.597882	149277	Node	186	0.033635	673646	Node	769	0.013309	256155	Node	1051	20.623458	74824431	Node	907	8.568778	69045163	Node	430	0.031741	105257	Node	938	11.358849	68522802	Node	619	54.905957	6946043	Branch	963	17.549159	66603894	Node	1001	5.470720	75594619	Node	940	41.015868	67819890	Node	28	0.001784	1392	Node	1074	0.032547	16492	Node	523	0.009988	18678	Node	724	0.011967	37845	Node	885	18.267126	72236388	Branch	770	0.001709	75517600	Branch	192	8.661974	52792511	Node	756	9.980791	75689499	Node	397	7.183443	75697300	Branch	745	0.504481	75698735	Branch	224	21.077814	69879740	Node	772	19.182616	67996195	Node	531	0.007177	60556	Node	976	48.420457	11289	Branch	804	18.746916	73520036	Branch	561	25.157578	73187471	Node	75	15.317894	75102238	Node	910	8.270923	75568134	Node	434	0.427806	75699562	Branch	401	7.689169	74829208	Node	436	12.076868	75280078	Node	889	75.028201	11252	Branch	746	21.434600	54968673	Branch	274	31.594792	58664443	Branch	649	10.281128	67746845	Node	1100	0.017043	69895954	Node	562	42.830923	74226322	Branch	253	0.159940	75388621	Node	473	0.000095	10989	Node	76	0.000063	3677	Node	944	41.564419	7039	Branch	388	39.550174	17681	Branch	83	16.233327	20066	Node	317	0.001480	36939	Node	117	18.630377	98934	Node	403	0.000597	504956	Node	952	0.009463	540623	Node	1025	0.000367	2723928	Node	157	18.605631	8817857	Branch	662	0.000781	21519283	Node	625	0.002288	24645866	Node	808	0.002589	30968898	Node	913	0.003971	40488032	Node	871	19.021206	60901472	Node	271	0.005488	65797266	Node	102	15.290207	70489564	Branch	425	10.684775	73150911	Branch	631	7.648581	74520904	Node	793	20.285042	74561542	Branch	47	20.146512	74839832	Branch	1061	25.105566	75251044	Node	741	44.854853	75361700	Branch	203	9.415620	75606690	Node	398	13.448019	75637763	Branch	703	4.115666	75698649	Branch	868	1.293399	75699229	Node	27	7.207205	75693272	Branch	443	0.948829	75693528	Branch	112	0.923586	75688958	Branch	676	0.011708	75673774	Branch	1021	9.409622	75592846	Branch	931	63.533070	8536	Branch	968	58.262087	11676	Branch	1043	30.124004	11956	Branch	89	39.745051	20224	Branch	196	18.614344	25468	Branch	668	46.032542	32795	Branch	290	28.727595	40586	Branch	1010	82.312422	46805	Branch	194	86.115653	49780	Branch	402	16.309153	71995	Branch	666	51.916344	100475	Branch	437	24.838107	103424	Branch	1067	27.359106	201232	Branch	989	19.261839	249240	Branch	773	42.494459	263338	Branch	760	10.915924	267244	Branch	349	33.989405	275283	Branch	185	72.636648	319078	Branch	378	16.582881	375862	Branch	121	14.151505	386113	Branch	1070	18.440601	425314	Branch	977	24.882447	438789	Branch	575	49.516443	456979	Branch	33	28.886848	521040	Branch	371	15.983110	629464	Branch	207	56.040881	632625	Branch	647	38.444575	21931575	Branch	115	27.178483	26216750	Branch	1031	13.464041	31946170	Branch	347	37.697294	32964306	Branch	327	16.671960	34579242	Branch	904	38.843757	35640061	Branch	928	77.560695	36053666	Branch	853	15.258216	45564475	Branch	841	21.595969	45954932	Branch	313	20.519614	46098889	Branch	524	11.763061	47448668	Branch	847	15.179667	48942426	Branch	201	7.962552	52449955	Branch	396	25.389174	53383569	Branch	747	25.745199	54965315	Branch	611	9.581275	57148577	Branch	851	10.830167	59995183	Branch	261	46.832825	63282579	Branch	229	51.965039	66072937	Branch	873	17.014052	67789711	Branch	962	12.574261	68248113	Branch	150	21.888784	70584135	Branch	152	28.371709	70587114	Branch	151	28.933005	70600285	Branch	1052	24.612100	70673820	Branch	766	19.346624	71618518	Branch	872	17.253181	71877451	Branch	1063	11.449751	72368035	Branch	1056	28.788810	72618658	Branch	807	10.518975	72910478	Branch	701	33.979832	73411966	Branch	606	18.246330	73663535	Branch	161	15.840527	73743927	Branch	174	21.587544	73957958	Branch	933	28.385440	74033665	Branch	751	73.890026	74188881	Branch	80	34.074163	74642605	Branch	265	8.844992	74859763	Branch	1050	25.586610	75085165	Branch	344	9.331200	75191508	Branch	925	4.778361	75250946	Branch	1034	24.701823	75317275	Branch	553	11.827901	75450985	Branch	709	20.062543	75456174	Branch	64	12.759153	75477128	Branch	365	20.863289	75500375	Branch	305	48.364570	75558742	Branch	475	10.811417	75593127	Branch	908	11.976882	75600519	Branch	997	21.470344	75622757	Branch	734	29.170441	75628633	Branch	1000	21.911750	75632126	Branch	307	14.378824	75637594	Branch	119	40.607133	75640639	Branch	12	11.568516	75641387	Branch	530	6.973247	75646830	Branch	351	7.520957	75653281	Branch	284	22.527657	75665207	Branch	820	12.772480	75669324	Branch	57	18.352402	75672604	Branch	332	1.486927	75682541	Branch	413	26.113935	75683419	Branch	226	9.552566	75683550	Branch	505	12.668445	75684119	Branch	85	3.292086	75687805	Branch	336	4.171417	75691870	Branch	170	0.171682	75692351	Branch	1004	1.771057	75693062	Branch	245	1.513530	75693080	Branch	826	8.635699	75693620	Branch	927	0.042436	75694440	Branch	429	3.838823	75694539	Branch	692	9.802157	75694816	Branch	992	0.033990	75696526	Branch	535	7.599543	75697350	Branch	727	0.208954	75697558	Branch	426	0.300151	75697717	Branch	435	18.655630	75698075	Branch	462	1.351968	75698692	Branch	479	7.964949	75699176	Branch	815	1.910000	75699886	Branch	");
//	LoadGeoData(Opt, Trees, CRates, "75700000	-25329.349319	2.000000	367339.503201	-11.966487	-34.832468	-31.680364	-50.591629	-14.894821	-14.549716	-17.363838	-34.398062	-17.183122	-28.843546	-18.495423	-32.820251	-17.063146	-29.431175	-10.920504	8.497495	-5.620074	7.450334	-7.575806	3.237263	-21.638991	-26.765918	-41.554646	-42.630088	-1.210231	-43.417107	8.112981	60.056737	-1.305884	-43.766598	-1.563355	-43.861783	-1.794825	-43.993427	4.903950	-48.014199	0.145436	-39.666893	-13.147223	38.042379	-10.565881	45.982037	-17.096365	41.104325	-16.639788	38.917105	-13.754015	41.446705	75.830644	25.973014	2.109328	50.261060	6.676783	43.463927	4.369688	40.696056	-22.478645	49.600469	-45.176985	37.275525	52.664056	16.112719	100.536791	40.787011	99.794758	42.090567	101.580792	42.787817	101.611074	44.950771	101.877625	44.257520	112.083166	46.229684	112.047165	46.256672	112.959350	39.922541	97.320094	43.725675	-76.649289	61.952060	98.893859	41.867497	96.096063	42.979057	97.473731	43.042624	98.448718	42.721499	91.591709	40.232995	-53.826057	38.030534	-55.308633	37.288041	-53.288754	37.955449	11.398206	42.730520	-53.859283	37.910290	-53.637429	38.245897	-50.700196	40.063620	-59.399504	41.754889	-58.815387	41.842993	-58.145462	37.056335	10.321777	65.812341	-57.348201	38.305689	-57.997357	37.850426	-60.098264	39.454753	-79.131369	55.645064	-80.001567	55.911528	-1.693837	-46.792315	112.210512	37.307548	113.280348	35.887342	114.820630	37.544428	116.073432	38.445852	97.356718	48.719555	94.851225	49.027323	91.858171	47.631746	90.791439	46.558160	-74.325286	58.693924	114.822217	38.827681	116.683976	39.202400	116.690165	39.084230	110.105131	41.565832	110.042419	43.224378	96.985912	38.993330	97.637212	38.892246	97.504206	38.429187	96.982125	39.004340	122.304748	43.967894	-83.126838	54.253933	-77.836567	53.727468	-78.076771	53.557331	-77.511315	53.918300	-77.285782	54.106066	-75.717108	55.511833	-77.256708	55.466778	-74.564625	55.762200	-74.215418	49.390007	-80.548558	38.222256	92.128707	49.296637	-82.632280	37.208792	-79.333348	53.979041	-76.362939	57.252215	-76.827216	57.175850	-77.378084	57.246288	-76.414716	57.835912	-77.293367	56.993048	-77.619681	56.587814	-79.170231	56.293738	-79.879419	56.032538	-72.500513	50.233776	-70.539655	52.887678	-70.711602	41.281668	-71.212791	44.890312	-82.834493	43.198566	-77.436715	45.668142	-77.262231	49.252979	-80.098328	50.633686	-79.094785	53.416416	-81.738816	56.189524	-81.932172	51.186386	-80.233341	52.963971	-80.440998	54.041896	-46.673149	35.833202	-77.407934	53.757753	-67.587705	53.256071	-72.105654	48.691608	-66.981324	53.529029	-78.207535	56.491320	-81.841422	56.380246	-80.321284	55.096332	111.655151	58.263490	124.357840	57.884495	-41.021333	33.085972	-37.729671	27.009138	-37.148395	28.130609	-56.986733	-42.294347	-34.537144	31.639793	1.224689	51.286085	18.859864	30.937515	16.885154	34.510227	15.463568	31.389241	-30.376735	32.963756	-58.214114	33.738282	-29.760465	33.296728	-32.161730	37.245020	-36.946475	41.199164	-41.265729	40.532119	-49.121297	38.574055	-60.010672	37.745820	-59.675391	38.420003	-63.277565	37.296364	-64.209205	36.502305	-61.813869	37.598419	-64.643790	40.264940	3.915117	38.243521	4.021557	43.750331	96.441772	44.846500	104.919866	40.680138	101.172727	42.214089	100.576682	43.760034	101.830427	42.450421	119.294454	41.423074	117.242025	40.570243	111.074311	37.622910	113.034964	40.309195	101.697270	42.296695	101.580999	47.562744	111.102081	49.862460	113.492967	50.580654	-58.508275	49.069658	-59.553757	50.218737	-53.660992	52.356537	-57.834096	51.810250	-66.760270	57.445022	-75.838406	54.170035	-77.038816	53.727854	-75.082572	55.702270	-72.925881	57.675390	-81.655410	54.105246	-79.397584	54.682901	-75.618302	58.229480	-77.712602	55.585651	-76.990762	55.045742	-57.604090	-39.719341	-84.888882	57.797327	-83.155586	57.764851	-85.111631	60.550376	-71.462754	56.476253	112.673960	56.034879	-79.498899	55.696502	-79.093222	55.721149	-53.312386	62.152949	-68.000512	43.866354	104.637775	37.233797	-51.935783	41.221998	-67.335919	50.829203	-9.486026	37.765055	-70.092375	55.523760	-65.696939	58.947882	-79.304558	52.807808	-74.075354	58.136986	112.606018	52.050264	-77.896987	60.087293	-76.535378	59.148591	-75.695174	57.518138	-76.797768	59.457049	-78.046053	56.074981	-73.722165	56.168276	-18.408459	-24.055849	-20.356592	-24.460835	-18.872660	-42.878458	-17.153380	-28.045578	-17.058799	-25.712601	-12.960607	-26.319917	-9.853430	-26.732357	-3.432882	-19.105040	-1.915063	-26.988180	-4.256344	-28.733140	-20.024007	-28.167485	-8.135037	-25.702965	-22.143372	-39.838242	-7.616756	-20.676046	-9.432752	-23.769695	-11.713247	-28.997518	-28.664828	-31.876483	-28.573401	-31.755145	-20.338380	-39.035595	-8.852277	-29.632994	-12.243134	-22.235232	-11.893693	-13.585732	-11.157527	-24.989138	-6.148198	-26.112610	-5.561643	-20.974488	-9.964653	-24.803887	-4.408162	-19.010210	-12.528053	-33.850905	-10.973168	-16.805331	-5.767620	-18.475309	-4.827520	-17.545198	5.469247	-32.601471	37.155693	-25.364254	43.531430	-26.372238	55.351623	-27.393913	37.682002	-27.520276	95.669045	24.709274	110.019119	27.612379	107.579977	27.099747	109.820096	30.211738	108.540077	29.915700	96.370835	33.773841	5.525366	31.786516	-0.529250	42.001472	-2.127285	40.231789	-6.773517	39.835819	-37.397284	33.035595	-46.272011	40.218618	-51.330004	35.589214	-27.078326	-41.455939	-18.973377	-37.327730	-46.928693	-45.975375	-52.953209	33.930568	-55.417583	32.934493	0.053339	-40.383454	-8.479502	-36.201567	-52.546061	35.367926	-54.113800	35.852518	-53.428586	35.027653	-51.426347	36.599968	-49.910743	37.154148	-51.097982	35.187667	-5.383041	36.636483	1.857466	32.568466	-39.858586	43.548436	1.945926	34.027669	-1.239908	-0.414523	2.114187	32.283842	-2.176282	23.519009	2.483826	24.212390	7.621965	27.949420	-0.949419	28.100219	8.237226	18.660332	0.266109	21.387915	-53.806291	33.985337	5.130929	36.674746	6.176094	40.628267	6.107674	38.638254	8.668996	24.668771	122.133760	12.956887	113.097810	15.960680	112.981249	16.155289	-44.293436	-41.836583	-44.228573	-47.151202	111.889624	-59.219866	-44.183693	-47.080414	-44.810938	-47.536566	-43.632425	-47.171799	-48.743292	-45.305340	-42.958457	-46.168529	-43.364961	-46.834789	-42.503549	-46.906196	-36.802973	-44.220773	-24.289963	-60.379913	-14.992964	-40.943826	-48.249652	-46.905178	-49.383993	-42.014233	-44.415345	-43.077958	-39.548274	-38.625359	-49.370606	-42.084597	-60.257643	-39.823487	-56.366086	-36.958688	-52.576077	-36.297117	-50.053386	-38.910021	-30.510564	-48.507695	-32.619352	-42.733866	-25.441108	-43.276965	-31.933523	16.454702	-33.346284	8.553583	-26.574351	-0.663951	-20.520709	7.016264	-20.622103	0.986985	-4.376084	10.549613	-27.324740	-10.875615	-21.724848	-8.132927	-8.451995	-13.198500	-2.482659	-14.290264	17.841308	-5.539626	-8.981808	-19.131686	-23.093510	25.494583	-18.464983	-25.015303	0.171090	-51.257190	-10.369530	-33.462454	-31.109774	-19.935822	-36.148105	-40.092781	-46.688447	-45.128292	48.461939	-25.523244	-45.086609	-45.296066	-53.953886	-44.782219	-45.687911	-44.984315	-23.290894	-5.236036	1.523191	52.004509	2.921347	54.228554	7.393125	38.582686	0.836906	24.338394	-18.788030	-34.778513	8.333263	39.143522	9.932444	39.884375	-10.640852	-9.955495	10.114901	40.642778	7.638662	39.847036	8.825841	40.282126	7.968140	37.489056	9.322671	40.899773	11.463141	43.116017	120.982589	54.212000	108.026542	24.823103	108.874160	25.507962	108.827997	25.512571	111.135749	37.836944	109.343727	25.208635	-37.881854	27.012763	-43.402327	33.979795	3.095330	41.875286	2.173019	44.052404	8.276709	39.343465	12.073917	47.267555	-34.820366	-37.843360	139.803304	40.934685	-7.049465	36.680287	-1.567971	38.043801	-37.679168	29.667850	-35.903309	36.792127	-39.644264	-49.148062	-45.141251	-45.991883	-45.606913	-46.252194	112.038867	56.487227	119.783186	47.760505	139.331343	52.420961	136.585334	64.864409	118.858492	57.702782	101.572107	49.455097	-21.338931	79.998780	8.154901	84.154227	40.970606	71.851364	68.000673	44.295591	-73.922900	48.815760	-74.129669	43.641871	-72.871835	46.184249	-72.652646	44.789020	-75.609278	56.689137	-80.844304	49.055742	-82.917675	52.532687	-83.094989	56.411251	-78.926864	56.791785	114.651335	48.500295	119.092106	29.507015	121.491412	46.529085	122.369852	45.633122	95.413293	44.591102	95.676936	44.384081	94.228164	42.174834	93.309411	42.089548	93.680867	42.918511	86.865884	45.426824	88.278990	42.698829	88.542575	47.196005	119.139430	43.367630	-13.500026	68.450372	3.031432	52.515496	21.342900	49.025717	123.602832	43.410867	90.499830	45.903916	-46.508487	-42.070584	-49.761671	-44.813081	88.105211	40.862019	123.750721	43.436761	124.090489	43.901492	123.957069	43.526450	94.462272	42.317797	95.158894	42.388217	94.224054	41.901858	95.724687	42.297100	97.054620	38.458385	-76.706838	57.458395	-76.869957	58.259136	91.844008	40.265449	92.197759	42.659528	91.278514	40.733086	90.887550	40.276997	90.858437	40.837015	90.411529	41.240474	102.172702	41.695366	102.199658	41.453062	101.550165	41.460071	101.476133	41.473803	101.000571	40.642326	96.817689	40.582975	-77.929236	40.247837	123.942646	43.252762	124.460832	42.334835	123.495950	42.019568	123.098701	42.145741	123.049414	42.093765	121.601714	44.093952	121.999979	44.358508	124.170361	43.160498	123.986661	43.691464	120.045289	42.585316	121.250665	43.413094	96.082143	39.256925	91.772682	40.044395	112.555881	42.914724	111.548793	42.697522	-42.524259	-46.532721	-45.373816	-47.518020	-47.646993	-44.807806	112.106068	42.809316	113.052596	42.703408	122.961113	43.685607	146.542508	85.888787	-136.597774	84.950290	10.266989	70.535991	59.892392	71.304676	80.159471	52.507207	77.230418	36.102159	120.396973	44.593894	120.209830	44.111244	119.992901	44.128938	120.463170	44.000480	120.604564	44.131535	121.101678	43.883663	121.382544	43.619456	122.513231	44.469273	122.516796	44.464924	120.084514	44.252449	119.776277	44.275479	120.442230	44.096508	118.654520	42.792969	120.841793	44.040424	120.058642	44.422926	120.653530	43.946678	121.665800	43.641913	122.024263	43.800942	121.133447	43.390074	119.453087	41.973241	118.089800	42.750321	116.858947	43.339247	117.338876	45.647777	113.255396	44.978827	88.440148	46.352193	120.680018	43.979685	119.346216	45.446157	121.072785	44.135419	122.182566	44.052330	122.044683	44.105390	120.890980	44.803906	120.470959	44.526653	104.743876	44.708480	104.046596	44.539473	101.742378	44.021253	-1.074514	53.616974	0.206498	46.803549	-6.502557	43.023278	-3.871977	48.882314	-67.363770	45.775371	-65.173145	43.987312	-65.227432	43.764045	");
	
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

	fflush(stdout);
	fflush(Opt->LogFile);


	StartT = GetSeconds();	
	for(Itters=1;;Itters++)
	{ 
 		CopyRates(NRates, CRates, Opt);

		MutateRates(Opt, NRates, Shed, Itters);

		if(Opt->NodeData == TRUE)
			SetTreeAsData(Opt, Trees, NRates->TreeNo);
		
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
					SetProgress(Env, Obj, Itters);
				#endif

				if(Opt->SaveModels == TRUE)
					SaveModelFile(SaveModelF, Opt, Trees, CRates);

				if(Opt->SaveTrees == TRUE)
					OutputTree(Opt, Trees, CRates, Itters, Opt->OutTrees);
			}

			if(Itters % Opt->Sample == 0)
			{
				// The schedule should be updated even when stones is running
				UpDateSchedule(Opt, Shed, CRates->RS);
				BlankSchedule(Shed);
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

