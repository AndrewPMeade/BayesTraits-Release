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

	if(Rates->NoEstData > 0 && Opt->DataType == CONTINUOUS)
		SetEstDataFromPrior(Rates);
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


	for(TNo=0;TNo<1000000;TNo++)
	{
		for(RIndex=0;RIndex<Rates->NoOfRates;RIndex++)
			Rates->Rates[RIndex] = RandFromPriorPosition(RIndex, Opt, Trees, Rates, RNG);

		SetDefMCMCParameters(Opt, Trees, Rates, RNG); 

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
	
//	printf("I:\t%lld\t%f\t%f\t", Itters, CRates->Lh + CRates->LhPrior, NRates->Lh + NRates->LhPrior);

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

void	MCMCTest(OPTIONS *Opt, TREES *Trees, RATES *Rates, SCHEDULE *Shed)
{
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
	double		StartT, EndT, TLh;
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
		
//	SetVarRatesFromStr(CRates, Opt, "64500000	-25475.320861	-26309.848926	199	NA	NA	NA	141	0.174389	3534	Node	442	0.352251	30941	Node	882	0.099346	722	Node	495	2.338966	53424250	Node	258	32.814691	5063946	Branch	538	0.034658	816289	Node	326	34.113059	40308102	Branch	334	39.578833	1181929	Branch	162	0.057030	12348	Node	1011	0.008253	142831	Node	181	20.988987	241571	Node	598	18.310070	2154078	Node	362	27.582164	1255752	Branch	206	32.883798	1033361	Branch	628	17.810821	1638	Branch	1066	8.128930	381439	Node	99	0.028381	5758	Node	545	28.746761	2095639	Branch	973	0.267223	52017370	Node	894	0.039756	12376244	Node	634	0.020298	206	Node	814	0.185443	2417834	Node	900	41.395113	35393690	Branch	67	0.040662	10191	Node	646	16.643959	1807632	Node	1042	18.918466	1206198	Node	951	0.035868	643862	Node	1083	13.268458	17615081	Node	212	0.010423	77727	Node	993	7.592800	51068404	Node	522	23.508201	52894758	Branch	830	21.952600	64491710	Branch	145	1.377479	64499192	Node	147	0.390768	64453437	Node	1089	11.379734	59251371	Node	857	0.018560	228445	Node	184	40.374413	8933620	Node	753	0.047487	1009909	Node	149	18.407831	58614576	Node	934	0.019203	42595778	Node	903	17.595728	3528309	Node	186	0.006977	2295164	Node	550	25.805436	1839789	Branch	769	0.007276	32091882	Node	907	11.156087	63492045	Node	619	48.271989	5268445	Branch	963	30.595260	58560691	Node	430	0.181241	64234561	Node	268	5.189382	64486351	Branch	1055	5.132223	62730709	Node	940	44.333492	62994154	Node	28	0.001132	12175	Node	724	0.010501	26721	Node	523	0.013930	1574952	Node	153	0.029219	59477815	Node	885	8.794267	61092622	Branch	304	11.819226	64297601	Node	529	3.041424	64482703	Node	772	25.549703	62205956	Node	224	9.188139	62140416	Node	531	0.001787	23261457	Node	976	48.808382	658418	Branch	1098	0.009902	201422	Node	292	0.047456	63542302	Node	13	1.494394	64495836	Branch	889	42.681207	698823	Branch	586	22.140366	2147035	Branch	739	9.323217	27085034	Branch	649	16.228929	47972211	Node	562	90.308495	59153607	Node	274	55.789812	62471144	Branch	825	6.239178	64480428	Node	79	19.361660	59210902	Node	808	0.020984	50672352	Node	944	37.873931	1915	Node	388	15.662692	4850	Branch	117	34.318299	18570	Node	83	37.514302	78808	Branch	1025	0.000197	138049	Node	76	0.001218	292943	Node	317	0.000517	435693	Node	741	0.000233	865497	Node	662	0.001399	1236786	Node	591	28.034581	2035871	Node	473	0.005704	32711324	Node	871	11.372988	36352275	Node	102	20.344051	58263136	Branch	157	46.825068	59599991	Branch	409	15.710331	59930029	Node	551	25.251304	61271690	Node	425	7.485539	62396365	Branch	398	10.886814	62405059	Branch	47	29.113278	64387676	Branch	793	10.518318	64441758	Node	460	3.264524	64485347	Node	554	0.089881	64491973	Node	293	10.044656	64496043	Node	197	0.071176	64497503	Branch	1061	6.181445	64498009	Node	845	0.758058	64490003	Branch	968	22.010831	1080	Branch	668	19.009054	2397	Branch	841	32.995659	5264	Branch	647	51.571943	6609	Branch	1043	51.401985	8801	Branch	207	42.663487	10250	Branch	1067	88.130415	13508	Branch	1010	56.188932	15269	Branch	396	26.932611	18355	Branch	575	49.767532	22253	Branch	773	64.548493	26867	Branch	89	19.462512	35074	Branch	1070	36.370604	110458	Branch	989	21.304852	141201	Branch	437	27.115959	216255	Branch	402	77.595733	362777	Branch	760	21.221893	1010018	Branch	378	16.139210	1221060	Branch	194	33.606830	1278367	Branch	313	33.557626	1338893	Branch	349	12.931661	1376510	Branch	347	12.091265	1377509	Branch	524	11.440514	1570432	Branch	558	24.908058	2042764	Branch	904	31.058639	4934168	Branch	928	38.542141	11829400	Branch	1056	27.956085	22454100	Branch	196	24.514365	23375840	Branch	1031	24.515751	25949474	Branch	229	60.363387	29714254	Branch	327	6.546362	29744235	Branch	931	43.353572	31804863	Branch	371	12.812288	41946283	Branch	666	16.658449	43785875	Branch	413	25.386549	46641117	Branch	1034	29.213369	46810378	Branch	407	11.372768	47174084	Branch	1000	10.544998	54343624	Branch	261	45.313607	57189403	Branch	997	17.372541	57535275	Branch	757	25.779158	57566617	Branch	807	13.097001	58150305	Branch	766	30.617780	58665166	Branch	121	28.781403	58703236	Branch	1063	48.091761	59327480	Branch	1004	18.125855	59758424	Branch	820	30.398937	60084996	Branch	611	10.069319	60386496	Branch	977	16.692000	61146146	Branch	1097	17.207482	61424859	Branch	80	22.668619	61687048	Branch	853	8.789114	61994568	Branch	420	23.909111	62343111	Branch	847	22.581812	62774236	Branch	115	52.530891	62943720	Branch	479	13.995481	62951739	Branch	873	26.492851	62970008	Branch	1052	56.087974	63023943	Branch	174	9.874775	63579937	Branch	64	13.426032	63615072	Branch	553	8.097071	63828448	Branch	734	10.062592	63844330	Branch	805	14.507245	63867628	Branch	826	3.540788	63968452	Branch	290	28.485825	64058522	Branch	481	4.326672	64120578	Branch	482	13.152312	64185338	Branch	851	7.854744	64198183	Branch	872	16.603379	64201464	Branch	505	30.142163	64253815	Branch	484	12.060096	64268023	Branch	33	47.670075	64324560	Branch	701	28.808396	64343665	Branch	1021	6.385329	64373937	Branch	1006	8.551196	64377480	Branch	815	8.770776	64393391	Branch	365	6.768640	64424115	Branch	962	17.905695	64424552	Branch	344	8.584431	64430343	Branch	632	12.736681	64438151	Branch	454	6.329011	64446859	Branch	925	11.851370	64453779	Branch	193	1.400707	64481772	Branch	351	10.069170	64485805	Branch	265	11.302300	64489245	Branch	490	0.434992	64489418	Branch	185	2.718866	64489469	Branch	1081	0.611101	64490354	Branch	778	7.784806	64492022	Branch	307	11.995525	64493175	Branch	57	9.687019	64493620	Branch	526	1.117921	64494710	Branch	267	0.334349	64495212	Branch	801	7.105603	64495575	Branch	541	3.902364	64495787	Branch	446	0.742743	64496628	Branch	114	0.336432	64496827	Branch	543	1.404922	64499636	Branch	527	0.452940	64499674	Branch	");
//	LoadGeoData(Opt, Trees, CRates, "64500000	-25475.320861	2.000000	476894.245712	-12.241430	-54.372313	-29.374763	-53.392316	-31.334839	-10.315918	-23.530626	-61.735250	-28.114094	-53.694947	-23.820705	-54.168004	1.021989	-52.100548	-2.206263	12.535661	-1.608834	14.459407	-4.302661	12.545568	-4.420746	-48.072172	-31.724486	-47.024754	-0.891350	-44.566302	-54.293652	50.872424	-1.325818	-44.483307	-1.419007	-44.473742	-1.297270	-44.416480	-6.737613	-28.580830	-9.523549	-18.427069	7.672958	16.344592	9.863560	25.900947	7.867131	32.911087	-49.870367	53.433011	-21.674372	45.463312	117.931232	19.799804	11.464746	43.389569	5.988624	37.199876	15.336518	52.190656	-47.915134	53.296391	-60.684210	42.349857	83.204407	46.318610	98.360462	41.263718	98.817109	41.354998	101.575021	40.005429	109.921719	44.226090	108.741927	44.518287	112.080048	46.207823	112.072375	46.221157	114.570236	49.024798	109.198514	43.229984	-79.235454	56.788108	107.340814	43.545189	105.346246	42.504104	104.027005	42.338974	104.727786	42.758098	94.801157	44.492894	-53.241945	34.715980	-53.501609	35.290647	-54.719204	36.301472	11.285120	42.686324	-51.110808	37.233261	-53.983467	37.611575	-53.979383	38.005905	-58.507055	41.941458	-59.926560	41.488940	-64.559652	38.071482	126.487389	43.109325	-62.158308	37.018600	-63.059413	39.934744	-62.916286	40.106081	-78.270786	55.218528	-79.467812	55.309942	-4.355717	-22.328095	125.147170	44.356423	115.983483	39.526680	111.233224	37.767544	117.825627	41.688626	99.177716	39.893173	96.846830	42.405929	91.355656	46.384545	91.371140	44.028781	-76.034987	55.559357	116.282032	41.799918	115.362331	42.278848	115.628402	41.106625	106.225696	43.255798	107.786229	43.031594	104.715454	37.982993	102.754258	37.532511	102.417766	38.410511	101.382468	38.014201	124.455955	41.926551	122.761867	45.698680	-74.025639	68.820582	-78.040869	59.005257	-75.695214	58.422174	-76.009595	57.792190	-77.933444	56.469321	-72.184406	54.710149	-72.304535	54.529816	106.631373	47.746217	109.443041	56.350721	89.988622	44.961861	-79.813561	47.874283	-81.131554	51.902196	-77.101694	56.803861	-77.508471	56.727208	-76.784674	56.918542	-76.633419	57.741747	-77.442358	56.998816	-77.791264	57.146058	-78.754358	56.779501	-78.366194	58.392937	-77.303262	43.537763	-75.323990	57.301913	-74.006076	43.044292	-77.843347	48.489784	-79.570515	46.004061	-73.073633	44.770631	-76.224000	44.835011	-76.742855	48.539620	-80.970093	47.078827	-78.693868	49.268991	-84.015783	51.229034	-78.890305	48.063134	-79.963810	51.738603	-36.769053	50.494055	-51.568818	44.051820	-47.614942	44.528276	-67.171729	47.583657	-70.181403	45.740287	-34.084986	45.226924	-35.362445	43.796333	-60.357672	56.400175	115.800556	52.321031	126.273746	48.280901	-37.280848	50.921778	-35.319233	44.622354	-30.182431	41.410338	-49.136702	-44.538674	-31.975548	39.608901	17.914937	27.490817	20.999243	30.971254	18.330421	35.067413	20.442339	31.692450	-35.124875	38.116802	-52.173945	32.853610	-32.674496	36.064607	-20.701976	37.558906	-50.347474	36.090128	-53.402178	36.055808	-53.306514	34.855839	-59.804106	38.168699	-59.829250	38.419042	-56.588582	41.333520	-47.801353	44.502021	-50.659036	43.107019	8.191991	47.064086	14.032516	41.970274	19.124608	34.619535	105.589790	48.167338	112.218461	46.629641	104.861952	47.309797	108.530282	45.166519	104.800865	43.195163	107.235534	37.811602	110.519776	39.380642	114.350327	41.989118	114.193012	37.826993	122.109118	38.910279	116.384332	40.931882	118.203042	39.527329	119.897676	44.083767	-57.426244	47.803410	-57.765559	45.128504	-63.229957	49.931486	-61.286702	53.378746	-76.349080	45.084536	-81.424488	51.101030	-78.949211	51.912147	-77.048935	52.886454	-81.220909	49.436458	-77.950059	54.833814	-79.285414	54.734795	-83.008427	58.098587	-81.946175	55.939674	-78.969822	55.385550	-49.612422	-39.872386	-77.300310	56.363084	-89.001582	57.905822	-84.308980	61.763195	-87.976429	57.827100	119.900799	52.396508	-84.535714	61.266711	-92.688215	55.378788	-56.363768	53.113295	-62.267757	50.012233	104.498843	30.262673	-69.483208	49.150982	-64.436785	52.793069	-0.658522	33.921926	-66.619832	53.202263	-81.293864	50.810347	-81.438368	49.401424	-73.116979	52.934071	116.811287	50.308601	-75.293460	57.466395	-74.939518	58.138787	-78.018478	57.154111	-72.208795	57.925871	-73.166119	58.424077	-74.324713	55.932198	-14.729569	-42.096936	-6.628741	-38.164771	-22.715657	-45.995486	-12.189739	-26.234935	-6.553696	-17.081645	-9.740867	-25.492125	-7.360107	-26.532070	-11.962598	-32.834133	-11.907448	-34.556855	-22.563069	-32.232078	-10.117580	-35.252763	-22.123027	-36.051757	-22.507051	-46.914969	-21.219765	-36.508127	-25.446903	-32.600112	-38.652545	-39.595416	-28.577271	-32.600513	-28.209319	-32.015351	-38.915889	-45.030536	-27.626605	-37.374830	-28.218543	-36.902038	-17.908039	-35.479652	-13.914346	-34.842636	-14.341979	-37.541959	-8.087776	-33.059295	-5.249140	-31.127079	3.469653	-36.909218	-22.516916	-42.061704	0.812134	-21.824458	-3.643922	-24.825354	3.400142	-31.530953	12.225971	-25.901893	33.646940	-25.562584	29.086933	-21.808228	34.918723	-26.822629	29.632528	-20.874425	-3.324419	34.280579	107.229454	26.078995	106.895877	26.639744	108.264742	30.004974	108.308230	29.305926	3.161926	31.101699	2.532747	30.848212	1.382522	29.022500	1.028279	30.206204	1.717944	32.147152	-56.024486	38.019969	-54.098501	39.418451	-53.787429	37.515374	-42.883306	-41.660090	-31.738819	-33.511964	-44.862971	-46.255894	-56.673186	35.431764	-58.072754	35.620179	-38.655707	19.388902	-42.018818	16.310043	-55.781283	36.745618	-52.791441	36.642657	-52.890329	35.254029	-54.847950	37.762164	-53.906567	36.407287	-55.471836	35.255316	1.003197	31.760854	1.232054	31.171443	-55.047113	32.806612	2.463042	30.472300	-6.698041	-23.277562	5.313645	31.337451	4.898307	30.895230	14.269930	28.601033	5.566176	26.653689	-0.509196	26.859393	11.807545	15.737272	20.279144	9.275757	-52.720077	46.318972	-9.143532	19.809020	-8.023651	13.415144	-2.227479	22.336440	-1.202863	39.731207	125.689454	31.802526	123.623912	32.508088	117.156211	33.176282	-37.123288	-54.906026	-43.694688	-45.107174	128.455630	-58.002019	-43.993427	-45.314961	-43.612136	-45.692385	-42.714110	-47.141979	-48.137915	-43.874833	-43.235771	-47.820500	-43.055107	-48.701403	-40.707692	-49.706207	-37.726444	-43.716800	-29.996543	-32.988102	-31.635897	-20.343771	-47.907232	-41.149117	-49.057972	-42.358962	-37.025331	-40.448852	-38.780988	-37.214421	-49.196320	-42.315908	-51.728746	-27.805367	-54.072081	-31.050007	-52.531416	-31.279070	-55.924184	-34.787853	-25.095995	-46.416338	-18.554177	-39.783580	-36.675843	-48.086661	-28.397381	-9.448268	-35.588422	-7.036870	-31.750556	-6.787735	-50.294655	-3.148709	-42.094067	-15.666340	-35.774094	-30.086812	8.969352	-17.360995	12.791630	-19.035628	-23.766145	-35.698835	9.054193	-27.286213	14.294157	-5.216139	-17.636639	-24.501574	-24.635846	17.733481	-27.592024	-30.581245	2.279867	22.968704	-42.772785	-38.924179	-19.057785	-27.373433	-26.152020	-36.074629	-46.533424	-45.427841	57.727544	-23.778935	-46.730798	-44.842577	-54.094146	-42.105197	-45.203055	-45.541631	26.496295	-32.577800	84.126456	35.631481	79.028185	40.669274	29.286472	45.422690	-29.305237	-40.500239	-29.304897	-40.431166	11.243719	43.897285	19.935230	32.462859	-14.036200	-7.739099	8.557910	38.162535	12.278853	42.797551	7.969308	42.878494	6.951825	39.834531	8.901396	39.105408	9.813555	39.526313	89.232192	45.214151	108.981179	24.688302	109.163832	25.532741	108.542067	25.759209	100.202369	32.515806	108.672868	24.909469	-46.806444	55.036178	-26.442127	50.772641	2.678326	45.782293	15.930648	46.302893	29.391303	50.102649	28.167028	44.526409	-57.750182	-32.657167	52.400360	58.538062	3.931317	25.633931	6.675960	28.333080	-13.184932	-13.209546	-16.696279	-24.084090	-16.689247	-20.599335	-42.883544	-46.818603	-45.915808	-46.364264	94.311315	41.864830	85.503992	40.002109	86.982928	53.029807	92.649657	57.296113	97.374688	58.108571	104.550057	60.695695	81.798605	40.773880	80.711390	44.698624	59.544305	44.664317	96.108367	48.631514	-71.799165	46.367413	-74.173873	43.036074	-71.297477	51.448975	-77.151160	50.362740	-77.204014	57.763566	-76.257800	52.044779	-76.159953	50.855669	-81.742385	58.435887	-74.276825	56.387295	96.430326	35.943030	36.734262	69.572319	77.806002	50.404921	93.535824	42.438248	90.506401	41.787942	90.494991	41.250914	97.310766	43.247135	95.392894	41.333276	96.084909	40.987880	93.138490	40.658745	90.779422	40.630597	96.545773	41.810290	119.664959	45.240630	20.274803	62.024273	33.924777	50.339294	35.843983	43.610201	118.947103	45.551024	93.905097	46.603010	-55.942037	-44.987455	-48.469341	-45.438172	89.480007	41.449681	119.688859	42.230679	122.511017	44.336419	123.016076	44.933614	92.577697	41.055327	92.522019	40.734607	89.809634	41.480110	85.839257	40.523498	79.893406	37.785707	-72.241082	58.437998	-79.743034	62.865890	91.101503	41.399355	91.891568	41.638150	90.332344	41.255453	90.375407	40.985314	89.758048	40.416351	90.093397	40.533564	95.572208	53.194543	99.464398	40.730366	99.635065	40.322652	99.436177	39.929495	100.111879	39.946951	97.347436	45.883575	-40.994489	76.258171	120.704740	43.092163	122.348664	44.699864	121.660245	42.992630	122.092092	43.226345	122.171482	43.101488	120.206861	41.788258	120.427498	41.901832	123.800542	44.209603	123.438879	44.245510	118.370452	38.985500	121.318925	39.804003	90.995443	41.359409	92.999094	38.492605	114.355401	43.983732	114.390308	44.872631	-46.447910	-46.791721	-46.630620	-46.932617	-47.705997	-46.923874	114.508777	44.672569	114.645089	45.444898	123.397571	44.655994	98.642325	52.181731	98.083313	54.833741	87.553629	43.228227	90.315380	42.060566	86.389714	46.165755	89.872723	42.884274	120.540175	44.761967	120.884964	44.755653	120.814251	44.671515	121.018320	44.805152	121.084459	44.878291	121.907195	44.219843	121.607164	44.268775	122.519141	44.471801	122.522009	44.469351	120.592298	44.983636	119.984237	44.438959	119.751699	44.424924	118.822422	42.993890	119.954372	44.418156	120.405625	44.795736	120.363543	44.538554	119.597114	43.838166	117.994624	44.837417	118.703079	44.067078	117.591399	43.604737	118.098624	43.445361	116.950552	44.451155	115.670587	44.324141	117.505826	44.978617	114.865731	47.453919	120.426788	44.840806	119.964153	44.918632	120.715697	44.839303	122.367712	43.846233	123.462963	44.216691	121.571967	44.351959	120.651891	44.887553	114.014902	45.361177	107.895500	47.309485	103.268166	46.297512	-57.850985	48.005929	-70.702869	38.084988	-68.935692	38.734778	-63.801139	39.089484	-69.789814	41.963082	-68.914745	42.707646	-68.449492	43.312862	");
	
//	CRates->Rates[0] = 0.948716989;
//	CRates->Rates[1] = -0.070138332;
//	CRates->Rates[2] = 0.005555377;
//	SetVarRatesFromStr(CRates, Opt, "250000000	-813.286755	-1372.696800	170	-1.980831	0.001681	1.100000	5	3.469827	240423064	Node	1943	2.600950	118283285	Node	3105	4.391442	231630559	Node	2690	3.904356	224265757	Node	1709	23.353375	249978771	Branch	946	1.333179	249995460	Node	1791	3.597936	247983697	Node	2521	1.746464	249987590	Node	244	3.758238	249942279	Node	1585	1.866354	249996573	Node	2812	1.882429	249994757	Node	3380	3.747226	244190618	Node	1955	3.025017	249681048	Node	1482	3.090089	249985549	Node	1247	2.289462	249995449	Node	1248	7.821833	249988633	Branch	3209	10.507939	249871146	Node	1722	6.372645	249991207	Node	2068	6.887873	249995295	Node	1326	1.167916	249996005	Node	461	0.113418	249994597	Node	3013	2.843083	249998544	Node	1796	16.990265	249905712	Node	2261	6.654649	249982903	Node	303	3.585946	249735993	Node	445	7.869760	249976715	Branch	949	2.749042	249847154	Node	1848	16.581278	249976078	Node	858	0.853010	249988029	Node	1765	4.117171	249998613	Branch	2661	3.722776	249999944	Node	3353	6.455643	249995620	Node	2354	1.137258	249988059	Branch	1020	3.172282	249997351	Node	2620	27.851215	249762025	Node	1249	3.699724	249987953	Branch	2849	22.063287	249997027	Branch	2325	1.695874	249990070	Node	1622	18.758698	249331403	Node	3547	12.890709	249472654	Node	2968	11.046752	249825064	Node	2884	7.009993	249953867	Branch	2310	6.187909	249975688	Node	1252	17.721129	249991376	Node	1492	30.556763	249973172	Branch	1219	6.591068	249991859	Node	2152	16.300132	249589723	Branch	2387	10.628437	249912158	Node	67	14.527080	249775050	Node	3279	16.156520	249995078	Node	1818	19.530555	249925573	Node	546	4.571831	249989833	Branch	2587	7.343901	249998091	Branch	2594	12.098696	249897336	Node	1778	16.742890	249993935	Branch	417	1.929979	249998751	Node	2652	0.567537	249994167	Node	1063	7.911753	249976982	Node	518	12.499728	249980179	Node	662	8.800268	249857948	Node	2281	9.992533	249975601	Node	463	15.508110	249978663	Branch	1569	3.764136	249990767	Node	346	20.241182	249978317	Node	2561	5.403297	249994758	Node	1197	6.200797	249992860	Node	2684	9.154263	249997816	Branch	114	11.108907	249997918	Node	2435	2.179189	249982919	Branch	3321	14.254226	249998356	Branch	1111	1.447878	249989641	Node	3330	2.197261	249999647	Node	1532	0.578284	249997503	Branch	2047	4.711460	249999996	Branch	1936	16.629118	249998284	Node	2534	4.355067	249993588	Node	3373	1.289376	249995090	Node	790	6.937615	249997100	Node	2256	7.990196	249975223	Branch	720	7.009595	249991241	Node	1404	8.036000	249688260	Node	1781	3.393456	249985059	Node	1229	43.433809	249956479	Node	2589	14.462304	249973449	Node	1084	6.601475	249328033	Node	447	8.809727	249742057	Node	196	8.491683	249910468	Node	1867	23.548116	249951544	Node	2301	19.732211	249975747	Node	213	2.426158	249984619	Node	2432	6.349351	249986281	Node	2794	7.022304	249989372	Node	237	13.757611	249980193	Branch	1052	3.045850	249990825	Node	158	5.570076	249991265	Node	3451	11.036735	249991582	Node	2475	6.977094	249996463	Node	572	12.005837	249997081	Node	44	11.214223	249990984	Branch	1980	14.698206	249986169	Branch	2089	2.007090	249993360	Branch	261	2.751004	249989080	Branch	663	1.876945	249999157	Node	2464	7.183595	249999177	Node	3181	1.221178	249998625	Branch	3267	5.007782	249987968	Branch	3153	3.107145	249993128	Branch	2202	7.569304	249994735	Branch	1701	1.086896	249991723	Branch	284	19.536830	249984692	Branch	3230	3.024693	249989111	Branch	1122	2.556125	249947178	Branch	1989	17.113265	248809434	Branch	1365	38.695545	248897612	Branch	380	26.498811	249000695	Branch	1225	7.942839	249582223	Branch	633	22.904522	249672209	Branch	1916	14.805290	249752398	Branch	1017	17.187397	249817534	Branch	3362	11.991721	249863064	Branch	1211	21.813305	249867606	Branch	714	25.223538	249874530	Branch	1240	10.493602	249877810	Branch	556	4.398262	249894258	Branch	1896	13.479238	249896593	Branch	333	10.511643	249919753	Branch	550	9.161837	249933939	Branch	1514	14.308321	249941538	Branch	1873	2.545131	249951753	Branch	3477	8.182326	249954453	Branch	2402	4.461188	249955799	Branch	1596	33.115344	249964202	Branch	1814	28.516644	249967460	Branch	2541	12.529041	249973992	Branch	2486	8.880906	249978609	Branch	1782	0.948882	249979819	Branch	1401	23.544714	249980674	Branch	956	7.504675	249982411	Branch	528	4.756283	249983316	Branch	2044	2.953979	249984401	Branch	3269	9.106695	249984739	Branch	116	8.056462	249985857	Branch	2593	13.568829	249987644	Branch	1113	0.944122	249989909	Branch	659	3.337960	249990220	Branch	3490	15.460627	249990673	Branch	310	13.994821	249991063	Branch	698	9.744439	249991588	Branch	3555	14.571446	249992880	Branch	1976	3.086000	249993417	Branch	3163	11.771334	249993465	Branch	2313	0.081916	249993535	Branch	1429	5.877593	249993591	Branch	282	24.611203	249993739	Branch	907	2.052023	249994111	Branch	833	0.019905	249995072	Branch	76	4.731401	249995313	Branch	3315	0.946128	249995385	Branch	2294	9.715803	249996301	Branch	2278	5.608676	249996380	Branch	1054	2.457148	249996426	Branch	2041	0.200671	249997022	Branch	1154	20.882203	249997079	Branch	932	18.337182	249997888	Branch	94	0.956474	249998139	Branch	496	22.342263	249998582	Branch	2789	4.918703	249999365	Branch	897	9.649782	249999454	Branch	3028	4.944496	249999534	Branch	2664	1.287332	249999665	Branch	");
	
	CRates->Lh	=	Likelihood(CRates, Trees, Opt);
	CalcPriors(CRates, Opt); 

//	printf("Lh:\t%f\n", CRates->Lh);fflush(stdout);
//	exit(0);
	
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

	MCMCTest(Opt, Trees, CRates, Shed);

	StartT = GetSeconds();	
	for(Itters=0;;Itters++)
	{ 
		SetCustomSchedule(Itters, Shed);

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
				BurntIn == TRUE
				&& StonesStarted(Opt->Stones, Itters) == FALSE)
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
				// The schedule should be updated even when stones are running
				UpDateSchedule(Opt, Shed, CRates->RS);
				BlankSchedule(Shed);
			}

			if(ExitMCMC(Opt, Itters) == TRUE)
			{


				if( (Opt->UseEqualTrees == FALSE) || 
					(CRates->TreeNo == Trees->NoTrees - 1))
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

