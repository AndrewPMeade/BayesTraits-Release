#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "typedef.h"
#include "praxis.h"
#include "genlib.h"
#include "likelihood.h"
#include "ml.h"

#ifdef NLOPT

typedef struct
{
	OPTIONS *Opt;
	TREES	*Trees;
	RATES	*Rates;
	ML_MAP	*MLMap;
	int		NoCalled;
} NLOPT_LH;



double	NLOptLh(unsigned N, const double *x, double *grad, void *Data)
{
	NLOPT_LH *NLOptLh;
	double	Lh;

	NLOptLh = (NLOPT_LH*)Data;
	
	memcpy(NLOptLh->MLMap->PVal, x, sizeof(double) * N);

	Lh = LikelihoodML(NLOptLh->MLMap, NLOptLh->Opt, NLOptLh->Trees, NLOptLh->Rates);
	
	NLOptLh->NoCalled++;

//	printf("%d\t%f\n", NLOptLh->NoCalled, Lh);fflush(stdout);

//	if(Lh == ERRLH)
//		Lh = -ERRLH;
	
	return Lh;
}


void	SetMinMax(nlopt_opt Opt, ML_MAP *MLMap)
{

	nlopt_set_lower_bounds(Opt, MLMap->PMin);
	nlopt_set_upper_bounds(Opt, MLMap->PMax);

/*
	N = PState->Rates->NoOfRates;

	Min = (double*)malloc(sizeof(double) * N);
	Max = (double*)malloc(sizeof(double) * N);
	if((Min == NULL) || (Max == NULL))
		MallocErr();

	for(Index=0;Index<N;Index++)
	{
		Min[Index] = MINRATE;
		Max[Index] = MAXRATE;
	}
*/
}

NLOPT_LH*	CreateNLOptLh(RATES *Rates, OPTIONS *Opt, TREES *Trees, ML_MAP *MLMap)
{
	NLOPT_LH *Ret;

	Ret = (NLOPT_LH*)SMalloc(sizeof(NLOPT_LH));

	Ret->NoCalled = 0;
	Ret->Rates = Rates;
	Ret->Opt = Opt;
	Ret->Trees = Trees;
	Ret->MLMap = MLMap;

	return Ret;
}


// ./Seq/Primates-25.trees ./Seq/Primates.txt <in.txt > sout.txt

double NLOptBT(RATES *Rates, OPTIONS *Opt, TREES *Trees, ML_MAP *MLMap)
{
	nlopt_opt NLOpt;
	double	*TRates;
	double	Lh;
	NLOPT_LH *OStruct;
	

	TRates = (double*)CloneMem(sizeof(double) * MLMap->NoP, MLMap->PVal);

	nlopt_srand(RandUSLong(Rates->RS));


	OStruct = CreateNLOptLh(Rates, Opt, Trees, MLMap);

	Lh 	= Likelihood(Rates, Trees, Opt);


/* Good ones */
	NLOpt = nlopt_create(NLOPT_LN_BOBYQA, MLMap->NoP);
//	NLOpt = nlopt_create(NLOPT_LN_NEWUOA, MLMap->NoP);
//	NLOpt = nlopt_create(NLOPT_LN_NELDERMEAD, MLMap->NoP);
//	NLOpt = nlopt_create(NLOPT_LN_PRAXIS, MLMap->NoP);

//	NLOpt = nlopt_create(NLOPT_LN_PRAXIS, MLMap->NoP);
//	NLOpt = nlopt_create(NLOPT_LN_COBYLA, MLMap->NoP);
//	NLOpt = nlopt_create(NLOPT_LN_NEWUOA, MLMap->NoP);
//	NLOpt = nlopt_create(NLOPT_LN_NEWUOA_BOUND, MLMap->NoP);
//	NLOpt = nlopt_create(NLOPT_LN_NELDERMEAD, MLMap->NoP);
//	NLOpt = nlopt_create(NLOPT_LN_SBPLX, MLMap->NoP);
//	NLOpt = nlopt_create(NLOPT_LN_AUGLAG, MLMap->NoP);
//	NLOpt = nlopt_create(NLOPT_LN_AUGLAG_EQ, MLMap->NoP);
//	NLOpt = nlopt_create(NLOPT_LN_BOBYQA, MLMap->NoP);

//	NLOpt = nlopt_create(NLOPT_GN_DIRECT, MLMap->NoP);
//	NLOpt = nlopt_create(NLOPT_GN_DIRECT_L, MLMap->NoP);
//	NLOpt = nlopt_create(NLOPT_GN_DIRECT_L_RAND, MLMap->NoP);
//	NLOpt = nlopt_create(NLOPT_GN_DIRECT_NOSCAL, MLMap->NoP);
//	NLOpt = nlopt_create(NLOPT_GN_DIRECT_L_NOSCAL, MLMap->NoP);
//	NLOpt = nlopt_create(NLOPT_GN_DIRECT_L_RAND_NOSCAL, MLMap->NoP);
//	NLOpt = nlopt_create(NLOPT_GN_ORIG_DIRECT, MLMap->NoP);
//	NLOpt = nlopt_create(NLOPT_GN_ORIG_DIRECT_L, MLMap->NoP);
//	NLOpt = nlopt_create(NLOPT_GN_CRS2_LM, MLMap->NoP);
//	NLOpt = nlopt_create(NLOPT_GN_MLSL, MLMap->NoP);
//	NLOpt = nlopt_create(NLOPT_GN_MLSL_LDS, MLMap->NoP);
//	NLOpt = nlopt_create(NLOPT_GN_ISRES, MLMap->NoP);
//	NLOpt = nlopt_create(NLOPT_GN_ESCH, MLMap->NoP);
			
	SetMinMax(NLOpt, MLMap);
	
	nlopt_set_xtol_rel(NLOpt, 0.000001);
	nlopt_set_maxeval(NLOpt, 20000);

	nlopt_set_max_objective(NLOpt, NLOptLh, (void*)OStruct);
		
	nlopt_optimize(NLOpt, TRates, &Lh);

	memcpy(MLMap->PVal, TRates, sizeof(double) * MLMap->NoP);

//	MLMapToRates(MLMap, Opt, Rates);
//	Rates->Lh = Likelihood(Rates, Trees, Opt);

	nlopt_destroy(NLOpt);

	free(TRates);
	free(OStruct);

	return Rates->Lh;
}
#endif
