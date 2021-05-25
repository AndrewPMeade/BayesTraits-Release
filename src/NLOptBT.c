#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "typedef.h"
#include "praxis.h"
#include "genlib.h"
#include "likelihood.h"

#ifdef NLOPT_BT
//	#include "nlopt.h"
	#include <nlopt.h>


double	NLOptLh(unsigned N, const double *x, double *grad, void *Data)
{
	PRAXSTATE *PState;
	double	Lh;

	PState = (PRAXSTATE*)Data;

	Lh = LhPraxis(PState, (double*)x);

//	if(Lh == ERRLH)
//		Lh = -ERRLH;
	
	return Lh;
}


void	SetMinMax(nlopt_opt Opt, PRAXSTATE *PState)
{

	nlopt_set_lower_bounds1(Opt, MINRATE);
	nlopt_set_upper_bounds1(Opt, MAXRATE);

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


double NLOptBT(double *X, PRAXSTATE *PState)
{
	nlopt_opt Opt;
	RATES	*Rates;
	double	*InitX, InitLh;
	int		Index;
	double	Lh;

	Rates = PState->Rates;	

	memcpy(Rates->Rates, X, Rates->NoOfRates);
//	Lh 	= Likelihood(Rates, PState->Trees, PState->Opt);
	

	Opt = nlopt_create(NLOPT_LN_COBYLA, Rates->NoOfRates); /* algorithm and dimensionality */
//	Opt = nlopt_create(NLOPT_LN_BOBYQA, Rates->NoOfRates); // Not very good on primaites full. 
//	Opt = nlopt_create(NLOPT_LN_NEWUOA, Rates->NoOfRates); /* algorithm and dimensionality */
//	Opt = nlopt_create(NLOPT_LN_PRAXIS, Rates->NoOfRates); /* algorithm and dimensionality */
//	Opt = nlopt_create(NLOPT_LN_NELDERMEAD, Rates->NoOfRates); /* algorithm and dimensionality */
	
	SetMinMax(Opt, PState);
	
//	nlopt_set_xtol_rel(Opt, 1e-4);	//	0.0001
	nlopt_set_xtol_rel(Opt, 0.0001);
	
	nlopt_set_min_objective(Opt, NLOptLh, (void*)PState);

//	nlopt_set_maxeval(Opt, 100000);
	
	nlopt_optimize(Opt, X, &Lh);

/*	printf("Lh:\t%f\t%d\t", Lh, NoOpt);

	for(Index=0;Index<Rates->NoOfRates;Index++)
		printf("%f\t", X[Index]);
	printf("\n");
*/
	memcpy(Rates->Rates, X, sizeof(double) * Rates->NoOfRates);

	Rates->Lh = Likelihood(Rates, PState->Trees, PState->Opt);

	return Rates->Lh;
}
#endif
