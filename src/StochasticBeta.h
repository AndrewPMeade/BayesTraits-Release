#ifndef STOCHASTIC_BETA_H
#define STOCHASTIC_BETA_H

#include <stdio.h>
#include <stdlib.h>

#include "GenLib.h"
#include "TypeDef.h"

void		InitStochasticBeta(OPTIONS *Opt, TREES *Trees, RATES *Rates);
void		FreeStochasticBeta(RATES *Rates);

void		MapStochasticBeta(TREES *Trees, RATES *Rates);
void		CopyStochasticBeta(RATES *RatesA, RATES *RatesB);

//void		ChangeStochasticBeta(RATES *Rates, SCHEDULE *Sched);
void		ChangeStochasticBeta(TREES *Trees, RATES *Rates, SCHEDULE *Sched);

void		CaclStochasticBetaRJ(RATES *Rates);


double		CaclStochasticBetaPrior(TREES *Trees, RATES *Rates);

int			GetNoStochasticBetaType(RATES *Rates, STOCHASTIC_BETA_TYPE SB_Type);

void		LogStochasticBetaResults(OPTIONS *Opt, TREES *Trees, RATES *Rates, long long It);


#endif