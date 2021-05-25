#if !defined (PHYLOPLASTYHEADDER)
#define PHYLOPLASTYHEADDER

#include "typedef.h"

int		UseNonParametricMethods(OPTIONS *Opt);

PLASTY*	CreatVarRates(RATES *Rates, TREES *Trees, OPTIONS *Opt);
void	FreeVarRates(PLASTY* Plasty);

void	VarRatesAddRemove(RATES *Rates, TREES *Trees, OPTIONS *Opt, SCHEDULE *Shed, long long It);
void	ChangeVarRatesScale(RATES *Rates, TREES *Trees, OPTIONS *Opt);
void	VarRatesMoveNode(RATES *Rates, TREES *Trees, OPTIONS *Opt);


void	VarRatesCopy(RATES *R1, RATES *R2);
void	VarRatesTree(OPTIONS *Opt, TREES *Trees, RATES *Rates, int Normalise);

void	InitVarRatesFiles(OPTIONS *Opt, TREES *Trees, RATES *Rates);
void	PrintVarRatesOutput(OPTIONS *Opt, TREES *Trees, RATES *Rates, long long It);

double	CalcVarRatesPriors(RATES *Rates, OPTIONS *Opt);
void	ChangeVarRatesHyperPrior(RATES *Rates, OPTIONS *Opt);

#endif
