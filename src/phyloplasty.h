#if !defined (PHYLOPLASTYHEADDER)
#define PHYLOPLASTYHEADDER

#include "typedef.h"

PLASTY*	CreatPlasty(RATES *Rates, TREES *Trees, OPTIONS *Opt);
void	FreePlasty(PLASTY* Plasty);

void	PPAddRemove(RATES *Rates, TREES *Trees, OPTIONS *Opt, long long It);
void	PPChangeScale(RATES *Rates, TREES *Trees, OPTIONS *Opt);
void	PPMoveNode(RATES *Rates, TREES *Trees, OPTIONS *Opt);


void	PlastyCopy(RATES *R1, RATES *R2);
void	Plasty(OPTIONS *Opt, TREES *Trees, RATES *Rates, int Normalise);

void	InitPPFiles(OPTIONS *Opt, TREES *Trees, RATES *Rates);
void	PrintPPOutput(OPTIONS *Opt, TREES *Trees, RATES *Rates, long long It);

double	CalcPPPriors(RATES *Rates, OPTIONS *Opt);
void	ChangePPHyperPrior(RATES *Rates, OPTIONS *Opt);

void	CheckPlasty(RATES *Rates, TREES *Trees, OPTIONS *Opt);

#endif
