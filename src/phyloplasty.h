#if !defined (PHYLOPLASTYHEADDER)
#define PHYLOPLASTYHEADDER

#include "typedef.h"

PLASTY*	CreatPlasty(RATES *Rates, TREES *Trees, OPTIONS *Opt);
void	FreePlasty(PLASTY* Plasty);

void	PlastyMove(RATES *Rates, TREES *Trees, OPTIONS *Opt);
void	PlastyCopy(RATES *R1, RATES *R2);
void	Plasty(RATES *Rates, TREES *Trees, OPTIONS *Opt);

#endif
