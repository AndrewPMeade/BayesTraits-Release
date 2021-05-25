#if !defined (LIKLIHEADDER)
#define LIKLIHEADDER

#include "praxis.h"

double	Likelihood(RATES* Rates, TREES *Trees, OPTIONS *Opt);
void	AllocLHInfo(TREES *Trees, OPTIONS *Opt);
void	FreeInvInfo(INVINFO* InvInfo);

double	LhPraxis(void* PState, double *List);

int		IsNum(double n);

#endif
