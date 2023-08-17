#ifndef RJMCMCHEADDER
#define RJMCMCHEADDER

#include "TypeDef.h"

void		MapRJRates(OPTIONS *Opt, RATES *Rates);
int			NoOfPramGroups(RATES* Rates, int *GroupID, int *GroupSize);
MAPINFO*	CreatMapping(RATES* Rates);
void		FreeMapInfo(MAPINFO *MapInfo);


int			RJSplit(RATES* Rates, OPTIONS* Opt);
int			RJMerge(RATES* Rates, OPTIONS* Opt);

int			RJAugment(RATES* Rates, OPTIONS* Opt);
int			RJReduce(RATES* Rates, OPTIONS* Opt);

#endif
