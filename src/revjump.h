#if !defined (RJMCMCHEADDER)
#define RJMCMCHEADDER

#include	"typedef.h"

void		MapRJRates(OPTIONS *Opt, RATES *Rates);
int			NoOfPramGroups(RATES* Rates, int *GroupID, int *GroupSize);
void		RJSplit(RATES* Rates, OPTIONS* Opt);
void		RJMerge(RATES* Rates, OPTIONS* Opt);
MAPINFO*	CreatMapping(RATES* Rates);
void		FreeMapInfo(MAPINFO *MapInfo);

void		RJAugment(RATES* Rates, OPTIONS* Opt);
void		RJReduce(RATES* Rates, OPTIONS* Opt);

#endif
