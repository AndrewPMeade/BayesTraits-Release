#ifndef CON_GLOBAL_TREND_H
#define CON_GLOBAL_TREND_H

#include "TypeDef.h"
#include "Trees.h"

void TestContrastGlobalTrend(OPTIONS* Opt, TREES* Trees, RATES* Rates);
void SetGlobalTrend(OPTIONS* Opt, TREES* Trees, RATES* Rates);

void	ChangeGlobalTrend(RATES* Rates, SCHEDULE* Shed);
double	CalcGlobalTrendPrior(RATES *Rates);


#endif