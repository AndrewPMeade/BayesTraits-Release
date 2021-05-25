#if !defined	CONTRASTS
#define CONTRASTS

#include "typedef.h"

void	InitContrast(OPTIONS *Opt, TREES* Trees);
void	FreeContrast(OPTIONS *Opt, TREES* Trees);

void	CalcContrast(TREES* Trees, RATES* Rates);
double	CalcContrastLh(TREES* Trees, RATES* Rates);

#endif