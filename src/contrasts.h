#ifndef CONTRASTS
#define CONTRASTS

#include "typedef.h"

void	InitContrastAll(OPTIONS *Opt, TREES* Trees);
void	InitContrastTree(OPTIONS *Opt, TREES* Trees, int TNo);
void	FreeAllContrast(OPTIONS *Opt, TREES* Trees);
void	FreeContrast(OPTIONS *Opt, TREES* Trees, int TreeNo);

void	CalcContrast(TREES* Trees, RATES* Rates);
double	CalcContrastLh(OPTIONS *Opt, TREES* Trees, RATES* Rates);

CONTRASTR*	AllocContrastRates(OPTIONS *Opt, RATES *Rates);

void	MLContrast(OPTIONS *Opt, TREES *Trees, RATES *Rates);

void	MutateContrastRates(OPTIONS *Opt, TREES* Trees, RATES* Rates);
void	CopyContrastRates(RATES* R1, RATES* R2, int NoSites);

#endif
