#ifndef CONTRASTS
#define CONTRASTS

#include "typedef.h"

void		InitContrastAll(OPTIONS *Opt, TREES* Trees);
void		InitContrastTree(OPTIONS *Opt, TREES* Trees, int TNo, int NoSites);
void		FreeAllContrast(OPTIONS *Opt, TREES* Trees);
void		FreeContrast(OPTIONS *Opt, TREES* Trees, int TreeNo);


void		CalcContrast(TREES* Trees, RATES* Rates);
double		CalcContrastLh(OPTIONS *Opt, TREES* Trees, RATES* Rates);

CONTRASTR*	CreatContrastRates(OPTIONS *Opt, RATES *Rates);
void		FreeContrastRates(RATES *Rates);

void		MLContrast(OPTIONS *Opt, TREES *Trees, RATES *Rates);

void		MapRatesToConVals(OPTIONS *Opt, RATES *Rates, CONTRASTR* ConR);

void		MutateContrastRates(OPTIONS *Opt, TREES* Trees, RATES* Rates, SCHEDULE*	Shed);
void		CopyContrastRates(OPTIONS *Opt, RATES* R1, RATES* R2, int NoSites);

//void		RecIntNode(NODE N, int SiteNo, double *Alpha, double *Sigma, double *Lh);
void		RecIntNode(NODE N, int SiteNo, double *Alpha);

void		NormaliseReg(OPTIONS *Opt, TREES *Trees, RATES *Rates);
#endif
