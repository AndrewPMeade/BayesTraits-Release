#ifndef CONTRASTS_H
#define CONTRASTS_H

#include "TypeDef.h"

REG_BETA_SPACE*	InitRegBetaSpace(int NoSites, int NoCont);
void			FreeRegBetaSpace(REG_BETA_SPACE* RSpace);

void		InitContrastAll(OPTIONS *Opt, TREES* Trees);
void		InitContrastTree(OPTIONS *Opt, TREES* Trees, int TNo, int NoSites);
void		FreeAllContrast(OPTIONS *Opt, TREES* Trees);
void		FreeContrast(OPTIONS *Opt, TREES* Trees, int TreeNo);


void		CalcContrast(OPTIONS *Opt, TREES* Trees, RATES* Rates);
double		CalcContrastLh(OPTIONS *Opt, TREES* Trees, RATES* Rates);

CONTRASTR*	CreatContrastRates(OPTIONS *Opt, RATES *Rates, TREES *Trees);
void		FreeContrastRates(RATES *Rates);

void		MapRatesToConVals(OPTIONS *Opt, RATES *Rates, CONTRASTR* ConR, TREES *Trees);

void		MutateContrastRates(OPTIONS *Opt, TREES* Trees, RATES* Rates, SCHEDULE*	Shed);
void		CopyContrastRates(OPTIONS *Opt, RATES* R1, RATES* R2, int NoSites);

void		RecIntNode(NODE N, int SiteNo, double *Alpha);

double		DataToZScore(double X, double Mean, double SD);

void		PrintContrast(RATES *Rates, TREES *Trees);

void		RetSetConTraitData(TREE *Tree, int NoSites);


void		CalcNodeContrast(NODE N, int NoSites);

#endif
