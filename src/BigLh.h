#ifndef BIGLH_HEADDER
#define BIGLH_HEADDER

#include "typedef.h"

void	InitTreeBigLh(OPTIONS *Opt, TREES *Trees);
void	FreeTreeBigLh(OPTIONS *Opt, TREES *Trees);

void	LhBigLh(NODE N, TREES *Trees, int Pre, int SiteNo, int *Err);

double 	CombineBigLh(RATES* Rates, TREES *Trees, OPTIONS *Opt, int SiteNo, int NOS);

void	SetBigLhNodeRec(NODE N, int NOS, int NoOfSites, RATES *Rates, OPTIONS *Opt);

#endif
