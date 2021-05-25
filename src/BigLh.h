#ifndef BIGLH_HEADDER
#define BIGLH_HEADDER

#include "typedef.h"

void	InitTreeBigLh(OPTIONS *Opt, TREES *Trees);
void	FreeTreeBigLh(OPTIONS *Opt, TREES *Trees);

void	LhBigLh(NODE N, TREES *Trees, int Pre, int SiteNo);

double 	CombineBigLh(RATES* Rates, TREES *Trees, OPTIONS *Opt, int SiteNo, int NOS);

void	SetBigLhNodeRec(NODE N, int NOS, int NoSites, RATES *Rates, OPTIONS *Opt);

void	FossilLhBig(NODE N, TREES *Trees, int *Mask, int SiteNo);

#endif
