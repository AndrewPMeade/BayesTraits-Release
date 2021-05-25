#ifndef LANDSCAPE_RJ_GROUPS_H
#define LANDSCAPE_RJ_GROUPS_H

#include "GenLib.h"
#include "TypeDef.h"

LAND_RATE_GROUPS*	CreateLandRateGroups(int NoRates);
void				FreeLandRateGroups(LAND_RATE_GROUPS* LandRateG);

void				CopyLandRateGroups(LAND_RATE_GROUPS* A, LAND_RATE_GROUPS* B);

void				MapLandRateGroups(OPTIONS *Opt, TREES *Trees, RATES *Rates);

#endif