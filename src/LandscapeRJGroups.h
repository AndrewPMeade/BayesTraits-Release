#ifndef LANDSCAPE_RJ_GROUPS_H
#define LANDSCAPE_RJ_GROUPS_H

#include "GenLib.h"
#include "TypeDef.h"

LAND_RATE_GROUPS*	CreateLandRateGroups(RATES *Rates, TREES *Trees);
void				FreeLandRateGroups(LAND_RATE_GROUPS* LandRateG);

void				CopyLandRateGroups(LAND_RATE_GROUPS* A, LAND_RATE_GROUPS* B);

void				MapLandRateGroups(OPTIONS *Opt, TREES *Trees, RATES *Rates);

void				SwapRateGroupParam(RATES *Rates);

void				InishLandRateGroupsOutput(OPTIONS *Opt, TREES *Trees, RATES *Rates);
void				OutputLandRateGroupsOutput(long long Itter, OPTIONS *Opt, TREES *Trees, RATES *Rates);

double				CalcPriorLandRateGoup(RATES *Rates);

void				ChangeLandRateGoup(RATES *Rates,  SCHEDULE* Shed);

void				SetNoFixedGroups(RATES *Rates, LAND_RATE_GROUPS* LandRGroup, int NoGroup);


#endif