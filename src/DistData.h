#ifndef DISTDATA_HEADDER
#define DISTDATA_HEADDER

#include "RandLib.h"

DIST_DATA*	LoadDistData(TREES *Trees, char *FName);
void		FreeDistData(DIST_DATA *DistData);
void		PrintDistData(FILE *Out, DIST_DATA *DistData);


DIST_DATE_RATES*	CreateDistDataRates(DIST_DATA* DistData, RANDSTATES *RS);
void	FreeDistDataRates(DIST_DATE_RATES* DistRates);
#endif