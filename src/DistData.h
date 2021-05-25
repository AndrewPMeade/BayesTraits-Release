#ifndef DISTDATA_HEADDER
#define DISTDATA_HEADDER

#include "RandLib.h"

DIST_DATA*	LoadDistData(OPTIONS *Opt, TREES *Trees, char *FName);
void		FreeDistData(DIST_DATA *DistData);
void		PrintDistData(FILE *Out, DIST_DATA *DistData);


DIST_DATE_RATES*	CreateDistDataRates(DIST_DATA* DistData, RANDSTATES *RS);
void				FreeDistDataRates(DIST_DATE_RATES* DistRates);

void				CopyDistDataRates(DIST_DATE_RATES* A, DIST_DATE_RATES* B);

void				SetTreeDistData(RATES *Rates, OPTIONS *Opt, TREES *Trees);

void				ChangeTreeDistData(OPTIONS *Opt, RATES *Rates);

void				OutputDataDistHeadder(FILE *Str, OPTIONS *Opt);
void				OutputDataDist(FILE *Str, RATES *Rates, OPTIONS *Opt);
#endif