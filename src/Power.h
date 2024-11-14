#ifndef POWER_H
#define POWER_H

#include "TypeDef.h"

int				GetNoPowerSites(OPTIONS *Opt);

void			PrintPowerOpt(FILE* Str, OPTIONS *Opt);

SITE_POWER*		CrateSitePowers(OPTIONS *Opt);
void			FreeSitePowers(SITE_POWER* SitePower);

void			CopySitePowers(SITE_POWER* A, SITE_POWER* B);

void			SetSitePower(RATES* Rates, TREES *Trees, OPTIONS *Opt);

void			ChangeSitePowers(RATES *Rates, SCHEDULE* Shed);
double			CalcSitePowersPrior(RATES *Rates);

#endif
