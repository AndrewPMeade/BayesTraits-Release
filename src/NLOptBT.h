#ifndef NLOPTBT_H
#define NLOPTBT_H

#include "typedef.h"
#include "praxis.h"
#include "ml.h"

int		ValidMLAlgName(char *Name);
void	PrintAlgNames(void);

double NLOptBT(RATES *Rates, OPTIONS *Opt, TREES *Trees, ML_MAP *MLMap);

#endif