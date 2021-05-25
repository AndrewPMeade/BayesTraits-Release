#ifndef NLOPTBT_H
#define NLOPTBT_H

#include "TypeDef.h"
#include "Praxis.h"
#include "ML.h"

int		ValidMLAlgName(char *Name);
void	PrintAlgNames(void);

double NLOptBT(RATES *Rates, OPTIONS *Opt, TREES *Trees, ML_MAP *MLMap);

#endif