#ifndef RJ_LOCACL_S_HEADDER
#define RJ_LOCACL_S_HEADDER

#include "typedef.h"

RJ_VARRATE_TYPE			NameToRJLocalType(char *Name);

int						UseRJLocalScalars(OPTIONS *Opt);

PRIORS*					GetPriorFromRJRatesScalar(OPTIONS *Opt, RJ_VARRATE_TYPE Type);

#endif