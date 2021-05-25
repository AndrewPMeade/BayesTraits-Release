#ifndef RJ_LOCACL_S_HEADDER
#define RJ_LOCACL_S_HEADDER

#include "typedef.h"

TRANSFORM_TYPE			NameToRJLocalType(char *Name);

int						UseRJLocalScalars(OPTIONS *Opt);

PRIORS*					GetPriorFromRJRatesScalar(OPTIONS *Opt, TRANSFORM_TYPE Type);

#endif