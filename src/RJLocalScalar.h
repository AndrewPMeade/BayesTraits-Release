#ifndef RJ_LOCACL_S_HEADDER
#define RJ_LOCACL_S_HEADDER

#include "TypeDef.h"

TRANSFORM_TYPE			NameToRJLocalType(char *Name, int *Err);

int						UseRJLocalScalars(OPTIONS *Opt);

PRIOR*					GetPriorFromRJRatesScalar(OPTIONS *Opt, TRANSFORM_TYPE Type);

#endif