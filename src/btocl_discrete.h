#ifndef BTOCL_DISCRETE_H
#define BTOCL_DISCRETE_H

#ifdef BTOCL
#include "btocl_runtime.h"
#include "btocl_kernels_bayestraits.h"
#include "typedef.h"

void btocl_AllocPMatrixInfo(TREES* Trees);
int	 btocl_SetAllPMatrix(RATES *Rates, TREES *Trees, OPTIONS *Opt, double RateMult);
void btocl_FreePMatrixInfo(TREES* Trees);

void btocl_AllocLhInfo(TREES* Trees);
int btocl_computePartialLh(RATES* Rates, TREES *Trees, OPTIONS *Opt);
void btocl_FreeLhInfo(TREES* Trees);


#endif

#endif
