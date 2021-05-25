#ifndef BTOCL_DISCRETE_H
#define BTOCL_DISCRETE_H

#ifdef BTOCL_DSC
#include "btocl_runtime.h"
#include "btocl_runtime_kernels.h"
#include "typedef.h"

void btocl_AllocLhInfo(TREES* Trees);
void btocl_FreeLhInfo(TREES* Trees);


int	 btocl_SetAllPMatrix(RATES *Rates, TREES *Trees, OPTIONS *Opt, double RateMult, double Kappa);

int btocl_computePartialLh(RATES* Rates, TREES *Trees, OPTIONS *Opt);

#endif

#endif
