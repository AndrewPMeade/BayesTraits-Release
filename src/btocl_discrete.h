#ifndef BTOCL_DISCRETE_H
#define BTOCL_DISCRETE_H

#ifdef BTOCL_DSC
#include "btocl_runtime.h"
#include "typedef.h"

int		SetAllPMatrixOpenCL(RATES *Rates, TREES *Trees, OPTIONS *Opt, double RateMult, double Kappa);

#endif

#endif
