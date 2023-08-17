#ifndef MCMC_H
#define MCMC_H

#include "TypeDef.h"

// Number of attempts to find a valid set of starting parameters
#define NO_RAND_START_TRIES 1000000


#ifdef	 JNIRUN
	void	MCMC(OPTIONS *Opt, TREES *Trees, JNIEnv *Env, jobject Obj);
#else
	void	MCMC(OPTIONS *Opt, TREES *Trees);
#endif

	double	ValidMCMCParameters(OPTIONS *Opt, TREES *Trees, RATES *Rates);
#endif
