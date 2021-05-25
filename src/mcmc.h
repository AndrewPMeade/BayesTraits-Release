#if !defined (MCMCHEADDER)
#define MCMCHEADDER

#include "typedef.h"

#ifdef	 JNIRUN
	void	MCMC(OPTIONS *Opt, TREES *Trees, JNIEnv *Env, jobject Obj);
#else
	void	MCMC(OPTIONS *Opt, TREES *Trees);
#endif

	double	ValidMCMCParameters(OPTIONS *Opt, TREES *Trees, RATES *Rates);
#endif
