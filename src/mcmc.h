#if !defined (MCMCHEADDER)
#define MCMCHEADDER

#include "typedef.h"

#ifdef	 JNIRUN
	void	MCMC(OPTIONS *Opt, TREES *Trees, JNIEnv *Env, jobject Obj);
#else
	void	MCMC(OPTIONS *Opt, TREES *Trees);
#endif

	void	LhOverAllModels(OPTIONS *Opt, TREES *Trees);
#endif
