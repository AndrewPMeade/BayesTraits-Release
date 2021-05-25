#if !defined (MLHEADDER)
#define MLHEADDER

#include "typedef.h"

#ifdef	 JNIRUN
	void	FindML(OPTIONS *Opt, TREES *Trees, JNIEnv *Env, jobject Obj);
#else
	void	FindML(OPTIONS *Opt, TREES *Trees);
#endif

double	PraxisGo(OPTIONS *Opt, RATES *Rates, TREES *Trees);
void	FindValidStartSet(double *Vec, RATES *Rates, TREES *Trees, OPTIONS *Opt);
#endif
