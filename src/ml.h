#if !defined (MLHEADDER)
#define MLHEADDER

#include "typedef.h"

#ifdef	 JNIRUN
	void	FindML(OPTIONS *Opt, TREES *Trees, JNIEnv *Env, jobject Obj);
#else
	void	FindML(OPTIONS *Opt, TREES *Trees);
#endif

void	MLTree(OPTIONS *Opt, TREES *Trees, RATES *Rates);
void	FindValidStartSet(OPTIONS *Opt, TREES *Trees, RATES *Rates);

#endif
