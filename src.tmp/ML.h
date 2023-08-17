#ifndef MLHEADDER
#define MLHEADDER

#include "TypeDef.h"

typedef enum
{
	ML_P_TYPE_NONE,
	ML_P_TYPE_RATE_S
} ML_P_TYPE;

typedef struct
{
	int NoP;

	double	*PVal;
	double	*PMin;
	double	*PMax;
	double	*PDef;
	ML_P_TYPE	*PType;

} ML_MAP;

#ifdef	 JNIRUN
	void	FindML(OPTIONS *Opt, TREES *Trees, JNIEnv *Env, jobject Obj);
#else
	void	FindML(OPTIONS *Opt, TREES *Trees);
#endif

void	MLTree(OPTIONS *Opt, TREES *Trees, RATES *Rates);
double	LikelihoodML(ML_MAP* MLMap, OPTIONS *Opt, TREES *Trees, RATES *Rates);
void	MLMapToRates(ML_MAP* MLMap, OPTIONS *Opt, RATES *Rates);


ML_MAP*	AllocMLMap(void);
void	FreeMLMap(ML_MAP *MLMap);

void	BuildMLMap(ML_MAP*	MLMap, OPTIONS *Opt, TREES *Trees, RATES *Rates);
ML_MAP*	MLMapTreeTry(OPTIONS *Opt, TREES *Trees, RATES *Rates, ML_MAP *Init);


void Opt1D(ML_MAP* Map, OPTIONS *Opt, TREES *Trees, RATES *Rates);
void CopyMLMap(ML_MAP *A, ML_MAP *B);

#endif
