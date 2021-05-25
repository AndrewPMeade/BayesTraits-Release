#if !defined (MLHEADDER)
#define MLHEADDER

#include "typedef.h"

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
#endif
