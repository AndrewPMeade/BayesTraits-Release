#if !defined (PHYLOPLASTYHEADDER)
#define PHYLOPLASTYHEADDER

#include "typedef.h"

int				UseNonParametricMethods(OPTIONS *Opt);

TRANSFORM_TYPE	StrToVarRatesType(char *Str);
char*			VarRatesTypeToStr(TRANSFORM_TYPE Type);

VARRATES*	CreatVarRates(RATES *Rates, TREES *Trees, OPTIONS *Opt);
void		FreeVarRates(VARRATES* Plasty);

void	VarRatesNode(TREES *Trees, TREE *Tree, NODE N, double Scale, TRANSFORM_TYPE Type);

void	VarRatesAddRemove(RATES *Rates, TREES *Trees, OPTIONS *Opt, SCHEDULE *Shed, long long It);
void	ChangeVarRatesScale(RATES *Rates, TREES *Trees, OPTIONS *Opt, SCHEDULE* Shed);
void	VarRatesMoveNode(RATES *Rates, TREES *Trees, OPTIONS *Opt);


void	VarRatesCopy(RATES *R1, RATES *R2);
void	VarRatesTree(OPTIONS *Opt, TREES *Trees, RATES *Rates, int Normalise);

void	InitVarRatesFiles(OPTIONS *Opt, TREES *Trees, RATES *Rates);
void	FinishVarRatesFiles(OPTIONS *Opt);
void	PrintVarRatesOutput(OPTIONS *Opt, TREES *Trees, RATES *Rates, long long It);

double	CaclVRPrior(double X, TRANSFORM_TYPE Type, RATES *Rates);
double	CalcVarRatesPriors(RATES *Rates, OPTIONS *Opt);
void	ChangeVarRatesHyperPrior(RATES *Rates, OPTIONS *Opt);

void	SetVarRatesFromStr(char *Str, RATES *Rates, OPTIONS *Opt);

double	CalcNormalHasting(double x, double SD);

#endif
