#ifndef VAR_RATES_HEADDER
#define VAR_RATES_HEADDER

#include "TypeDef.h"


#define NO_FABRIC_HOMO_P	2		
#define FABRIC_HOMO_PRIOR_COST	-2		

int				UseRJLocalScalar(OPTIONS *Opt);
int				UseNonParametricMethods(OPTIONS *Opt);

TRANSFORM_TYPE	StrToVarRatesType(char *Str);
char*			VarRatesTypeToStr(TRANSFORM_TYPE Type);

VARRATES*	CreatVarRates(RATES *Rates, TREES *Trees, OPTIONS *Opt);
void		FreeVarRates(VARRATES* Plasty);

void	VarRatesNode(TREES *Trees, TREE *Tree, NODE N, double Scale, TRANSFORM_TYPE Type);

void	VarRatesAddRemove(RATES *Rates, TREES *Trees, OPTIONS *Opt, SCHEDULE *Shed, size_t It);
void	ChangeVarRatesScale(RATES *Rates, TREES *Trees, OPTIONS *Opt, SCHEDULE* Shed);
void	VarRatesMoveNode(RATES *Rates, TREES *Trees, OPTIONS *Opt);


void	VarRatesCopy(TREES *Trees, RATES *R1, RATES *R2);
void	VarRatesTree(OPTIONS *Opt, TREES *Trees, RATES *Rates, int Normalise);

void	InitVarRatesFiles(OPTIONS *Opt, TREES *Trees, RATES *Rates);
void	FinishVarRatesFiles(OPTIONS *Opt);
void	PrintVarRatesOutput(OPTIONS *Opt, TREES *Trees, RATES *Rates, size_t It);

double	CaclVRPrior(double X, TRANSFORM_TYPE Type, RATES *Rates);
double	CalcVarRatesPriors(RATES *Rates, OPTIONS *Opt);
void	ChangeVarRatesHyperPrior(RATES *Rates, OPTIONS *Opt);

void	SetVarRatesFromStr(RATES *Rates, OPTIONS *Opt, char *Str, TREES *Trees);

double	CalcNormalHasting(double x, double SD);

NODE	GetVRNode(TREES *Trees, int TreeNo, VAR_RATES_NODE *VR_Node);

void	OutputVarRatesType(FILE *Out, TRANSFORM_TYPE Type);
void	DumpVarRates(FILE *Out, TREES *Trees, VARRATES* VarRates);

void	SetFabricAncStates(OPTIONS *Opt, TREES *Trees, RATES *Rates);

void	SetHomoFabric(RATES *Rates);
void	ChangeHomoFabric(RATES *Rates, SCHEDULE* Shed);


#endif
