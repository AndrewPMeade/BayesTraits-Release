#ifndef LIKELIHOOD_HEADDER
#define LIKELIHOOD_HEADDER

#include "praxis.h"

void	FossilLh(NODE N, OPTIONS *Opt, TREES *Trees, int SiteNo);

double		Likelihood(RATES* Rates, TREES *Trees, OPTIONS *Opt);
void		AllocLHInfo(TREES *Trees, OPTIONS *Opt);
void		FreeInvInfo(INVINFO* InvInfo);
INVINFO*	AllocInvInfo(int NOS);

int			ValidLh(double LH, MODEL_TYPE MT);

double		LhPraxis(void* PState, double *List);

int			IsNum(double n);

void		SetUpAMatrix(MODEL Model, RATES *Rates, TREES *Trees, int NOS, INVINFO *InvInfo, double *RateP, double *Pi);

int			SetAllPMatrix(RATES* Rates, TREES *Trees, OPTIONS *Opt, double Gamma);
int			SetUpAllAMatrix(RATES* Rates, TREES *Trees, OPTIONS *Opt);


void		LhTransformTree(RATES* Rates, TREES *Trees, OPTIONS *Opt);

#endif
