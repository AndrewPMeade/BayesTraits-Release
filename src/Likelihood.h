#ifndef LIKELIHOOD_HEADDER
#define LIKELIHOOD_HEADDER

#include "praxis.h"

void	FossilLh(NODE N, TREES *Trees, int SiteNo);

double		Likelihood(RATES* Rates, TREES *Trees, OPTIONS *Opt);
void		AllocLHInfo(TREES *Trees, OPTIONS *Opt);
void		FreeInvInfo(INVINFO* InvInfo);
INVINFO*	AllocInvInfo(int NOS);

double		LhPraxis(void* PState, double *List);

int			IsNum(double n);

void		CheckRatesVec(double *Vec, int Size);

void		MapConParams(OPTIONS *Opt, RATES *Rates, double *List);


int		SetUpAMatrix(RATES* Rates, TREES *Trees, OPTIONS *Opt);
int		SetAllPMatrix(RATES* Rates, TREES *Trees, OPTIONS *Opt, double Gamma, double Kappa);

#endif
