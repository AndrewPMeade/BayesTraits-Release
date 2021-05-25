#if !defined (LIKLIHEADDER)
#define LIKLIHEADDER


double	Likelihood(RATES* Rates, TREES *Trees, OPTIONS *Opt);
void	AllocLHInfo(TREES *Trees, OPTIONS *Opt);
void	FreeInvInfo(INVINFO* InvInfo);


void	SetUpPrarix(RATES *Rates, TREES *Trees, OPTIONS *Opt);

double	LhPraxisCon(double *List);
int		IsNum(double n);

#endif