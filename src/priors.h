#if !defined (PRIOSHEADDER)
#define PRIOSHEADDER

void		CreatPriors(OPTIONS *Opt, RATES* Rates);

void		CalcPriors(RATES* Rates, OPTIONS* Opt);
void		FreePriors(RATES* Rates);
void		FreePrior(PRIORS* P);


double		RateToBetaLh(double Rate, int NoOfCats, double* Prams);

void		CopyPriors(PRIORS **APriosList, PRIORS **BPriosList, int NoOfPriors);
void		MutatePriors(PRIORS **PriosList, int NoOfPriors);
void		MutatePriorsNormal(PRIORS **PriosList, int NoOfPriors, double Dev);
void		CopyRatePriors(PRIORS**APriosList, PRIORS **BPriosList, int NoOfPriors);
void		CopyMutPriors(PRIORS **APriosList, PRIORS **BPriosList, int NoOfPriors, double Dev);
void		SetRatesToPriors(OPTIONS *Opt, RATES* Rates);

void		CopyPrior(PRIORS *A, PRIORS *B);

double		LnExpP(double x, double Mean, double CatSize);


#endif