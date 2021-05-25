#ifndef PRIOSHEADDER
#define PRIOSHEADDER

void		CreatPriors(OPTIONS *Opt, RATES* Rates);

void		CalcPriors(RATES* Rates, OPTIONS* Opt);
void		FreePriors(RATES* Rates);
void		FreePrior(PRIORS* P);


void		MutatePriors(RATES *Rates, PRIORS **PriosList, int NoOfPriors);
/* void		MutatePriors(PRIORS **PriosList, int NoOfPriors); */
void		MutatePriorsNormal(RATES *Rates, PRIORS **PriosList, int NoOfPriors, double Dev);
void		CopyRatePriors(PRIORS**APriosList, PRIORS **BPriosList, int NoOfPriors);
/*void		CopyMutPriors(PRIORS **APriosList, PRIORS **BPriosList, int NoOfPriors, double Dev); */
void		CopyMutPriors(RANDSTATES *RandStates, PRIORS **APriosList, PRIORS **BPriosList, int NoOfPriors, double Dev);

void		SetRatesToPriors(OPTIONS *Opt, RATES* Rates);

void		CopyPrior(PRIORS *A, PRIORS *B);

double		LnExpP(double x, double Mean, double CatSize);


PRIORS*		CreateBetaPrior(double Mean, double Var);
PRIORS*		CreateGammaPrior(double Mean, double Var);
PRIORS*		CreateUniformPrior(double Min, double Max);
PRIORS*		CreateChiPrior(double Mean, double Var);
PRIORS*		CreateExpPrior(double Mean);
PRIORS*		CreateInvGammaPrior(double Alpha, double Beta);

PRIORDIST	StrToPriorDist(char* Str);

PRIORS*		CreatePrior(int Tokes, char **Passed);

double		RandFromPrior(RANDSTATES *RS, PRIORS* Prior);

#endif
