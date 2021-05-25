#ifndef PRIOR_H
#define PRIOR_H

#include <gsl/gsl_randist.h>

void		AddPriorToOpt(OPTIONS *Opt, PRIOR *Prior);
void		RemovePriorFormOpt(char *Name, OPTIONS *Opt);
void		ReplacePrior(OPTIONS *Opt, PRIOR *Prior);

PRIOR*		GetPriorFromName(char *Name, PRIOR** PList, int NoPrior);
PRIOR*		ClonePrior(PRIOR *Prior);
PRIOR**		ClonePriors(PRIOR **PList, int NoPriors);

double		RandFromPrior(gsl_rng *RNG, PRIOR *Prior);

double		CalcLhPriorP(double X, PRIOR *Prior);
void		CalcPriors(RATES* Rates, OPTIONS* Opt);

void		CrateRatePriors(OPTIONS* Opt, RATES* Rates);
void		FreePriors(RATES* Rates);

void		MutatePriors(RATES *Rates, PRIOR **PriosList, int NoOfPriors);
void		MutatePriorsNormal(RATES *Rates, PRIOR **PriosList, int NoOfPriors, double Dev);

void		CopyPrior(PRIOR *A, PRIOR *B);

PRIOR*		CreateGammaPrior(char *Name, double Shape, double Scale);
PRIOR*		CreateUniformPrior(char *Name, double Min, double Max);
PRIOR*		CreateChiPrior(char *Name, double Mean);
PRIOR*		CreateExpPrior(char *Name, double Alpha);
PRIOR*		CreateSGammaPrior(char *Name, double Alpha, double Beta);
PRIOR*		CreateLogNormalPrior(char *Name, double Location, double Scale);
PRIOR*		CreateNormalPrior(char *Name, double Mean, double SD);


PRIOR*		CreatePrior(char *Name, PRIORDIST PDist, double *PVal);
PRIOR*		CreateHyperPrior(char *Name, PRIORDIST PDist, double *PVal);

void		FreePrior(PRIOR* P);

PRIORDIST	StrToPriorDist(char* Str);

PRIOR*		CreatePriorFromStr(char *Name, int Tokes, char **Passed);
PRIOR*		CreateHyperPriorFromStr(char *Name, int Tokes, char **Passed);

#endif
