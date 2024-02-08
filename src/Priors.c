/*
*  BayesTriats 4.0
*
*  copyright 2022
*
*  Andrew Meade
*  School of Biological Sciences
*  University of Reading
*  Reading
*  Berkshire
*  RG6 6BX
*
* BayesTriats is free software: you can redistribute it and/or modify
* it under the terms of the GNU General Public License as published by
* the Free Software Foundation, either version 3 of the License, or
* (at your option) any later version.
* 
* This program is distributed in the hope that it will be useful,
* but WITHOUT ANY WARRANTY; without even the implied warranty of
* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
* GNU General Public License for more details.
* 
* You should have received a copy of the GNU General Public License
* along with this program.  If not, see <http://www.gnu.org/licenses/>
*
*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "TypeDef.h"
#include "GenLib.h"
#include "Priors.h"
#include "RandLib.h"
#include "RevJump.h"
#include "Likelihood.h"
#include "VarRates.h"
#include "RandDists.h"
#include "RandLib.h"
#include "Prob.h"
#include "LocalTransform.h"
#include "GlobalTrend.h"


#include <gsl/gsl_cdf.h>
#include <gsl/gsl_randist.h>



/*
extern double beta(double a, double b);
extern double incbet(double aa, double bb, double xx );

extern double igam(double a, double x);
extern double igamc ( double, double );
extern double igami ( double, double );
extern double gamma(double a);
*/

//extern double chdtr(double df, double x);

int		ValidPriorLh(double LH)
{
	if(LH == LH+1 || LH != LH || LH == ERRLH)
		return FALSE;

	return TRUE;
}


void	FreePrior(PRIOR* P)
{
	if(P->DistVals != NULL)
		free(P->DistVals);

	if(P->HP != NULL)
		free(P->HP);

	if(P->Name != NULL)
		free(P->Name);

	free(P);
}

void		CrateRatePriors(OPTIONS* Opt, RATES* Rates)
{
	Rates->Priors = ClonePriors(Opt->AllPriors, Opt->NoAllPriors);
	Rates->NoPriors = Opt->NoAllPriors;
}

void	FreePriors(RATES *Rates)
{
	int	Index;

	for(Index=0;Index<Rates->NoPriors;Index++)
		FreePrior(Rates->Priors[Index]);

	free(Rates->Priors);
}

PRIOR*	AllocBlankPrior(int NoP)
{
	PRIOR* Ret;

	Ret = (PRIOR*)SMalloc(sizeof(PRIOR));

	Ret->DistVals = NULL;
	if(NoP > 0)
		Ret->DistVals = (double*)SMalloc(sizeof(double) * NoP);
	
	Ret->HP			= NULL;
	Ret->UseHP		= FALSE;
	Ret->Name		= NULL;

	Ret->Discretised= FALSE;
	Ret->Width		= -1.0;

	return Ret;
}


PRIOR*		CreatePrior(char *Name, PRIORDIST PDist, double *PVal)
{
	PRIOR *Ret;
	int NoPram;

	NoPram = DISTPRAMS[PDist];
	Ret = AllocBlankPrior(NoPram);

	memcpy(Ret->DistVals, PVal, sizeof(double) * NoPram);
	
	Ret->Dist = PDist;
	Ret->Name = StrMake(Name);

	return Ret;
}

PRIOR*		CreateGammaPrior(char *Name, double Shape, double Scale)
{
	PRIOR* Ret;

	Ret = AllocBlankPrior(2);
	
	Ret->Dist			= PDIST_GAMMA;
	Ret->DistVals[0]	= Shape;
	Ret->DistVals[1]	= Scale;

	Ret->Name = StrMake(Name);
	
	return Ret;
}

PRIOR*		CreateUniformPrior(char *Name, double Min, double Max)
{
	PRIOR* Ret;

	Ret = AllocBlankPrior(2);
	
	Ret->Dist			= PDIST_UNIFORM;
	Ret->DistVals[0]	= Min;
	Ret->DistVals[1]	= Max;

	Ret->Name = StrMake(Name);
	
	return Ret;
}

PRIOR*		CreateChiPrior(char *Name, double Mean)
{
	PRIOR* Ret;

	Ret = AllocBlankPrior(1);
	
	Ret->Dist			= PDIST_CHI;
	Ret->DistVals[0]	= Mean;
	
	Ret->Name = StrMake(Name);
	
	return Ret;
}

PRIOR*		CreateExpPrior(char *Name, double Alpha)
{
	PRIOR* Ret;

	Ret = AllocBlankPrior(1);
	
	Ret->Dist			= PDIST_EXP;
	Ret->DistVals[0]	= Alpha;

	Ret->Name = StrMake(Name);
		
	return Ret;
}


PRIOR*		CreateSGammaPrior(char *Name, double Alpha, double Beta)
{
	int		NoP;
	PRIOR	*Ret;

	NoP = DISTPRAMS[PDIST_SGAMMA];
	Ret = AllocBlankPrior(NoP);

	Ret->Dist = PDIST_SGAMMA;
	Ret->DistVals[0] = Alpha;
	Ret->DistVals[1] = Beta;
	
	Ret->Name = StrMake(Name);

	return Ret;
}

PRIOR*		CreateLogNormalPrior(char *Name, double Location, double Scale)
{
	int		NoP;
	PRIOR	*Ret;

	NoP = DISTPRAMS[PDIST_LOGNORMAL];
	Ret = AllocBlankPrior(NoP);

	Ret->Dist = PDIST_LOGNORMAL;
	Ret->DistVals[0] = Location;
	Ret->DistVals[1] = Scale;

	Ret->Name = StrMake(Name);

	return Ret;
}

PRIOR*		CreateNormalPrior(char *Name, double Mean, double SD)
{
	int		NoP;
	PRIOR	*Ret;

	if(SD <= 0)
	{
		printf("prior %s, normal distribution SD must be >0", Name);
		exit(0);
	}

	NoP = DISTPRAMS[PDIST_NORMAL];
	Ret = AllocBlankPrior(NoP);

	Ret->Dist = PDIST_NORMAL;
	Ret->DistVals[0] = Mean;
	Ret->DistVals[1] = SD;

	Ret->Name = StrMake(Name);

	return Ret;
}

PRIOR* CreateWeibullPrior(char* Name, double Scale, double Exponent)
{
	int		NoP;
	PRIOR* Ret;

	if(Scale <= 0 || Exponent <= 0)
	{
		printf("prior %s, Weibull distribution Scale and Exponent must be >0", Name);
		exit(0);
	}

	NoP = DISTPRAMS[PDIST_WEIBULL];
	Ret = AllocBlankPrior(NoP);

	Ret->Dist = PDIST_WEIBULL;
	Ret->DistVals[0] = Scale;
	Ret->DistVals[1] = Exponent;

	Ret->Name = StrMake(Name);

	return Ret;
}


PRIOR*		CrateUndefinedPrior(char *Name)
{
	PRIOR* Ret;

	Ret = AllocBlankPrior(0);
	Ret->Dist = PDIST_UNDEFINED;
	Ret->Name = StrMake(Name);

	return Ret;
}

void		SetHPDistParam(int Pos, PRIOR* Prior)
{
	int OS;
	double Min, Max;

	OS = Pos * 2;

	Min = Prior->HP[OS];
	Max = Prior->HP[OS+1];

	if(Min > Max)
	{
		printf("Hyper prior %s min is larger than max %f", Prior->Name, Max);
		exit(0);
	}

	Prior->DistVals[Pos] = Min + ((Max - Min) * 0.5);
}

PRIOR*		CreateHyperPrior(char *Name, PRIORDIST PDist, double *PVal)
{
	PRIOR* Ret;
	int		Index, NoParam, NoHPParam;

	NoParam = DISTPRAMS[PDist];
	Ret = AllocBlankPrior(NoParam);

	Ret->Name = StrMake(Name);
	Ret->Dist = PDist;
	Ret->UseHP = TRUE;	
		
	NoHPParam = NoParam * 2;

	Ret->HP = (double*)SMalloc(sizeof(double) * NoHPParam);
	memcpy(Ret->HP, PVal, sizeof(double) * NoHPParam);

	for(Index=0;Index<NoParam;Index++)
		SetHPDistParam(Index, Ret);

	return Ret;	
}

PRIOR*		ClonePrior(PRIOR *Prior)
{
	PRIOR *Ret;

	if(Prior->UseHP == FALSE)
		Ret = CreatePrior(Prior->Name, Prior->Dist, Prior->DistVals);
	else
		Ret = CreateHyperPrior(Prior->Name, Prior->Dist, Prior->HP);

	Ret->Discretised = Prior->Discretised;
	Ret->Width = Prior->Width;

	return Ret;
}

PRIOR**	ClonePriors(PRIOR **PList, int NoPriors)
{
	int Index;
	PRIOR **Ret;

	if(NoPriors == 0)
		return NULL;
		

	Ret = (PRIOR**)SMalloc(sizeof(PRIOR*) * NoPriors);
	
	for(Index=0;Index<NoPriors;Index++)
		Ret[Index] = ClonePrior(PList[Index]);

	return Ret;
}

int		FindCat(double x, double CatWidth, int NoOfCats)
{
	double Ret;

	Ret = x / CatWidth;

	Ret = floor(Ret);

	if(Ret>=NoOfCats)
		Ret = NoOfCats-1;

	return (int)Ret;
}


/*
double	FindBetaBeta(double Mean, double Var, double Scale)
{
	double	Ret=0;

	Ret = Mean - Scale;
	Ret = Ret * Ret;
	Ret = Mean * Ret;
	Ret = Ret / Var;
	Ret = (Ret + Mean) - Scale;
	Ret = Ret / Scale;

	return Ret;
}

double FindBetaAlpha(double Mean, double Scale, double Beta)
{
	double Ret=0;

	Ret = Mean * Beta;
	Ret = Ret / (Scale - Mean);

	return Ret;
}

*/

/*
double	RateToBetaLh(double Rate, int NoOfCats, double* Prams)
{
	double	Alpha;
	double	Beta;
	double	Scale=100;
	double	CatWidth;
	int		Cat;
	double	Ret;
	double	Mean;
	double	Var;

	if(Rate >= Scale)
		return 0;

	Mean	= Prams[0];
	Var		= Prams[1];

	Beta	=	FindBetaBeta(Mean, Var, Scale);
	Alpha	=	FindBetaAlpha(Mean, Scale,  Beta);

	if((Beta < 0) || (Alpha < 0))
	{
		return 0.0;
	}

	Rate = Rate / Scale;

	CatWidth=(double)1/NoOfCats;
	Cat = FindCat(Rate, CatWidth, NoOfCats);

	Ret = BetaProb(Cat, Alpha, Beta, CatWidth);

	return Ret;
}
*/


/*
double	FindBetaBeta(double Mue, double Sigma)
{
	double	Ret=0;

	Ret = Mue * ((Mue - 1) * (Mue - 1));
	Ret = Ret / Sigma;
	Ret = Ret + (Mue - 1);

	return Ret;
}

double FindBetaAlpha(double Beta, double Mue)
{
	double Ret=0;

	Ret = Mue * Beta;
	Ret = Ret / (1-Mue);

	return Ret;
}


double	RateToBetaLh(double Rate, int NoOfCats, double* Prams)
{
	double	Alpha;
	double	Beta;
	double	Scale;
	double	CatWidth;
	int		Cat;
	double	Ret;
	double	Mean;
	double	Var;
	double	Max;

	Mean	= Prams[0];
	Var		= Prams[1];

	Max		= (sqrt(Var) * 4) + Mean;

	Mean	=	Mean / Max;
	Var		=	Var / (Max * Max);

	Beta	=	FindBetaBeta(Mean, Var);
	Alpha	=	FindBetaAlpha(Beta, Mean);

	if((Beta < 0) || (Alpha < 0))
	{
		return 0;
	}

	Scale = Alpha / (Alpha + Beta);
	Scale = Mean / Scale;

	Rate = Rate / Max;
	Rate = Rate / Scale;

	CatWidth=(double)1/NoOfCats;
	Cat = FindCat(Rate, CatWidth, NoOfCats);

	Ret = BetaProb(Cat, Alpha, Beta, CatWidth);
	return Ret;
}
*/


double	LogLogNormalP(double X, PRIOR *Prior)
{
	double Ret, A, B;

	if(X < 0.0)
		return ERRLH;

	if(Prior->Discretised == FALSE)
		Ret = gsl_ran_lognormal_pdf(X, Prior->DistVals[0], Prior->DistVals[1]);
	else
	{
		A = gsl_cdf_lognormal_P (X, Prior->DistVals[0], Prior->DistVals[1]);
		B = gsl_cdf_lognormal_P (X+Prior->Width, Prior->DistVals[0], Prior->DistVals[1]);
		Ret = B - A;
	}

	Ret = log(Ret);

	if(ValidPriorLh(Ret) == FALSE)
		return ERRLH;

	return Ret;
}


double	LogNormalP(double X, PRIOR *Prior)
{
	double Ret, A, B;
	
	X = X - Prior->DistVals[0];
	
	if(Prior->Discretised == FALSE)
		Ret = gsl_ran_gaussian_pdf(X, Prior->DistVals[1]);
	else
	{
		A = gsl_cdf_gaussian_P(X, Prior->DistVals[1]);
		B = gsl_cdf_gaussian_P(X + Prior->Width, Prior->DistVals[1]);
		Ret = B - A;
	}
	Ret = log(Ret);

	if(ValidPriorLh(Ret) == FALSE)
		return ERRLH;

	return Ret;
}

void	TestGammaPrior(PRIOR *Prior)
{
	double X, LogLh;

	Prior->DistVals[0] = 9.0;
	Prior->DistVals[1] = 1.0;

	for(X=0;X<20;X+=0.01)
	{
		LogLh = gsl_cdf_gamma_P(X, Prior->DistVals[0], Prior->DistVals[1]);
		LogLh = gsl_ran_gamma_pdf(X, Prior->DistVals[0], Prior->DistVals[1]);
		printf("%f\t%f\n", X, LogLh);
	}

	exit(0);
}

double	LogGammaP(double X, PRIOR *Prior)
{
	double Ret, A, B;

	if(X < 0.0)
		return ERRLH;


	if(Prior->Discretised == FALSE)
		Ret = gsl_ran_gamma_pdf(X, Prior->DistVals[0], Prior->DistVals[1]);
	else
	{
		A = gsl_cdf_gamma_P(X, Prior->DistVals[0], Prior->DistVals[1]);
		B = gsl_cdf_gamma_P(X+Prior->Width, Prior->DistVals[0], Prior->DistVals[1]);
		Ret = B - A;
	}

	return log(Ret);
}

double	LogWeibullP(double X, PRIOR* Prior)
{
	double Ret;

	if(X < 0.0)
		return ERRLH;
 
	Ret = 1.0;
	if(Prior->Discretised == FALSE)
		Ret = gsl_ran_weibull_pdf(X, Prior->DistVals[0], Prior->DistVals[1]);
	else
	{
		printf("Disceretised weiball is not supported.\n");
		exit(1);
	}

	return log(Ret);
}

double LogUniP(double X, PRIOR *Prior)
{
	double Ret;
	
	if(X < Prior->DistVals[0] || X > Prior->DistVals[1])
		return ERRLH;

	Ret = 1.0 / (Prior->DistVals[1] - Prior->DistVals[0]);

	return log(Ret);
}


double LogChiSquaredP(double X, PRIOR *Prior)
{
	double Ret, A, B;

	if(X < 0.0)
		return ERRLH;

	if(Prior->Discretised == FALSE)
		Ret = gsl_ran_chisq_pdf(X, Prior->DistVals[0]);
	else
	{
		A = gsl_cdf_chisq_P(X, Prior->DistVals[0]);
		B = gsl_cdf_chisq_P(X+Prior->Width, Prior->DistVals[0]);
		Ret = B - A;
	}

	return log(Ret);
}

void	PLhTest(PRIOR *Prior)
{
	double X, Lh;

	for(X=0;X<=10;X+=0.001)
	{
		Lh  = PDFSGamma(X, Prior->DistVals[0], Prior->DistVals[1]);

		printf("%f\t%f\n", X, Lh);
	}

	exit(0);
}

double LogSGammaP(double X, PRIOR *Prior)
{
	double Ret;


//	PLhTest(Prior);

	if(X < 0.0)
		return ERRLH;

	if(Prior->Discretised == TRUE)
	{
		printf("Discretised scaled gamma is not supoorted\n");
		exit(1);
	}

	Ret = PDFSGamma(X, Prior->DistVals[0], Prior->DistVals[1]);

	return log(Ret);
}

double	LogExpContinuous(double X, PRIOR *Prior)
{
	double A, B, Alpha;

	Alpha = 1.0 / Prior->DistVals[0];

	A = -Alpha * X;
	B = log(Alpha);
	return A + B;
}

double LogExpP(double X, PRIOR *Prior)
{
	double Ret, A, B, Alpha;

	if(X < 0.0)
		return ERRLH;

	Alpha = Prior->DistVals[0];

	if(Prior->Discretised == FALSE)
//		Ret = gsl_ran_exponential_pdf(X, Alpha);
		return LogExpContinuous(X, Prior);
	else
	{
		A = gsl_cdf_exponential_P(X, Alpha);
		B = gsl_cdf_exponential_P(X+Prior->Width, Alpha);
		Ret = B - A;
	}

	return log(Ret);
}



void	ChiSTest(void)
{
	PRIOR *Prior;
	double X, P;

	Prior = CreateExpPrior("SG", 1.0);
	Prior = CreateChiPrior("SG", 9);

	for(X=0.000001;X<8;X+=0.01)
	{
		P = LogChiSquaredP(X, Prior);
		P = exp(P);
		printf("%f\t%f\n", X, P);
	}

	exit(0);
}

void	PriorTest(void)
{
	PRIOR* Prior;
	double X, P;

	Prior = CreateLogNormalPrior("Prior", -0.6, 0.55);

	for(X=0.00001;X<5;X+=0.01)
	{
		P = LogLogNormalP(X, Prior);
		printf("%f\t%f\n", X, P);
	}

	exit(0);
}

double	CalcLhPriorP(double X, PRIOR *Prior)
{
	double Ret;

//	PriorTest();
	
	switch(Prior->Dist)
	{
		case PDIST_GAMMA:
			Ret = LogGammaP(X, Prior);
		break;

		case PDIST_UNIFORM:
			Ret = LogUniP(X, Prior);
		break;

		case PDIST_EXP:
			Ret = LogExpP(X, Prior);
		break;

		case PDIST_CHI:
			Ret = LogChiSquaredP(X, Prior);
		break;

		case PDIST_SGAMMA:
			Ret = LogSGammaP(X, Prior);
		break;

		case PDIST_LOGNORMAL:
			Ret = LogLogNormalP(X, Prior);
		break;

		case PDIST_NORMAL:
			Ret = LogNormalP(X, Prior);
		break;

		case PDIST_WEIBULL:
			Ret = LogWeibullP(X, Prior);

	}

	return Ret;
}

double		RandFromPrior(gsl_rng *RNG, PRIOR *Prior)
{
	switch(Prior->Dist)
	{
		case PDIST_GAMMA:
			return gsl_ran_gamma(RNG, Prior->DistVals[0], Prior->DistVals[1]);
	
		case PDIST_UNIFORM:
			return gsl_ran_flat(RNG,  Prior->DistVals[0], Prior->DistVals[1]);
		
		case PDIST_EXP:
			return gsl_ran_exponential(RNG, Prior->DistVals[0]);

		case PDIST_CHI:
			return gsl_ran_chisq(RNG, Prior->DistVals[0]);

		case PDIST_SGAMMA:
			return gsl_ran_gamma(RNG, Prior->DistVals[0], Prior->DistVals[1]);

		case PDIST_LOGNORMAL:
			return gsl_ran_lognormal(RNG, Prior->DistVals[0], Prior->DistVals[1]);

		case PDIST_NORMAL:
			return gsl_ran_gaussian_ziggurat(RNG, Prior->DistVals[1]) + Prior->DistVals[0];

		case PDIST_WEIBULL:
			return gsl_ran_weibull(RNG, Prior->DistVals[0], Prior->DistVals[1]);
	}

	printf("Prior dist not found for %s.\n", Prior->Name);
	exit(1);

	return -1.0;
}



double	CalcTreeTransPrior(RATES *Rates, OPTIONS *Opt)
{
	double PLh, Ret;
	PRIOR	*Prior;

	Ret = 0;

	if(Opt->EstKappa == TRUE)
	{
		Prior = GetPriorFromName("Kappa", Rates->Priors, Rates->NoPriors);
		PLh = CalcLhPriorP(Rates->Kappa, Prior);
		if(PLh == ERRLH)
			return ERRLH;
		Ret += PLh;
	}

	if(Opt->EstLambda == TRUE)
	{
		Prior = GetPriorFromName("Lambda", Rates->Priors, Rates->NoPriors);
		PLh = CalcLhPriorP(Rates->Lambda, Prior);
		if(PLh == ERRLH)
			return ERRLH;
		Ret += PLh;
	}

	if(Opt->EstDelta == TRUE)
	{
		Prior = GetPriorFromName("Delta", Rates->Priors, Rates->NoPriors);
		PLh = CalcLhPriorP(Rates->Delta, Prior);
		if(PLh == ERRLH)
			return ERRLH;
		Ret += PLh;
	}

	if(Opt->EstOU == TRUE)
	{
		Prior = GetPriorFromName("OU", Rates->Priors, Rates->NoPriors);
		PLh = CalcLhPriorP(Rates->OU, Prior);
		if(PLh == ERRLH)
			return ERRLH;
		Ret += PLh;
	}

	if(Opt->EstGamma == TRUE)
	{
		Prior = GetPriorFromName("Gamma", Rates->Priors, Rates->NoPriors);
		PLh = CalcLhPriorP(Rates->Gamma, Prior);
		if(PLh == ERRLH)
			return ERRLH;
		Ret += PLh;
	}
	
	return Ret;
}



double CalcRJDummyPriors(OPTIONS *Opt, RATES* Rates)
{
	RJDUMMY *RJDummy;
	DUMMYCODE *DC;
	int		Index;
	double	Ret, P;

	Ret = 0;
	RJDummy = Rates->RJDummy;

	for(Index=0;Index<RJDummy->NoDummyCode;Index++)
	{
		DC = RJDummy->DummyList[Index];
		
		P = gsl_ran_gaussian_pdf(DC->Beta[0], 1.0);

		if(DC->Type == RJDUMMY_INTER_SLOPE)
			P *= gsl_ran_gaussian_pdf(DC->Beta[1], 1.0);
				
		Ret += log(P);
	}

	return Ret;
}

double	CalcRatePrior(RATES* Rates, OPTIONS* Opt)
{
	int NoRatePriors, Index;
	double PLh, R, Ret;
	PRIOR *Prior;

	Prior = NULL;
	NoRatePriors = Rates->NoOfRates;
	if(Opt->UseRJMCMC == TRUE)
	{
		Prior = GetPriorFromName("RJRates", Rates->Priors, Rates->NoPriors);
		NoRatePriors = Rates->NoOfRJRates;
	}
	
	Ret = 0;

	for(Index=0;Index<NoRatePriors;Index++)
	{
		R = Rates->Rates[Index];
		if(Opt->UseRJMCMC == FALSE)
			Prior = GetPriorFromName(Rates->RateNames[Index], Rates->Priors, Rates->NoPriors);
	
		PLh = CalcLhPriorP(R, Prior);

		if(PLh == ERRLH || IsNum(PLh) == FALSE)
			return ERRLH;

		Ret += PLh;
	}

	return Ret;
}

double	CalcRegVarPrior(RATES* Rates)
{
	double Ret, Var;
	PRIOR *Prior;

	Var = Rates->Contrast->GlobalVar;

	Prior = GetPriorFromName("Var", Rates->Priors, Rates->NoPriors);

	if(Prior == NULL)
		return 0;

	Ret = CalcLhPriorP(Var, Prior);

	return Ret;
}

double	CalcTimeSlicePriors(RATES *Rates)
{
	int Index;
	PRIOR *PTime, *PScale;
	TIME_SLICES *TSlices;
	TIME_SLICE	*TS;
	double PVal, Ret;


	Ret = 0;

	TSlices = Rates->TimeSlices;

	PTime = GetPriorFromName("TimeSlice-Time", Rates->Priors, Rates->NoPriors);
	PScale = GetPriorFromName("TimeSlice-Scale", Rates->Priors, Rates->NoPriors);

	for(Index=0;Index<TSlices->NoTimeSlices;Index++)
	{
		TS = TSlices->TimeSlices[Index];
		PVal = 0;

		if(TS->FixedTime == FALSE)
		{
			PVal = CalcLhPriorP(TS->Time, PTime);
			if(PVal == ERRLH)
				return PVal;
		}

		if(TS->FixedScale == FALSE)
		{
			PVal += CalcLhPriorP(TS->Scale, PScale);
			if(PVal == ERRLH)
				return PVal;
		}
				
		Ret += PVal;
	}
	return Ret;
}

double	CaclAnsStatePriors(RATES *Rates, OPTIONS *Opt)
{
	int Index, SiteNo;
	PRIOR *Prior;
	double Ret, PLh;

	if(Opt->DataType == DISCRETE)
		return 0.0;

	Ret = 0;
	for(Index=0;Index<Rates->NoEstData;Index++)
	{
		SiteNo = Rates->EstDataSiteNo[Index];
		Prior = GetAnsStatePrior(Opt->Model, SiteNo, Rates->Priors, Rates->NoPriors);
		
		PLh = CalcLhPriorP(Rates->EstData[Index], Prior);
		
		if(PLh == ERRLH)
			return ERRLH;		

		Ret += PLh;
	}

	return Ret;
}

double	CaclNormQPrior(RATES* Rates, OPTIONS* Opt)
{
	double Ret;
	PRIOR *Prior;

	Prior = GetPriorFromName("GlobalRate", Rates->Priors, Rates->NoPriors);
	Ret = CalcLhPriorP(Rates->GlobablRate, Prior);

	return Ret;
}

double CaclCVPrior(OPTIONS *Opt,RATES *Rates)
{
	double Ret,POnToOff;
	PRIOR *Prior;

	Prior = GetPriorFromName("CVSwichRate", Rates->Priors, Rates->NoPriors);

	Ret = CalcLhPriorP(Rates->OffToOn, Prior);
	if(Ret == ERRLH)
		return ERRLH;

	if(Rates->OffToOn == Rates->OnToOff)
		return Ret;

	POnToOff = CalcLhPriorP(Rates->OnToOff, Prior);
	if(POnToOff == ERRLH)
		return ERRLH;

	return POnToOff + Ret;
}

void	CalcPriors(RATES* Rates, OPTIONS* Opt)
{
	double	PLh, Ret;
		
	Ret = 0;
	PLh = 0;	

	Rates->LhPrior = ERRLH;

	if(Opt->LoadModels == FALSE)
		PLh = CalcRatePrior(Rates, Opt);

	if(PLh == ERRLH)
		return;

	Ret += PLh;
	
	PLh = CalcTreeTransPrior(Rates, Opt);
	if(PLh == ERRLH)
		return;
	Ret += PLh;

	PLh = CaclAnsStatePriors(Rates, Opt);
	if(PLh == ERRLH)
		return;
	Ret += PLh;
	
	if(Opt->Model == M_CONTRAST_REG && Opt->NoLh == FALSE)
	{
		PLh = CalcRegVarPrior(Rates);
		if(PLh == ERRLH)
			return;
		Ret += PLh;
	}

	if(Opt->RJDummy == TRUE)
	{
		printf("Must check RJ dummy priors.\n");
		exit(0);
		PLh = CalcRJDummyPriors(Opt, Rates);
		if(PLh == ERRLH)
			return;
		Ret += PLh;
	}
	
	if(UseNonParametricMethods(Opt) == TRUE)
	{
		PLh = CalcVarRatesPriors(Rates, Opt);

		if(PLh == ERRLH)
			return;

		Ret += PLh;
	}

	if(Rates->NoLocalTransforms > 0)
	{
		PLh = CaclLocalTransformsPrior(Rates);
		
		if(PLh == ERRLH)
			return;

		Ret += PLh;
	}


	if(Rates->TimeSlices != NULL)
	{
		PLh = CalcTimeSlicePriors(Rates);
		
		if(PLh == ERRLH)
			return;

		Ret += PLh;
	}

	if(Opt->NormQMat == TRUE)
	{
		PLh = CaclNormQPrior(Rates, Opt);

		if(PLh == ERRLH)
			return;

		Ret += PLh;
	}

	if(Opt->UseGlobalTrend == TRUE)
	{
		PLh = CalcGlobalTrendPrior(Rates);
		
		if(PLh == ERRLH)
			return;

		Ret += PLh;
	}

	if(Opt->UseCovarion == TRUE)
	{
		PLh = CaclCVPrior(Opt, Rates);

		if(PLh == ERRLH)
			return;

		Ret += PLh;
	}

	Rates->LhPrior = Ret;
}

void	MutatePriors(RATES *Rates, PRIOR **PriosList, int NoOfPriors)
{
	int PIndex;

	for(PIndex=0;PIndex<NoOfPriors;PIndex++)
	{
		PriosList[PIndex]->DistVals[0] = RandDouble(Rates->RS) * 100;
		PriosList[PIndex]->DistVals[1] = RandDouble(Rates->RS) * 100;
	}
}

double	ChangePriorNorm(RATES *Rates, double Val, double Dev, double Min, double Max)
{
	double	Ret=0;

	do
	{
		Ret = RandNormal(Rates->RS, Val, Dev);
	} while(Ret > Max || Ret < Min);

	return Ret;
}

void	MutatePriorsNormal(RATES *Rates, PRIOR **PriosList, int NoOfPriors, double Dev)
{
	int		PIndex;
	int		RIndex;
	PRIOR	*Prior;
	double	Min;
	double	Max;

	for(PIndex=0;PIndex<NoOfPriors;PIndex++)
	{
		Prior = PriosList[PIndex];

		if(Prior->UseHP == TRUE)
		{
			for(RIndex=0;RIndex<DISTPRAMS[Prior->Dist];RIndex++)
			{
				Min = Prior->HP[RIndex*2];
				Max = Prior->HP[(RIndex*2)+1];


				Prior->DistVals[RIndex] = ChangePriorNorm(Rates, Prior->DistVals[RIndex], Dev, Min, Max);
			//	Prior->DistVals[RIndex] = RandUniDouble(Rates->RS, Min, Max);
			}
		}
	}
}

double ChangePrior(RANDSTATES *RandStates, double Rate, double dev)
{
	double	Ret;

	Ret = (RandDouble(RandStates) * dev) - (dev / 2.0); 

	Ret += Rate;

	if(Ret > 100)
		Ret = 100 - (Ret - 100);

	if(Ret < 0)
		Ret = -Ret;

	return Ret;
}


void	CopyPrior(PRIOR *A, PRIOR *B)
{
	A->Dist		= B->Dist;
	A->UseHP	= B->UseHP;

	memcpy(A->DistVals, B->DistVals, sizeof(double) * DISTPRAMS[A->Dist]);

	if(B->UseHP == TRUE)
		memcpy(A->HP, B->HP, sizeof(double) * DISTPRAMS[A->Dist] * 2);

	A->Discretised = B->Discretised;
	A->Width = B->Width;
}

PRIORDIST	StrToPriorDist(char* Str, int *Err)
{
	int			Index;

	MakeLower(Str);

	for(Index=0;Index<NO_PRIOR_DIST;Index++)
	{
		if(strcmp(Str, DISTNAMES[Index])==0)
		{
			*Err = FALSE;
			return (PRIORDIST)(Index);
		}
	}

	*Err = TRUE;
	return (PRIORDIST)(0);
}

int			CheckPriorDistVals(PRIORDIST PDist, int Tokes, char **Passed)
{
	int Index;
	double P;

	for(Index=0;Index<Tokes;Index++)
	{
		if(IsValidDouble(Passed[Index]) == FALSE)
		{
			printf("Cannot convert %s to a valid prior parameter.\n", Passed[Index]);
			exit(1);
		}

		P = atof(Passed[Index]);
		
		if(!(PDist == PDIST_UNIFORM || PDist == PDIST_NORMAL || PDist == PDIST_LOGNORMAL))
		{
			if(P < 0)
			{
				printf("Prior parameters values must be greater than 0, value %f is invalid.\n", P);
				exit(1);
			}
		}
	}

	return TRUE;
}

double*		MakePriorParam(int Tokes, char **Passed)
{
	double *Ret;
	int Index;

	Ret = (double*)SMalloc(sizeof(double) * Tokes);

	for(Index=0;Index<Tokes;Index++)
		Ret[Index] = atof(Passed[Index]);

	return Ret;
}

PRIOR*		CreatePriorFromStr(char *Name, int Tokes, char **Passed)
{
	PRIORDIST	PD;
	double		*PVal;
	PRIOR		*Ret;
	int			Err;

	Ret = NULL;

	if(Tokes != 2 && Tokes != 3)
	{
		printf("Prior requires a distribution name, (beta, gamma, uniform, chi, exp, invgamma, normal) and distribution parameters.\n");
		exit(1);
	}

	StrToPriorDist(Passed[0], &Err);
	
	if(Err == TRUE)
	{
		printf("Invalid prior distribution name. Valid names are beta, gamma, uniform, chi, exp, invgamma, normal.\n");
		exit(1);
	}

	PD = StrToPriorDist(Passed[0], &Err);

	if(Tokes -1 != DISTPRAMS[PD])
	{
		printf("Prior %s (%s) requires %d parameters.\n", Name, DISTNAMES[PD], DISTPRAMS[PD]);
		exit(0);
	}

	if(CheckPriorDistVals(PD, Tokes-1, &Passed[1]) == FALSE)
		exit(0);

	PVal = MakePriorParam(Tokes-1, &Passed[1]);

	if(PD == PDIST_GAMMA)
		Ret = CreateGammaPrior(Name, PVal[0], PVal[1]);

	if(PD == PDIST_UNIFORM)
		Ret = CreateUniformPrior(Name, PVal[0], PVal[1]);

	if(PD == PDIST_CHI)
		Ret = CreateChiPrior(Name, PVal[0]);

	if(PD == PDIST_EXP)
		Ret = CreateExpPrior(Name, PVal[0]);

	if(PD == PDIST_LOGNORMAL)
		Ret = CreateLogNormalPrior(Name, PVal[0], PVal[1]);

	if(PD == PDIST_NORMAL)
		Ret = CreateNormalPrior(Name, PVal[0], PVal[1]);

	if(PD == PDIST_SGAMMA)
		Ret = CreateSGammaPrior(Name, PVal[0], PVal[1]);

	if(PD == PDIST_WEIBULL)
		Ret = CreateWeibullPrior(Name, PVal[0], PVal[1]);

	free(PVal);

	return Ret;
}

PRIOR*		CreateHyperPriorFromStr(char *Name, int Tokes, char **Passed)
{
	PRIOR		*Ret;
	PRIORDIST	PDist;
	int			NoParam;
	double		*PVal;
	int			Err;

	if(Tokes < 3)
	{
		printf("A hyper prior require a distribution and min max values for each parameters");
		exit(0);
	}

	PDist = StrToPriorDist(Passed[0], &Err);

	if(Err == TRUE)
	{
		printf("Invalid prior distribution name,. valid names are beta, gamma, uniform, chi, exp, invgamma.\n");
		exit(1);
	}

	NoParam = DISTPRAMS[PDist] * 2;

	if(NoParam != Tokes - 1)
	{
		printf("The %s hyper prior require %d distribution parameter %d supplied.\n", Passed[0], NoParam, Tokes - 1);
		exit(0);
	}
	
	if(CheckPriorDistVals(PDist, Tokes-1, &Passed[1]) == FALSE)
		return NULL;

	PVal = MakePriorParam(Tokes-1, &Passed[1]);
		
	Ret = CreateHyperPrior(Name, PDist, PVal);

	free(PVal);

	return Ret;
}

PRIOR*		GetPriorFromName(char *Name, PRIOR** PList, int NoPrior)
{
	int Index;

	for(Index=0;Index<NoPrior;Index++)
		if(StrICmp(Name, PList[Index]->Name) == 0)
			return PList[Index];

	return NULL;
}

PRIOR*		GetAnsStatePrior(MODEL Model, int SiteNo, PRIOR** PList, int NoPrior)
{
	char *Buffer;
	PRIOR *Ret;

	Buffer = (char*)SMalloc(sizeof(char) * 128);
	if(Model == M_CONTINUOUS_REG && SiteNo == -1)
		sprintf(Buffer, "AncState-Dep");
	else
		sprintf(Buffer, "AncState-%d", SiteNo+1);

	Ret = GetPriorFromName(Buffer, PList, NoPrior);

	free(Buffer);

	return Ret;
}

int		GetPriorPosFromName(char *Name, PRIOR** PList, int NoPrior)
{
	int Index;

	for(Index=0;Index<NoPrior;Index++)
		if(StrICmp(Name, PList[Index]->Name) == 0)
			return Index;

	return -1;
}

void		AddPriorToOpt(OPTIONS *Opt, PRIOR *Prior)
{
	if(GetPriorFromName(Prior->Name, Opt->AllPriors, Opt->NoAllPriors) != NULL)
	{
		printf("prior %s allready exists", Prior->Name);
		exit(1);
	}

	Opt->AllPriors = (PRIOR**)AddToList(&Opt->NoAllPriors, (void**)Opt->AllPriors, (void*)Prior);
}

void		RemovePriorFormOpt(char *Name, OPTIONS *Opt)
{
	PRIOR **NPList, *Prior;
	int Index, Pos;

	Prior = GetPriorFromName(Name, Opt->AllPriors, Opt->NoAllPriors);
	if(Prior == NULL)
		return;

	NPList = (PRIOR**)SMalloc(sizeof(PRIOR*) * Opt->NoAllPriors);

	Pos = 0;
	for(Index=0;Index<Opt->NoAllPriors;Index++)
	{
		if(strcmp(Name, Opt->AllPriors[Index]->Name) != 0)
			NPList[Pos++] = Opt->AllPriors[Index];
		else
			Prior = Opt->AllPriors[Index];
	}

	FreePrior(Prior);

	free(Opt->AllPriors);

	if(Pos == 0)
	{
		free(NPList);
		Opt->AllPriors = NULL;
	}
	else
		Opt->AllPriors = NPList;
	
	Opt->NoAllPriors = Pos;
}

void	ReplacePrior(OPTIONS *Opt, PRIOR *Prior)
{
	int Pos;

	Pos = GetPriorPosFromName(Prior->Name, Opt->AllPriors, Opt->NoAllPriors);
	if(Pos != -1)
		Opt->AllPriors[Pos] = Prior;
	else
	{
		printf("Cannot find prior name %s.\n", Prior->Name);
		exit(0);
	}
}

double	CalcNormalHasting(double x, double SD)
{
	double Ret;

	Ret = gsl_cdf_gaussian_P(x, SD);
	//	Ret = ndtr(x/SD);

	return log(Ret);
}

void TestPrior(PRIOR *Prior, size_t NoSamples)
{
	gsl_rng *rng;
	size_t Index;
	
	rng = gsl_rng_alloc(gsl_rng_mt19937);
	printf("Sample form prior: %s\n", Prior->Name);
	
	for(Index=0;Index<NoSamples;Index++)
		printf("%zu\t%f\n", Index, RandFromPrior(rng, Prior));	

	gsl_rng_free(rng);
}

