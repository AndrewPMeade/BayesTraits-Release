#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "typedef.h"
#include "genlib.h"
#include "priors.h"
#include "RandLib.h"
#include "revjump.h"
#include "likelihood.h"
#include "VarRates.h"
#include "randdists.h"
#include "RandLib.h"
#include "Prob.h"
#include "gsl\gsl_cdf.h"



/*
extern double beta(double a, double b);
extern double incbet(double aa, double bb, double xx );

extern double igam(double a, double x);
extern double igamc ( double, double );
extern double igami ( double, double );
extern double gamma(double a);
*/

extern double chdtr(double df, double x);

int		ValidPriorLh(double LH)
{
	if(LH == LH+1 || LH != LH || LH == ERRLH)
		return FALSE;

	return TRUE;
}

double		CalcPriorPDFLh(double X, PRIOR *Prior)
{
	double Ret;

	Ret = 0.0;

	if(Prior->Dist == UNIFORM)
	{
		if(X < Prior->DistVals[0] || X > Prior->DistVals[1])
			return ERRLH;

		Ret = 1.0 / (Prior->DistVals[1]  - Prior->DistVals[0]);
		return log(Ret);
	}
	
	if(X < 0.0)
		return ERRLH;
	
	if(Prior->Dist == GAMMA)
		Ret = PDFGamma(X, Prior->DistVals[0], Prior->DistVals[1]);

	if(Prior->Dist == CHI)
		Ret = PDFChi(X, Prior->DistVals[0], Prior->DistVals[1]);
	
	if(Prior->Dist == EXP)
		Ret = PDFExp(X, Prior->DistVals[0]);
	
	if(Prior->Dist == INVGAMMA)
		Ret = PDFInvGamma(X, Prior->DistVals[0], Prior->DistVals[1]);

	if(Prior->Dist == SGAMMA)
		Ret = PDFSGamma(X, Prior->DistVals[0], Prior->DistVals[1]);

	Ret = log(Ret);

	if(ValidPriorLh(Ret) == FALSE)
		return ERRLH;

	return Ret;	
}

void	FreePrior(PRIOR* P)
{
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
}

PRIOR*			CreatPrior(PRIOR* P, int RateNo)
{
	PRIOR* Ret;
	int		NoOfDistVals;
	
	Ret = (PRIOR*)SMalloc(sizeof(PRIOR));

	Ret->RateNo		= RateNo;

	Ret->Dist		= P->Dist;
	Ret->NoOfCats	= -1;

	NoOfDistVals	= DISTPRAMS[P->Dist];

	Ret->DistVals = (double*)SMalloc(sizeof(double) * NoOfDistVals);

	memcpy(Ret->DistVals, P->DistVals, sizeof(double) * NoOfDistVals);

	Ret->UseHP		= P->UseHP;

	if(Ret->UseHP == TRUE)
	{
		Ret->HP = (double*)SMalloc(sizeof(double) * (NoOfDistVals * 2));
		memcpy(Ret->HP, P->HP, sizeof(double) * (NoOfDistVals * 2));
	}
	else
		Ret->HP = NULL;

	Ret->Name = NULL;

	return Ret;
}

PRIOR*	AllocBlankPrior(int NoP)
{
	PRIOR* Ret;

	Ret = (PRIOR*)SMalloc(sizeof(PRIOR));

	Ret->DistVals = (double*)SMalloc(sizeof(double) * NoP);
	
	Ret->RateNo		= -1;
	Ret->HP			= NULL;
	Ret->UseHP		= FALSE;
	Ret->OffSet		= 0;
	Ret->Name		= NULL;

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

PRIOR*		CreateGammaPrior(char *Name, double Mean, double Var)
{
	PRIOR* Ret;

	Ret = AllocBlankPrior(2);
	
	Ret->Dist			= GAMMA;
	Ret->DistVals[0]	= Mean;
	Ret->DistVals[1]	= Var;

	Ret->Name = StrMake(Name);
	
	return Ret;
}

PRIOR*		CreateUniformPrior(char *Name, double Min, double Max)
{
	PRIOR* Ret;

	Ret = AllocBlankPrior(2);
	
	Ret->Dist			= UNIFORM;
	Ret->DistVals[0]	= Min;
	Ret->DistVals[1]	= Max;

	Ret->Name = StrMake(Name);
	
	return Ret;
}

PRIOR*		CreateChiPrior(char *Name, double Mean, double Var)
{
	PRIOR* Ret;

	Ret = AllocBlankPrior(2);
	
	Ret->Dist			= CHI;
	Ret->DistVals[0]	= Mean;
	Ret->DistVals[1]	= Var;

	Ret->Name = StrMake(Name);
	
	return Ret;
}

PRIOR*		CreateExpPrior(char *Name, double Mean)
{
	PRIOR* Ret;

	Ret = AllocBlankPrior(1);
	
	Ret->Dist			= EXP;
	Ret->DistVals[0]	= Mean;

	Ret->Name = StrMake(Name);
		
	return Ret;
}

PRIOR*		CreateInvGammaPrior(char *Name, double Alpha, double Beta)
{
	PRIOR* Ret;

	Ret = AllocBlankPrior(2);
	
	Ret->Dist			= INVGAMMA;
	Ret->DistVals[0]	= Alpha;
	Ret->DistVals[1]	= Beta;

	Ret->Name = StrMake(Name);
			
	return Ret;
}

PRIOR*		CreateSGammaPrior(char *Name, double Alpha, double Beta)
{
	int		NoP;
	PRIOR	*Ret;

	NoP = DISTPRAMS[SGAMMA];
	Ret = AllocBlankPrior(NoP);

	Ret->Dist = SGAMMA;
	Ret->DistVals[0] = Alpha;
	Ret->DistVals[1] = Beta;
	
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

double	PBetaWidth(double X, double Alpha, double Beta, double Width)
{
	double	P1, P2;
	double	X1, X2;
	double	Scale;
//	int		Cat;

	Scale = 100;
	
	X1 = (X+Width) / Scale;
	X2 = X / Scale;
		
	if(X1 > 1.0 || X2 > 1.0)
		return 0.0;

	P1 = CDFBeta(X1, Alpha, Beta);
	P2 = CDFBeta(X2, Alpha, Beta);

	return P1 - P2;
}

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

double	ChiToZ(double Rate, double Mue, double Sig)
{
	double Ret;

	Ret = (Rate - Mue) / Sig;
	Ret = Ret * Ret;

	return Ret;
}

double	FindChiP(double Z, double DF, double Width)
{
	double	P1, P2;

	P1 = chdtr(DF, Z+Width);
	P2 = chdtr(DF, Z);
	

	return P1 - P2;
}

double	PChiWidth(double Rate, double Mean, double Sig, double Width)
{
	double Z, Ret;

	Z = ChiToZ(Rate, Mean, Sig);
	Ret = FindChiP(Z, 1, Width);

	return Ret;
}

/*
double	CalcChiPriors(RATES* NRates, RATES* CRates, OPTIONS* Opt)
{
	PRIORS	*Prior;
	int		PIndex=0;
	double	NSum;
	double	CSum;
	double	Ret=0;
	int		NoOfChi;

	NoOfChi = 0;
	CSum = 0;
	NSum = 0;
	for(PIndex=0;PIndex<CRates->NoOfRates;PIndex++)
	{
		Prior = CRates->Prios[PIndex];

		if(Prior->Dist == CHI)
		{
			CSum += ChiToZ(CRates->FullRates[Prior->RateNo], CRates->Prios[PIndex]->DistVals[0], CRates->Prios[PIndex]->DistVals[1]);
			NSum += ChiToZ(NRates->FullRates[Prior->RateNo], NRates->Prios[PIndex]->DistVals[0], NRates->Prios[PIndex]->DistVals[1]);

			NoOfChi++;
		}

		if(NRates->FullRates[Prior->RateNo] > 100)
			return -99999;
	}

	if(NoOfChi == 0)
		return 0;

	CRates->LhPrior = log(FindChiP(CSum, NoOfChi, 1.0 / Opt->PriorCats));

	if(NSum > 1000)
		return -9999;

	NRates->LhPrior = log(FindChiP(NSum, NoOfChi, 1.0 / Opt->PriorCats));
	

    return NRates->LhPrior - CRates->LhPrior;
}
*/

double PUni(double X, double Min, double Max)
{
	double Ret;

	if(X < Min || X > Max)
		return 0.0;
	
	Ret = 1.0 / (Max - Min);

	return Ret;	
}

double	PExpWidth(double X, double Alpha, double Width)
{
	double P1, P2;

	P1 = CDFExp(X+Width, Alpha);
	P2 = CDFExp(X, Alpha);

	return P1 - P2;
}

double	PGamaWidth(double X, double Alpha, double Beta, double Width)
{
	double	P1,P2;
	
	P1 = CDFGamma(X+Width, Alpha, Beta);
	P2 = CDFGamma(X, Alpha, Beta);

	return P1 - P2;
}

double PInvGammaWidth(double X, double Alpha, double Beta, double Width)
{
	double P1, P2;

	P1 = CDFInvGamma(X+Width, Alpha, Beta);
	P2 = CDFInvGamma(X, Alpha, Beta);

	return P1 - P2;
}

void InverseGammaTest(void)
{
	double X, P;
	double Alpha, Beta;

	Alpha = 5;
	Beta = 1;
	
//	P = InvGammaCDF(0.5, 3.88, 8.5);
	
	for(X=0.001;X<100;X+=0.1)
	{
//		P = PInvGamma(X, 2, 0.130435);
	//	P = PGamaWidth(X, 1.0, 2.0, 0.01);
	//	P = PExpWidth(X, 1.0, 0.01);
		
		P = PBetaWidth(X, Alpha, Beta, 1.0 / 100);

		P = log(P);
		printf("%f\t%f\n", X, P);
	}

	exit(0);
}

void	PriorLhTest(PRIOR *Prior, double Width)
{
	double *List;
	int Index, No;
	double P;

	No = 0;
	List = LoadDouble("./Seq/Bug/LHTest.txt", &No);

	for(Index=0;Index<No;Index++)
	{
	//	P = PDFInvGamma(List[Index], Prior->DistVals[0], Prior->DistVals[1]);
		P = PInvGammaWidth(List[Index], Prior->DistVals[0], Prior->DistVals[1], Width);
		printf("%f\t%f\n", log(P), List[Index]);
	}
}

double TestPT(void)
{
	double S, Lh;


	for(S=0;S<500;S+=0.01)
	{
		Lh = PDFInvGamma(S, 2.00, 0.13);

		printf("%f\t%f\n", S, log(Lh));
	}
//	PDFInvGamma(Val, Prior->DistVals[0], Prior->DistVals[1]);
	exit(0);
}

double	CalcGPriorGSL(double X, PRIOR *P, double Width)
{
	double A, B;

	A = gsl_cdf_gamma_P(X, P->DistVals[0], P->DistVals[1]);
	B = gsl_cdf_gamma_P(X+Width, P->DistVals[0], P->DistVals[1]);

	A = B - A;
	A = log(A);
	printf("A:\t%f\n", A);

	return A;
}

void GammaPTest(void)
{
	double GLh, Lh, X, Width; 
	PRIOR *Prior;

	Prior = CreateGammaPrior("Hello", 1.0, 1.0);


	Width = 1.0 / 100;

	X = 0.5;
	Lh = PGamaWidth(X, Prior->DistVals[0], Prior->DistVals[1], Width);

	CalcGPriorGSL(X, Prior, Width);
	Lh = log(Lh);
	printf("Lh:\t%f\n", Lh);
	exit(0);

}

double	CaclPriorCost(double X, PRIOR *Prior, int NoCats)
{
	double Ret, Width, Alpha;

	Width = 1.0 / NoCats;

	GammaPTest();

	switch(Prior->Dist)
	{
		case GAMMA:
			Ret = PGamaWidth(X, Prior->DistVals[0], Prior->DistVals[1], Width);
		break;

		case UNIFORM:
			Ret = PUni(X, Prior->DistVals[0], Prior->DistVals[1]);
		break;

		case EXP:
			Alpha = 1.0 / Prior->DistVals[0];
			Ret = PExpWidth(X, Alpha, Width);
		break;

		case CHI:
			printf("Chi is not currenlty implmented. \n");
			exit(0);
			Ret = PChiWidth(X, Prior->DistVals[0], Prior->DistVals[1], Width);
		break;

		case INVGAMMA:
			Ret = PDFInvGamma(X, Prior->DistVals[0], Prior->DistVals[1]);
		break;
	}

	Ret = log(Ret);

	return Ret;
}


double	CalcTreeTransPrior(RATES *Rates, OPTIONS *Opt)
{
	double PLh, Ret;
	PRIOR	*Prior;

	Ret = 0;

	if(Opt->EstKappa == TRUE)
	{
		Prior = GetPriorFromName("Kappa", Rates->Priors, Rates->NoPriors);
		PLh = CaclPriorCost(Rates->Kappa, Prior, Opt->PriorCats);
		if(PLh == ERRLH)
			return ERRLH;
		Ret += PLh;
	}

	if(Opt->EstLambda == TRUE)
	{
		Prior = GetPriorFromName("Lambda", Rates->Priors, Rates->NoPriors);
		PLh = CaclPriorCost(Rates->Lambda, Prior, Opt->PriorCats);
		if(PLh == ERRLH)
			return ERRLH;
		Ret += PLh;
	}

	if(Opt->EstDelta == TRUE)
	{
		Prior = GetPriorFromName("Delta", Rates->Priors, Rates->NoPriors);
		PLh = CaclPriorCost(Rates->Delta, Prior, Opt->PriorCats);
		if(PLh == ERRLH)
			return ERRLH;
		Ret += PLh;
	}

	if(Opt->EstOU == TRUE)
	{
		Prior = GetPriorFromName("OU", Rates->Priors, Rates->NoPriors);
		PLh = CaclPriorCost(Rates->OU, Prior, Opt->PriorCats);
		if(PLh == ERRLH)
			return ERRLH;
		Ret += PLh;
	}

	if(Opt->EstGamma == TRUE)
	{
		Prior = GetPriorFromName("Gamma", Rates->Priors, Rates->NoPriors);
		PLh = CaclPriorCost(Rates->Gamma, Prior, Opt->PriorCats);
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
		
		P = PDFNorm(DC->Beta[0],  0.0, 1.0);

		if(DC->Type == RJDUMMY_INTER_SLOPE)
			P *= PDFNorm(DC->Beta[1],  0.0, 1.0);
				
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
	
		PLh = CaclPriorCost(R, Prior, Opt->PriorCats);

		if(PLh == ERRLH || IsNum(PLh) == FALSE)
			return ERRLH;

		Ret += PLh;
	}

	return Ret;
}

void	CalcPriors(RATES* Rates, OPTIONS* Opt)
{
	double	PLh, Ret;
		
	Ret = 0;

	Rates->LhPrior = ERRLH;

	if(Opt->LoadModels == TRUE)
	{
		Rates->LhPrior = 0;
		return;
	}	

	PLh = CalcRatePrior(Rates, Opt);
	if(PLh == ERRLH)
		return;

	Ret += PLh;
	
	PLh = CalcTreeTransPrior(Rates, Opt);
	if(PLh == ERRLH)
		return;

	Ret += PLh;

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
	} while((Ret > Max) || (Ret < Min));

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
		for(RIndex=0;RIndex<DISTPRAMS[Prior->Dist];RIndex++)
		{
			Min = Prior->HP[RIndex*2];
			Max = Prior->HP[(RIndex*2)+1];
			Prior->DistVals[RIndex] = ChangePriorNorm(Rates, Prior->DistVals[RIndex], Dev, Min, Max);
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


void	CopyMutPriors(RANDSTATES *RandStates, PRIOR **APriosList, PRIOR **BPriosList, int NoOfPriors, double Dev)
{
	int		PIndex;
	PRIOR *APrios=NULL;
	PRIOR *BPrios=NULL;

	int		RIndex;

	for(PIndex=0;PIndex<NoOfPriors;PIndex++)
	{
		APrios = APriosList[PIndex];
		BPrios = BPriosList[PIndex];

		APrios->Dist	= BPrios->Dist;
		APrios->NoOfCats= BPrios->NoOfCats;
		APrios->RateNo	= BPrios->RateNo;

		for(RIndex=0;RIndex<DISTPRAMS[APrios->Dist];RIndex++)
			APrios->DistVals[RIndex] = ChangePrior(RandStates, BPrios->DistVals[RIndex], Dev);
	}
}

void	CopyPrior(PRIOR *A, PRIOR *B)
{
	A->Dist		= B->Dist;
	A->NoOfCats	= B->NoOfCats;
	A->RateNo	= B->RateNo;
	A->UseHP	= B->UseHP;

	memcpy(A->DistVals, B->DistVals, sizeof(double) * DISTPRAMS[A->Dist]);

	if(B->UseHP == TRUE)
		memcpy(A->HP, B->HP, sizeof(double) * DISTPRAMS[A->Dist] * 2);
}

PRIORDIST	StrToPriorDist(char* Str)
{
	int			Index;

	MakeLower(Str);

	for(Index=0;Index<NO_PRIOR_DIST;Index++)
	{
		if(strcmp(Str, DISTNAMES[Index])==0)
			return (PRIORDIST)(Index);
	}

	return (PRIORDIST)-1;
}

int			CheckPriorDistVals(int Tokes, char **Passed)
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
		if(P < 0)
		{
			printf("Prior parameters values must be greater than 0, value %f is invalid.\n", P);
			exit(1);
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

	Ret = NULL;

	if(Tokes != 2 && Tokes != 3)
	{
		printf("Prior requires a distribution name, (beta, gamma, uniform, chi, exp, invgamma) and distribution parameters.\n");
		exit(1);
	}

	if(StrToPriorDist(Passed[0]) == -1)
	{
		printf("Invalid prior distribution name,. valid names are beta, gamma, uniform, chi, exp, invgamma.\n");
		exit(1);
	}

	PD = StrToPriorDist(Passed[0]);

	if(Tokes -1 != DISTPRAMS[PD])
	{
		printf("Prior %s (%s) requires %d parameters.\n", Name, DISTNAMES[PD], DISTPRAMS[PD]);
		exit(0);
	}

	if(CheckPriorDistVals(Tokes-1, &Passed[1]) == FALSE)
		exit(0);

	PVal = MakePriorParam(Tokes-1, &Passed[1]);


	if(PD == GAMMA)
		Ret = CreateGammaPrior(Name, PVal[0], PVal[1]);

	if(PD == UNIFORM)
		Ret = CreateUniformPrior(Name, PVal[0], PVal[1]);

	if(PD == CHI)
		Ret = CreateChiPrior(Name, PVal[0], PVal[1]);

	if(PD == EXP)
		Ret = CreateExpPrior(Name, PVal[0]);

	if(PD == INVGAMMA)
		Ret = CreateInvGammaPrior(Name, PVal[0], PVal[1]);
	
	free(PVal);

	return Ret;
}


PRIOR*		CreateHyperPriorFromStr(char *Name, int Tokes, char **Passed)
{
	PRIOR		*Ret;
	PRIORDIST	PDist;
	int			NoParam;
	double		*PVal;

	if(Tokes < 3)
	{
		printf("A hyper prior require a distribution and min max values for each parameters");
		exit(0);
	}

	PDist = StrToPriorDist(Passed[0]);

	if(PDist == -1)
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
	
	if(CheckPriorDistVals(Tokes-1, &Passed[1]) == FALSE)
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
	PRIOR **NPList;
	int Index, Pos;

	NPList = (PRIOR**)SMalloc(sizeof(PRIOR*) * Opt->NoAllPriors);

	Pos = 0;
	for(Index=0;Index<Opt->NoAllPriors;Index++)
	{
		if(strcmp(Name, Opt->AllPriors[Index]->Name) != 0)
			NPList[Pos++] = Opt->AllPriors[Index];
	}

	free(Opt->AllPriors);
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