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

/*
extern double beta(double a, double b);
extern double incbet(double aa, double bb, double xx );

extern double igam(double a, double x);
extern double igamc ( double, double );
extern double igami ( double, double );
extern double gamma(double a);
*/

extern double chdtr(double df, double x);

void	FreePrior(PRIORS* P)
{
	free(P->DistVals);

	if(P->HP != NULL)
		free(P->HP);

	if(P->RateName != NULL)
		free(P->RateName);

	free(P);
}

void	FreePriors(RATES *Rates)
{
	int	PIndex;
	
	if(Rates->Prios == NULL)
		return;
	
	for(PIndex=0;PIndex<Rates->NoOfPriors;PIndex++)
		FreePrior(Rates->Prios[PIndex]);

	free(Rates->Prios);

	if(Rates->PriorGamma != NULL)
		FreePrior(Rates->PriorGamma);

	if(Rates->PriorOU != NULL)
		FreePrior(Rates->PriorOU);

	if(Rates->PriorDelta != NULL)
		FreePrior(Rates->PriorDelta);

	if(Rates->PriorKappa != NULL)
		FreePrior(Rates->PriorKappa);

	if(Rates->PriorLambda != NULL)
		FreePrior(Rates->PriorLambda);
}

PRIORS*			CreatPrior(PRIORS* P, int RateNo)
{
	PRIORS* Ret;
	int		NoOfDistVals;
	
	Ret = (PRIORS*) malloc(sizeof(PRIORS));
	if(Ret == NULL)
		MallocErr();

	Ret->RateNo		= RateNo;

	Ret->Dist		= P->Dist;
	Ret->NoOfCats	= -1;

	NoOfDistVals	= DISTPRAMS[P->Dist];

	Ret->DistVals = (double*)malloc(sizeof(double) * NoOfDistVals);
	if(Ret==NULL)
		MallocErr();

	memcpy(Ret->DistVals, P->DistVals, sizeof(double) * NoOfDistVals);

	Ret->UseHP		= P->UseHP;

	if(Ret->UseHP == TRUE)
	{
		Ret->HP = (double*)malloc(sizeof(double) * (NoOfDistVals * 2));
		if(Ret->HP == NULL)
			MallocErr();
		memcpy(Ret->HP, P->HP, sizeof(double) * (NoOfDistVals * 2));
	}
	else
		Ret->HP = NULL;

	Ret->RateName = NULL;

	return Ret;
}

PRIORS*	AllocBlankPrior(int NoP)
{
	PRIORS* Ret;

	Ret = (PRIORS*)malloc(sizeof(PRIORS));
	if(Ret == NULL)
		MallocErr();

//	Ret->Dist = UNIFORM;
	Ret->DistVals = (double*)malloc(sizeof(double) * NoP);
	if(Ret->DistVals == NULL)
		MallocErr();
	
	Ret->RateNo		= -1;
	Ret->HP			= NULL;
	Ret->UseHP		= FALSE;
	Ret->OffSet		= 0;
	Ret->RateName	= NULL;

	return Ret;
}

PRIORS*		CreateBetaPrior(double Mean, double Var)
{
	PRIORS* Ret;

	Ret = AllocBlankPrior(2);
	
	Ret->Dist			= BETA;
	Ret->DistVals[0]	= Mean;
	Ret->DistVals[1]	= Var;
	
	return Ret;
}

PRIORS*		CreateGammaPrior(double Mean, double Var)
{
	PRIORS* Ret;

	Ret = AllocBlankPrior(2);
	
	Ret->Dist			= GAMMA;
	Ret->DistVals[0]	= Mean;
	Ret->DistVals[1]	= Var;
	
	return Ret;
}

PRIORS*		CreateUniformPrior(double Min, double Max)
{
	PRIORS* Ret;

	Ret = AllocBlankPrior(2);
	
	Ret->Dist			= UNIFORM;
	Ret->DistVals[0]	= Min;
	Ret->DistVals[1]	= Max;
	
	return Ret;
}

PRIORS*		CreateChiPrior(double Mean, double Var)
{
	PRIORS* Ret;

	Ret = AllocBlankPrior(2);
	
	Ret->Dist			= CHI;
	Ret->DistVals[0]	= Mean;
	Ret->DistVals[1]	= Var;
	
	return Ret;
}

PRIORS*		CreateExpPrior(double Mean)
{
	PRIORS* Ret;

	Ret = AllocBlankPrior(1);
	
	Ret->Dist			= EXP;
	Ret->DistVals[0]	= Mean;
		
	return Ret;
}

PRIORS*		CreateInvGammaPrior(double Alpha, double Beta)
{
	PRIORS* Ret;

	Ret = AllocBlankPrior(2);
	
	Ret->Dist			= INVGAMMA;
	Ret->DistVals[0]	= Alpha;
	Ret->DistVals[1]	= Beta;
			
	return Ret;
}

void		CreatTreeTransformPriors(OPTIONS *Opt, RATES *Rates)
{
	if(Opt->UseGamma == TRUE)
		Rates->PriorGamma = CreatPrior(Opt->PriorGamma, -1);
	else
		Rates->PriorGamma = NULL;

	if(Opt->EstKappa == TRUE)
		Rates->PriorKappa = CreatPrior(Opt->PriorKappa, -1);
	else
		Rates->PriorKappa = NULL;


	if(Opt->EstDelta == TRUE)
		Rates->PriorDelta = CreatPrior(Opt->PriorDelta, -1);
	else
		Rates->PriorDelta = NULL;


	if(Opt->EstLambda == TRUE)
		Rates->PriorLambda = CreatPrior(Opt->PriorLambda, -1);
	else
		Rates->PriorLambda = NULL;

	if(Opt->EstOU == TRUE)
		Rates->PriorOU = CreatPrior(Opt->PriorOU, -1);
	else
		Rates->PriorOU = NULL;
}

void		CreatPriors(OPTIONS *Opt, RATES* Rates)
{
	PRIORS**	Ret=NULL;
	int			NoOfPriors;
	int			Index;
	int			RIndex;

	CreatTreeTransformPriors(Opt, Rates);
	
	if(Opt->UseRJMCMC == TRUE)
	{
		Ret = (PRIORS**)malloc(sizeof(PRIORS*));
		if(Ret==NULL)
			MallocErr();

		Ret[0] = CreatPrior(Opt->RJPrior, -1);

		Rates->NoOfPriors = 1;
		Rates->Prios = Ret;
		return;
	}

	NoOfPriors = 0;
	
	for(Index=0;Index<Opt->NoOfRates;Index++)
	{
		if(Opt->ResTypes[Index] == RESNONE)
			NoOfPriors++;
	}
	fflush(stdout);
	
	Rates->NoOfPriors = NoOfPriors;
	
	if(NoOfPriors==0)
	{
		Rates->Prios = NULL;
	}

	Ret = (PRIORS**)malloc(sizeof(PRIORS*)*NoOfPriors);
	if(Ret==NULL)
		MallocErr();

	Index=0;
	for(RIndex=0;RIndex<Opt->NoOfRates;RIndex++)
	{
		if(Opt->ResTypes[RIndex] == RESNONE)
		{
			Ret[Index] = CreatPrior(Opt->Priors[RIndex], Index) ;
			Ret[Index]->NoOfCats = Opt->PriorCats;

			Index++;
		}
	}


	Rates->Prios = Ret;
}

void	SetRatesToPriors(OPTIONS *Opt, RATES* Rates)
{
	PRIORS	*Prior;
	int		PIndex=0;
	int		NoOfPriors;

	NoOfPriors = Rates->NoOfRates;
	if(Opt->UseRJMCMC == TRUE)
		NoOfPriors = 1;

	for(PIndex=0;PIndex<NoOfPriors;PIndex++)
	{
		Prior = Rates->Prios[PIndex];
		
		if(Prior->Dist == UNIFORM)
		{
			if(Opt->LoadModels == FALSE)
			{
				if(Opt->ModelType != MT_CONTRAST)
					Rates->Rates[PIndex] = Prior->DistVals[0]+((Prior->DistVals[1]-Prior->DistVals[0])/2.0);
			}
		}
	}	
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

void	PriorLhTest(PRIORS *Prior, double Width)
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

double	CaclPriorCost(double Val, PRIORS *Prior, int NoCats)
{
	double Ret, Width, Alpha;

	Width = 1.0 / NoCats;

//	InverseGammaTest();
//	FatTailTest(Prior, Width);

	switch(Prior->Dist)
	{
		case BETA:
			Ret = PBetaWidth(Val, Prior->DistVals[0], Prior->DistVals[1], Width);
		break;

		case GAMMA:
			Ret = PGamaWidth(Val, Prior->DistVals[0], Prior->DistVals[1], Width);
		break;

		case UNIFORM:
			Ret = PUni(Val, Prior->DistVals[0], Prior->DistVals[1]);
		break;

		case EXP:
			Alpha = 1.0 / Prior->DistVals[0];
			Ret = PExpWidth(Val, Alpha, Width);
		break;

		case CHI:
			printf("Chi is not currenlty implmented. \n");
			exit(0);
			Ret = PChiWidth(Val, Prior->DistVals[0], Prior->DistVals[1], Width);
		break;

		case INVGAMMA:
			Ret = PDFInvGamma(Val, Prior->DistVals[0], Prior->DistVals[1]);
//			InvGammaTest(Prior, Width);
//			Ret = PInvGammaWidth(Val, Prior->DistVals[0], Prior->DistVals[1], Width);
		break;
	}

	Ret = log(Ret);

	return Ret;
}


double	CalcTreeTransPrior(RATES *Rates, OPTIONS *Opt)
{
	double Ret;

	Ret = 0;

	if(Opt->EstKappa == TRUE)
		Ret += CaclPriorCost(Rates->Kappa, Rates->PriorKappa, Opt->PriorCats);

	if(Opt->EstLambda == TRUE)
		Ret += CaclPriorCost(Rates->Lambda, Rates->PriorLambda, Opt->PriorCats);

	if(Opt->EstDelta == TRUE)
		Ret += CaclPriorCost(Rates->Delta, Rates->PriorDelta, Opt->PriorCats);

	if(Opt->EstOU == TRUE)
		Ret += CaclPriorCost(Rates->OU, Rates->PriorOU, Opt->PriorCats);

	if(Opt->EstGamma == TRUE)
		Ret += CaclPriorCost(Rates->Gamma, Rates->PriorGamma, Opt->PriorCats);
	
	return Ret;
}



double CalcRJDummyPriors(OPTIONS *Opt, RATES* Rates)
{
	RJDUMMY *RJDummy;
	DUMMYCODE *DC;
	int		Index;
	double	Ret, P;

//	Testxx2();
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

void	CalcPriors(RATES* Rates, OPTIONS* Opt)
{
	PRIORS	*Prior;
	int		PIndex, Err;
	double	TLh;
	double	Rate;
	int		NoPRates;
	double	CalcP;

	CalcP = 0;
	Err = FALSE;

	if(Opt->LoadModels == TRUE)
		return;

	if(Opt->UseRJMCMC == TRUE)
		NoPRates = Rates->NoOfRJRates;
	else
		NoPRates = Rates->NoOfRates;

	for(PIndex=0;PIndex<NoPRates;PIndex++)
	{
		if(Opt->UseRJMCMC == FALSE)
		{
			Prior = Rates->Prios[PIndex];
			Rate = Rates->Rates[Prior->RateNo];
		}
		else
		{
			Prior = Rates->Prios[0];
			Rate = Rates->Rates[PIndex];
		}

		TLh = CaclPriorCost(Rate, Prior, Opt->PriorCats);

		if(TLh == ERRLH || IsNum(TLh) == FALSE)
		{
			Rates->LhPrior = ERRLH;
			return;
		}

		CalcP += TLh;
	}

	CalcP += CalcTreeTransPrior(Rates, Opt);

	if(Opt->RJDummy == TRUE)
		CalcP += CalcRJDummyPriors(Opt, Rates);

	if(UseNonParametricMethods(Opt) == TRUE)
		CalcP += CalcVarRatesPriors(Rates, Opt, &Err);
	
	if((CalcP == ERRLH) || (IsNum(CalcP) == FALSE) || (Err == TRUE))
	{
		Rates->LhPrior = ERRLH;
		return;
	}

	Rates->LhPrior = CalcP;
}

void	MutatePriors(RATES *Rates, PRIORS **PriosList, int NoOfPriors)
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

void	MutatePriorsNormal(RATES *Rates, PRIORS **PriosList, int NoOfPriors, double Dev)
{
	int		PIndex;
	int		RIndex;
	PRIORS	*Prior;
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


void	CopyMutPriors(RANDSTATES *RandStates, PRIORS **APriosList, PRIORS **BPriosList, int NoOfPriors, double Dev)
{
	int		PIndex;
	PRIORS *APrios=NULL;
	PRIORS *BPrios=NULL;

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

void	CopyPrior(PRIORS *A, PRIORS *B)
{
	A->Dist		= B->Dist;
	A->NoOfCats	= B->NoOfCats;
	A->RateNo	= B->RateNo;
	A->UseHP	= B->UseHP;

	memcpy(A->DistVals, B->DistVals, sizeof(double) * DISTPRAMS[A->Dist]);

	if(B->UseHP == TRUE)
		memcpy(A->HP, B->HP, sizeof(double) * DISTPRAMS[A->Dist] * 2);
}

void	CopyRatePriors(PRIORS**APriosList, PRIORS **BPriosList, int NoOfPriors)
{
	int		PIndex;

	for(PIndex=0;PIndex<NoOfPriors;PIndex++)
		CopyPrior(APriosList[PIndex], BPriosList[PIndex]);
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
			return FALSE;
		}

		P = atof(Passed[Index]);
		if(P < 0)
		{
			printf("Prior parameters values must be greater than 0, value %f is invalid.\n", P);
			return FALSE;
		}
	}

	return TRUE;
}

double*		MakePriorParam(int Tokes, char **Passed)
{
	double *Ret;
	int Index;

	Ret = (double*)malloc(sizeof(double) * Tokes);
	if(Ret == NULL)
		MallocErr();

	for(Index=0;Index<Tokes;Index++)
		Ret[Index] = atof(Passed[Index]);

	return Ret;
}


PRIORS*		CreatePrior(int Tokes, char **Passed)
{
	PRIORDIST	PD;
	double		*PVal;
	PRIORS		*Ret;

	Ret = NULL;
	
	if(Tokes != 2 && Tokes != 3)
	{
		printf("Prior requires a distribution name, (beta, gamma, uniform, chi, exp, invgamma) and distribution parameters.\n");
		return NULL;
	}

	if(StrToPriorDist(Passed[0]) == -1)
	{
		printf("Invalid prior distribution name,. valid names are beta, gamma, uniform, chi, exp, invgamma.\n");
		return NULL;
	}

	PD = StrToPriorDist(Passed[0]);

	if(Tokes -1 != DISTPRAMS[PD])
	{
		printf("Prior %s requires %d parameters.\n", DISTNAMES[PD], DISTPRAMS[PD]);
		return NULL;
	}

	if(CheckPriorDistVals(Tokes-1, &Passed[1]) == FALSE)
		return NULL;

	PVal = MakePriorParam(Tokes-1, &Passed[1]);

	if(PD == BETA)
		Ret = CreateBetaPrior(PVal[0], PVal[1]);

	if(PD == GAMMA)
		Ret = CreateGammaPrior(PVal[0], PVal[1]);

	if(PD == UNIFORM)
		Ret = CreateUniformPrior(PVal[0], PVal[1]);

	if(PD == CHI)
		Ret = CreateChiPrior(PVal[0], PVal[1]);

	if(PD == EXP)
		Ret = CreateExpPrior(PVal[0]);

	if(PD == INVGAMMA)
		Ret = CreateInvGammaPrior(PVal[0], PVal[1]);

	free(PVal);

	return Ret;
}

