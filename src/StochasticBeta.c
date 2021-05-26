#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

#include "StochasticBeta.h"
#include "Priors.h"


int			GetNoSB(TREES *Trees)
{
	return Trees->Tree[0]->NoNodes;
}

void		InitStochasticBeta(OPTIONS *Opt, TREES *Trees, RATES *Rates)
{	
	int Index;
		
	Rates->NoSB	= GetNoSB(Trees);
	
	Rates->SB_Type_Map = (STOCHASTIC_BETA_TYPE*)SMalloc(sizeof(STOCHASTIC_BETA_TYPE) * Rates->NoSB);
	Rates->SB_Vect = (double*)SMalloc(sizeof(double) * Rates->NoSB);

	for(Index=0;Index<Rates->NoSB;Index++)
	{
		Rates->SB_Type_Map[Index] = SB_NONE;
		Rates->SB_Vect[Index] = 0.0;
	}
}

void		FreeStochasticBeta(RATES *Rates)
{
	free(Rates->SB_Type_Map);
	free(Rates->SB_Vect);
}

void		MapStochasticBeta(TREES *Trees, RATES *Rates)
{
	int Index;
	NODE	Node;
	TREE *Tree;

	Tree = Trees->Tree[0];
		
	for(Index=0;Index<Rates->NoSB;Index++)
	{
		Node = Tree->NodeList[Index];

		// Add to the Landscape beta, something may be there.
		Node->LandscapeBeta += Rates->SB_Vect[Index];
	}
}

void		CopyStochasticBeta(RATES *RatesA, RATES *RatesB)
{
	int NoNodes;

	NoNodes = RatesA->NoSB;

	memcpy(RatesA->SB_Type_Map, RatesB->SB_Type_Map, sizeof(STOCHASTIC_BETA_TYPE) * NoNodes);
	memcpy(RatesA->SB_Vect, RatesB->SB_Vect, sizeof(double) * NoNodes);
}

double	CalcStocaticPrior(double BetaT, double T, double Sig2)
{
	double Z, Ret;

	Z = BetaT / sqrt(Sig2 * T);

	Ret = CaclNormalLogLh(Z, 1.0, 1.0);

	if(ValidPriorLh(Ret) == FALSE)
		return ERRLH;

	return Ret;
}


void	PriorP(void)
{
	double Beta, P, BL;


	BL = 1.0;
	for(Beta=-10;Beta<10;Beta+=0.01)
	{
		P = CalcStocaticPrior(Beta*BL, BL, 0.011807697213);
		printf("%f\t%f\n", Beta, P);
	}

	exit(0);
}

double CaclStochasticBetaPrior(TREES *Trees, RATES *Rates)
{
	double Ret, P, Sig2;
	int Index;
	NODE Node;
	PRIOR *Prior;
	TREE *Tree;
	
//	PriorP();
	Tree = Trees->Tree[0];

	Sig2 = Rates->Rates[1];

	Prior = GetPriorFromName("StochasticBeta", Rates->Priors, Rates->NoPriors);
//	PriorP(Prior);

	Ret = 0;
	for(Index=1;Index<Rates->NoSB;Index++)
	{
		Node = Tree->NodeList[Index];

		if(Rates->SB_Type_Map[Index] == SB_NONE)
			P = CalcStocaticPrior(Rates->SB_Vect[Index], Tree->NodeList[Index]->Length, Sig2);

		if(Rates->SB_Type_Map[Index] == SB_RJ)
		{
			if(Rates->SB_Vect[Index] < 0)
				P = CalcLhPriorP(-Rates->SB_Vect[Index], Prior);
			else
				P = CalcLhPriorP(Rates->SB_Vect[Index], Prior);
		}	
		
		if(P == ERRLH)
			return ERRLH;
		
		Ret += P;
	}

	return Ret;
}

int		GetCaclStochasticBetaPoint(RATES *Rates)
{
	int Point;

	Point = RandUSInt(Rates->RS) % (Rates->NoSB - 1);

	return Point + 1;
}

void	CaclStochasticBetaRJ(RATES *Rates)
{
	int Point;
	PRIOR *Prior;

	Point = GetCaclStochasticBetaPoint(Rates);

	Prior = GetPriorFromName("StochasticBeta", Rates->Priors, Rates->NoPriors);

	if(Rates->SB_Type_Map[Point] == SB_NONE)
	{
		Rates->SB_Type_Map[Point] = SB_RJ;
//		Rates->SB_Vect[Point] = RandFromPrior(Rates->RNG, Prior);
//		if(RandDouble(Rates->RS) < 0.5)
//			Rates->SB_Vect[Point] = Rates->SB_Vect[Point] * -1;
	}
	else
	{
		Rates->SB_Type_Map[Point] = SB_NONE;
//		Rates->SB_Vect[Point] = 0.0;
	}

//	Rates->SB_Vect[Point] = 0.0;
}

void	ChangeStochasticBeta(TREES *Trees, RATES *Rates, SCHEDULE *Sched)
{
	TREE *Tree;
	int Point;
	double Dev, Change;

	Point = GetCaclStochasticBetaPoint(Rates);

	if(Rates->SB_Type_Map[Point] == SB_NONE)
		Sched->CurrentAT = Sched->StochasticBeta;
	else
		Sched->CurrentAT = Sched->StochasticBetaPrior;
		
	Dev = Sched->CurrentAT->CDev;

	Change = RandNormal(Rates->RS, Rates->SB_Vect[Point], Dev);
	Change = Change - Rates->SB_Vect[Point];

//	printf("Change\t%f\n", Change);

	Rates->SB_Vect[Point] = Change;

	Tree = Trees->Tree[0];
	if(Tree->NodeList[Point]->Tip == TRUE)
		Rates->SB_Vect[Point] = 0.0;
}

int	GetNoStochasticBetaType(RATES *Rates, STOCHASTIC_BETA_TYPE SB_Type)
{
	int Ret, Index;

	Ret = 0;

	for(Index=0;Index<Rates->NoSB;Index++)
		if(Rates->SB_Type_Map[Index] == SB_Type)
			Ret++;

	return Ret;
}

void	LogStochasticBetaResults(OPTIONS* Opt, TREES* Trees, RATES* Rates, long long It)
{
	int Index;
	FILE* Out;

	Out = Opt->VarRatesLog;

	for (Index = 0; Index < Rates->NoSB; Index++)
	{
		fprintf(Out, "0\t");
		fprintf(Out, "%f\t", Rates->SB_Vect[Index]);
		fprintf(Out, "%d\t", Index);
		if (Rates->SB_Type_Map[Index] == SB_NONE)
			fprintf(Out, "None\t");
		else
			fprintf(Out, "RJ\t");
	}

}

