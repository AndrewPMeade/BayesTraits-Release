#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "TypeDef.h"
#include "Trees.h"
#include "Likelihood.h"
#include "GenLib.h"

double MLLambda(TREE *Tree)
{
	double SumBL;
	int Index;

	SumBL = 0;
	for(Index=1;Index<Tree->NoNodes;Index++)
		SumBL += Tree->NodeList[Index]->Length;
	
	return SumBL / (Tree->NoNodes - 1);
}

double	NormPLhVect(double *Vect, int NOS)
{
	int Index;
	double Ret;

	Ret = 0;
	for(Index=0;Index<NOS;Index++)
		Ret += Vect[Index];

	return Ret;
}

double	CaclNodeLhRate(NODE N, int NOS, double *Lambda)
{
	int Index;
	double *PartLh;

	double Ret, CLabda, Lh, NC;

	
//	PartLh = N->Ans->Partial[0];
	PartLh = N->Partial[0];

	NC = NormPLhVect(PartLh, NOS);

	Ret = 0;

	for(Index=0;Index<NOS;Index++)
	{
		CLabda = Lambda[Index];
		if(N->Tip == TRUE)
			Lh = exp(-(N->Length/CLabda));
		else
			Lh = (1.0/CLabda) * exp(-(N->Length/CLabda));

		Lh = Lh * (PartLh[Index] / NC);

		Ret += Lh;
	}

	return log(Ret);
}

double	LhStateSpeciationRate(OPTIONS *Opt, TREES *Trees, RATES *Rates, double *Lambda)
{
	int Index, NOS;
	double Ret;
	TREE *Tree;
	NODE N;

	NOS = Trees->NoStates;
	Ret = 0;

	Tree = Trees->Tree[Rates->TreeNo];

	for(Index=0;Index<Tree->NoNodes;Index++)
	{
		N = Tree->NodeList[Index];

		if(N != Tree->Root)
			Ret += CaclNodeLhRate(N, NOS, Lambda);
	}

	if(IsNum(Ret) == FALSE)
		return ERRLH;

	return Ret;
}


void	StateSpeciationRateTemp(OPTIONS *Opt, TREES *Trees, RATES *Rates)
{
	TREE *Tree;
	int Index, NOS, SIndex;
	NODE N;
	double NC;

	Tree = Trees->Tree[0];

	NOS = Trees->NoStates;

	for(Index=0;Index<Tree->NoNodes;Index++)
	{
		N = Tree->NodeList[Index];

		if(N->Ans != NULL)
		{
			NC = NormPLhVect(N->Partial[0], NOS);

			printf("%d\t%d\t%f\t", Index, N->Tip, N->Length);

			for(SIndex=0;SIndex<NOS;SIndex++)
			{
				printf("%f\t", N->Partial[0][SIndex]/NC);			
			}
			printf("\n");
		}
	}

	exit(0);
}

void	SetTreeRandBL(TREE *Tree)
{
	int Index;
	RANDSTATES *RS;

	RS = CreateSeededRandStates(1977);

	for(Index=0;Index<Tree->NoNodes;Index++)
	{
		Tree->NodeList[Index]->Length = RandExp(RS, 0.1);
	//	printf("%f\n", Tree->NodeList[Index]->Length);
	}
			
//	exit(0);
	FreeRandStates(RS);
}

double	CaclStateSpeciationRateLh(OPTIONS *Opt, TREES *Trees, RATES *Rates)
{
	double *Lambda;
	double Lh, L;
	int Index;
	TREE *Tree;

	Likelihood(Rates, Trees, Opt);
	
//	StateSpeciationRateTemp(Opt, Trees, Rates);
	
	Lambda = (double*)SMalloc(sizeof(double) * Trees->NoStates);
	
	Tree = Trees->Tree[Rates->TreeNo];

	SetTreeRandBL(Tree);
	
	Index=0;
	for(L=0.01;L<1.0;L+=0.001)
//	for(Index=0;Index<10000;Index++)
	{
//		Lambda[0] = RandDouble(Rates->RS)*2.5;
//		Lambda[1] = RandDouble(Rates->RS)*2.5;
	
		Lambda[0] = L;
		Lambda[1] = L;

		Lh = LhStateSpeciationRate(Opt, Trees, Rates, Lambda);

		printf("%d\t%f\t%f\t%f\n", Index, Lh, Lambda[0], Lambda[1]);
	}

	exit(0);

	return 0;
}