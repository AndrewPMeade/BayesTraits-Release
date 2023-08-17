#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

#include "GlobalTrend.h"
#include "GenLib.h"
#include "Likelihood.h"
#include "Contrasts.h"
#include "RandLib.h"
#include "Priors.h"

void	SetTipGlobalTrend(NODE Node, double GlobalBeta, double SumBL, int NoSites)
{
	int Index;
	CONDATA* Con;

	Con = Node->ConData;

	for(Index=0;Index<NoSites;Index++)
		Con->Contrast[0]->Data[Index] -= GlobalBeta * SumBL;
}

void	RecApplyGlobalTrend(NODE Node, double GlobalBeta, double SumBL, int NoSites)
{
	int Index;

	SumBL += Node->Length;

	if(Node->Tip == TRUE)
		SetTipGlobalTrend(Node, GlobalBeta, SumBL, NoSites);

	for(Index=0;Index<Node->NoNodes;Index++)
		RecApplyGlobalTrend(Node->NodeList[Index], GlobalBeta, SumBL, NoSites);

}

void SetGlobalTrend(OPTIONS* Opt, TREES* Trees, RATES* Rates)
{
	int Index;
	TREE *Tree;
	NODE Root;
	
	Tree = Trees->Tree[Rates->TreeNo];
	Root = Tree->Root;
	
	for(Index=0;Index<Root->NoNodes;Index++)
		RecApplyGlobalTrend(Root->NodeList[Index], Rates->GlobalTrend, 0.0, Trees->NoSites);
}

void	ChangeGlobalTrend(RATES *Rates, SCHEDULE *Shed)
{
	double Dev, Change;

 	Shed->CurrentAT = Shed->GlobalTrendAT;

	Dev = Shed->CurrentAT->CDev;
	
	Change = RandNormal(Rates->RS, 0, Dev);

	Rates->GlobalTrend += Change;
}

double	CalcGlobalTrendPrior(RATES* Rates)
{
	double Lh;
	PRIOR *Prior;
	
	Prior = GetPriorFromName("GlobalTrend", Rates->Priors, Rates->NoPriors);

	Lh = CalcLhPriorP(Rates->GlobalTrend, Prior);

	return Lh;
}