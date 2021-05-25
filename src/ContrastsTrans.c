#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "typedef.h"
#include "genlib.h"
#include "trees.h"
#include "ContrastsTrans.h"

void	SetFixedConTrans(OPTIONS *Opt, TREES *Trees);
void	SetConTrans(TREE *Tree, double Kappa, double Lambda, double Delta, double OU);

void	RecGetSumBL(NODE N, double *SumBL)
{
	int Index;

	*SumBL += N->Length;

	if(N->Tip == TRUE)
		return;

	for(Index=0;Index<N->NoNodes;Index++)
		RecGetSumBL(N->NodeList[Index], SumBL);
}

double	SumNodeBL(NODE N)
{
	double Ret;

	Ret = 0;
	RecGetSumBL(N, &Ret);

	Ret = Ret - N->Length;

	return Ret;
}

void	RecScaleSubTree(NODE N, double Scale)
{
	int Index;

	N->Length = N->Length * Scale;

	if(N->Tip == TRUE)
		return;

	for(Index=0;Index<N->NoNodes;Index++)
		RecScaleSubTree(N->NodeList[Index], Scale);
}

void	ScaleSubTree(NODE N, double Scale)
{
	int Index;

	if(N->Tip == TRUE)
	{
		printf("Sub Tree can only scale internal nodes.\n");
		exit(0);
	}

	for(Index=0;Index<N->NoNodes;Index++)
		RecScaleSubTree(N->NodeList[Index], Scale);
}

void	RecTransContNodeDelta(NODE N, double Delta, double PathLen)
{
	int Index;
	double TLen;

	TLen = N->Length + PathLen;

	N->Length = pow(PathLen+N->Length, Delta) - pow(PathLen, Delta) ;

	if(N->Tip == TRUE)
		return;

	for(Index=0;Index<N->NoNodes;Index++)
		RecTransContNodeDelta(N->NodeList[Index], Delta, TLen);
}

void	TransContNodeDelta(NODE N, double Delta, int Norm)
{
	double SumBL,Scale;
	int Index;

	if(Norm == TRUE)
		SumBL = SumNodeBL(N);

	for(Index=0;Index<N->NoNodes;Index++)
		RecTransContNodeDelta(N->NodeList[Index], Delta, 0);

	if(Norm == FALSE)
		return;

	Scale = SumBL / SumNodeBL(N);
	ScaleSubTree(N, Scale);
}

void	RecTransContNodeKappa(NODE N, double Kappa, double KappaPathLen)
{
	int Index;
	double TLen;

	TLen = pow(N->Length, Kappa) + KappaPathLen;

	N->Length = TLen - KappaPathLen;

	if(N->Tip == TRUE)
		return;

	for(Index=0;Index<N->NoNodes;Index++)
		RecTransContNodeKappa(N->NodeList[Index], Kappa, TLen);
}

void	TransContNodeKappa(NODE N, double Kappa, int Norm)
{
	double SumBL,Scale;
	int Index;

	if(Norm == TRUE)
		SumBL = SumNodeBL(N);

	for(Index=0;Index<N->NoNodes;Index++)
		RecTransContNodeKappa(N->NodeList[Index], Kappa, 0);

	if(Norm == FALSE)
		return;

	Scale = SumBL / SumNodeBL(N);
	ScaleSubTree(N, Scale);
}


void	SetContTrans(TREE *Tree, double Kappa, double Lambda, double Delta, double OU)
{
	
}



double	CaclOU(double PathLen, double OU, double T)
{
	double Ret;

	Ret = exp(-2.0 * OU * (T - PathLen));
	Ret *=  1.0 - exp(-2.0 * OU * PathLen);

	Ret *= 1.0 / (2.0 * OU);

	return Ret;
}

void	RecTransContNodeOU(NODE N, double OU, double T, double PathLen)
{
	int Index;
	double TLen;

	TLen = N->Length + PathLen;

//	N->Length = pow(PathLen+N->Length, OU) - pow(PathLen, Delta) ;
	
	N->Length = CaclOU(PathLen+N->Length, OU, T) - CaclOU(PathLen, OU, T);

	if(N->Tip == TRUE)
		return;

	for(Index=0;Index<N->NoNodes;Index++)
		RecTransContNodeOU(N->NodeList[Index], OU, T, TLen);
}

void FindOUT(NODE N, double *T)
{
	int Index;

	if(N->Tip == TRUE)
	{
		if(N->DistToRoot > *T)
			*T = N->DistToRoot;
		return;
	}

	for(Index=0;Index<N->NoNodes;Index++)
		FindOUT(N->NodeList[Index], T);
}

void	TransContNodeOU(NODE N, double OU, int Norm)
{
	double SumBL,Scale, T;
	int Index;

	if(Norm == TRUE)
		SumBL = SumNodeBL(N);

	T = -1;
	FindOUT(N, &T);

	for(Index=0;Index<N->NoNodes;Index++)
		RecTransContNodeOU(N->NodeList[Index], OU, T,  0);

	if(Norm == FALSE)
		return;

	Scale = SumBL / SumNodeBL(N);
	ScaleSubTree(N, Scale);
}

void	RecTransContNodeLambda(NODE N, double Lambda,  double PathLen)
{
	double	TLen;
	int		Index;
	
	if(N->Tip == TRUE)
	{
		N->Length = N->DistToRoot - PathLen;
		return;
	}


	N->Length = N->Length * Lambda;

	TLen = N->Length + PathLen;

	for(Index=0;Index<N->NoNodes;Index++)
		RecTransContNodeLambda(N->NodeList[Index], Lambda, TLen);
}

void	TransContNodeLambda(NODE N, double Lambda, int Norm)
{
	double SumBL, Scale;
	int Index;

	if(Norm == TRUE)
		SumBL = SumNodeBL(N);
	
	for(Index=0;Index<N->NoNodes;Index++)
		RecTransContNodeLambda(N->NodeList[Index], Lambda,  N->DistToRoot);

	

	if(Norm == FALSE)
		return;

	Scale = SumBL / SumNodeBL(N);
	ScaleSubTree(N, Scale);
}

int		NeedToReSetBL(OPTIONS *Opt)
{
	if(Opt->UseVarRates == TRUE)
		return TRUE;

	if(Opt->UseKappa  == TRUE)
		return TRUE;

	if(Opt->UseOU == TRUE)
		return TRUE;

	if(Opt->UseDelta == TRUE)
		return TRUE;
	
	if(Opt->UseLambda == TRUE)
		return TRUE;
	
	return FALSE;
}

void	TransformContrastTree(OPTIONS *Opt, TREES *Trees, RATES *Rates)
{
	TREE *Tree;
	NODE Root;
	
	if(NeedToReSetBL(Opt) == FALSE)
		return;

	Tree = Trees->Tree[Rates->TreeNo];
	Root = Tree->Root;
	

	if(Opt->UseKappa == TRUE)
	{
		if(Opt->EstKappa == TRUE)
			TransContNodeKappa(Root, Rates->Kappa, FALSE);
		else
			TransContNodeKappa(Root, Opt->FixKappa, FALSE);
	}

	if(Opt->UseOU == TRUE)
	{
		if(Opt->EstOU == TRUE)
			TransContNodeOU(Root, Rates->OU, FALSE);
		else
			TransContNodeOU(Root, Opt->FixOU, FALSE);
	}
	
	if(Opt->UseDelta == TRUE)
	{
		if(Opt->EstDelta == TRUE)
			TransContNodeDelta(Root, Rates->Delta, FALSE);
		else
			TransContNodeDelta(Root, Opt->FixDelta, FALSE);
	}

	
	if(Opt->UseLambda == TRUE)
	{
		SetTreeDistToRoot(Tree);
		if(Opt->EstLambda == TRUE)
			TransContNodeLambda(Root, Rates->Lambda, FALSE);
		else
			TransContNodeLambda(Root, Opt->FixLambda, FALSE);
	}
}

void	TransformContrastTreeFixed(OPTIONS *Opt, TREES *Trees)
{
	TREE *Tree;
	NODE Root;
	int		Index;

	for(Index=0;Index<Trees->NoOfTrees;Index++)
	{
		Tree = Trees->Tree[Index];
		Root = Tree->Root;
	
		if((Opt->UseKappa == TRUE) && (Opt->FixKappa != -1))
			TransContNodeKappa(Root, Opt->FixKappa, FALSE);
	
		if((Opt->UseOU == TRUE) && (Opt->FixOU != -1))
			TransContNodeOU(Root, Opt->FixOU, FALSE);
	
		if((Opt->UseDelta == TRUE) && (Opt->FixDelta != -1))
			TransContNodeDelta(Root, Opt->FixDelta, FALSE);

		if((Opt->UseLambda == TRUE) && (Opt->FixLambda != -1))
		{
			SetTreeDistToRoot(Tree);
			TransContNodeLambda(Root, Opt->FixLambda, FALSE);
		}

		SetAsUserBranchLength(Tree);
		SetTreeDistToRoot(Tree);
	}
}

