#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "typedef.h"
#include "genlib.h"
#include "trees.h"
#include "TransformTree.h"
#include "continuous.h"

void	SetFixedConTrans(OPTIONS *Opt, TREES *Trees);
void	SetConTrans(TREE *Tree, double Kappa, double Lambda, double Delta, double OU);

void	RecTransContNodeDelta(NODE N, double Delta, double PathLen)
{
	int Index;
	double TLen;

	TLen = N->Length + PathLen;

	N->Length = pow(PathLen+N->Length, Delta) - pow(PathLen, Delta);

	if(N->Tip == TRUE)
		return;

	for(Index=0;Index<N->NoNodes;Index++)
		RecTransContNodeDelta(N->NodeList[Index], Delta, TLen);
}

void	TransformTreeDelta(NODE N, double Delta, int Norm)
{
	int Index;
	double SumBL,Scale;

	if(Norm == TRUE)
		SumBL = SumNodeBL(N);

	for(Index=0;Index<N->NoNodes;Index++)
		RecTransContNodeDelta(N->NodeList[Index], Delta, 0);
//	RecTransContNodeDelta(N, Delta, 0);

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

void	TransformTreeKappa(NODE N, double Kappa, int Norm)
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



void	GetOUTTaxaDistRec(NODE N, int ID, double *Dist)
{
	int Index;

	if(N->Tip == TRUE)
	{
		if(N->TipID == ID)
			*Dist = N->DistToRoot;
		return;
	}

	for(Index=0;Index<N->NoNodes;Index++)
		GetOUTTaxaDistRec(N->NodeList[Index], ID, Dist);
}

double GetOUTaxaDist(TREES *Trees, NODE N, double PathLen)
{
	int XID, YID;
	double Ret, XDist, YDist;

	if(N->Tip == TRUE)
		return 0;	

	if(N->Ans == NULL)
		return N->Length;
	

	XID = Trees->Taxa[N->VPosX]->No;
	YID = Trees->Taxa[N->VPosY]->No;


	XDist = 0;
	YDist = 0;

	GetOUTTaxaDistRec(N, XID, &XDist);
	GetOUTTaxaDistRec(N, YID, &YDist);
	Ret = (XDist - PathLen) + (YDist - PathLen);
	return Ret;
}

double	CaclOU(double PathLen, double OU, double T, double D)
{
	double Ret;

//	Ret = exp(-2.0 * OU * (T - PathLen));
	Ret = exp(-OU * D);
	Ret *=  1.0 - exp(-2.0 * OU * PathLen);

	Ret *= 1.0 / (2.0 * OU);

	return Ret;
}

void	RecTransContNodeOU(TREES *Trees, NODE N, double OU, double T, double PathLen)
{
	int Index;
	double TLen;
	double Dist1, Dist2;

	if(N->VPosY == 0 && N->VPosX == 3)
		printf("ee\n");

	TLen = N->Length + PathLen;
	Dist1 = GetOUTaxaDist(Trees, N, TLen);
	Dist2 = GetOUTaxaDist(Trees, N->Ans, PathLen);
	
//	N->Length = CaclOU(PathLen+N->Length, OU, T, Dist) - CaclOU(PathLen, OU, T, Dist);
	N->Length = CaclOU(PathLen+N->Length, OU, T, Dist1) - CaclOU(PathLen, OU, T, Dist2);

	printf("Len:\t%f\n", N->Length);

	if(N->Tip == TRUE)
		return;

	for(Index=0;Index<N->NoNodes;Index++)
		RecTransContNodeOU(Trees, N->NodeList[Index], OU, T, TLen);
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

void	TestNOUT(NODE N)
{
	int Index;

	if(N->Tip == TRUE)
		return;

	for(Index=0;Index<N->NoNodes;Index++)
		TestNOUT(N->NodeList[Index]);
}

void	TransformTreeOU(TREES *Trees, NODE N, double OU, int Norm)
{
	double SumBL,Scale, T;
	int Index;
	
//	TestNOUT(N);	exit(0);

	if(Norm == TRUE)
		SumBL = SumNodeBL(N);

	T = -1;
	FindOUT(N, &T);

	for(Index=0;Index<N->NoNodes;Index++)
		RecTransContNodeOU(Trees, N->NodeList[Index], OU, T,  0);

	MATRIX *V;

	V = AllocMatrix(Trees->NoOfTaxa, Trees->NoOfTaxa);
	RecSetDistToRoot(N);
	TreeToV(Trees, Trees->Tree[0], V);
//	VToTree(V, Trees->Tree[0]);
	PrintMatrix(V, "V=", stdout);

	SaveTrees("testout.trees", Trees);
	exit(0);

//	VToTree(Trees, Tree);
	if(Norm == FALSE)
		return;

	Scale = SumBL / SumNodeBL(N);
	ScaleSubTree(N, Scale);

	
}


/*
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

void	TestNOUT(NODE N)
{
	int Index;

	if(N->Tip == TRUE)
		return;

	

	for(Index=0;Index<N->NoNodes;Index++)
		TestNOUT(N->NodeList[Index]);
}

void	TransformTreeOU(TREES *Trees, NODE N, double OU, int Norm)
{
	double SumBL,Scale, T;
	int Index;

	TestNOUT(N);
	exit(0);

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

	SaveTrees("testout.trees", Trees);
}

*/

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

void	TransformTreeLambda(NODE N, double Lambda, int Norm)
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

int		NeedToReSetBL(OPTIONS *Opt, RATES *Rates)
{
	if(Rates->UseLocalTransforms == TRUE)
		return TRUE;

	if(Rates->VarRates != NULL)
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

void	TransformTree(OPTIONS *Opt, TREES *Trees, RATES *Rates, int Norm)
{
	TREE *Tree;
	NODE Root;

	Tree = Trees->Tree[Rates->TreeNo];
	Root = Tree->Root;
	
	if(Opt->UseKappa == TRUE)
	{
		if(Opt->EstKappa == TRUE)
			TransformTreeKappa(Root, Rates->Kappa, Norm);
		else
			TransformTreeKappa(Root, Opt->FixKappa, Norm);
	}

	if(Opt->UseOU == TRUE)
	{
		SetTreeDistToRoot(Tree);

		if(Opt->EstOU == TRUE)
			TransformTreeOU(Trees, Root, Rates->OU, Norm);
		else
			TransformTreeOU(Trees, Root, Opt->FixOU, Norm);
	}

	if(Opt->UseDelta == TRUE)
	{
		if(Opt->EstDelta == TRUE)
			TransformTreeDelta(Root, Rates->Delta, Norm);
		else 
			TransformTreeDelta(Root, Opt->FixDelta, Norm);
	}
		
	if(Opt->UseLambda == TRUE)
	{
		SetTreeDistToRoot(Tree);

		if(Opt->EstLambda == TRUE)
			TransformTreeLambda(Root, Rates->Lambda, Norm);
		else
			TransformTreeLambda(Root, Opt->FixLambda, Norm);
	}
}

double	GetTransformDefValue(TRANSFORM_TYPE TranType)
{
	if(TranType == VR_OU)
		return MIN_OU;

	return 1.0;
}

