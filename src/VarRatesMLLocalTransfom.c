#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

#include "GenLib.h"
#include "TypeDef.h"
#include "Part.h"
#include "Likelihood.h"
#include "ML.h"

extern double GoldenOpt(int maxiter, int count, double ax, double bx, double cx, double (*f)(double), double tol, double *xmin);

void SetAllLTVals(RATES *Rates, double Val)
{
	int Index;

	for(Index=0;Index<Rates->NoLocalTransforms;Index++)
		Rates->LocalTransforms[Index]->Scale = Val;
}

void PrintNodeInfo(NODE Node, TREES *Trees)
{
	printf("%d\t%f\t", Node->ID, Node->UserLength);
	PrintPart(stdout, Trees, Node->Part);
	printf("\n");
}

void PrintAllVarRatesNode(VARRATES *VarRates, TREES *Trees)
{
	int Index;
	VAR_RATES_NODE *VRNode;
	NODE Node;

	for(Index=0;Index<VarRates->NoNodes;Index++)
	{
		VRNode = VarRates->NodeList[Index];

		if(VRNode->Type == VR_NODE)
		{
			Node = VRNode->NodeList[0];
			PrintNodeInfo(Node, Trees);
		}
	}
}

void PrintLocalTransfomrs(RATES *Rates, TREES *Trees)
{
	LOCAL_TRANSFORM *LT;
	int Index;
	NODE Node;

	for(Index=0;Index<Rates->NoLocalTransforms;Index++)
	{
		LT = Rates->LocalTransforms[Index];
		Node = LT->TagList[0]->NodeList[0];
		PrintNodeInfo(Node, Trees);

	}
}

VAR_RATES_NODE* GetVarRateOnLocalTransfomr(LOCAL_TRANSFORM *LT,  RATES *Rates)
{
	VARRATES *VarRates;
	VAR_RATES_NODE *VRNode;
	int Index;

	VarRates = Rates->VarRates;

	for(Index=0;Index<VarRates->NoNodes;Index++)
	{
		VRNode = VarRates->NodeList[Index];
		if(VRNode->NodeList[0] == LT->TagList[0]->NodeList[0] && VRNode->Type == VR_NODE)
			return VRNode;
	}

	return NULL;
}

RATES*		GRates;
TREES*		GTrees;
OPTIONS*	GOpt;
LOCAL_TRANSFORM* GLT;

double	Opt1DFuncLT(double Pram)
{
	double	Ret;
		
	GLT->Scale = Pram;

	Ret = Likelihood(GRates, GTrees, GOpt);
	
	return -Ret;
}

void	OptLocalTransform(OPTIONS *Opt, TREES *Trees, RATES *Rates, LOCAL_TRANSFORM *LT)
{
	VAR_RATES_NODE *VRNode;
	double VRNodeVal,InitLh,  Lh, Scale;

	VRNode = GetVarRateOnLocalTransfomr(LT, Rates);

	VRNodeVal = 1.0;
	if(VRNode != NULL)
	{
		VRNodeVal = VRNode->Scale;
		VRNode->Scale = 1.0;
	}
	
	InitLh = Likelihood(Rates, Trees, Opt);

	GRates = Rates;
	GTrees= Trees;
	GOpt = Opt;
	GLT = LT;

	Lh = -GoldenOpt(500, 1, 0., 1.0, 100, Opt1DFuncLT, 0.00001, &Scale);
	printf("%f\t%f\t%f\t%f\t%f\t", InitLh, Lh, Scale, VRNodeVal, Lh-InitLh);


	if(VRNode != NULL)
		VRNode->Scale = VRNodeVal;
}

void	VarRatesMLLocalTransfom(OPTIONS *Opt, TREES *Trees, RATES *Rates)
{
	int Index;
	
//	PrintAllVarRatesNode(Rates->VarRates, Trees);
//	PrintLocalTransfomrs(Rates, Trees);
//	exit(0);

//	Rates->GlobalTrend = 0.004889909;
//	Rates->GlobalTrend = 0.007261267;
//	Rates->GlobalTrend = 0.006129202;
//	Rates->GlobalTrend = 0.005071155;
//	Rates->GlobalTrend = 0.007770396;
	Rates->GlobalTrend = 0.00573875;


	printf("Md5sum\tLh Scale == 1.0\tLh Opt\tOpt Scale\tInitScale\tLhDiff\n");
	for(Index=0;Index<Rates->NoLocalTransforms;Index++)
	{
		SetAllLTVals(Rates, 1.0);

		printf("%s\t", Rates->LocalTransforms[Index]->Name);
		OptLocalTransform(Opt, Trees, Rates, Rates->LocalTransforms[Index]);


		printf("\n");

		fflush(stdout);
	}

	exit(0);
}