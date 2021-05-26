/*
*  BayesTriats 3.0
*
*  copyright 2017
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
#include <math.h>
#include <string.h>

#include "GenLib.h"
#include "TypeDef.h"
#include "Trees.h"
#include "Rates.h"
#include "ML.h"
#include "Likelihood.h"
#include "Praxis.h"
#include "Priors.h"
#include "Continuous.h"
#include "Threaded.h"
#include "Options.h"
#include "Data.h"
#include "TimeSlices.h"
#include "NLOptBT.h"
#include "Landscape.h"
#include "Part.h"
#include "StateSpeciationRate.h"
#include "LocalTransformMLAllNodes.h"

#include <gsl/gsl_matrix.h>

#define MAX_ML_FREE_P 1048576
#define NO_RAND_TRIES 10000


void Opt1D(ML_MAP* Map, OPTIONS *Opt, TREES *Trees, RATES *Rates);


ML_MAP*	AllocMLMap(void)
{
	ML_MAP*	Ret;
	

	Ret = (ML_MAP*)SMalloc(sizeof(ML_MAP));
	
	Ret->PVal = (double*)SMalloc(sizeof(double) * MAX_ML_FREE_P);
	Ret->PMin = (double*)SMalloc(sizeof(double) * MAX_ML_FREE_P);
	Ret->PMax = (double*)SMalloc(sizeof(double) * MAX_ML_FREE_P);
	Ret->PDef = (double*)SMalloc(sizeof(double) * MAX_ML_FREE_P);
	Ret->PType = (ML_P_TYPE*)SMalloc(sizeof(ML_P_TYPE) * MAX_ML_FREE_P);
		
	Ret->NoP = 0;

	return Ret;
}

void	FreeMLMap(ML_MAP *MLMap)
{
	free(MLMap->PVal);
	free(MLMap->PDef);
	free(MLMap->PMin);
	free(MLMap->PMax);
	free(MLMap->PType);
	
	free(MLMap);
}

void	CopyMLMap(ML_MAP *A, ML_MAP *B)
{
	A->NoP = B->NoP;

	memcpy(A->PVal, B->PVal, sizeof(double) * A->NoP);
	memcpy(A->PMin, B->PMin, sizeof(double) * A->NoP);
	memcpy(A->PMax, B->PMax, sizeof(double) * A->NoP);
	memcpy(A->PDef, B->PDef, sizeof(double) * A->NoP);
	memcpy(A->PType, B->PType, sizeof(ML_P_TYPE) * A->NoP);
}

void	AddTypePToMLMap(ML_MAP*	MLMap, double DefV, double MinV, double MaxV, ML_P_TYPE Type)
{
	MLMap->PDef[MLMap->NoP] = DefV;
	MLMap->PMin[MLMap->NoP] = MinV;
	MLMap->PMax[MLMap->NoP] = MaxV;
	MLMap->PType[MLMap->NoP] = Type;
	MLMap->NoP++;
}

void	AddPToMLMap(ML_MAP*	MLMap, double DefV, double MinV, double MaxV)
{
	AddTypePToMLMap(MLMap, DefV, MinV, MaxV, ML_P_TYPE_NONE);
}

void	BuildMLMap(ML_MAP*	MLMap, OPTIONS *Opt, TREES *Trees, RATES *Rates)
{
	int Index;
	double MidPoint;
	TIME_SLICE *TS;
	LOCAL_TRANSFORM *LT;

	for(Index=0;Index<Rates->NoOfRates;Index++)
	{
		if(Opt->RateMax > 1.0)
			AddPToMLMap(MLMap, 1.0, Opt->RateMin, Opt->RateMax);
		else
		{
			MidPoint = (Opt->RateMax + Opt->RateMin) * .5;
			AddPToMLMap(MLMap, MidPoint, Opt->RateMin, Opt->RateMax);
		}

	}

	if(Opt->EstKappa == TRUE)
		AddPToMLMap(MLMap, 1.0, MIN_KAPPA, MAX_KAPPA);

	if(Opt->EstLambda == TRUE)
		AddPToMLMap(MLMap, 1.0, MIN_LAMBDA, MAX_LAMBDA);

	if(Opt->EstDelta == TRUE)
		AddPToMLMap(MLMap, 1.0, MIN_DELTA, MAX_DELTA);

	if(Opt->EstOU == TRUE)
		AddPToMLMap(MLMap, 1.0, MIN_OU, MAX_OU);
	//	AddPToMLMap(MLMap, MIN_OU, MIN_OU, MAX_OU);

	if(Opt->EstGamma == TRUE)
		AddPToMLMap(MLMap, 1.0, MIN_GAMMA, MAX_GAMMA);

	if(Opt->UseGlobalTrend == TRUE)
		AddPToMLMap(MLMap, 0.0, MIN_GLOBAL_TREND, MAX_GLOBAL_TREND);

	for(Index=0;Index<Rates->NoLocalTransforms;Index++)
	{
		LT = Rates->LocalTransforms[Index];
		if(LT->Est == TRUE)
		{
			if(LT->Type == VR_OU)
				AddPToMLMap(MLMap, MIN_OU, MIN_LOCAL_RATE, MAX_LOCAL_RATE);
			else
			{
				if(LT->Type == VR_NODE || LT->Type == VR_BL)
					AddTypePToMLMap(MLMap, 1.0, MIN_LOCAL_RATE, MAX_LOCAL_RATE, ML_P_TYPE_RATE_S);
				else
				{
					if(LT->Type == VR_LS_BL)
						AddPToMLMap(MLMap, 1.0, -MAX_LOCAL_RATE, MAX_LOCAL_RATE);
					else
						AddPToMLMap(MLMap, 1.0, MIN_LOCAL_RATE, MAX_LOCAL_RATE);
				}
			}
		}
	}

	if(Rates->TimeSlices != NULL)
	{
		for(Index=0;Index<Opt->TimeSlices->NoTimeSlices;Index++)
		{
			TS = Opt->TimeSlices->TimeSlices[Index];
			if(TS->FixedTime == FALSE)
				AddPToMLMap(MLMap, RandDouble(Rates->RS), 0.0, 1.0);

			if(TS->FixedScale == FALSE)
				AddTypePToMLMap(MLMap, 1.0, MIN_LOCAL_RATE, MAX_LOCAL_RATE, ML_P_TYPE_RATE_S);
		}
	}
}

void	CheckMLMapVals(ML_MAP* MLMap)
{
	int Index;

	for(Index=0;Index<MLMap->NoP;Index++)
	{
		if(MLMap->PVal[Index] > MLMap->PMax[Index])
			MLMap->PVal[Index] = MLMap->PMax[Index];

		if(MLMap->PVal[Index] < MLMap->PMin[Index])
			MLMap->PVal[Index] = MLMap->PMin[Index];
	}
}

void	MLMapToRatesTimeSlices(ML_MAP* MLMap, OPTIONS *Opt, RATES *Rates, int *Pos)
{
	int Index;
	TIME_SLICE *TS;

	for(Index=0;Index<Opt->TimeSlices->NoTimeSlices;Index++)
	{
		TS = GetTimeSlice(Rates->TimeSlices, Opt->TimeSlices->TimeSlices[Index]->Name);
		if(TS->FixedTime == FALSE)
			TS->Time = MLMap->PVal[(*Pos)++];

		if(TS->FixedScale == FALSE)
			TS->Scale = MLMap->PVal[(*Pos)++];
	}
}

void	MLMapToRates(ML_MAP* MLMap, OPTIONS *Opt, RATES *Rates)
{
	int Index, Pos;

	Pos = 0;

	CheckMLMapVals(MLMap);

	for(Index=0;Index<Rates->NoOfRates;Index++)
		Rates->Rates[Index] = MLMap->PVal[Pos++];

	if(Opt->EstKappa == TRUE)
		Rates->Kappa = MLMap->PVal[Pos++];

	if(Opt->EstLambda == TRUE)
		Rates->Lambda = MLMap->PVal[Pos++];

	if(Opt->EstDelta == TRUE)
		Rates->Delta = MLMap->PVal[Pos++];

	if(Opt->EstOU == TRUE)
		Rates->OU = MLMap->PVal[Pos++];

	if(Opt->EstGamma == TRUE)
		Rates->Gamma = MLMap->PVal[Pos++];

	if(Opt->UseGlobalTrend == TRUE)
		Rates->GlobalTrend = MLMap->PVal[Pos++];

	for(Index=0;Index<Rates->NoLocalTransforms;Index++)
		if(Rates->LocalTransforms[Index]->Est == TRUE)
			Rates->LocalTransforms[Index]->Scale = MLMap->PVal[Pos++];

	MLMapToRatesTimeSlices(MLMap, Opt, Rates, &Pos);


}

double	LikelihoodML(ML_MAP* MLMap, OPTIONS *Opt, TREES *Trees, RATES *Rates)
{
	MLMapToRates(MLMap, Opt, Rates);

 	return Likelihood(Rates, Trees, Opt);
}

void	MLMapSetDefVals(ML_MAP* MLMap)
{
	memcpy(MLMap->PVal, MLMap->PDef, sizeof(double) * MLMap->NoP);
}

double	GetRandRateScale(RANDSTATES *RS, double Min, double Max)
{
	if(RandDouble(RS) < 0.5)
		return RandUniDouble(RS, Min, 1.0);

	return RandUniDouble(RS, 1.0, Max);
}

void	MLMapSetRandVals(ML_MAP* MLMap, RANDSTATES *RS)
{
	int Index;

	for(Index=0;Index<MLMap->NoP;Index++)
	{
		if(MLMap->PType[Index] == ML_P_TYPE_NONE)
			MLMap->PVal[Index] = RandUniDouble(RS, MLMap->PMin[Index], MLMap->PMax[Index]);
		
		if(MLMap->PType[Index] == ML_P_TYPE_RATE_S)
			MLMap->PVal[Index] = GetRandRateScale(RS, MLMap->PMin[Index], MLMap->PMax[Index]);
	}
}

void	MLMapSetRatesFixedVals(ML_MAP *MLMap, RATES *Rates, double Val)
{
	int Index;

	MLMapSetDefVals(MLMap);

	for(Index=0;Index<Rates->NoOfRates;Index++)
		MLMap->PVal[Index] = Val;
}

void	FindValidMLStartSet(ML_MAP *MLMap, OPTIONS *Opt, TREES *Trees, RATES *Rates)
{
	int Index;
	double Lh;
	
	for(Index=0;Index<NO_RAND_TRIES;Index++)
	{
		MLMapSetRandVals(MLMap, Rates->RS);

		Lh = LikelihoodML(MLMap, Opt, Trees, Rates);

		if(Lh != ERRLH)
			return;
	}

	MLMapSetDefVals(MLMap);
	Lh = LikelihoodML(MLMap, Opt, Trees, Rates);
	if(Lh != ERRLH)
		return;

	printf("Cannot find a valid starting set of parameters.\n");
	exit(0);
}

double	LhPraxis(void *P, double *List)
{
	double		Ret;
	PRAXSTATE	*PState;
	ML_MAP		*MLMap;

	PState = (PRAXSTATE*)P;
	MLMap = (ML_MAP*)PState->Pt;

	memcpy(MLMap->PVal, List, sizeof(double) * MLMap->NoP);
	
	Ret = LikelihoodML(MLMap ,  PState->Opt, PState->Trees, PState->Rates);

	PState->NoLhCalls++;

	if(Ret == ERRLH)
		return -ERRLH;
	
	return -Ret;
}

double*	MLMapClonePVect(ML_MAP*	MLMap)
{
	return (double*)CloneMem(sizeof(double) * MLMap->NoP, (void*)MLMap->PVal);
}

ML_MAP*	MLMapTreeTry(OPTIONS *Opt, TREES *Trees, RATES *Rates, ML_MAP *Init)
{
	ML_MAP		*Ret;
	PRAXSTATE	*PState;
	double		*TVect, Lh;

	Ret = AllocMLMap();

	BuildMLMap(Ret, Opt, Trees, Rates);
	
	FindValidMLStartSet(Ret, Opt, Trees, Rates);	

	if(Init != NULL)
		CopyMLMap(Ret, Init);


	Lh = LikelihoodML(Ret, Opt, Trees, Rates);

	TVect = MLMapClonePVect(Ret);

	PState = NULL;

#ifndef NLOPT
	if(Ret->NoP > 1)
		PState = IntiPraxis(LhPraxis, TVect, Ret->NoP, 0, 1, 4, 50000);
	else
		PState = IntiPraxis(LhPraxis, TVect, Ret->NoP, 0, 1, 1, 250);

	PState->Opt		= Opt;
	PState->Trees	= Trees;
	PState->Rates	= Rates;
	PState->Pt		= (void*)Ret;

	Lh = praxis(PState);
	memcpy(Ret->PVal, TVect, sizeof(double) * Ret->NoP);
	FreePracxStates(PState);
#else
	NLOptBT(Rates, Opt, Trees, Ret);
#endif
	
	Lh = LikelihoodML(Ret, Opt, Trees, Rates);

	
	free(TVect);
	return Ret;
}


void	MLTree(OPTIONS *Opt, TREES *Trees, RATES *Rates)
{
	ML_MAP*	CMap, *BMap;
	double CLh, BLh;
	int Index;
			
	BMap = AllocMLMap();

	BuildMLMap(BMap, Opt, Trees, Rates);
	BLh = LikelihoodML(BMap, Opt, Trees, Rates);

	if(ValidLh(BLh, Opt->Model) == FALSE)
	{
		FindValidMLStartSet(BMap, Opt, Trees, Rates);	
		BLh = LikelihoodML(BMap, Opt, Trees, Rates);
	}

	if(BMap->NoP != 0)
	{
#ifndef NLOPT
		if(BMap->NoP == 1)
			Opt1D(BMap, Opt, Trees, Rates);
		else
#endif
		{
			for(Index=0;Index<Opt->MLTries;Index++)
			{
				if(Index!=0)
					CMap = MLMapTreeTry(Opt, Trees, Rates, NULL);
				else
					CMap = MLMapTreeTry(Opt, Trees, Rates, BMap);

				
				
				CLh = LikelihoodML(CMap, Opt, Trees, Rates);
			
				if(CLh > BLh)
				{
					CopyMLMap(BMap, CMap);
					BLh = CLh;
				}

				FreeMLMap(CMap);
			}
		}
	}
	
	Rates->Lh = LikelihoodML(BMap, Opt, Trees, Rates);
	FreeMLMap(BMap);
}


/*
void	CalclAllNodeBeta(OPTIONS *Opt, TREES *Trees, RATES *Rates)
{
	TREE *Tree;
	int Index;
	NODE Node;
	double ILh, OLh, BL;

	Tree = Trees->Tree[0];

	
	BlankLandscape(Rates->Landscape);
	AddLandscapeFromPart(Rates->Landscape, Tree->Root->Part, Trees, 0.0);

	for(Index=1;Index<Tree->NoNodes;Index++)
	{
		Node = Tree->NodeList[Index];
		
		Rates->Landscape->NodeList[0]->Beta = 0.0;
		Rates->Landscape->NodeList[0]->Part = Node->Part;
		Rates->Landscape->NodeList[0]->NodeList[0] = NULL;

		ILh = Likelihood(Rates, Trees, Opt);
		
				
		MLTree(Opt, Trees, Rates);

		OLh = Likelihood(Rates, Trees, Opt);
		printf("Lh:\t%f\t%f\t%f\n", ILh, OLh, Rates->Landscape->NodeList[0]->Beta, Node->Length);
		fflush(stdout);
	}
}
*/

/*
void	CalclAllNodeBeta(OPTIONS *Opt, TREES *Trees, RATES *Rates)
{
	TREE *Tree;
	int Index;
	NODE Node;
	double ILh, OLh, BL;

	Tree = Trees->Tree[0];
	
	for(Index=1;Index<Tree->NoNodes;Index++)
	{
		Node = Tree->NodeList[Index];

		BlankLandscape(Rates->Landscape);

		AddLandscapeFromPart(Rates->Landscape, Node->Part, Trees, 0.0);

		ILh = Likelihood(Rates, Trees, Opt);
		
		MLTree(Opt, Trees, Rates);

		OLh = Likelihood(Rates, Trees, Opt);
		printf("Lh:\t%f\t%f\t", ILh, OLh);

		printf("%f\t%f\t", Rates->Landscape->NodeList[0]->Beta, Node->Length);

		PrintPart(stdout, Trees, Node->Part);

		printf("\n");
		fflush(stdout);
	}
}
*/

void	PrintPartData(FILE *Str, TREES *Trees, PART *Part)
{
	int Index, ID;
	TAXA *T;

	fprintf(Str, "%zu\t%d\t%f\t", Part->PartID, Part->Freq, Part->Prob);
	fprintf(Str, "\t%d\t", Part->NoTaxa);

	for(Index=0;Index<Part->NoTaxa;Index++)
	{
		ID = Part->Taxa[Index];
		T = Trees->Taxa[ID];
		fprintf(Str, "%f\t", T->ConData[0]);
	}
}

void	CalcAllNodeTransfroms(OPTIONS *Opt, TREES *Trees, RATES *Rates)
{
	TREE *Tree;
	int Index;
	NODE Node;
	double ILh, OLh, Height;
	
	Tree = Trees->Tree[0];
	
	for(Index=1;Index<Tree->NoNodes;Index++)
	{
		Node = Tree->NodeList[Index];

		if(Node->Tip == FALSE)
		{
			Rates->LocalTransforms[0]->TagList[0]->NodeList[0] = Node;
			Rates->LocalTransforms[0]->Scale = 1.0;
			
			ILh = Likelihood(Rates, Trees, Opt);

			MLTree(Opt, Trees, Rates);

			OLh = Likelihood(Rates, Trees, Opt);
			fprintf(Opt->LogFile, "%s\tLh:\t%f\t%f\t", Opt->DataFN, ILh, OLh);

			fprintf(Opt->LogFile, "%.12f\t%f\t", Rates->LocalTransforms[0]->Scale, Node->UserLength);

			Height = GetNodeHeight(Node);
			fprintf(Opt->LogFile, "%f\t", Height);
							
			PrintPart(Opt->LogFile, Trees, Node->Part);
	//		PrintPartData(stdout, Trees, Node->Part);
	//		RecPrintNode(Node);
			fprintf(Opt->LogFile, "\n");
			fflush(Opt->LogFile);
		}
	}

	exit(0);
}

void	AddTipVal(TREE *Tree, double Scale)
{
	int Index;
	NODE N;

	for(Index=0;Index<Tree->NoNodes;Index++)
	{
		N = Tree->NodeList[Index];
		if(N->Tip == TRUE)
		{
			N->Length = N->UserLength + (1.0 + Scale);
			//N->Length = N->UserLength + Scale;
//			N->UserLength = N->Length;
		}
	}
}

void	TestErrModel(OPTIONS *Opt, RATES *Rates, TREES *Trees)
{
	TREE *Tree;
	double	Scale, Lh;
//	return;
	Tree = Trees->Tree[0];


//	AddTipVal(Tree, 18.254);
//	SaveTrees("sout.trees", Trees);
//	exit(0);

	for(Scale=0;Scale<25;Scale+=0.01)
	{
		AddTipVal(Tree, Scale);
		Lh = Likelihood(Rates, Trees, Opt);
		printf("%f\t%f\n", Scale, Lh);
	}

	exit(0);
}

void	MLTest(OPTIONS *Opt, TREES *Trees, RATES *Rates)
{
	double X, Lh;
	return;
	Rates->Rates[0] = 9.998893;
	Rates->Rates[1] = 0.233734;
	Rates->Rates[2] = 8.434539;
	Rates->Rates[3] = 2.276221;
	Rates->Rates[4] = 0;
	Rates->Rates[5] = 3.056001;
	Rates->Rates[6] = 0.055404;
	Rates->Rates[7] = 0.205653;

//	Rates->OffToOn = 0.205653;
//	Rates->OnToOff = 0.205653;


	for(X=0.00000001;X<1000; X+=0.01)
	{
		Rates->Rates[0] = X;

		Lh = Likelihood(Rates, Trees, Opt);

		printf("%f\t%f\n", X, Lh);
	}

	printf("Lh:\t%f\n", Lh);

	exit(0);
}

void	PrintMLHeader()
{
	printf("Tree No\tLh\tLh Elapsed Seconds");
	printf("\n");
	fflush(stdout);
}

void PrintMLStdOutHeader(void)
{
	printf("Tree No\tLh\tElapsed Seconds\n");

}

void PrintMLTree(int TreeNo, double Lh, double Sec)
{
	printf("%d\t%f\t%f", TreeNo, Lh, Sec);
	printf("\n");
	fflush(stdout);
}


void	FindML(OPTIONS *Opt, TREES *Trees)
{
	RATES *Rates; 
	int Index;
	double Lh;
	double	TStart, TEnd;

	Rates = CreatRates(Opt);

	PrintOptions(Opt->LogFile, Opt);
	PrintRatesHeadder(Opt->LogFile, Opt);

	PrintOptions(stdout, Opt);
	PrintMLStdOutHeader();
//	PrintRatesHeadder(stdout, Opt);


	fflush(stdout);
	fflush(Opt->LogFile);

	TStart = GetSeconds();

	if(Opt->UseMLLandscape == TRUE)
		LocalTransformMLAllNodes(Opt, Trees, Rates, TRUE, FALSE, TRUE);
		

//	CalcAllNodeTransfroms(Opt, Trees, Rates); return;
//	MLTest(Opt, Trees, Rates);
//	TestErrModel(Opt, Rates, Trees);
//	TestContrastGlobalTrend(Opt, Trees, Rates);
	
	for(Index=0;Index<Trees->NoTrees;Index++)
	{
		Rates->TreeNo = Index;

		if(Opt->ModelType == MT_CONTINUOUS)
			InitContinusTree(Opt, Trees, Rates->TreeNo);

		if(Opt->NodeData == TRUE || Opt->NodeBLData == TRUE)
			SetTreeAsData(Opt, Trees, Rates->TreeNo);


		if(Opt->UseMLLandscape == TRUE)
			MLLandscape(Opt, Trees, Rates);

		MLTree(Opt, Trees, Rates);

		Lh = Likelihood(Rates, Trees, Opt);
		
		PrintRates(Opt->LogFile, Rates, Opt, NULL);
		fprintf(Opt->LogFile, "\n");
		fflush(Opt->LogFile);

		PrintMLTree(Index, Rates->Lh, GetSeconds() - TStart);

		if(Opt->SaveTrees == TRUE)
			OutputTree(Opt, Trees, Rates, Index+1, Opt->OutTrees);
		
		if(Opt->ModelType == MT_CONTINUOUS)
		{
			FreeConVar(Trees->Tree[Rates->TreeNo]->ConVars, Trees->NoTaxa);
			Trees->Tree[Rates->TreeNo]->ConVars = NULL;
		}
	}

	TEnd = GetSeconds();
	printf("Sec:\t%f\n", TEnd - TStart);

//	CaclStateSpeciationRateLh(Opt, Trees, Rates);

	FreeRates(Rates, Trees);
}

#define SHFT2(a,b,c) (a)=(b);(b)=(c);
#define SHFT3(a,b,c,d) (a)=(b);(b)=(c);(c)=(d);

double GoldenOpt(int maxiter, int count, double ax, double bx, double cx, double (*f)(double), double tol, double *xmin)
{
	double f1,f2,x0,x1,x2,x3;
	double R, C;

	R =  0.61803399 ;
	C = 1.0 - R;

	if (maxiter<=0)
	{
		*xmin = bx;
		count = 0;
		return 0;
	}
	else if (maxiter==1)
	{
		*xmin = bx;
		count = 1;
		return (*f)(bx);
	}

	x0 = ax; 
	x3 = cx;
	if (fabs(cx-bx) > fabs(bx-ax)) 
	{ 
		x1 = bx;
		x2 = bx+C*(cx-bx); 
	} 
	else 
	{
		x2 = bx;
		x1 = bx-C*(bx-ax);
	}
	f1 = (*f)(x1); 

	f2 = (*f)(x2);
	count = 2;
	while (count<maxiter && fabs(x3-x0)>tol*(fabs(x1)+fabs(x2))) 
	{
		if (f2 < f1) 
		{ 
			SHFT3(x0,x1,x2,R*x1+C*x3) 
				SHFT2(f1,f2,(*f)(x2)) 
		} 
		else 
		{
			SHFT3(x3,x2,x1,R*x2+C*x0)
				SHFT2(f2,f1,(*f)(x1)) 
		}
		count++;
	} 
	if (f1 < f2) 
	{ 
		*xmin=x1;
		return f1;
	} 
	else 
	{
		*xmin=x2;
		return f2;
	}
}

RATES*		FRates;
TREES*		FTrees;
OPTIONS*	FOpt;
ML_MAP*		FMap;

double	Opt1DFunc(double Pram)
{
	double	Ret;
	
	FMap->PVal[0] = Pram;

	Ret = LikelihoodML(FMap, FOpt, FTrees, FRates);

	return -Ret;
}

void Opt1DAcc(double *BLh, double *BVal, double *CLh, double *CVal, MODEL_TYPE MT)
{
	if(ValidLh(*CLh, MT) == FALSE)
		return;

	if(*CLh > *BLh)
	{
		*BLh = *CLh;
		*BVal = *CVal;
	}
}

void Opt1D(ML_MAP* Map, OPTIONS *Opt, TREES *Trees, RATES *Rates)
{
	double BLh, BVal;
	double CLh, CVal;
	double Tol, Val;
	int Index;

	Tol = 0.0000001;

	FRates = Rates;
	FTrees = Trees;
	FOpt = Opt;
	FMap = Map;
	
	BLh = LikelihoodML(Map, Opt, Trees, Rates);
	BVal = Map->PVal[0];

	// Test Min
	CLh = -GoldenOpt(500, 1, Map->PMin[0], Map->PMin[0], Map->PMax[0], Opt1DFunc, Tol, &CVal);
	Opt1DAcc(&BLh, &BVal, &CLh, &CVal, Opt->ModelType);

	// Test Max
	CLh = -GoldenOpt(500, 1, Map->PMin[0], Map->PMax[0], Map->PMax[0], Opt1DFunc, Tol, &CVal);
	Opt1DAcc(&BLh, &BVal, &CLh, &CVal, Opt->ModelType);

	// Def Min
	CLh = -GoldenOpt(500, 1, Map->PMin[0], Map->PDef[0], Map->PMax[0], Opt1DFunc, Tol, &CVal);
	Opt1DAcc(&BLh, &BVal, &CLh, &CVal, Opt->ModelType);

	// A set of Random
	for(Index=0;Index<Opt->MLTries;Index++)
	{
		Val = Map->PMin[0] + (RandDouble(Rates->RS) * (Map->PMax[0] - Map->PMin[0]));
		CLh = -GoldenOpt(500, 1, Map->PMin[0], Val, Map->PMax[0], Opt1DFunc, Tol, &CVal);
		Opt1DAcc(&BLh, &BVal, &CLh, &CVal, Opt->ModelType);
	}

	Map->PVal[0] = BVal;
	BLh = LikelihoodML(Map, Opt, Trees, Rates);
}
