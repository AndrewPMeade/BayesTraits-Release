#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

#include "genlib.h"
#include "typedef.h"
#include "trees.h"
#include "Rates.h"
#include "ml.h"
#include "likelihood.h"
#include "praxis.h"
#include "priors.h"
#include "continuous.h"
#include "Threaded.h"
#include "options.h"
#include "data.h"

#define MAX_ML_FREE_P 1048576
#define NO_RAND_TRIES 10000


typedef struct
{
	int NoP;

	double	*PVal;
	double	*PMin;
	double	*PMax;
	double	*PDef;
	PRIOR	**Priors;

} ML_MAP;

void Opt1D(ML_MAP* Map, OPTIONS *Opt, TREES *Trees, RATES *Rates);


ML_MAP*	AllocMLMap(void)
{
	ML_MAP*	Ret;

	Ret = (ML_MAP*)malloc(sizeof(ML_MAP));
	if(Ret == NULL)
		MallocErr();

	Ret->PVal = (double*)malloc(sizeof(double) * MAX_ML_FREE_P);
	Ret->PMin = (double*)malloc(sizeof(double) * MAX_ML_FREE_P);
	Ret->PMax = (double*)malloc(sizeof(double) * MAX_ML_FREE_P);
	Ret->PDef = (double*)malloc(sizeof(double) * MAX_ML_FREE_P);

	if( Ret->PVal == NULL ||
		Ret->PMin == NULL ||
		Ret->PMax == NULL ||
		Ret->PDef == NULL)
		MallocErr();

	Ret->NoP = 0;

	return Ret;
}

void	FreeMLMap(ML_MAP *MLMap)
{
	free(MLMap->PVal);
	free(MLMap->PDef);
	free(MLMap->PMin);
	free(MLMap->PMax);
	
	free(MLMap);
}

void	CopyMLMap(ML_MAP *A, ML_MAP *B)
{
	A->NoP = B->NoP;

	memcpy(A->PVal, B->PVal, sizeof(double) * A->NoP);
	memcpy(A->PMin, B->PMin, sizeof(double) * A->NoP);
	memcpy(A->PMax, B->PMax, sizeof(double) * A->NoP);
	memcpy(A->PDef, B->PDef, sizeof(double) * A->NoP);
}

void	AddPToMLMap(ML_MAP*	MLMap, double DefV, double MinV, double MaxV)
{
	MLMap->PDef[MLMap->NoP] = DefV;
	MLMap->PMin[MLMap->NoP] = MinV;
	MLMap->PMax[MLMap->NoP] = MaxV;
	MLMap->NoP++;
}

void	BuildMLMap(ML_MAP*	MLMap, OPTIONS *Opt, TREES *Trees, RATES *Rates)
{
	int Index;

	for(Index=0;Index<Rates->NoOfRates;Index++)
		AddPToMLMap(MLMap, 1.0, MINRATE, MAXRATE);

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

	for(Index=0;Index<Rates->NoLocalTransforms;Index++)
		if(Rates->LocalTransforms[Index]->Est == TRUE)
			if(Rates->LocalTransforms[Index]->Type == VR_OU)
				AddPToMLMap(MLMap, MIN_OU, MIN_LOCAL_RATE, MAX_LOCAL_RATE);
			else
				AddPToMLMap(MLMap, 1.0, MIN_LOCAL_RATE, MAX_LOCAL_RATE);
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

	for(Index=0;Index<Rates->NoLocalTransforms;Index++)
		if(Rates->LocalTransforms[Index]->Est == TRUE)
			Rates->LocalTransforms[Index]->Scale = MLMap->PVal[Pos++];
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

void	MLMapSetRandVals(ML_MAP* MLMap, RANDSTATES *RS)
{
	int Index;

	for(Index=0;Index<MLMap->NoP;Index++)
	{
		MLMap->PVal[Index] = RandDouble(RS) * (MLMap->PMax[Index] - MLMap->PMin[Index]);
		MLMap->PVal[Index] += MLMap->PMin[Index];
	}
}

void	MLMapSetRatesFixedVals(ML_MAP* MLMap, RATES *Rates, double Val)
{
	int Index;

	MLMapSetDefVals(MLMap);

	for(Index=0;Index<Rates->NoOfRates;Index++)
		MLMap->PVal[Index] = Val;
}

void	FindValidMLStartSet(ML_MAP* MLMap, OPTIONS *Opt, TREES *Trees, RATES *Rates)
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

double	LhPraxis(void* P, double *List)
{
	double		Ret;
	PRAXSTATE	*PState;
	ML_MAP		*MLMap;

	PState = (PRAXSTATE*)P;
	MLMap = (ML_MAP*)PState->Pt;

	memcpy(MLMap->PVal, List, sizeof(double) * MLMap->NoP);
	
	Ret = LikelihoodML(MLMap ,  PState->Opt, PState->Trees, PState->Rates);

	PState->NoLhCalls++;

	if(ValidLh(Ret) == FALSE)
		return -ERRLH;
	
	return -Ret;
}

double*	MLMapClonePVect(ML_MAP*	MLMap)
{
	double *Ret;

	Ret = (double*)malloc(sizeof(double) * MLMap->NoP);
	if(Ret == NULL)
		MallocErr();

	memcpy(Ret, MLMap->PVal, sizeof(double) * MLMap->NoP);

	return Ret;
}

ML_MAP*	MLMapTreeTry(OPTIONS *Opt, TREES *Trees, RATES *Rates)
{
	ML_MAP		*Ret;
	PRAXSTATE	*PState;
	double		*TVect, Lh;

	Ret = AllocMLMap();

	BuildMLMap(Ret, Opt, Trees, Rates);
	
	FindValidMLStartSet(Ret, Opt, Trees, Rates);	

	Lh = LikelihoodML(Ret, Opt, Trees, Rates);

	TVect = MLMapClonePVect(Ret);

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

	Lh = LikelihoodML(Ret, Opt, Trees, Rates);
	


	FreePracxStates(PState);

	free(TVect);
	return Ret;
}

void	MLTest(OPTIONS *Opt, TREES *Trees, RATES *Rates)
{
	double OU, Lh;

	for(OU=0.0001;OU<10.0;OU+=0.01)
	{
		Rates->OU = OU;
		Lh = Likelihood(Rates, Trees, Opt);

		printf("%f\t%f\n", OU, Lh);
	}

	exit(0);
}

//./Seq/AASacra.trees ./Seq/AASacra.txt < in.txt > sout.txt
//./Seq/a67.1000.trees ./Seq/descent.txt < in.txt > sout.txt

// ./Seq/MLTest/Tree-00000717.trees ./Seq/MLTest/Data-00000717.txt < in.txt > sout.txt
// ./Seq/Testing/PrimatesBody.trees ./Seq/Testing/PrimatesBody.txt < ./Seq/Testing/in.txt > ./Seq/Testing/sout.txt
void	MLTree(OPTIONS *Opt, TREES *Trees, RATES *Rates)
{
	ML_MAP*	CMap, *BMap;
	double CLh, BLh;
	int Index;

//	MLTest(Opt, Trees, Rates);
		
	BMap = AllocMLMap();

	BuildMLMap(BMap, Opt, Trees, Rates);
	FindValidMLStartSet(BMap, Opt, Trees, Rates);	
	BLh = LikelihoodML(BMap, Opt, Trees, Rates);

	if(BMap->NoP != 0)
	{
		if(BMap->NoP == 1)
			Opt1D(BMap, Opt, Trees, Rates);
		else
		{
			for(Index=0;Index<Opt->MLTries;Index++)
			{
				CMap = MLMapTreeTry(Opt, Trees, Rates);
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


void	FindML(OPTIONS *Opt, TREES *Trees)
{
	RATES *Rates; 
	int Index;
	double Lh;
	double	TStart, TEnd;

	Rates = CreatRates(Opt);
	
	PrintOptions(stdout, Opt);
	PrintRatesHeadder(stdout, Opt);

	PrintOptions(Opt->LogFile, Opt);
	PrintRatesHeadder(Opt->LogFile, Opt);

	fflush(stdout);
	fflush(Opt->LogFile);

	TStart = GetSeconds();

	for(Index=0;Index<Trees->NoOfTrees;Index++)
	{
		Rates->TreeNo = Index;

		if(Opt->ModelType == MT_CONTINUOUS)
			InitContinusTree(Opt, Trees, Rates->TreeNo);

		if((Opt->NodeData == TRUE) || (Opt->NodeBLData == TRUE))
			SetTreeAsData(Opt, Trees, Rates->TreeNo);

		MLTree(Opt, Trees, Rates);

		Lh = Likelihood(Rates, Trees, Opt);


		PrintRates(Opt->LogFile, Rates, Opt, NULL);
		fprintf(Opt->LogFile, "\n");
		fflush(Opt->LogFile);

		PrintRates(stdout, Rates, Opt, NULL);
		printf("\n");
		fflush(stdout);
		
		if(Opt->ModelType == MT_CONTINUOUS)
		{
			FreeConVar(Trees->Tree[Rates->TreeNo]->ConVars, Trees->NoOfTaxa);
			Trees->Tree[Rates->TreeNo]->ConVars = NULL;
		}
	}

	TEnd = GetSeconds();
	printf("Sec:\t%f\n", TEnd - TStart);

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

void Opt1DAcc(double *BLh, double *BVal, double *CLh, double *CVal)
{
	if(ValidLh(*CLh) == FALSE)
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
	Opt1DAcc(&BLh, &BVal, &CLh, &CVal);

	// Test Max
	CLh = -GoldenOpt(500, 1, Map->PMin[0], Map->PMax[0], Map->PMax[0], Opt1DFunc, Tol, &CVal);
	Opt1DAcc(&BLh, &BVal, &CLh, &CVal);

	// Def Min
	CLh = -GoldenOpt(500, 1, Map->PMin[0], Map->PDef[0], Map->PMax[0], Opt1DFunc, Tol, &CVal);
	Opt1DAcc(&BLh, &BVal, &CLh, &CVal);

	// A set of Random
	for(Index=0;Index<Opt->MLTries;Index++)
	{
		Val = Map->PMin[0] + (RandDouble(Rates->RS) * (Map->PMax[0] - Map->PMin[0]));
		CLh = -GoldenOpt(500, 1, Map->PMin[0], Val, Map->PMax[0], Opt1DFunc, Tol, &CVal);
		Opt1DAcc(&BLh, &BVal, &CLh, &CVal);
	}

	Map->PVal[0] = BVal;
	BLh = LikelihoodML(Map, Opt, Trees, Rates);
}
