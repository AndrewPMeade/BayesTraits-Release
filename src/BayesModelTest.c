#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "typedef.h"
#include "trees.h"
#include "data.h"
#include "options.h"
#include "rates.h"
#include "likelihood.h"
#include "rand.h"
#include "priors.h"
#include "mcmc.h"
#include "praxis.h"
#include "ml.h"
#include "genlib.h"
#include "continuous.h"


void	GenDataFormModel(OPTIONS *Opt, TREES *Trees, double *Model, double *Error)
{
	int		NIndex;
	TAXA	*Taxa;
	TREE	*Tree;
	double	Y, Err;

	Tree = &Trees->Tree[0];

/*	printf("{");
	for(NIndex=0;NIndex<Trees->NoOfTaxa;NIndex++)
	{
		Taxa = &Trees->Taxa[NIndex];
		printf("\"%s\",", Taxa->Name);
	}
	printf("}\n");

	printf("{");
	for(NIndex=0;NIndex<Trees->NoOfTaxa;NIndex++)
	{
		Taxa = &Trees->Taxa[NIndex];
		printf("%f,", Taxa->ConData[0]);
	}
	printf("}\n");

	exit(0);
*/
	for(NIndex=0;NIndex<Trees->NoOfTaxa;NIndex++)
	{
		Taxa = &Trees->Taxa[NIndex];
		printf("%f\t", Taxa->ConData[0]);
		

		Y = Model[0] + (Taxa->ConData[0] * Model[1]);
		Err = Error[NIndex]; 
		Taxa->Dependant = Y + Err;
		
		
	}


}

void	SetData(OPTIONS *Opt, TREES *Trees, RATES *Rates, double *Data)
{
	int		NIndex;
	TAXA	*Taxa;
	TREE	*Tree;
	
	Tree = &Trees->Tree[0];

	Rates->Rates[0] = Data[0];
	Rates->Rates[1] = Data[1];	
	
	for(NIndex=2;NIndex<Trees->NoOfTaxa;NIndex++)
	{
		Taxa = &Trees->Taxa[NIndex];
		Taxa->Dependant = Data[NIndex];
	}	
}

void	BayesModeTest(OPTIONS *Opt, TREES *Trees)
{
		RATES	*Rates;
		NUMFILE	*NumFile;
/*		NUMFILE	*Err; */
		int		Index;

		Rates = CreatRates(Opt);

		InitContinusTree(Opt, Trees, 0);

		NumFile = LoadNumFile("./Seq/Res.txt");
/*		Err		= LoadNumFile("./Seq/MNErr.txt");

		NumFile = LoadNumFile("Models.txt");
		Err		= LoadNumFile("cMNErr.txt");
*/
		for(Index=0;Index<NumFile->NoOfLines;Index++)
		{
		/*	GenDataFormModel(Opt, Trees, NumFile->Data[Index], Err->Data[rand() % Err->NoOfLines]); */
			SetData(Opt, Trees, Rates, NumFile->Data[Index]);
			CalcZ(Trees, &Trees->Tree[0], Opt);

			Rates->Lh = Likelihood(Rates, Trees, Opt);
			printf("%d\t%f\n", Index, Rates->Lh);
		}

		exit(0);
}
