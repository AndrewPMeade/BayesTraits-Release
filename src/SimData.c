/*
*  BayesTriats 4.0
*
*  copyright 2022
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
#include <string.h>
#include <math.h>

#include "GenLib.h"
#include "TypeDef.h"
#include "Trees.h"
#include "RandLib.h"
#include "Likelihood.h"

FILE*	OpenRandSimOutFile(int NOS)
{
	FILE *Ret;
	char *Buffer;

	Buffer = (char*)malloc(sizeof(char) * 1024);
	if(Buffer == NULL)
		MallocErr();

	sprintf(Buffer, "RandNOS-%03d.txt", NOS);

	Ret = OpenWrite(Buffer);

	free(Buffer);

	return Ret;
}
	

void	BuildStateDS(int NOS, TREES *Trees)
{
	char *SymList = "0123456789abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ";
	FILE *Out;
	int Index, Len;
		
	Len = (int)strlen(SymList);
	Out = OpenRandSimOutFile(NOS);

	for(Index=0;Index<Trees->NoTaxa;Index++)
	{
		fprintf(Out, "%s\t%c\n", Trees->Taxa[Index]->Name, SymList[Index%NOS]);
	}

	fclose(Out);
}

void	BuildAllSateDS(TREES *Trees)
{
	int Index;

	for(Index=2;Index<63;Index++)
		BuildStateDS(Index, Trees);

	exit(0);
}


int		GetRootRates(TREES *Trees, gsl_rng	*RNG)
{
	return gsl_rng_uniform_int(RNG, Trees->NoStates);
}

int		EndState(int StartS, MATRIX *P, gsl_rng	*RNG)
{
	double Point, Sum;
	int Index, NOS;

	NOS = P->NoOfCols;

	Point = gsl_rng_uniform(RNG); 
	Sum = 0;

	for(Index=0;Index<NOS-1;Index++)
	{
		if((Point > Sum) && (Point <= Sum + P->me[StartS][Index]))
			return Index;

		Sum += P->me[StartS][Index];
	}

	return Index;
}

void	RecSimData(int StartS, NODE N, TREES *Trees, gsl_rng *RNG, FILE *Out, int *NoChange)
{
	MATRIX *P;
	int Index, EndS;

	P = Trees->PList[N->ID];
	EndS = EndState(StartS, P, RNG);

	if(EndS != StartS)
		(*NoChange)++;

	if(N->Tip == TRUE)
	{
		fprintf(Out, "%s\t%c\n", N->Taxa->Name, Trees->SymbolList[EndS]);
		return;
	}

	for(Index=0;Index<N->NoNodes;Index++)
		RecSimData(EndS, N->NodeList[Index], Trees, RNG, Out, NoChange);
}

FILE*	OpenSimOutFile(char *BaseFN)
{
	FILE *Ret;
	char *Buffer;

	Buffer = (char*)SMalloc(sizeof(char) * (strlen(BaseFN)+64));

	sprintf(Buffer, "%s.Sim.txt", BaseFN);

	Ret = OpenWrite(Buffer);

	free(Buffer);

	return Ret;
}

int PrintNoChanges(NODE *NList, int NoNodes)
{
	int Index, Ret;
	NODE N;

	Ret = 0;

	for(Index=0;Index<NoNodes;Index++)
	{
		N = NList[Index];

		if(N->Ans != NULL)
		{
		
		}
	}

	return Ret;
}
// ./Seq/CommunalSongSubSet.trees ./Seq/Sim.data.txt < in.txt > sout.txt 
// ./Seq/CommunalSongSubSet.trees ./Seq/MattingSystem.txt < in.txt > sout.txt
// ./Seq/BirdTreeSubSet.trees ./Seq/MediumCommunalTerritory.txt < in.txt > sout.txt

// ./Seq/Tobias/MatingSystem.trees ./Seq/Tobias/MatingSystem.txt < ./Seq/Tobias/in.txt > ./Seq/Tobias/sout.txt

void	SimData(OPTIONS *Opt, TREES *Trees, RATES *Rates)
{
	double	Lh;
	int		Err;
	gsl_rng	*RNG;
	
	int		Index, RootS, SNo, NoChanges, No0;
	NODE	Root;
	FILE	*OutF;

	return;
//	BuildAllSateDS(Trees);

	OutF = OpenWriteWithExt(Opt->BaseOutputFN, OUTPUT_EXT_SIM);
	
//	RS = CreateSeededRandStates(Opt->Seed);


	RNG	= gsl_rng_alloc(RNG_TYPE);
	gsl_rng_set(RNG, Opt->Seed);


	No0 = 0;
	for(SNo=0;SNo<1;SNo++)
	{
		Rates->Rates[0] = gsl_rng_uniform_pos(RNG) * 0.0001;
	//	Rates->Rates[0] = 0.00002169;
		Rates->Rates[0] = 0.00001;
		Rates->Rates[0] = 0.000009;

		Rates->Rates[0] = 1.0;
		Rates->Rates[0] = 0.01;

	//	for(Index=0;Index<Rates->NoOfRates;Index++)
	//		Rates->Rates[Index] = 1.0;
		Rates->TreeNo = 0;

		Lh = Likelihood(Rates, Trees, Opt);
	//	printf("SimTreeLh:\t%f\n", Lh);

		Err = SetUpAllAMatrix(Rates, Trees, Opt);

		if(Err > 0)
		{
			printf("A matrix err\n");
			exit(0);
		}

		Err = SetAllPMatrix(Rates, Trees, Opt, 1.0);
		if(Err == TRUE)
		{
			printf("P matrix err\n");
			exit(0);
		}
	
		RootS = GetRootRates(Trees, RNG);
		RootS = 0;

		Root = Trees->Tree[0]->Root;

	//	printf("Root State = %d\t%c\n", RootS, Trees->SymbolList[RootS]);
		

		NoChanges = 0;
		for(Index=0;Index<Root->NoNodes;Index++)
		{
			RecSimData(RootS, Root->NodeList[Index], Trees, RNG, OutF, &NoChanges);
		}

		printf("Sim:\t%d\t%.012f\t%d\n", SNo, Rates->Rates[0], NoChanges);
		fflush(stdout);

		if(NoChanges == 0)
			No0++;
	}

	printf("No0\t%d\n", No0);

	fclose(OutF);

	gsl_rng_free(RNG);

	exit(0);
}