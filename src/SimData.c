#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "genlib.h"
#include "typedef.h"
#include "trees.h"
#include "RandLib.h"
#include "likelihood.h"

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
	int Index, Len, Pos;
		
	Len = strlen(SymList);
	Out = OpenRandSimOutFile(NOS);

	for(Index=0;Index<Trees->NoOfTaxa;Index++)
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


int		GetRootRates(TREES *Trees, RANDSTATES *RS)
{
	int Ret;

	Ret = (int)(RandUSInt(RS) % Trees->NoOfStates);

	return Ret;
}

int		EndState(int StartS, MATRIX *P, RANDSTATES *RS)
{
	double Point, Sum;
	int Index, NOS;

	NOS = P->NoOfCols;

	Point = RandDouble(RS);
	Sum = 0;

	for(Index=0;Index<NOS-1;Index++)
	{
		if((Point > Sum) && (Point <= Sum + P->me[StartS][Index]))
			return Index;

		Sum += P->me[StartS][Index];
	}

	return Index;
}

void	RecSimData(int StartS, NODE N, TREES *Trees, RANDSTATES *RS, FILE *Out)
{
	MATRIX *P;
	int Index, EndS;

	P = Trees->PList[N->ID];
	EndS = EndState(StartS, P, RS);

	if(N->Tip == TRUE)
	{
		fprintf(Out, "%s\t%c\n", N->Taxa->Name, Trees->SymbolList[EndS]);
		return;
	}

	for(Index=0;Index<N->NoNodes;Index++)
		RecSimData(EndS, N->NodeList[Index], Trees, RS, Out);
}

FILE*	OpenSimOutFile(char *BaseFN)
{
	FILE *Ret;
	char *Buffer;

	Buffer = (char*)malloc(sizeof(char) * (strlen(BaseFN)+64));
	if(Buffer == NULL)
		MallocErr();

	sprintf(Buffer, "%s.sim.txt", BaseFN);

	Ret = OpenWrite(Buffer);

	free(Buffer);

	return Ret;
}

void	SimData(OPTIONS *Opt, TREES *Trees, RATES *Rates)
{
	double	Lh;
	int		Err;
	RANDSTATES *RS;
	int		Index, RootS;
	NODE	Root;
	FILE	*OutF;
	return;
//	BuildAllSateDS(Trees);

	OutF = OpenSimOutFile(Opt->LogFN);

	
	RS = CreateSeededRandStates(Opt->Seed);

//	Rates->Rates[0] = 0.000001;
//	Rates->Rates[0] = 1.503441;

//	for(Index=0;Index<Rates->NoOfRates;Index++)
//		Rates->Rates[Index] = RandDouble(RS) * 10;

	Rates->TreeNo = 0;

	Lh = Likelihood(Rates, Trees, Opt);
	printf("SimTreeLh:\t%f\n", Lh);

	Err = SetUpAMatrix(Rates, Trees, Opt);

	if(Err > 0)
	{
		printf("A matrix err\n");
		exit(0);
	}

	Err = SetAllPMatrix(Rates, Trees, Opt, 1.0, 1.0);
	if(Err == TRUE)
	{
		printf("P matrix err\n");
		exit(0);
	}
	
	RootS = GetRootRates(Trees, RS);

	Root = Trees->Tree[0]->Root;

	printf("Root State = %d\t%c\n", RootS, Trees->SymbolList[RootS]);

	for(Index=0;Index<Root->NoNodes;Index++)
	{
		RecSimData(RootS, Root->NodeList[Index], Trees, RS, OutF);
	}

	fclose(OutF);

	exit(0);
}