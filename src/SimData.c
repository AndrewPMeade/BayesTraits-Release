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


int		GetRootRates(TREES *Trees, RANDSTATES *RS)
{
	int Ret;

	Ret = (int)(RandUSInt(RS) % Trees->NoStates);

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

void	RecSimData(int StartS, NODE N, TREES *Trees, RANDSTATES *RS, FILE *Out, int *NoChange)
{
	MATRIX *P;
	int Index, EndS;

	P = Trees->PList[N->ID];
	EndS = EndState(StartS, P, RS);

	if(EndS != StartS)
		(*NoChange)++;

	if(N->Tip == TRUE)
	{
		fprintf(Out, "%s\t%c\n", N->Taxa->Name, Trees->SymbolList[EndS]);
		return;
	}

	for(Index=0;Index<N->NoNodes;Index++)
		RecSimData(EndS, N->NodeList[Index], Trees, RS, Out, NoChange);
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
	RANDSTATES *RS;
	int		Index, RootS, SNo, NoChanges, No0;
	NODE	Root;
	FILE	*OutF;

	return;
//	BuildAllSateDS(Trees);

	OutF = OpenWriteWithExt(Opt->BaseOutputFN, OUTPUT_EXT_SIM);
	
	RS = CreateSeededRandStates(Opt->Seed);

	No0 = 0;
	for(SNo=0;SNo<100000;SNo++)
	{
		Rates->Rates[0] = RandDouble(RS) * 0.0001;
	//	Rates->Rates[0] = 0.00002169;
		Rates->Rates[0] = 0.00001;
		Rates->Rates[0] = 0.000009;

	//	Rates->Rates[1] = 1;

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
	
		RootS = GetRootRates(Trees, RS);
		RootS = 0;

		Root = Trees->Tree[0]->Root;

	//	printf("Root State = %d\t%c\n", RootS, Trees->SymbolList[RootS]);
		

		NoChanges = 0;
		for(Index=0;Index<Root->NoNodes;Index++)
		{
			RecSimData(RootS, Root->NodeList[Index], Trees, RS, OutF, &NoChanges);
		}

		printf("Sim:\t%d\t%.012f\t%d\n", SNo, Rates->Rates[0], NoChanges);
		fflush(stdout);

		if(NoChanges == 0)
			No0++;
	}

	printf("No0\t%d\n", No0);

	fclose(OutF);

	exit(0);
}