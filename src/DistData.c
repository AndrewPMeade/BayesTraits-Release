#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

#include "typedef.h"
#include "genlib.h"
#include "DistData.h"
#include "trees.h"


int	GetLinked(char *Taxa, char *Str)
{
	MakeLower(Str);

	if(strcmp(Str, "linked") == 0)
		return TRUE;

	if(strcmp(Str, "unlinked") == 0)
		return FALSE;

	printf("Cannot convert %s to linked or unliked for taxa %s.\n", Str, Taxa);
	exit(0);
}

double*	PassDataDist(char *Taxa, char *Str, int *NoPoint)
{
	size_t s;
	int Index, Tokes;
	char **Passed;
	double *Ret;

	s = strlen(Str) + 1;

	ReplaceChar(',', ' ', Str);

	Passed = (char**)SMalloc(sizeof(char*) * s);
	
	Tokes = MakeArgv(Str, Passed, (int)s);
	
	*NoPoint = Tokes;
	Ret = (double*)SMalloc(sizeof(double) * Tokes);
	
	for(Index=0;Index<Tokes;Index++)
	{
		if(IsValidDouble(Passed[Index]) == FALSE)
		{
			printf("Cannot convert %s to a valid double for taxa %s.\n", Passed[Index], Taxa);
			exit(0);
		}

		Ret[Index] = atof(Passed[Index]);
	}

	return Ret;
}

DIST_DATA_TAXA*	MakeDistDataTaxa(TREES *Trees, char **Passed, size_t MaxSize)
{
	DIST_DATA_TAXA* Ret;
	int	Index;

	Ret = (DIST_DATA_TAXA*)malloc(sizeof(DIST_DATA_TAXA));
	if(Ret == NULL)
		MallocErr();

	
	Ret->Taxa = GetTaxaFromName(Passed[0], Trees->Taxa, Trees->NoOfTaxa);

	if(Ret->Taxa == NULL)
	{
		printf("Cannot find valid taxa %s", Passed[0]);
		exit(0);
	}

	Ret->Linked = GetLinked(Passed[0], Passed[1]);

	Ret->NoSites = (int*)malloc(sizeof(int) * Trees->NoOfSites);
	Ret->Data = (double**)malloc(sizeof(double*) * Trees->NoOfSites);

	for(Index=0;Index<Trees->NoOfSites;Index++)
		Ret->Data[Index] = PassDataDist(Passed[0], Passed[Index+2], &Ret->NoSites[Index]);

	if(Ret->Linked == TRUE)
	{
		for(Index=1;Index<Trees->NoOfSites;Index++)
			if(Ret->NoSites[0] != Ret->NoSites[Index])
			{
				printf("Dist Data taxa (%s) linked site %d has %d point expecting %d", Passed[0], Index+1, Ret->NoSites[Index], Ret->NoSites[0]);
				exit(0);
			}
	}

	return Ret;
}

void		ProcDistDataFile(TREES *Trees, DIST_DATA *DistData, TEXTFILE *TF)
{
	char *Buffer, **Passed;
	int Index, Tokes;
	DIST_DATA_TAXA	*DTaxa;

	Buffer = (char*)malloc(sizeof(char) * (TF->MaxLine + 1));
	Passed = (char**)malloc(sizeof(char*) * (TF->MaxLine + 1));

	if((Buffer == NULL) || (Passed == NULL))
		MallocErr();

	for(Index=0;Index<TF->NoOfLines;Index++)
	{
		strcpy(Buffer, TF->Data[Index]);
		Tokes = MakeArgv(Buffer, Passed, TF->MaxLine + 1);
		if(Tokes != Trees->NoOfSites + 2 && Tokes > 0)
		{
			printf("Dist Data file (%s) Line %d, expecting %d tokens found %d.\n", DistData->FName, Index, Trees->NoOfSites+2, Tokes);
			exit(0);
		}

		if(Tokes > 0)
		{
			DTaxa = MakeDistDataTaxa(Trees, Passed, TF->MaxLine);
			DistData->DistDataTaxa[DistData->NoTaxa++] = DTaxa;
		}
	}

	free(Buffer);
	free(Passed);
}

DIST_DATA*	LoadDistData(TREES *Trees, char *FName)
{
	DIST_DATA* Ret;
	TEXTFILE *TF;

	Ret = (DIST_DATA*)malloc(sizeof(DIST_DATA));
	if(Ret == NULL)
		MallocErr();

	Ret->NoSites = Trees->NoOfSites;

	Ret->FName = StrMake(FName);

	TF = LoadTextFile(FName, FALSE);
	
	Ret->NoTaxa = 0;
	Ret->DistDataTaxa = (DIST_DATA_TAXA**)SMalloc(sizeof(DIST_DATA_TAXA*) * TF->NoOfLines);

	ProcDistDataFile(Trees, Ret, TF);

//	exit(0);
	
	return Ret;
}

void		FreeDataDistTaxa(DIST_DATA_TAXA *DTaxa, int NoSites)
{
	int Index;

	free(DTaxa->NoSites);

	for(Index=0;Index<NoSites;Index++)
		free(DTaxa->Data[Index]);

	free(DTaxa->Data);

	free(DTaxa);
}

void		FreeDistData(DIST_DATA *DistData)
{
	int Index;	

	for(Index=0;Index<DistData->NoTaxa;Index++)
		FreeDataDistTaxa(DistData->DistDataTaxa[Index], DistData->NoSites);
	
	free(DistData->DistDataTaxa);
	free(DistData->FName);
}

void		PrintDistDataTaxa(FILE *Out, DIST_DATA_TAXA *TData, int NoSites)
{
	int Index;

	fprintf(Out, "             %s ", TData->Taxa->Name);

	if(TData->Linked == TRUE)
	{
		fprintf(Out, "%d Linked sites\n", TData->NoSites[0]);
		return;
	}

	for(Index=0;Index<NoSites-1;Index++)
		fprintf(Out, "%d,", TData->NoSites[Index]);
	fprintf(Out, "%d ", TData->NoSites[Index]);

	fprintf(Out, "Unlinked sites\n");
}

void		PrintDistData(FILE *Out, DIST_DATA *DistData)
{
	int Index;

	fprintf(Out, "Distribution Data:\n");
	fprintf(Out, "             No Taxa %d\n", DistData->NoTaxa);

	for(Index=0;Index<DistData->NoTaxa;Index++)
		PrintDistDataTaxa(Out, DistData->DistDataTaxa[Index], DistData->NoSites);
}

void SetRandDistDataRates(DIST_DATE_RATES *DistRates, DIST_DATA* DistData, RANDSTATES *RS)
{
	int TIndex, SIndex;
	DIST_DATA_TAXA *DistT;

	for(TIndex=0;TIndex<DistData->NoTaxa;TIndex++)
	{
		DistT = DistData->DistDataTaxa[TIndex];

		for(SIndex=0;SIndex<DistData->NoSites;SIndex++)
		{
			if(DistT->Linked == TRUE && SIndex != 0)
				DistRates->SiteMap[TIndex][SIndex] = DistRates->SiteMap[TIndex][0];
			else
				DistRates->SiteMap[TIndex][SIndex] = RandUSInt(RS) % DistT->NoSites[SIndex];
		}
	}
}

DIST_DATE_RATES*	CreateDistDataRates(DIST_DATA* DistData, RANDSTATES *RS)
{
	DIST_DATE_RATES* Ret;
	int Index;
	
	Ret = (DIST_DATE_RATES*)SMalloc(sizeof(DIST_DATE_RATES));
	
	Ret->NoSites = DistData->NoSites;
	Ret->NoTaxa = DistData->NoTaxa;

	Ret->SiteMap = (int**)SMalloc(sizeof(int*) * DistData->NoTaxa);
	
	for(Index=0;Index<DistData->NoTaxa;Index++)
		Ret->SiteMap[Index] = (int*)SMalloc(sizeof(int) * DistData->NoSites);
	
	SetRandDistDataRates(Ret, DistData, RS);

	return Ret;
}

void	FreeDistDataRates(DIST_DATE_RATES* DistRates)
{
	int Index;

	for(Index=0;Index<DistRates->NoTaxa;Index++)
		free(DistRates->SiteMap[Index]);

	free(DistRates);
}