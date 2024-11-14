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
#include <math.h>
#include <string.h>

#include "TypeDef.h"
#include "GenLib.h"
#include "DistData.h"
#include "Trees.h"
#include "Geo.h"

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

	free(Passed);

	return Ret;
}

void ProcGeoDistDataTaxa(DIST_DATA_TAXA* DDTaxa, TREES *Trees, char **Passed)
{
	double *Lat, *Long, X, Y, Z;
	int Index, NoLat, NoLong, NoData;
			
	Long = PassDataDist(DDTaxa->Taxa->Name, Passed[2], &NoLong);
	Lat = PassDataDist(DDTaxa->Taxa->Name, Passed[3], &NoLat);

	if(NoLong != NoLat)
	{
		printf("Data distribution , number of long and lat must be the same for taxa %s.\n", DDTaxa->Taxa->Name);
		exit(1);
	}

	NoData = NoLong;

	for(Index=0;Index<Trees->NoSites;Index++)
		DDTaxa->NoSites[Index] = NoData;
	
	for(Index=0;Index<Trees->NoSites;Index++)
		DDTaxa->Data[Index] = (double*)SMalloc(sizeof(double) * NoData);

	for(Index=0;Index<NoData;Index++)
	{
		LongLatToXYZ(Long[Index], Lat[Index], &X, &Y, &Z);
		DDTaxa->Data[0][Index] = X;
		DDTaxa->Data[1][Index] = Y;
		DDTaxa->Data[2][Index] = Z;
	}

	free(Long);
	free(Lat);
}

DIST_DATA_TAXA*	MakeDistDataTaxa(OPTIONS *Opt, TREES *Trees, char **Passed, size_t MaxSize)
{
	DIST_DATA_TAXA* Ret;
	int	Index;

	Ret = (DIST_DATA_TAXA*)SMalloc(sizeof(DIST_DATA_TAXA));
	
	Ret->Taxa = GetTaxaFromName(Passed[0], Trees->Taxa, Trees->NoTaxa);

	if(Ret->Taxa == NULL)
	{
		printf("Cannot find valid taxa %s", Passed[0]);
		exit(0);
	}

	Ret->Linked = GetLinked(Passed[0], Passed[1]);

	if(Ret->Linked == FALSE && Opt->Model == M_GEO)
	{
		printf("Distribution data for geo models must be linked.");
		exit(1);	
	}

	Ret->NoSites = (int*)SMalloc(sizeof(int) * Trees->NoSites);
	Ret->Data = (double**)SMalloc(sizeof(double*) * Trees->NoSites);

	if(Opt->Model == M_GEO)
	{
		ProcGeoDistDataTaxa(Ret, Trees, Passed);
	}
	else
	{
		for(Index=0;Index<Trees->NoSites;Index++)
			Ret->Data[Index] = PassDataDist(Passed[0], Passed[Index+2], &Ret->NoSites[Index]);
	}

	if(Ret->Linked == TRUE)
	{
		for(Index=1;Index<Trees->NoSites;Index++)
			if(Ret->NoSites[0] != Ret->NoSites[Index])
			{
				printf("Dist Data taxa (%s) linked site %d has %d point expecting %d", Passed[0], Index+1, Ret->NoSites[Index], Ret->NoSites[0]);
				exit(1);
			}
	}

	return Ret;
}

int			ValidDistDataLine(OPTIONS *Opt, TREES *Trees, int Tokes)
{
	if(Tokes == 0)
		return TRUE;

	if(Opt->Model == M_GEO)
	{
		if(Tokes != Trees->NoSites + 1)
			return FALSE;
		else
			return TRUE;
	}

	if(Tokes != Trees->NoSites + 2)
		return FALSE;

	return TRUE;
}

void		ProcDistDataFile(OPTIONS *Opt, TREES *Trees, DIST_DATA *DistData, TEXTFILE *TF)
{
	char *Buffer, **Passed;
	int Index, Tokes;
	DIST_DATA_TAXA	*DTaxa;

	Buffer = (char*)SMalloc(sizeof(char) * (TF->MaxLine + 1));
	Passed = (char**)SMalloc(sizeof(char*) * (TF->MaxLine + 1));
	
	for(Index=0;Index<TF->NoOfLines;Index++)
	{
		strcpy(Buffer, TF->Data[Index]);
		Tokes = MakeArgv(Buffer, Passed, TF->MaxLine + 1);
		if(ValidDistDataLine(Opt, Trees, Tokes) == FALSE)
		{
			printf("Dist Data file (%s) Line %d, expecting %d tokens found %d.\n", DistData->FName, Index, Trees->NoSites+2, Tokes);
			exit(0);
		}

		if(Tokes > 0)
		{
			DTaxa = MakeDistDataTaxa(Opt, Trees, Passed, TF->MaxLine);
			DistData->DistDataTaxa[DistData->NoTaxa++] = DTaxa;
		}
	}

	free(Buffer);
	free(Passed);
}

DIST_DATA*	LoadDistData(OPTIONS *Opt, TREES *Trees, char *FName)
{
	DIST_DATA* Ret;
	TEXTFILE *TF;

	Ret = (DIST_DATA*)SMalloc(sizeof(DIST_DATA));

	Ret->NoSites = Trees->NoSites;

	Ret->FName = StrMake(FName);

	TF = LoadTextFile(FName, FALSE);
	
	Ret->NoTaxa = 0;
	Ret->DistDataTaxa = (DIST_DATA_TAXA**)SMalloc(sizeof(DIST_DATA_TAXA*) * TF->NoOfLines);

	ProcDistDataFile(Opt, Trees, Ret, TF);

	FreeTextFile(TF);

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
	free(DistData);
}

void		PrintDistDataTaxa(FILE *Out, DIST_DATA_TAXA *TData, int NoSites)
{
	int Index;

	fprintf(Out, "             %s ", TData->Taxa->Name);

	if(TData->Linked == TRUE)
	{
		fprintf(Out, "%d Linked samples\n", TData->NoSites[0]);
		return;
	}

	for(Index=0;Index<NoSites-1;Index++)
		fprintf(Out, "%d,", TData->NoSites[Index]);
	fprintf(Out, "%d ", TData->NoSites[Index]);

	fprintf(Out, "Unlinked samples\n");
}

void		PrintDistData(FILE *Out, DIST_DATA *DistData)
{
	int Index;

	fprintf(Out, "Distribution Data:\n");
	fprintf(Out, "             No Taxa %d\n", DistData->NoTaxa);

	for(Index=0;Index<DistData->NoTaxa;Index++)
		PrintDistDataTaxa(Out, DistData->DistDataTaxa[Index], DistData->NoSites);
}

void SetRandDistDataRates(DIST_DATE_RATES *DistRates, DIST_DATA* DistData, gsl_rng *RNG)
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
				DistRates->SiteMap[TIndex][SIndex] = (int)gsl_rng_uniform_int(RNG, DistT->NoSites[SIndex]);
		}
	}
}

int**	AllocDistDataRateMem(int NoTaxa, int NoSites)
{
	int **Ret;
	int Index;

	Ret = (int**)SMalloc(sizeof(int*) * NoTaxa);
	Ret[0] = (int*)SMalloc(sizeof(int) * NoTaxa * NoSites);

	for(Index=1;Index<NoTaxa;Index++)
		Ret[Index] = Ret[0] + Index * NoSites;

	return Ret;
}

DIST_DATE_RATES*	CreateDistDataRates(DIST_DATA* DistData, gsl_rng *RNG)
{
	DIST_DATE_RATES* Ret;
	
	Ret = (DIST_DATE_RATES*)SMalloc(sizeof(DIST_DATE_RATES));
	
	Ret->NoSites = DistData->NoSites;
	Ret->NoTaxa = DistData->NoTaxa;
	Ret->SiteMap = AllocDistDataRateMem(DistData->NoTaxa, DistData->NoSites);
	
	SetRandDistDataRates(Ret, DistData, RNG);

	return Ret;
}

void	FreeDistDataRates(DIST_DATE_RATES* DistRates)
{
	free(DistRates->SiteMap[0]);
	free(DistRates->SiteMap);

	free(DistRates);
}

void	CopyDistDataRates(DIST_DATE_RATES* A, DIST_DATE_RATES* B)
{
	memcpy(A->SiteMap[0], B->SiteMap[0], sizeof(int) * A->NoSites * A->NoTaxa);
}

void	SetDistDataTaxaCon(NODE N, DIST_DATA_TAXA *DTaxa, int NoSites, int *PosList)
{
	int SIndex;

	CONTRAST *Con;

	Con = N->ConData->Contrast[0];

	for(SIndex=0;SIndex<NoSites;SIndex++)
		Con->Data[SIndex] = DTaxa->Data[SIndex][PosList[SIndex]];
}

void	SetDistDataTaxaGeo(NODE N, DIST_DATA_TAXA *DTaxa, int NoSites, int *PosList)
{
	int SIndex;
	FATTAILNODE *FTN;

	FTN = N->FatTailNode;
	
	for(SIndex=0;SIndex<NoSites;SIndex++)
		FTN->Ans[SIndex] = DTaxa->Data[SIndex][PosList[SIndex]];
}

void	SetDistDataTaxa(DIST_DATA *DD, DIST_DATA_TAXA *DTaxa, TREE *Tree, OPTIONS *Opt, int *PosList)
{
	NODE N;

	N = GetTreeTaxaNode(Tree, DTaxa->Taxa->No);

	if(Opt->ModelType == MT_CONTRAST)
		SetDistDataTaxaCon(N, DTaxa, DD->NoSites, PosList);

	if(Opt->ModelType == MT_FATTAIL)	
		SetDistDataTaxaGeo(N, DTaxa, DD->NoSites, PosList);
}


void	SetTreeDistData(RATES *Rates, OPTIONS *Opt, TREES *Trees)
{
	int TIndex;
	DIST_DATE_RATES	*DDRates;
	DIST_DATA		*DD;
	TREE			*Tree;
	DIST_DATA_TAXA	*DDTaxa;

	DD = Opt->DistData;
	DDRates = Rates->DistDataRates;
	
	Tree = Trees->Tree[Rates->TreeNo];

	for(TIndex=0;TIndex<DD->NoTaxa;TIndex++)
	{
		DDTaxa = DD->DistDataTaxa[TIndex];
		SetDistDataTaxa(DD, DDTaxa, Tree, Opt, DDRates->SiteMap[TIndex]);
	}
}

void	ChangeTreeDistData(OPTIONS *Opt, RATES *Rates)
{
	int Index, TaxaNo, NPos;
	DIST_DATA		*DD;
	DIST_DATE_RATES	*DDRates;
	DIST_DATA_TAXA	*DDTaxa;

	DD = Opt->DistData;
	DDRates = Rates->DistDataRates;

	TaxaNo = (int)gsl_rng_uniform_int(Rates->RNG, DDRates->NoTaxa);

	DDTaxa = DD->DistDataTaxa[TaxaNo];
	
	if(DDTaxa->Linked == TRUE)
	{
		NPos = (int)gsl_rng_uniform_int(Rates->RNG, DDTaxa->NoSites[0]);
		for(Index=0;Index<DD->NoSites;Index++)
			DDRates->SiteMap[TaxaNo][Index] = NPos;
	}
	else
	{
		Index = (int)gsl_rng_uniform_int(Rates->RNG, DD->NoSites);
		DDRates->SiteMap[TaxaNo][Index] = (int)gsl_rng_uniform_int(Rates->RNG, DDTaxa->NoSites[Index]);
	}
}

void OutputDataDistHeadder(FILE *Str, OPTIONS *Opt)
{
	int				TIndex, SIndex;
	DIST_DATA		*DData;
	DIST_DATA_TAXA	*DDTaxa;

	DData = Opt->DistData;

	for(TIndex=0;TIndex<DData->NoTaxa;TIndex++)
	{
		DDTaxa = DData->DistDataTaxa[TIndex];
		if(Opt->Model == M_GEO)
			fprintf(Str, "%s Long\t%s Lat\t", DDTaxa->Taxa->Name, DDTaxa->Taxa->Name);
		else
			for(SIndex=0;SIndex<DData->NoSites;SIndex++)	
				fprintf(Str, "%s-%d\t", DDTaxa->Taxa->Name, SIndex+1);
	}
}

void OutputDataDistGeo(FILE *Str, RATES *Rates, OPTIONS *Opt, int *Map, DIST_DATA_TAXA *DDTaxa)	
{
	double Long, Lat, X, Y, Z;
		
	X = DDTaxa->Data[0][Map[0]];
	Y = DDTaxa->Data[1][Map[1]];
	Z = DDTaxa->Data[2][Map[2]];

	XYZToLongLat(X, Y, Z, &Long, &Lat);

	fprintf(Str, "%f\t%f\t", Long, Lat);
}

void OutputDataDist(FILE *Str, RATES *Rates, OPTIONS *Opt)
{
	int				Pos, TIndex, SIndex;
	double			Val;
	DIST_DATA		*DData;
	DIST_DATE_RATES	*DDRates;
	DIST_DATA_TAXA	*DDTaxa;

	DData = Opt->DistData;
	DDRates = Rates->DistDataRates;

	for(TIndex=0;TIndex<DData->NoTaxa;TIndex++)
	{
		DDTaxa = DData->DistDataTaxa[TIndex];
		if(Opt->Model == M_GEO)
			OutputDataDistGeo(Str, Rates, Opt, DDRates->SiteMap[TIndex], DDTaxa);
		else
		{
			for(SIndex=0;SIndex<DData->NoSites;SIndex++)	
			{
			
				Pos = DDRates->SiteMap[TIndex][SIndex];
				Val = DData->DistDataTaxa[TIndex]->Data[SIndex][Pos];

				fprintf(Str, "%f\t", Val);
			}
		}
	}
}


