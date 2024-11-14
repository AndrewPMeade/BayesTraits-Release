#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

#include "TypeDef.h"
#include "GenLib.h"
#include "RestrictionMap.h"
#include "BinaryCompressedResMaps.h"

//4499328773

BIN_MAP_INFO*	CreateBinMapInfo(char *FName, double AgeMin, double AgeMax)
{
	BIN_MAP_INFO* Ret;

	Ret = (BIN_MAP_INFO*)SMalloc(sizeof(BIN_MAP_INFO));

	Ret->FileName = StrMake(FName);
	Ret->AgeMin = AgeMin;
	Ret->AgeMax = AgeMax;

	return Ret;
}

void	FreeBinMapInfo(BIN_MAP_INFO *BinMapInfo)
{
	free(BinMapInfo->FileName);
	free(BinMapInfo);
}

void	FreeBinMapInfoList(BIN_MAP_INFO **BinMapInfo, int Size)
{
	int Index;

	for(Index=0;Index<Size;Index++)
		FreeBinMapInfo(BinMapInfo[Index]);

	free(BinMapInfo);
}

int BinMapComp(const void *AP, const void *BP)
{
	BIN_MAP_INFO *A, *B;

	A = *(BIN_MAP_INFO**)AP;
	B = *(BIN_MAP_INFO**)BP;

	if(A->AgeMin < B->AgeMax)
		return -1;

	if(A->AgeMin == B->AgeMin)
		return 0;

	return 1;
}

BIN_MAP_INFO**	LoadBinMapInfo(TEXTFILE* TF, int *Size)
{
	
	BIN_MAP_INFO** Ret;
	int Index, Tokes;
	char **Passed;
	
	

	*Size = 0;
	Ret = (BIN_MAP_INFO**)SMalloc(sizeof(BIN_MAP_INFO*) * TF->NoOfLines);
	Passed = (char**)SMalloc(sizeof(char*)* TF->MaxLine);

	for(Index=1;Index<TF->NoOfLines;Index++)
	{
		Tokes = MakeArgvChar(TF->Data[Index], Passed, TF->MaxLine, '\t');
		if(Tokes == 3)
		{
			Ret[*Size] = CreateBinMapInfo(Passed[0], atof(Passed[1]), atof(Passed[2]));
			(*Size)++;
		}
	}
	
	free(Passed);

	qsort(Ret, *Size, sizeof(BIN_MAP_INFO*), BinMapComp);

	return Ret;
}

void	GetMapGroupInfo(MAP_GROUP* MapGroup, char* Line)
{
	char **Passed;
	char *Copy;
	int Tokes;

	Copy = StrMake(Line);
	Passed = (char**)SMalloc(sizeof(char*) * (strlen(Copy)+1));
		
	Tokes = MakeArgv(Copy, Passed, (int)strlen(Copy));

	if(Tokes != 3)
	{
		printf("%s\nshould have three tokens, separated by tabs, a filename, number of longitude boxes, number of latitude boxes.\n", Copy);
		exit(1);
	}

	MapGroup->FileName = StrMake(Passed[0]);
	MapGroup->NoLongBoxes = atoi(Passed[1]);
	MapGroup->NoLatBoxes = atoi(Passed[2]);

	if(MapGroup->NoLongBoxes < 2 || MapGroup->NoLatBoxes < 2)
	{
		printf("Map File Error, number of boxes must be greater than 1.\n");
		printf("error with paramters %s %s when building map files.\n", Passed[1], Passed[2]);
		exit(1);
	}

	free(Copy);
	free(Passed);
}

MAP_GROUP* LoadMapGroup(char* FName)
{
	MAP_GROUP* MapGroup;
	TEXTFILE* TF;

	MapGroup = (MAP_GROUP*)SMalloc(sizeof(MAP_GROUP));

	TF = LoadTextFile(FName, FALSE);

	GetMapGroupInfo(MapGroup, TF->Data[0]);

	MapGroup->BinMaps = LoadBinMapInfo(TF, &MapGroup->NoMaps);

	FreeTextFile(TF);

	return MapGroup;
}

void FreeMapGroup(MAP_GROUP* MapGroup)
{
	free(MapGroup->FileName);
	FreeBinMapInfoList(MapGroup->BinMaps, MapGroup->NoMaps);
	free(MapGroup);
}

void	PrintBinMapInfo(char *FName, BIN_MAP_INFO** BinMapInfo, int Size)
{
	int Index;

	printf("Input (%s) contains.\n", FName);
	for(Index=0;Index<Size;Index++)
		printf("%s\t%f\t%f\n", BinMapInfo[Index]->FileName, BinMapInfo[Index]->AgeMin, BinMapInfo[Index]->AgeMax);

	fflush(stdout);
}

void	WriteResMapInfo(FILE *File, RESTRICTION_MAP *ResMap)
{
	size_t FNSize;

	fwrite(&ResMap->AgeMin, sizeof(double), 1, File);
	fwrite(&ResMap->AgeMax, sizeof(double), 1, File);
	FNSize = strlen(ResMap->FileName);

	fwrite(&FNSize, sizeof(size_t), 1, File);
	fwrite(ResMap->FileName, sizeof(char), FNSize, File);
}

void	WriteResMapData(FILE *File, RESTRICTION_POINT *GeoPointListConMem, size_t NoResPoints)
{
	double *CMem;
	int *ICMem;
	size_t Index;

	fwrite(&NoResPoints, sizeof(size_t), 1, File);
//	Writing the strucutre has padding, non machen indepdent and intel inspector does not like it. 
//	fwrite(GeoPointListConMem, sizeof(RESTRICTION_POINT), NoResPoints, File);

	CMem = (double*)SMalloc(sizeof(double) * NoResPoints);

	for(Index=0;Index<NoResPoints;Index++)
		CMem[Index] = GeoPointListConMem[Index].Long;
	fwrite(CMem, sizeof(double), NoResPoints, File);

	for(Index=0;Index<NoResPoints;Index++)
		CMem[Index] = GeoPointListConMem[Index].Lat;
	fwrite(CMem, sizeof(double), NoResPoints, File);

	ICMem = (int*)CMem;
	for(Index=0;Index<NoResPoints;Index++)
		ICMem[Index] = GeoPointListConMem[Index].Value;
	fwrite(CMem, sizeof(int), NoResPoints, File);

	free(CMem);
}

void	FindNoDiffPoints(RESTRICTION_MAP *ResMapA, RESTRICTION_MAP *ResMapB, size_t *NoDiff, RESTRICTION_POINT *SavePoints)
{
	size_t Index;

	for(Index=0;Index<ResMapA->NoResPoint;Index++)
	{
		if(ResMapA->GeoPointListConMem[Index].Value != ResMapB->GeoPointListConMem[Index].Value)
		{
			
			if(SavePoints != NULL)
			{
				SavePoints[*NoDiff].Lat = ResMapB->GeoPointListConMem[Index].Lat;
				SavePoints[*NoDiff].Long = ResMapB->GeoPointListConMem[Index].Long;
				SavePoints[*NoDiff].Value = ResMapB->GeoPointListConMem[Index].Value;
			}

			(*NoDiff)++;
		}
	}
}


RESTRICTION_POINT*	GetDiffPoints(RESTRICTION_MAP *ResMapA, RESTRICTION_MAP *ResMapB, size_t *NoDiff)
{
	RESTRICTION_POINT* Ret;

	*NoDiff = 0;
	
	FindNoDiffPoints(ResMapA, ResMapB, NoDiff, NULL);
	Ret = (RESTRICTION_POINT*)SMalloc(sizeof(RESTRICTION_POINT) * *NoDiff);

	*NoDiff = 0;
	FindNoDiffPoints(ResMapA, ResMapB, NoDiff, Ret);
		
	return Ret;
}

RESTRICTION_MAP**	LoadBinCompResMaps(char *FName, int *Size);

void	PrintBuildBinResMapsHeader(void)
{
	printf("\n\n\n");
	printf("BuildBinaryMaps (%s) - Creates a binary map file from csv files.\n", __DATE__);
	printf("The fist line, should be the output file name, number of longitude bins, the number of latitude bins.\n");
	printf("The other lines should be, csv file defining the map, start heigh (from the root), end high (from the root).\n");
	printf("See manual for details.\n");
	printf("Building the maps file can take a while and require considerable memory, depending on the number of maps and number of points per map.\n");
	printf("The program will exit once the maps have been build.\n");
	printf("\n\n\n");
	fflush(stdout);
}

void	BuildBinResMaps(char *MapFileName)
{
	FILE *File;
	MAP_GROUP*	MapGroup;
	BIN_MAP_INFO**	BinMapInfo;
	int Index;
	RESTRICTION_MAP *ResMapA, *ResMapB;
	size_t NoDiff;
	RESTRICTION_POINT* DiffPoints;

	PrintBuildBinResMapsHeader();

	MapGroup = LoadMapGroup(MapFileName);
	BinMapInfo = MapGroup->BinMaps;

	PrintBinMapInfo(MapFileName, BinMapInfo, MapGroup->NoMaps);

	File = fopen(MapGroup->FileName, "wb");

	fwrite(&MapGroup->NoMaps, sizeof(Index), 1, File);
	fwrite(&MapGroup->NoLongBoxes, sizeof(Index), 1, File);
	fwrite(&MapGroup->NoLatBoxes, sizeof(Index), 1, File);


	ResMapA = LoadResMap(BinMapInfo[0]->FileName, BinMapInfo[0]->AgeMin, BinMapInfo[0]->AgeMax);
	WriteResMapInfo(File, ResMapA);
	WriteResMapData(File, ResMapA->GeoPointListConMem, ResMapA->NoResPoint);
	printf("Processed Map %s\n", ResMapA->FileName);
	ResMapB = NULL;
	for(Index=0;Index<MapGroup->NoMaps-1;Index++)
	{
		ResMapB = LoadResMap(BinMapInfo[Index+1]->FileName, BinMapInfo[Index+1]->AgeMin, BinMapInfo[Index+1]->AgeMax);
				
		WriteResMapInfo(File, ResMapB);
		DiffPoints = GetDiffPoints(ResMapA, ResMapB, &NoDiff);
		WriteResMapData(File, DiffPoints, NoDiff);
		free(DiffPoints);

		FreeResMap(ResMapA);
		ResMapA = ResMapB;
		printf("Processed Map %s\n", ResMapA->FileName);fflush(stdout);
	}

	if(ResMapB == NULL && ResMapA != NULL)
		FreeResMap(ResMapA);

	if(ResMapB != NULL)
		FreeResMap(ResMapB);

	printf("\nFile (%s) created successfully.\n", MapGroup->FileName);

	FreeMapGroup(MapGroup);	

	fclose(File);

	exit(0);
}

void				LoadBinCompResMapData(FILE *File, size_t NoResPoints, RESTRICTION_POINT* ResPoints)
{
	double *CMem;
	int *ICMem;
	int Index;

	//	fread(Ret->GeoPointListConMem, sizeof(RESTRICTION_POINT), Ret->NoResPoint, File);	

	CMem = (double*)SMalloc(sizeof(double) * NoResPoints);
	
	fread(CMem, sizeof(double), NoResPoints, File);
	for(Index=0;Index<NoResPoints;Index++)
		ResPoints[Index].Long = CMem[Index];

	fread(CMem, sizeof(double), NoResPoints, File);
	for(Index=0;Index<NoResPoints;Index++)
		ResPoints[Index].Lat= CMem[Index];

	ICMem = (int*)CMem;
	fread(ICMem, sizeof(int), NoResPoints, File);
	for(Index=0;Index<NoResPoints;Index++)
		ResPoints[Index].Value= ICMem[Index];

	free(CMem);
}

RESTRICTION_MAP*	LoadBinCompResMap(FILE *File)
{
	RESTRICTION_MAP* Ret;
	size_t StrSize;

	Ret = (RESTRICTION_MAP*)SMalloc(sizeof(RESTRICTION_MAP));

	fread(&Ret->AgeMin, sizeof(double), 1, File);
	fread(&Ret->AgeMax, sizeof(double), 1, File);
	fread(&StrSize, sizeof(size_t), 1, File);

	Ret->FileName = (char*)SMalloc(sizeof(char) * (StrSize+1));
	fread(Ret->FileName, sizeof(char), StrSize, File);
	Ret->FileName[StrSize] = '\0';

	fread(&Ret->NoResPoint, sizeof(size_t), 1, File);

	Ret->GeoPointListConMem = (RESTRICTION_POINT*)SBMalloc(sizeof(RESTRICTION_POINT) * Ret->NoResPoint);
	
	LoadBinCompResMapData(File, Ret->NoResPoint, Ret->GeoPointListConMem );

	Ret->Grid = NULL;

//	printf("%s\t%f\t%f\t%zu\n", Ret->FileName, Ret->AgeMin, Ret->AgeMax, Ret->NoResPoint);

	return Ret;
}

void	LinkFistMap(RESTRICTION_MAP* Map)
{
	size_t Index;

	Map->GeoPointList = (RESTRICTION_POINT**)SMalloc(sizeof(RESTRICTION_POINT*) * Map->NoResPoint);

	for(Index=0;Index<Map->NoResPoint;Index++)
		Map->GeoPointList[Index] = &Map->GeoPointListConMem[Index];
}

void	CombinePointLists(RESTRICTION_MAP* Full, RESTRICTION_MAP* Current, size_t MaxDataPoints)
{
	size_t Index, CIndex;
	int Done;

	Done = FALSE;
	CIndex = 0;
	for(Index=0;Index<MaxDataPoints;Index++)
	{
		if( Done == FALSE &&
			Full->GeoPointList[Index]->Lat == Current->GeoPointListConMem[CIndex].Lat &&
			Full->GeoPointList[Index]->Long == Current->GeoPointListConMem[CIndex].Long)
		{
			Current->GeoPointList[Index] = &Current->GeoPointListConMem[CIndex];
			CIndex++;

			// Used to stop CIndex going out of bounds. 
			if(CIndex == Current->NoResPoint)
			{
				CIndex = 0;
				Done = TRUE;
			}
		}
		else
			Current->GeoPointList[Index] = Full->GeoPointList[Index];
	}

}

void	ReLinkMaps(RESTRICTION_MAP** Maps, int NoMaps, int NoLongBoxes, int NoLatBoxes)
{
	int Index;
	size_t MaxDataPoints;

	LinkFistMap(Maps[0]);

	MaxDataPoints = Maps[0]->NoResPoint;
	CreateGeoGrid(Maps[0], NoLongBoxes, NoLatBoxes);

	for(Index=1;Index<NoMaps;Index++)
	{
		Maps[Index]->GeoPointList = (RESTRICTION_POINT**)SMalloc(sizeof(RESTRICTION_POINT*) * MaxDataPoints);
//		memcpy(Maps[Index]->GeoPointList, Maps[0]->GeoPointList, sizeof(RESTRICTION_POINT*) * MaxDataPoints);

		CombinePointLists(Maps[Index-1], Maps[Index], MaxDataPoints);

		Maps[Index]->NoResPoint = MaxDataPoints;

		CreateGeoGrid(Maps[Index], NoLongBoxes, NoLatBoxes);
	}

}

RESTRICTION_MAP**	LoadBinCompResMaps(char *FName, int *Size)
{
	FILE *File;
	RESTRICTION_MAP** Ret;
	int Index;
	int NoLongBoxes, NoLatBoxes;

	File = fopen(FName, "rb");

	if(File == NULL)
	{
		printf("Cannot load binary maps from file %s\n", FName);
		exit(1);
	}

	fread(Size, sizeof(int), 1, File);
	fread(&NoLongBoxes, sizeof(int), 1, File);
	fread(&NoLatBoxes, sizeof(int), 1, File);

	Ret = (RESTRICTION_MAP**)SMalloc(sizeof(RESTRICTION_MAP*) * *Size);

	for(Index=0;Index<*Size;Index++)
		Ret[Index] = LoadBinCompResMap(File);

	ReLinkMaps(Ret, *Size, NoLongBoxes, NoLatBoxes);
	
	fflush(stdout);
	fclose(File);	

	return Ret;
}
