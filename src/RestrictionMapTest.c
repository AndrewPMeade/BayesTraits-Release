#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "RestrictionMapTest.h"
#include "RestrictionMap.h"

void TestLine(RESTRICTION_MAP** ResMapList, int NoMaps, int Index, double Height, double Long, double Lat)
{
	RESTRICTION_POINT*	Point;
	RESTRICTION_MAP* GeoMap;
	

	GeoMap = GetMapFromHeight(Height, ResMapList, NoMaps);

	Point = GetClosestPointBox(GeoMap, Long, Lat);

	printf("%d\t%s\t%f\t%f\t%f\t", Index, GeoMap->FileName, Height, Long, Lat);

	printf("%f\t%f\t%d\t", Point->Long, Point->Lat, Point->Value);

	printf("\n");
}

void TestRestrictionMap(int argc, char **argv)
{
	int NoResMaps, Index, Tokes;
	RESTRICTION_MAP** ResMapList;
	TEXTFILE *TF;
	char **Passed;

	if(argc != 3)
	{
		printf("TestRestrictionMap takes a binary map file and three col file for heights, longitude and latitude, tab separated.\n");
		exit(1);
	}

	ResMapList = LoadBinCompResMaps(argv[1], &NoResMaps);


	TF = LoadTextFile(argv[2], FALSE);
	Passed = (char**)SMalloc(sizeof(char*) * TF->MaxLine);

	for(Index=1;Index<TF->NoOfLines;Index++)
	{
		ReplaceChar(',', '\t', TF->Data[Index]);
		Tokes = MakeArgv(TF->Data[Index], Passed, TF->MaxLine);

		if(Tokes == 3)
			TestLine(ResMapList, NoResMaps, Index, atof(Passed[0]), atof(Passed[1]), atof(Passed[2]));
	}


	FreeTextFile(TF);
	free(Passed);
	for(Index=0;Index<NoResMaps;Index++)
		FreeResMap(ResMapList[Index]);

	free(ResMapList);
}