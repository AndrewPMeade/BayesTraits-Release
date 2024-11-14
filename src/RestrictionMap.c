#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "GenLib.h"
#include "RestrictionMap.h"
#include "Trees.h"
#include "Part.h"

#define BOX_CHUNK_SIZE 1024
#define MAX_LINE_SIZE 1048576	//Max line length 1MB

RESTRICTION_MAP*	InitGeoMap(void)
{
	RESTRICTION_MAP* Ret;

	Ret = (RESTRICTION_MAP*)SMalloc(sizeof(RESTRICTION_MAP));

	Ret->FileName = NULL;
	Ret->GeoPointList = NULL;
	Ret->NoResPoint = 0;
	Ret->AgeMin = -1.0;
	Ret->AgeMax = -1.0;
	Ret->GeoPointListConMem = NULL;

	return Ret;
}

void PrintResMap(RESTRICTION_MAP* GeoMap)
{
	int Index;
	for(Index=0;Index<GeoMap->NoResPoint;Index++)
		printf("%d\t%f\t%f\t%d\n", Index, GeoMap->GeoPointList[Index]->Long, GeoMap->GeoPointList[Index]->Lat, GeoMap->GeoPointList[Index]->Value);
}

RESTRICTION_POINT*	GeoPointFromLine(char **Passed)
{
	RESTRICTION_POINT* Ret;

	Ret = (RESTRICTION_POINT*)SMalloc(sizeof(RESTRICTION_POINT));
	Ret->Long = atof(Passed[0]);
	Ret->Lat = atof(Passed[1]);
	Ret->Value = atoi(Passed[2]);

	return Ret;
}

RESTRICTION_BOX* CrateGeoBox(void)
{
	RESTRICTION_BOX* Ret;

	Ret = (RESTRICTION_BOX*)SMalloc(sizeof(RESTRICTION_BOX));
	Ret->NoResPoints = 0;
	Ret->PointList = (RESTRICTION_POINT**)SMalloc(sizeof(RESTRICTION_POINT*) * BOX_CHUNK_SIZE);

	return Ret;
}

void	FreeGeoBox(RESTRICTION_BOX *Box)
{
	free(Box->PointList);
	free(Box);
}

RESTRICTION_BOX**	CreateGridRow(size_t NoLong)
{
	RESTRICTION_BOX** Ret;
	size_t Index;

	Ret = (RESTRICTION_BOX**)SMalloc(sizeof(RESTRICTION_BOX*) * NoLong);

	for(Index=0;Index<NoLong;Index++)
		Ret[Index] = CrateGeoBox();

	return Ret;
}

void	InitGeoGrid(RESTRICTION_MAP* GeoMap, size_t NoLong, size_t NoLat)
{
	size_t IndexY;

	GeoMap->Grid = (RESTRICTION_BOX***)SMalloc(sizeof(RESTRICTION_BOX**)*NoLat);

	for(IndexY=0;IndexY<NoLat;IndexY++)
		GeoMap->Grid[IndexY] = CreateGridRow(NoLong);
}

void	GetPointLocation(RESTRICTION_MAP* GeoMap, double Long, double Lat, size_t *LongI, size_t *LatI)
{
	*LongI = (size_t)((Long + 180) * (GeoMap->NoGridLong/360.0));
	*LatI = (size_t)((Lat + 90) * (GeoMap->NoGridLat/180.0));
}

void	PlacePointInBox(RESTRICTION_BOX *Box, RESTRICTION_POINT *Point)
{
	RESTRICTION_POINT **NList;

	if((Box->NoResPoints + 1) % BOX_CHUNK_SIZE == 0)
	{
		NList = (RESTRICTION_POINT**)SMalloc(sizeof(RESTRICTION_POINT*) * (Box->NoResPoints + BOX_CHUNK_SIZE));
		memcpy(NList, Box->PointList, sizeof(RESTRICTION_POINT*) * Box->NoResPoints);
		free(Box->PointList);
		Box->PointList = NList;
	}

	Box->PointList[Box->NoResPoints++] = Point;
}

void	IndexToLongLat(RESTRICTION_MAP* GeoMap, size_t LongI, size_t LatI, double *Long, double *Lat)
{
	*Long = (((LongI + 0.5) / GeoMap->NoGridLong) * 360.0) - 180; 
	*Lat = (((LatI + 0.5) / GeoMap->NoGridLat) * 180.0) - 90; 
}

void PrintBox(RESTRICTION_BOX *Box)
{
	size_t Index;

	printf("Box:\t%zu\n", Box->NoResPoints);

	for(Index=0;Index<Box->NoResPoints;Index++)
		printf("%zu\t%f\t%f\n", Index, Box->PointList[Index]->Long, Box->PointList[Index]->Lat);

	exit(0);
}

void FillEmptyBox(RESTRICTION_MAP *GeoMap, RESTRICTION_BOX *Box, size_t LongI, size_t LatI)
{
	double Long, Lat;
	RESTRICTION_POINT*	Point;


	IndexToLongLat(GeoMap, LongI, LatI, &Long, &Lat);
	Point = GetClosestPoint(GeoMap, Long, Lat);

	PlacePointInBox(GeoMap->Grid[LatI][LongI], Point);

}

void FillZeroMap(RESTRICTION_MAP* GeoMap)
{
	size_t LatI, LongI;


	for(LatI=0;LatI<GeoMap->NoGridLat;LatI++)
		for(LongI=0;LongI<GeoMap->NoGridLong;LongI++)
			if(GeoMap->Grid[LatI][LongI]->NoResPoints == 0)
				FillEmptyBox(GeoMap, GeoMap->Grid[LatI][LongI], LongI, LatI);

}

void PlaceGeoGrid(RESTRICTION_MAP* GeoMap)
{
	size_t Index;
	size_t LatI, LongI;
	RESTRICTION_POINT *Point;


	LatI = LongI = 0;
	for(Index=0;Index<GeoMap->NoResPoint;Index++)
	{
		Point = GeoMap->GeoPointList[Index];

		GetPointLocation(GeoMap, Point->Long, Point->Lat, &LongI, &LatI);

		if(LongI >= GeoMap->NoGridLong || LatI >= GeoMap->NoGridLat)
		{
			printf("point (%zu\t%f\t%f) in file (%s) falls outside range.\n", Index, Point->Long, Point->Lat, GeoMap->FileName);
			exit(1);
		}
		else
			PlacePointInBox(GeoMap->Grid[LatI][LongI], GeoMap->GeoPointList[Index]);
	}

	FillZeroMap(GeoMap);

	return;

	for(LatI=0;LatI<GeoMap->NoGridLat;LatI++)
	{
		for(LongI=0;LongI<GeoMap->NoGridLong;LongI++)
			printf("%zu\t", GeoMap->Grid[LatI][LongI]->NoResPoints);
		printf("\n");
	}

	exit(0);
}

void	CompressGeoGridBox(RESTRICTION_BOX *Box)
{
	size_t Index;

	Box->No0 = Box->No1 = 0;
	for(Index=0;Index<Box->NoResPoints;Index++)
	{
		if(Box->PointList[Index]->Value == 0)
			Box->No0++;

		if(Box->PointList[Index]->Value == 1)
			Box->No1++;
	}

	if(Box->No0 == Box->NoResPoints || Box->No1 == Box->NoResPoints)
		Box->NoResPoints = 1;
}

void	OrderGeoGridBox(RESTRICTION_BOX *Box)
{

}

void CreateGeoGrid(RESTRICTION_MAP* GeoMap, size_t NoLong, size_t NoLat)
{
	InitGeoGrid(GeoMap, NoLong, NoLat);
	GeoMap->NoGridLong = NoLong;
	GeoMap->NoGridLat = NoLat;

	PlaceGeoGrid(GeoMap);

	size_t LatI, LongI;

	for(LatI=0;LatI<GeoMap->NoGridLat;LatI++)
	{
		for(LongI=0;LongI<GeoMap->NoGridLong;LongI++)
		{
			CompressGeoGridBox(GeoMap->Grid[LatI][LongI]);
			OrderGeoGridBox(GeoMap->Grid[LatI][LongI]);
		}
	}
}

size_t	GetNoLines(FILE *InFile)
{
	char *Buffer;
	size_t NoLines;

	NoLines = 0;	
	Buffer = (char*)SMalloc(sizeof(char) * MAX_LINE_SIZE);

	while(fgets(Buffer,MAX_LINE_SIZE,InFile) != NULL)
		NoLines++;
	
	(void)fseek(InFile, 0L, SEEK_SET);

	free(Buffer);
	return NoLines;
}


void	SetPointList(RESTRICTION_MAP* ResMap, size_t MaxDataPoints)
{
	int Index;

	ResMap->GeoPointList = (RESTRICTION_POINT**)SMalloc(sizeof(RESTRICTION_POINT*) * MaxDataPoints);

	for(Index=0;Index<MaxDataPoints;Index++)
	{
		ResMap->GeoPointList[Index] = &ResMap->GeoPointListConMem[Index];
	}

}

void	AllocDataPoints(RESTRICTION_MAP* ResMap, size_t MaxDataPoints)
{
//	size_t Index;


	ResMap->GeoPointListConMem = (RESTRICTION_POINT*)SMalloc(sizeof(RESTRICTION_POINT) * MaxDataPoints);

	SetPointList(ResMap, MaxDataPoints);
}

void	ReadPointData(FILE *FIn, RESTRICTION_MAP* ResMap)
{
	char *Buffer;
	char **Passed;
	int	Tokes;
	size_t Index;
	RESTRICTION_POINT *Point;


	Buffer = (char*)SMalloc(sizeof(char) * (MAX_LINE_SIZE+1));
	Passed = (char**)SMalloc(sizeof(char*) * MAX_LINE_SIZE);

	// disgard the header line. 
	fgets(Buffer,MAX_LINE_SIZE, FIn);
	
	Index = 0;
	while(fgets(Buffer,MAX_LINE_SIZE, FIn) != NULL)
	{
		ReplaceChar('\t', ',', Buffer);
		Tokes = MakeArgvChar(Buffer, Passed, MAX_LINE_SIZE, ',');
		if(Tokes == 3)
		{
			Point = ResMap->GeoPointList[Index];
				
			Point->Long = atof(Passed[0]);
			Point->Lat = atof(Passed[1]);
			Point->Value = atoi(Passed[2]);

			Index++;
		}
	}

	ResMap->NoResPoint = Index;

	free(Buffer);
	free(Passed);
}



RESTRICTION_MAP*	LoadResMap(char *FileName, double AgeMin, double AgeMax)
{
	RESTRICTION_MAP* Ret;
	size_t NoLines;
	FILE *FIn;

	FIn = fopen(FileName, "rb");
	if(FIn == NULL)
	{
		fprintf(stdout, "Could not open file %s for reading\n", FileName);
		exit(1);
	}

	NoLines = GetNoLines(FIn);

	Ret = InitGeoMap();
	Ret->FileName = StrMake(FileName);
	Ret->AgeMin = AgeMin;
	Ret->AgeMax = AgeMax;
	Ret->NoResPoint = 0;
	

	AllocDataPoints(Ret, NoLines);
	ReadPointData(FIn, Ret);

	CreateGeoGrid(Ret, NO_GEO_MAP_BOX_SIZE_LONG, NO_GEO_MAP_BOX_SIZE_LAT);

	fclose(FIn);

	return Ret;
}

/*
RESTRICTION_MAP*	LoadResMap(char *FileName, double AgeMin, double AgeMax)
{
	RESTRICTION_MAP* Ret;
	TEXTFILE* TF;
	int Index, Tokes;
	char **Passed;


	NewLoadResMap(FileName, AgeMin, AgeMax);

	Ret = InitGeoMap();
	Ret->FileName = StrMake(FileName);
	Ret->AgeMin = AgeMin;
	Ret->AgeMax = AgeMax;
	
	TF = LoadTextFile(FileName, FALSE);

	Ret->GeoPointList = (RESTRICTION_POINT**)SMalloc(sizeof(RESTRICTION_POINT*)*TF->NoOfLines);
	Passed = (char**)SMalloc(sizeof(char*) * TF->MaxLine);

	for(Index=0;Index<TF->NoOfLines;Index++)
	{
		Tokes = MakeArgvChar(TF->Data[Index], Passed, TF->NoOfLines, ',');
		if(Tokes == 3)
			Ret->GeoPointList[Ret->NoResPoint++] = GeoPointFromLine(Passed);
		
	}

//	PrintGeoMap(Ret);
	
	CreateGeoGrid(Ret, 250, 150);

	FreeTextFile(TF);
	free(Passed);
	return Ret;
}
*/

void FreeGeoGrid(RESTRICTION_MAP* GeoMap)
{
	size_t LongI, LatI;

	if(GeoMap->Grid == NULL)
		return;

	for(LatI=0;LatI<GeoMap->NoGridLat;LatI++)
	{
		for(LongI=0;LongI<GeoMap->NoGridLong;LongI++)
			FreeGeoBox(GeoMap->Grid[LatI][LongI]);

		free(GeoMap->Grid[LatI]);
	}

	free(GeoMap->Grid);
}

void FreeResMap(RESTRICTION_MAP* GeoMap)
{
	FreeGeoGrid(GeoMap);

	free(GeoMap->GeoPointList);
	free(GeoMap->GeoPointListConMem);

	free(GeoMap->FileName);
	free(GeoMap);
}

double	GeoDist(double Long1, double Lat1, double Long2, double Lat2)
{
	double Long, Lat;

	Long = (Long1-Long2) * (Long1-Long2);
	Lat = (Lat1-Lat2) * (Lat1-Lat2);
	return sqrt(Long + Lat);
}

/*
double	GeoDist(double Long1, double Lat1, double Long2, double Lat2)
{
	double nDLat, nDLong, nA, nC;

	nDLat = (Lat1 - Lat2) *  0.017453293;
	nDLong = (Long1 - Long2) *  0.017453293;

	Lat1 = Lat1 * 0.017453293;
	Lat2 = Lat2 * 0.017453293;

	nA = (sin(nDLat/2) * sin(nDLat/2)) + cos(Lat1) * cos(Lat2) * (sin(nDLong/2) * sin(nDLong/2));
	nC = 2 * atan2(sqrt(nA), sqrt(1 - nA));
	return nC * 6372.797;
}
*/

double GeoPointDist(RESTRICTION_POINT* Point, double Long, double Lat)
{
	return GeoDist(Point->Long, Point->Lat, Long, Lat);
}

RESTRICTION_POINT*	GetClosestPointList(RESTRICTION_POINT **PointList, size_t Size, double Long, double Lat)
{
	RESTRICTION_POINT* Ret;
	double MinDist, Dist;
	size_t Index;

	if(Size == 0)
		return NULL;

	Ret = PointList[0];
	MinDist = GeoPointDist(Ret, Long, Lat);

	for(Index=1;Index<Size;Index++)
	{
		Dist = GeoPointDist(PointList[Index], Long, Lat);
		if(Dist < MinDist)
		{
			MinDist = Dist;
			Ret = PointList[Index];
		}
	}

	return Ret;
}


RESTRICTION_POINT*	GetClosestPoint(RESTRICTION_MAP* GeoMap, double Long, double Lat)
{
	return GetClosestPointList(GeoMap->GeoPointList, GeoMap->NoResPoint, Long, Lat);
}

RESTRICTION_POINT*	GetClosestPointBox(RESTRICTION_MAP* GeoMap, double Long, double Lat)
{
	size_t	LongI, LatI;
	RESTRICTION_BOX *Box;

	GetPointLocation(GeoMap, Long, Lat, &LongI, &LatI);

	Box = GeoMap->Grid[LatI][LongI];

	if(Box->NoResPoints == 1)
		return Box->PointList[0];
		
	return GetClosestPointList(Box->PointList, Box->NoResPoints, Long, Lat);
}


int	ValidResPoint(RESTRICTION_MAP* Map, double Long, double Lat)
{
	RESTRICTION_POINT *CPoint;

	if(Map == NULL)
		return TRUE;

	CPoint = GetClosestPointBox(Map, Long, Lat);
	
	if(CPoint->Value == 1)
		return TRUE;

	return FALSE;
}

int ValidNodeResPoint(NODE_RES_MAP *NodeResMap, double Long, double Lat)
{
	int Valid;

	Valid = ValidResPoint(NodeResMap->ResMap, Long, Lat);

	if(NodeResMap->Flipped == TRUE)
	{
		if(Valid == TRUE)
			return FALSE;
		else
			return TRUE;
	}

	return Valid;
}

int	GetNoGlobal(RESTRICTION_MAP** Maps, int NoMaps)
{
	int Index, Ret;

	Ret = 0;
	for(Index=0;Index<NoMaps;Index++)
		if(Maps[Index]->AgeMin == GLOBAL_AGE_VAL)
			Ret++;

	return Ret;
}

int HeightInMapRange(double Height, RESTRICTION_MAP* Map)
{
	if(Height >= Map->AgeMin && Height < Map->AgeMax)
		return TRUE;

	return FALSE;
}

int	OverLap(RESTRICTION_MAP *MapX, RESTRICTION_MAP *MapY)
{
	if(MapX == MapY)
		return FALSE;
	
	if(fmax(MapX->AgeMin, MapY->AgeMin) < fmin(MapX->AgeMax, MapY->AgeMax))
		return TRUE;

	return FALSE;
}

void CheckMapOverLaps(RESTRICTION_MAP** Maps, int NoMaps)
{
	int X, Y;
	RESTRICTION_MAP *MapX, *MapY;
	
	if(GetNoGlobal(Maps, NoMaps) == 1)
		return;

	for(X=0;X<NoMaps;X++)
	{
		MapX = Maps[X];
		for(Y=X+1;Y<NoMaps;Y++)
		{			
			MapY = Maps[Y];

			if(OverLap(MapX, MapY) == TRUE)
			{
				printf("Maps overlap\n%s\t%f\t%f\n%s\t%f\t%f\n", MapX->FileName, MapX->AgeMin, MapX->AgeMax, MapY->FileName, MapY->AgeMin, MapY->AgeMax);
				exit(1);
			}
		}
	}
}

void CheckResMaps(RESTRICTION_MAP** Maps, int NoMaps)
{
	if(NoMaps == 0)
		return;

	if(NoMaps > 1 && GetNoGlobal(Maps, NoMaps) >= 1)
	{
		printf("A global map is being used with either multiple global maps or age restricted maps.");
		exit(1);
	}

	CheckMapOverLaps(Maps, NoMaps);
}



RESTRICTION_MAP* GetMapFromHeight(double Height, RESTRICTION_MAP** Maps, int NoMaps)
{
	int Index;
	RESTRICTION_MAP* Map;

	for(Index=0;Index<NoMaps;Index++)
	{
		Map = Maps[Index];

		if(Map->AgeMax == GLOBAL_AGE_VAL || HeightInMapRange(Height, Map) == TRUE)
			return Map;
		
	}

	return NULL;
}

void PrintResMaps(TREES *Trees)
{
	TREE *Tree;
	int TIndex, NIndex;
	NODE Node;


	for(TIndex=0;TIndex<Trees->NoTrees;TIndex++)
	{
		Tree = Trees->Tree[TIndex];
		for(NIndex=0;NIndex<Tree->NoNodes;NIndex++)
		{
			Node = Tree->NodeList[NIndex];

			printf("%f\t", Node->Height);
			if(Node->NodeResMap!= NULL)
				printf("%s\t", Node->NodeResMap->ResMap->FileName);
			else
				printf("NA\t");
			PrintPartTaxaOnly(stdout, Trees, Tree->NodeList[NIndex]->Part);
			printf("\n"); 
		}
	}
}

NODE_RES_MAP*	CrateNodeResMap(RESTRICTION_MAP *ResMap)
{
	NODE_RES_MAP* NodeResMap;

	NodeResMap = SMalloc(sizeof(NODE_RES_MAP));

	NodeResMap->Flipped = FALSE;
	NodeResMap->ResMap = ResMap;


	return NodeResMap;

}
void	FreeNodeResMap(NODE_RES_MAP* NodeResMap)
{
	free(NodeResMap);
}

void	FlipTag(NODE Node)
{
	int Index;

	if(Node->NodeResMap->Flipped == TRUE)
		Node->NodeResMap->Flipped = FALSE;
	else
		Node->NodeResMap->Flipped = TRUE;

	for(Index=0;Index<Node->NoNodes;Index++)
		FlipTag(Node->NodeList[Index]);
}


void	SetFlippedMaps(OPTIONS *Opt, TREES *Trees)
{
	int Index;
	TAG *Tag;

	for(Index=0;Index<Opt->NoFlippedNodes;Index++)
	{
		Tag = Opt->FlippedNodes[Index];
		FlipTag(Tag->NodeList[0]);
	}
}

void MapResMaps(OPTIONS *Opt, TREES *Trees, RESTRICTION_MAP** Maps, int NoMaps)
{
	TREE *Tree;
	int TIndex, NIndex;
	NODE Node;
	RESTRICTION_MAP *RMap;

	if(NoMaps == 0)
		return;

	for(TIndex=0;TIndex<Trees->NoTrees;TIndex++)
	{
		Tree = Trees->Tree[TIndex];
		for(NIndex=0;NIndex<Tree->NoNodes;NIndex++)
		{
			Node = Tree->NodeList[NIndex];
			RMap = GetMapFromHeight(Node->Height, Maps, NoMaps);
			Node->NodeResMap = CrateNodeResMap(RMap);
		}
	}

	SetFlippedMaps(Opt, Trees);
//	PrintResMaps(Trees);
//	exit(0);
}

size_t FindDiff(RESTRICTION_MAP* Map1, RESTRICTION_MAP* Map2)
{
	size_t Index, Diff;
	RESTRICTION_POINT *P1, *P2;
	
	Diff = 0;
	Index = 0;

	for(Index=0;Index<Map1->NoResPoint;Index++)
	{
		P1 = Map1->GeoPointList[Index];
		P2 = Map2->GeoPointList[Index];

		if(P1->Lat != P2->Lat || P1->Long != P2->Long)
		{
			printf("Diff points\n");
			exit(1);
		}

		if(P1->Value != P2->Value)
			Diff++;
	}

	return Diff;
}



