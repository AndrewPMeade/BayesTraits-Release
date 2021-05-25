#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "typedef.h"
#include "Geo.h"
#include "genlib.h"
#include "threaded.h"
#include "FatTail.h"

#ifndef M_PI
	#define M_PI 3.14159265358979323846
#endif

#define EARTH_KM 6371.0
#define EARTH_KM_EARTH_KM 40589641

int		ValidGeoData(TREES *Trees)
{
	return FALSE;
}

void	LongLatToXYZ(double Long, double Lat, double *X, double *Y, double *Z)
{
	Long = Long * M_PI / 180.0;
	Lat = Lat * M_PI / 180.0;

	*X = EARTH_KM * cos(Lat) * cos(Long);
	*Y = EARTH_KM * cos(Lat) * sin(Long);
	*Z = EARTH_KM * sin(Lat);
}

void	XYZToLongLat(double X, double Y, double Z, double *Long, double *Lat)
{
	*Long = atan2(Y, X);
	*Lat = atan2(Z, sqrt(X * X + Y * Y));

	*Lat = *Lat * (180.0 / M_PI);
	*Long = *Long * (180.0 / M_PI);
}

void	GetRandLongLat(RANDSTATES *RS, double Long, double Lat, double *RLong, double *RLat, double Radius)
{
	double Dist, Bearing;

	Long = Long * M_PI / 180.0;
	Lat = Lat * M_PI / 180.0;

	Dist = RandDouble(RS) * Radius;

	Bearing	= RandDouble(RS) * 360;
	Bearing = Bearing * M_PI / 180.0;
	
	Bearing = fmod(Bearing + (2.0 * M_PI), (2.0 * M_PI));
	
	*RLat = sin(Lat) * cos(Dist/EARTH_KM) + cos(Lat) * sin(Dist/EARTH_KM) * cos(Bearing);
	*RLat = asin(*RLat);

	*RLong = Long + atan2(sin(Bearing) * sin(Dist/EARTH_KM) * cos(Lat), cos(Dist/EARTH_KM) - sin(Lat) * sin(*RLat));
	
	*RLat = *RLat * (180.0 / M_PI);
		
	if(*RLong * (180.0 / M_PI) > 180)
		*RLong = *RLong * (180.0 / M_PI) - 360.0;
	else
		*RLong = *RLong * (180.0 / M_PI);
}



void	GetRandXYZPoint(RANDSTATES *RS, double SX, double SY, double SZ, double *RX, double *RY, double *RZ, double Radius)
{
	double Long, Lat, RLong, RLat;

	XYZToLongLat(SX, SY, SZ, &Long, &Lat);
	
	GetRandLongLat(RS, Long, Lat, &RLong, &RLat, Radius);

//	printf("CLongLat:\t\t%f\t%f\t%f\t%f\t", Long, Lat, RLong, RLat);
	
	LongLatToXYZ(RLong, RLat, RX, RY, RZ);
}

double	GeoCalcAnsStateLh(double X, double Y, double Z, NODE N, STABLEDIST *SD)
{
	double Ret, Val;
	int Index;

	Ret = 0;

	for(Index=0;Index<N->NoNodes;Index++)
	{
		Val = X - N->NodeList[Index]->FatTailNode->Ans[0];
		Ret += StableDistTPDF(SD, Val, N->NodeList[Index]->Length);

		Val = Y - N->NodeList[Index]->FatTailNode->Ans[1];
		Ret += StableDistTPDF(SD, Val, N->NodeList[Index]->Length);

		Val = Z - N->NodeList[Index]->FatTailNode->Ans[2];
		Ret += StableDistTPDF(SD, Val, N->NodeList[Index]->Length);
	}

	if(N->Ans != NULL)
	{
		Val = X - N->Ans->FatTailNode->Ans[0];
		Ret += StableDistTPDF(SD, Val, N->Length);

		Val = Y - N->Ans->FatTailNode->Ans[1];
		Ret += StableDistTPDF(SD, Val, N->Length);

		Val = Z - N->Ans->FatTailNode->Ans[2];
		Ret += StableDistTPDF(SD, Val, N->Length);
	}

	return Ret;
}

void	NodeToXYZ(NODE N, double *X, double *Y, double *Z)
{
	*X = N->FatTailNode->Ans[0];
	*Y = N->FatTailNode->Ans[1];
	*Z = N->FatTailNode->Ans[2];
}

void	XYZToNode(NODE N, double X, double Y, double Z)
{
	N->FatTailNode->Ans[0] = X;
	N->FatTailNode->Ans[1] = Y;
	N->FatTailNode->Ans[2] = Z;
}

void	NodeToLongLat(NODE N, double *Long, double *Lat)
{
	double X, Y, Z;

	NodeToXYZ(N, &X, &Y, &Z);
	XYZToLongLat(X, Y, Z, Long, Lat);
}

void	GeoUpDateNode(NODE N, RATES *Rates)
{
	FATTAILRATES *FTR;
	double	NLh, CLh, X, Y, Z, RX, RY, RZ;
	STABLEDIST *SD;
	int Index;
	
	FTR = Rates->FatTailRates;
	SD = FTR->SDList[0];

	NodeToXYZ(N, &X, &Y, &Z);
	
	CLh = GeoCalcAnsStateLh(X, Y, Z, N, SD);
		

	for(Index=0;Index<10000;Index++)
	{
		GetRandXYZPoint(Rates->RS, X, Y, Z, &RX, &RY, &RZ, 1000);
		
		NLh = GeoCalcAnsStateLh(RX, RY, RZ, N, SD);

		if(log(RandDouble(Rates->RS)) < (NLh - CLh))
		{
			XYZToNode(N, RX, RY, RZ);
			return;
		}

//		printf("%d\t%f\t%f\t%f\n", Index, CLh-NLh, NLh, CLh);
	//	printf("%f\t%f\t%f\t%f\t%f\t%f\t\n", X, Y, Z, RX, RY, RZ);
	}

	//exit(0);
}

void	GeoUpDateAllAnsStates(OPTIONS *Opt, TREES *Trees, RATES *Rates)
{
	int NIndex;
	TREE *Tree;
	FATTAILRATES *FTR;
	NODE N;
	
	Tree = Trees->Tree[Rates->TreeNo];

	FTR = Rates->FatTailRates;

	FatTailSetAnsSates(Tree, Trees->NoOfSites, FTR);
//	CorrectIntGeoNodes(Tree);
//	MapRatesToFatTailRate(Rates, FTR);
	SetStableDist(FTR->SDList[0], FTR->Alpha[0], FTR->Scale[0]);

	for(NIndex=0;NIndex<Tree->NoNodes;NIndex++)
	{
		N = Tree->NodeList[NIndex];
		if(N->Tip == FALSE)
			GeoUpDateNode(N, Rates);
		
	}
//	exit(0);
	

	FatTailGetAnsSates(Tree, Trees->NoOfSites, FTR);
}


void	CorrectIntGeoNodes(TREE *Tree)
{
	NODE			N;
	int				NIndex;
	double			X,Y,Z, Long, Lat;
	
	for(NIndex=0;NIndex<Tree->NoNodes;NIndex++)
	{
		N = Tree->NodeList[NIndex];
		if(N->Tip == FALSE)
		{
			NodeToXYZ(N, &X, &Y, &Z);

			XYZToLongLat(X, Y, Z, &Long, &Lat);
			LongLatToXYZ(Long, Lat, &X, &Y, &Z);

			XYZToNode(N, X, Y, Z);
		}
	}
}


void	NewGeoTest(RANDSTATES *RS)
{
	double Lat, Long, RLong, RLat;
	double x, y, z;
	int i;

//	GetRandLongLat(RS, -18, 65, &RLong, &RLat, 1000);
//	exit(0);
	Lat = 0;
	Long = 0;
	for(i=0;i<10000;i++)
	{
		GetRandLongLat(RS, 0.0, 0.0, &RLong, &RLat, 1000);
//		GetRandLongLat(RS, Long, Lat, &RLong, &RLat, EARTH_KM*10);
//		GetRandLongLat(RS, Long, Lat, &RLong, &RLat, 100);
		LongLatToXYZ(RLong, RLat, &x, &y, &z);
		printf("%d\t%f\t%f\t%f\t%f\t%f\n", i, RLong, RLat, x, y, z);

		Long = RLong;
		Lat = RLat;
	}


	exit(0);
}

void GeoTest(OPTIONS *Opt, TREES *Trees, RATES *Rates)
{
	double Long, Lat, x, y, z;
	double S;
	int i;

	NewGeoTest(Rates->RS);

		Lat = RandIntBetween(Rates->RS, -89, 89);

		Long = RandIntBetween(Rates->RS, -180, 180);

//		XYZToLatLong(2252.488651, 3901.424788, 4504.977303, &Lat, &Long);

	
	S = GetSeconds();
	for(i=0;i<100000;i++)
	{

		Lat = (double)RandIntBetween(Rates->RS, -90, 90);
		Long = (double)RandIntBetween(Rates->RS, -180, 180);

		LongLatToXYZ(Long, Lat, &x, &y, &z);

		printf("%d\t%f\t%f\t%f\t%f\t%f\t", i, Lat, Long, x, y, z);

		XYZToLongLat(x, y, z, &Long, &Lat);
		printf("%f\t%f\n", Lat, Long);

	}
	S = GetSeconds() - S;

	printf("Sec:\t%f\t%lu\n", S, i);

	printf("%f\t%f\t%f\n", x, y, z);

	exit(0);
}