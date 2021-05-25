#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "typedef.h"
#include "Geo.h"
#include "genlib.h"
#include "Threaded.h"
#include "FatTail.h"
#include "likelihood.h"
#include "Rates.h"


#ifndef M_PI
	#define M_PI 3.14159265358979323846
#endif

#define EARTH_KM 6371.0
//#define EARTH_KM 1.0

#define EARTH_KM_EARTH_KM 40589641

void	ValidGeoData(TREES *Trees)
{
	int Index;
	double Long, Lat;
	TAXA *T;

	if(Trees->NoOfSites != 2)
	{
		printf("Geo Data must be long lat\n");
		exit(1);
	}

	for(Index=0;Index<Trees->NoOfTaxa;Index++)
	{
		T = Trees->Taxa[Index];

		if(T->Exclude == FALSE)
		{
			Long = T->ConData[0];
			Lat = T->ConData[1];

			if(Long < -180 || Long > 180)
			{
				printf("Invalid longitude (%f) for taxa %s.\n", Long, T->Name);
				exit(1);
			}

			if(Lat < -90 || Lat > 90)
			{
				printf("Invalid latitude (%f) for taxa %s.\n", Lat, T->Name);
				exit(1);
			}
		}
	}
}

void	SetGeoMissingData(TAXA *Taxa)
{
	char *New;

	New = (char*)SMalloc(sizeof(char) * 3);

	if(Taxa->EstDataP[0] == TRUE || Taxa->EstDataP[1] == TRUE)
		New[0] = New[1] = New[2] = TRUE;
	else
		New[0] = New[1] = New[2] = FALSE;

	free(Taxa->EstDataP);
	Taxa->EstDataP = New;
}

void	PreProcessGeoData(TREES *Trees)
{
	int TIndex;
	TAXA *Taxa;
	double X, Y, Z;

	for(TIndex=0;TIndex<Trees->NoOfTaxa;TIndex++)
	{
		Taxa = Trees->Taxa[TIndex];
		if(Taxa->Exclude == FALSE)
		{
			free(Taxa->DesDataChar[0]);
			free(Taxa->DesDataChar[1]);
			free(Taxa->DesDataChar);
			Taxa->DesDataChar = NULL;


			LongLatToXYZ(Taxa->ConData[0], Taxa->ConData[1], &X, &Y, &Z);
			free(Taxa->ConData);

			Taxa->ConData = (double*)SMalloc(sizeof(double) * 3);
			Taxa->ConData[0] = X;
			Taxa->ConData[1] = Y;
			Taxa->ConData[2] = Z;

			SetGeoMissingData(Taxa);
		}
	}
	Trees->NoOfSites = 3;
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

void	LongLatToNode(NODE N, double Long, double Lat)
{
	double X, Y, Z;

	LongLatToXYZ(Long, Lat, &X, &Y, &Z);
	XYZToNode(N, X, Y, Z);
}

/* Self MCMC
void	GeoUpDateNode(NODE N, RATES *Rates, RANDSTATES *RS)
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
		GetRandXYZPoint(RS, X, Y, Z, &RX, &RY, &RZ, 1000);
		
		NLh = GeoCalcAnsStateLh(RX, RY, RZ, N, SD);

		if(log(RandDouble(RS)) < (NLh - CLh))
		{
			XYZToNode(N, RX, RY, RZ);
			return;
		}
	}
//	printf("Fail.\n");
}
*/

void	GeoUpDateNode(NODE N, RATES *Rates, RANDSTATES *RS)
{
	FATTAILRATES *FTR;
	double NLh, CLh, X, Y, Z, RX, RY, RZ;
	STABLEDIST *SD;
	

	FTR = Rates->FatTailRates;
	SD = FTR->SDList[0];

	NodeToXYZ(N, &X, &Y, &Z);
	
	CLh = GeoCalcAnsStateLh(X, Y, Z, N, SD);

	GetRandXYZPoint(RS, X, Y, Z, &RX, &RY, &RZ, 2000);
		
	NLh = GeoCalcAnsStateLh(RX, RY, RZ, N, SD);
	
	XYZToNode(N, RX, RY, RZ);

	Rates->Lh = Rates->Lh + (NLh - CLh);
}


void	GeoForceUpDateNode(NODE N, RATES *Rates, RANDSTATES *RS)
{
	FATTAILRATES *FTR;
	double NLh, CLh, X, Y, Z, RX, RY, RZ;
	STABLEDIST *SD;
	int Changed;


	FTR = Rates->FatTailRates;
	SD = FTR->SDList[0];

	NodeToXYZ(N, &X, &Y, &Z);

	CLh = GeoCalcAnsStateLh(X, Y, Z, N, SD);
	Changed = FALSE;
	do
	{
		GetRandXYZPoint(RS, X, Y, Z, &RX, &RY, &RZ, 2000);
		NLh = GeoCalcAnsStateLh(RX, RY, RZ, N, SD);

		if(log(RandDouble(Rates->RS)) < (NLh - CLh))
			Changed = TRUE;

	} while(Changed == FALSE);

	XYZToNode(N, RX, RY, RZ);

	Rates->Lh = Rates->Lh + (NLh - CLh);
}

/* // Serial
void	GeoUpDateAllAnsStates(OPTIONS *Opt, TREES *Trees, RATES *Rates)
{
	int NIndex;
	TREE *Tree;
	FATTAILRATES *FTR;
	NODE N;
	
	Tree = Trees->Tree[Rates->TreeNo];

	FTR = Rates->FatTailRates;

	FatTailSetAnsSates(Tree, Trees->NoOfSites, FTR);
	SetStableDist(FTR->SDList[0], FTR->Alpha[0], FTR->Scale[0]);
	
	for(NIndex=0;NIndex<Tree->NoNodes;NIndex++)
	{
		N = Tree->NodeList[NIndex];
		if(N->Tip == FALSE)
			GeoUpDateNode(N, Rates);
		
	}

	FatTailGetAnsSates(Tree, Trees->NoOfSites, FTR);
}
*/

void	GeoUpDateAllAnsStates(OPTIONS *Opt, TREES *Trees, RATES *Rates)
{
	int NIndex;
	
	TREE *Tree;
	FATTAILRATES *FTR;
	NODE N;
	
	LhTransformTree(Rates, Trees, Opt);

	Tree = Trees->Tree[Rates->TreeNo];

	FTR = Rates->FatTailRates;

	MapRatesToFatTailRate(Rates, FTR);

	FatTailSetAnsSates(Tree, Trees->NoOfSites, FTR);

	SetStableDist(FTR->SDList[0], FTR->Alpha[0], FTR->Scale[0]);
	
	for(NIndex=0;NIndex<Tree->NoNodes;NIndex++)
	{
		N = Tree->NodeList[NIndex];
		if(N->Tip == FALSE)
			GeoForceUpDateNode(N, Rates, Rates->RS);
	}
			
	FatTailGetAnsSates(Tree, Trees->NoOfSites, FTR);
	
	Rates->Lh = Likelihood(Rates, Trees, Opt);

	Rates->AutoAccept = TRUE;
}

/*
void	GeoUpDateAnsStates(OPTIONS *Opt, TREES *Trees, RATES *Rates)
{
	int NIndex;

	TREE *Tree;
	FATTAILRATES *FTR;
	NODE N;

	Tree = Trees->Tree[Rates->TreeNo];

	FTR = Rates->FatTailRates;

	MapRatesToFatTailRate(Rates, FTR);

	FatTailSetAnsSates(Tree, Trees->NoOfSites, FTR);

	SetStableDist(FTR->SDList[0], FTR->Alpha[0], FTR->Scale[0]);

	do
	{
		NIndex = RandUSInt(Rates->RS) % Tree->NoNodes;
		N = Tree->NodeList[NIndex];
	} while(N->Tip == TRUE);

	GeoUpDateNode(N, Rates, Rates->RS);
	
	FatTailGetAnsSates(Tree, Trees->NoOfSites, FTR);
}
*/

void	GeoUpDateAnsStates(OPTIONS *Opt, TREES *Trees, RATES *Rates)
{
	int NIndex;
	TREE *Tree;
	FATTAILRATES *FTR;
	NODE N;

	LhTransformTree(Rates, Trees, Opt);

	Tree = Trees->Tree[Rates->TreeNo];
	
	FTR = Rates->FatTailRates;

	MapRatesToFatTailRate(Rates, FTR);

	FatTailSetAnsSates(Tree, Trees->NoOfSites, FTR);

	SetStableDist(FTR->SDList[0], FTR->Alpha[0], FTR->Scale[0]);

	do
	{
		NIndex = RandUSInt(Rates->RS) % Tree->NoNodes;
		N = Tree->NodeList[NIndex];
	} while(N->Tip == TRUE);

	GeoUpDateNode(N, Rates, Rates->RS);

	FatTailGetAnsSates(Tree, Trees->NoOfSites, FTR);
	Rates->CalcLh = FALSE;
}

/*
void	GeoUpDateAllAnsStates(OPTIONS *Opt, TREES *Trees, RATES *Rates)
{
	int PIndex, NIndex, TNo;
	TREE *Tree;
	FATTAILRATES *FTR;
	NODE N;
	
	Tree = Trees->Tree[Rates->TreeNo];

	FTR = Rates->FatTailRates;

	FatTailSetAnsSates(Tree, Trees->NoOfSites, FTR);
	SetStableDist(FTR->SDList[0], FTR->Alpha[0], FTR->Scale[0]);
		
	for(PIndex=0;PIndex<Tree->NoFGroups;PIndex++)
	{
//		#pragma omp parallel for num_threads(Opt->Cores) private(TNo, N) schedule(dynamic, 1)
		for(NIndex=0;NIndex<Tree->NoFNodes[PIndex];NIndex++)
		{
			TNo = GetThreadNo();
			N = Tree->FNodes[PIndex][NIndex];

			GeoUpDateNode(N, Rates, Rates->RSList[TNo]);
		}

	}

	FatTailGetAnsSates(Tree, Trees->NoOfSites, FTR);
}
*/

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

	printf("Sec:\t%f\t%d\n", S, i);

	printf("%f\t%f\t%f\n", x, y, z);

	exit(0);
}

void	SetLoadGeoData(char **AnsState, TREES *Trees, TREE *Tree)
{
	int Index, Pos;
	double Long, Lat;
	NODE N;
	int NoInt;

	Pos = 0;

	NoInt = 0;
	for(Index=0;Index<Tree->NoNodes;Index++)
	{
		N = Tree->NodeList[Index];
		if(N->Tip == FALSE)
		{
			NoInt++;
			Long = atof(AnsState[Pos++]);
			Lat = atof(AnsState[Pos++]);

			LongLatToNode(N, Long, Lat);
		}
	}
}


void	LoadGeoDataTest(OPTIONS *Opt, TREES *Trees, RATES *CRates)
{
	double Alpha, Lh;


	for(Alpha=0;Alpha<2.0;Alpha+=0.001)
	{
		CRates->Rates[0] = Alpha;
		Lh = Likelihood(CRates, Trees, Opt);
		printf("%f\t%f\n", Alpha, Lh);
	}

	exit(0);
}

void	LoadGeoData(OPTIONS *Opt, TREES *Trees, RATES *CRates, char *Str)
{
	TREE *Tree;
	FATTAILRATES *FTR;
	char *S, **Passed;
	int Tokes;

	S = StrMake(Str);
	
	Tree = Trees->Tree[CRates->TreeNo];

	FTR = CRates->FatTailRates;

	Passed = (char**)SMalloc(sizeof(char*) * (strlen(S) + 1));

	Tokes = MakeArgv(S, Passed, (int)(strlen(S) + 1));

	FTR->Alpha[0] = atof(Passed[2]);
	FTR->Scale[0] = atof(Passed[3]);
		
	SetLoadGeoData(&Passed[4], Trees, Tree);

	MapFatTailRateToRates(CRates, FTR);

	FatTailGetAnsSates(Tree, Trees->NoOfSites, FTR);

	CRates->Lh = Likelihood(CRates, Trees, Opt);

	printf("checkpoint lh:\t%f\n", CRates->Lh);
	

	free(Passed);
	free(S);
}