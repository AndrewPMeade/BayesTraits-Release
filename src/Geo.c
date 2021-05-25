#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "typedef.h"
#include "Geo.h"
#include "genlib.h"
#include "threaded.h"
#include "FatTail.h"
#include "likelihood.h"
#include "rates.h"

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
	double	NLh, CLh, X, Y, Z, RX, RY, RZ;
	STABLEDIST *SD;
	

	FTR = Rates->FatTailRates;
	SD = FTR->SDList[0];

	NodeToXYZ(N, &X, &Y, &Z);
	
	CLh = GeoCalcAnsStateLh(X, Y, Z, N, SD);

//	if(CLh > 0)
//		printf("Err\n");
		
//	GetRandXYZPoint(RS, X, Y, Z, &RX, &RY, &RZ, 1000);
	GetRandXYZPoint(RS, X, Y, Z, &RX, &RY, &RZ, 100);
		
	NLh = GeoCalcAnsStateLh(RX, RY, RZ, N, SD);
	

	XYZToNode(N, RX, RY, RZ);

	Rates->Lh = Rates->Lh + (NLh - CLh);
	
//	printf("%f\t%f\t%f\n", CLh, NLh, (NLh - CLh));


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

void	SetLoadGeoData(double *Vect, TREES *Trees, TREE *Tree)
{
	int Index, SIndex, Pos;
	NODE N;

	Pos = 0;

	for(Index=0;Index<Tree->NoNodes;Index++)
	{
		N = Tree->NodeList[Index];
		if(N->Tip == FALSE)
		{
			for(SIndex=0;SIndex<Trees->NoOfSites;SIndex++)
				N->FatTailNode->Ans[SIndex] = Vect[Pos++];
		}
	}
}

void	LoadGeoDataTest(OPTIONS *Opt, TREES *Trees, RATES *CRates, RATES *NRates)
{
	int Index;

	
	for(Index=0;Index<10000;Index++)
	{
		CopyRates(NRates, CRates, Opt);
		GeoUpDateAllAnsStates(Opt, Trees, NRates);
	}

	exit(0);
}

void	LoadGeoData(OPTIONS *Opt, TREES *Trees, RATES *CRates, RATES *NRates)
{
	int Index;
	TREE *Tree;
	//double Vect[] = {0.200161,0.053138,753.169564,-3436.161322,5311.795551,753.169982,-3436.161073,5311.795653,753.169522,-3436.161057,5311.795728,753.177014,-3436.158419,5311.796373,784.523995,-3268.579972,5412.074285,753.179424,-3436.161622,5311.793959,753.119873,-3436.202718,5311.775818,753.119997,-3436.200966,5311.776934,753.306416,-3436.556721,5311.520343,617.180851,-4439.383553,4527.7591,806.044111,-3266.771262,5410.003643,617.126734,-4439.396095,4527.754179,617.120046,-4439.397077,4527.754128,617.08903,-4439.437188,4527.719027,617.332621,-4439.364255,4527.757331,617.372624,-4439.358731,4527.757292,617.372944,-4439.407372,4527.709557,632.608183,-3712.834353,5138.512329,632.716539,-3712.766949,5138.54769,629.624153,-3712.571715,5139.068562,629.600841,-3712.572437,5139.070897,629.600857,-3712.572319,5139.07098};
	FATTAILRATES *FTR;
	double Vect[] = {0.259795,0.066530,751.378629,-3437.069454,5311.461637,753.258104,-3436.195650,5311.760789,753.264434,-3436.198048,5311.758340,753.254442,-3436.195880,5311.761160,784.523499,-3268.580864,5412.073818,753.249728,-3436.196573,5311.761380,753.247693,-3436.199677,5311.759660,684.119337,-3445.857367,5314.855477,685.966751,-3449.193539,5312.452781,1174.404203,-4308.824310,4543.616272,912.632312,-3308.485778,5367.556718,1174.440167,-4308.836315,4543.595592,1174.438958,-4308.836987,4543.595267,616.839183,-4439.468506,4527.722364,1174.325102,-4308.200488,4544.228219,629.148487,-3194.496934,5476.221555,857.173841,-3352.341239,5349.458124,629.094699,-3194.502862,5476.224277,629.012338,-3187.840835,5480.114532,629.569706,-3712.553591,5139.088326,629.550511,-3712.552243,5139.091651,629.549949,-3712.582582,5139.069802};
	return;
	Tree = Trees->Tree[CRates->TreeNo];

	FTR = CRates->FatTailRates;

//	Rates->Rates[0] = Vect[0];
//	Rates->Rates[1] = Vect[1];

	FTR->Alpha[0] = Vect[0];
	FTR->Scale[0] = Vect[1];

//	memcpy(Tree->FatTailTree->AnsVect, &Vect[2], sizeof(double) * Tree->NoNodes * 3);


	SetLoadGeoData(&Vect[2], Trees, Tree);
	MapFatTailRateToRates(CRates, FTR);

	FatTailGetAnsSates(Tree, 3, FTR);

	CRates->Lh = Likelihood(CRates, Trees, Opt);
	
	printf("Init Lh:\t%f\n", CRates->Lh);

	LoadGeoDataTest(Opt, Trees, CRates, NRates);

	exit(0);
}