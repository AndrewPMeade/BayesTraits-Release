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
#include <string.h>
#include <math.h>

#include "TypeDef.h"
#include "Geo.h"
#include "GenLib.h"
#include "Threaded.h"
#include "FatTail.h"
#include "Likelihood.h"
#include "Rates.h"
#include "DistData.h"
#include "RestrictionMap.h"
#include "IntraNode.h"
#include "Trees.h"


#ifndef M_PI
	#define M_PI 3.14159265358979323846
#endif

#define EARTH_KM 6371.0
//#define EARTH_KM 1.0

#define EARTH_KM_EARTH_KM 40589641

int		ValidLongLat(double Long, double Lat)
{
	if(Long < -180 || Long >= 180)
		return FALSE;

	if(Lat < -90 || Lat >= 90)
		return FALSE;

	return TRUE;
}

int		ValidGeoNodeXYZ(NODE Node, double RX, double RY, double RZ)
{
	double Long, Lat;

	if(Node->NodeResMap == NULL)
		return TRUE;

	XYZToLongLat(RX, RY, RZ, &Long, &Lat);

	return ValidNodeResPoint(Node->NodeResMap, Long, Lat);
}

int		ValidGeoNodeLongLat(NODE Node, double Long, double Lat)
{
	if(ValidLongLat(Long, Lat) == FALSE)
		return FALSE;

	if(Node->NodeResMap == NULL)
		return TRUE;


	return ValidNodeResPoint(Node->NodeResMap, Long, Lat);
}

void	ValidGeoData(TREES *Trees)
{
	int Index;
	double Long, Lat;
	TAXA *T;

	if(Trees->NoSites != 2)
	{
		printf("Geo Data must be long lat\n");
		exit(1);
	}

	for(Index=0;Index<Trees->NoTaxa;Index++)
	{
		T = Trees->Taxa[Index];

		if(T->Exclude == FALSE)
		{
			Long = T->ConData[0];
			Lat = T->ConData[1];

			if(ValidLongLat(Long, Lat) == FALSE)
			{
				printf("Invalid longitude (%f) or latitude (%f) for taxa %s.\n", Long, Lat, T->Name);
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

	for(TIndex=0;TIndex<Trees->NoTaxa;TIndex++)
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
	Trees->NoSites = 3;
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

void	GetRandLongLatCore(gsl_rng *RNG, double Long, double Lat, double *RLong, double *RLat, double Radius)
{
	double Dist, Bearing;

	Long = Long * M_PI / 180.0;
	Lat = Lat * M_PI / 180.0;

	Dist = gsl_rng_uniform(RNG) * Radius;

	Bearing	= gsl_rng_uniform(RNG) * 360;
	Bearing = Bearing * M_PI / 180.0;
	
	Bearing = fmod(Bearing + (2.0 * M_PI), (2.0 * M_PI));
	
	*RLat = sin(Lat) * cos(Dist/EARTH_KM) + cos(Lat) * sin(Dist/EARTH_KM) * cos(Bearing);
	*RLat = asin(*RLat);

	*RLong = Long + atan2(sin(Bearing) * sin(Dist/EARTH_KM) * cos(Lat), cos(Dist/EARTH_KM) - sin(Lat) * sin(*RLat));
	
	*RLat = *RLat * (180.0 / M_PI);
		
	if(*RLong * (180.0 / M_PI) > 180)
		*RLong = *RLong * (180.0 / M_PI) - 360.0;
	else
	{
		*RLong = *RLong * (180.0 / M_PI);
		if (*RLong < -180)
			*RLong = *RLong + 360;
	}
}

void	GetRandLongLat(gsl_rng* RNG, double Long, double Lat, double* RLong, double* RLat, double Radius)
{
	do
		GetRandLongLatCore(RNG, Long, Lat, RLong, RLat, Radius);
	while(ValidLongLat(*RLong, *RLat) == FALSE);
}

void	GetGloablRandLongLat(gsl_rng *RNG, double *RLong, double *RLat)
{
	*RLong = gsl_ran_flat(RNG, -179.999999, 180.0);
	*RLat = gsl_ran_flat(RNG, -89.999999, 90.0);
}

void RRest(NODE Node, double Long, double Lat, double Radius)
{
	int Index;
	double RLong, RLat;

	printf("Long:\tLat:\n");
	printf("%f\t%f\n\n\n", Long, Lat);

	for(Index=0;Index<10000;Index++)
	{
		GetRandLongLat(Node->RNG, Long, Lat, &RLong, &RLat, Radius);
		printf("%f\t%f\n", RLong, RLat);
	}

	exit(0);
}

void	GetRandXYZPoint(NODE Node, double SX, double SY, double SZ, double *RX, double *RY, double *RZ, double Radius)
{
	double Long, Lat, RLong, RLat;
	int Tried;

	XYZToLongLat(SX, SY, SZ, &Long, &Lat);
	
	Tried = 0;
	do
	{
		if(Tried < 100)
			GetRandLongLat(Node->RNG, Long, Lat, &RLong, &RLat, Radius);
		else
			GetGloablRandLongLat(Node->RNG, &RLong, &RLat);

		Tried++;
	} while(ValidGeoNodeLongLat(Node, RLong, RLat) == FALSE);
	
	LongLatToXYZ(RLong, RLat, RX, RY, RZ);
}

double	GeoCalcAnsStateLh(double X, double Y, double Z, NODE N, double Scale)
{
	double Ret, Val;
	int Index;

	Ret = 0;

	for(Index=0;Index<N->NoNodes;Index++)
	{
		Val = X - N->NodeList[Index]->FatTailNode->Ans[0];
		Ret += StableDistTPDF(Scale, Val, N->NodeList[Index]->Length);

		Val = Y - N->NodeList[Index]->FatTailNode->Ans[1];
		Ret += StableDistTPDF(Scale, Val, N->NodeList[Index]->Length);

		Val = Z - N->NodeList[Index]->FatTailNode->Ans[2];
		Ret += StableDistTPDF(Scale, Val, N->NodeList[Index]->Length);
	}

	if(N->Ans != NULL)
	{
		Val = X - N->Ans->FatTailNode->Ans[0];
		Ret += StableDistTPDF(Scale, Val, N->Length);

		Val = Y - N->Ans->FatTailNode->Ans[1];
		Ret += StableDistTPDF(Scale, Val, N->Length);

		Val = Z - N->Ans->FatTailNode->Ans[2];
		Ret += StableDistTPDF(Scale, Val, N->Length);
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




void	GeoUpDateNode(NODE N, RATES *Rates)
{
	FATTAILRATES *FTR;
	double NLh, CLh, X, Y, Z, RX, RY, RZ;
	double Scale;
	

	FTR = Rates->FatTailRates;
	Scale = Rates->Rates[0];	

	NodeToXYZ(N, &X, &Y, &Z);
	
	CLh = GeoCalcAnsStateLh(X, Y, Z, N, Scale );

	do
	{
		GetRandXYZPoint(N, X, Y, Z, &RX, &RY, &RZ, 2000);
	} while(ValidGeoNodeXYZ(N, RX, RY, RZ) == FALSE);
		
	NLh = GeoCalcAnsStateLh(RX, RY, RZ, N, Scale);
	
	XYZToNode(N, RX, RY, RZ);

	Rates->Lh = Rates->Lh + (NLh - CLh);
}

void PrintGeoNodeInfo(NODE N, double Lh)
{
	double X,Y,Z;

	NodeToXYZ(N, &X, &Y, &Z);


	printf("%f\t%f\t%f\t%f\t|\t", Lh, X, Y, Z);
}


void	GeoForceUpDateNode(NODE N, RATES *Rates)
{
	FATTAILRATES *FTR;
	double NLh, CLh, X, Y, Z, RX, RY, RZ;
	int Changed, Tried;
	double Scale;

	FTR = Rates->FatTailRates;
	Scale = Rates->Rates[0];

	NodeToXYZ(N, &X, &Y, &Z);

	CLh = GeoCalcAnsStateLh(X, Y, Z, N, Scale);

	Changed = FALSE;
	Tried = 0;
	do
	{
		GetRandXYZPoint(N, X, Y, Z, &RX, &RY, &RZ, 2000);

		NLh = GeoCalcAnsStateLh(RX, RY, RZ, N, Scale);
		
		if(log(gsl_rng_uniform_pos(N->RNG)) < (NLh - CLh))
			Changed = TRUE;
		else
			Tried++;

	} while(Changed == FALSE && Tried < 1);

	
	if(Changed == TRUE)
		XYZToNode(N, RX, RY, RZ);

//	Rates->Lh = Rates->Lh + (NLh - CLh);
}

void	PrintGeoData(OPTIONS *Opt, TREES *Trees, RATES *Rates)
{
	int Index;
	TREE *Tree;
	NODE N;
	double X,Y,Z;

	Tree = Trees->Tree[0];

	for(Index=0;Index<Tree->NoNodes;Index++)
	{
		N = Tree->NodeList[Index];

		NodeToXYZ(N, &X, &Y, &Z);
		printf("%d\t%d\t%f\t%f\t%f\t%f\n", Index, N->Part->NoTaxa, N->Length, X, Y, Z);
	}

}

void	GeoUpDateAllAnsStates(OPTIONS *Opt, TREES *Trees, RATES *Rates)
{
	int NIndex, GIndex;
	FATTAILTREE *FTT;
	size_t	GroupSize;
	NODE	*NodeList;
	TREE *Tree;
	FATTAILRATES *FTR;

	CheckRestictedMaps(Trees, Rates->FatTailRates);

	LhTransformTree(Rates, Trees, Opt);


	Tree = Trees->Tree[Rates->TreeNo];
	FTT = Tree->FatTailTree;

	FTR = Rates->FatTailRates;

	FatTailSetAnsSates(Tree, Trees->NoSites, FTR);
	
	if(Opt->UseDistData == TRUE)
		SetTreeDistData(Rates, Opt, Trees);


	if(Opt->UseIntraNode == TRUE)
		ChangeAllIntraNodeLh(Opt, Trees, Rates);
	else
	{
		for(GIndex=0;GIndex<FTT->NoParallelGroups;GIndex++)
		{
			GroupSize = FTT->ParallelNodeListLength[GIndex];
			NodeList = FTT->ParallelNodeList[GIndex];

	#ifdef OPENMP_THR
	// Dynamic slows things down oddly and can leve long long running process. 
	//	#pragma omp parallel for num_threads(Opt->Cores) schedule(dynamic)
		#pragma omp parallel for num_threads(Opt->Cores) schedule(static)
	#endif
			for(NIndex=0;NIndex<GroupSize;NIndex++)
				GeoForceUpDateNode(NodeList[NIndex], Rates);
		}
	}

	FatTailGetAnsSates(Tree, Trees->NoSites, FTR);

	CheckRestictedMaps(Trees, FTR);

	Rates->AutoAccept = TRUE;
}


void	GeoUpDateAnsStates(OPTIONS *Opt, TREES *Trees, RATES *Rates)
{
	int NIndex;
	TREE *Tree;
	FATTAILRATES *FTR;
	NODE N;

	LhTransformTree(Rates, Trees, Opt);

	Tree = Trees->Tree[Rates->TreeNo];
	
	FTR = Rates->FatTailRates;
		
	FatTailSetAnsSates(Tree, Trees->NoSites, FTR);

	do
	{
		NIndex = (int)gsl_rng_uniform_int(Rates->RNG, Tree->NoNodes);
		N = Tree->NodeList[NIndex];
	} while(N->Tip == TRUE);

	GeoUpDateNode(N, Rates);

	FatTailGetAnsSates(Tree, Trees->NoSites, FTR);
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

void	SetValidNode(NODE Node)
{
	double			X,Y,Z, Long, Lat;
	double			RX, RY, RZ;

	NodeToXYZ(Node, &X, &Y, &Z);
	XYZToLongLat(X, Y, Z, &Long, &Lat);
	
	while(ValidGeoNodeLongLat(Node, Long, Lat) == FALSE)
	{
		GetRandXYZPoint(Node, X, Y, Z, &RX, &RY, &RZ, 2000);
		XYZToLongLat(RX, RY, RZ, &Long, &Lat);
	}

	LongLatToNode(Node, Long, Lat);
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

			SetValidNode(N);
		}
	}
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
		
	SetLoadGeoData(&Passed[4], Trees, Tree);
		
	FatTailGetAnsSates(Tree, Trees->NoSites, FTR);

	CRates->Lh = Likelihood(CRates, Trees, Opt);

	printf("checkpoint lh:\t%f\n", CRates->Lh);
	

	free(Passed);
	free(S);
}


void	CheckRestictedMaps(TREES *Trees, FATTAILRATES *FTR)
{
	TREE *Tree;
	int Index, Valid;
	NODE Node;

	double X,Y,Z, Long, Lat;

	Tree = Trees->Tree[0];

	FatTailSetAnsSates(Tree, Trees->NoSites, FTR);

	for(Index=0;Index<Tree->NoNodes;Index++)
	{
		Node = Tree->NodeList[Index];
		if(Node->Tip == FALSE)
		{
			NodeToXYZ(Node, &X, &Y, &Z);

			XYZToLongLat(X, Y, Z, &Long, &Lat);

			Valid = ValidGeoNodeLongLat(Node, Long, Lat);


			if(Valid == FALSE)
			{
				printf("%d\t%d\t%d\t%f\t%f\n", Index, Valid, Node->Tip, Long, Lat);
				printf("Invalid point.\n");
				exit(1);
			}
		}
	}
}

void	RecSimGeoData(NODE Node, RATES *Rates, double Scale)
{
	double X, Y, Z;
	double Long, Lat;
	int Index;

	X = Node->Ans->FatTailNode->Ans[0];
	Y = Node->Ans->FatTailNode->Ans[1];
	Z = Node->Ans->FatTailNode->Ans[2];

	X += gsl_ran_gaussian(Rates->RNG, sqrt(Scale * Node->Length));
	Y += gsl_ran_gaussian(Rates->RNG, sqrt(Scale * Node->Length));
	Z += gsl_ran_gaussian(Rates->RNG, sqrt(Scale * Node->Length));

	XYZToLongLat(X, Y, Z, &Long, &Lat);
	LongLatToXYZ(Long, Lat, &X, &Y, &Z);

	Node->FatTailNode->Ans[0] = X;
	Node->FatTailNode->Ans[1] = Y;
	Node->FatTailNode->Ans[2] = Z;

	for(Index=0;Index<Node->NoNodes;Index++)
		RecSimGeoData(Node->NodeList[Index], Rates, Scale);
}

void	PrintSimGeoData(TREE *Tree)
{
	int Index;
	double Long, Lat;
	NODE Node;

	printf("\n\n\nAncestral states\n\n\n");

	for(Index=0;Index<Tree->NoNodes;Index++)
	{
		Node = Tree->NodeList[Index];

		XYZToLongLat(Node->FatTailNode->Ans[0], Node->FatTailNode->Ans[1], Node->FatTailNode->Ans[2], &Long, &Lat);

		printf("%f\t%f\t", Long, Lat);
		RecPRintNodeTaxa(Node, ',');
		printf("\n");
	}

	printf("\n\n\n\nTrait Data.\n\n\n");

	for(Index=0;Index<Tree->NoNodes;Index++)
	{
		Node = Tree->NodeList[Index];

		XYZToLongLat(Node->FatTailNode->Ans[0], Node->FatTailNode->Ans[1], Node->FatTailNode->Ans[2], &Long, &Lat);

		if(Node->Tip == TRUE)
		{
			printf("%s\t%f\t%f\t", Node->Taxa->Name, Long, Lat);
			printf("\n");
		}
	}
}

void	SimGeoData(OPTIONS *Opt, TREES *Trees, RATES *Rates)
{
	double X, Y, Z;
	double Scale;
	TREE *Tree;
	NODE Root;
	FATTAILRATES *FTR;
	int Index;
	
//	return;
	Tree = Trees->Tree[0];
	FTR = Rates->FatTailRates;
	Scale = 500000;
//	Scale = 3;

	FatTailSetAnsSates(Tree, Trees->NoSites, FTR);

	Root = Tree->Root;

	LongLatToXYZ(0.0, 0.0, &X, &Y, &Z);

	Root->FatTailNode->Ans[0] = X;
	Root->FatTailNode->Ans[1] = Y;
	Root->FatTailNode->Ans[2] = Z;

	for(Index=0;Index<Root->NoNodes;Index++)
		RecSimGeoData(Root->NodeList[Index], Rates, Scale);

	PrintSimGeoData(Tree);
	exit(0);
//	GetGloablRandLongLat(Rates->RNG, &RLong, &RLat);
//	XYZToLongLat(X, Y, Z, &Long, &Lat);
//	LongLatToXYZ(Long, Lat, &X, &Y, &Z);

}