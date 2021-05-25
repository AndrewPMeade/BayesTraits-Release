#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "typedef.h"
#include "rates.h"
#include "genlib.h"
#include "RandLib.h"
#include "trees.h"
#include "continuous.h"
#include "revjump.h"

void	RJSplitPropRatio(RATES* Rates, OPTIONS* Opt, MAPINFO *MapInfo, int G0Size, int G1Size, double Rate, int OGroup);


int		FindGroupsGrater(int No, MAPINFO *MapInfo)
{
	int Ret;
	int Index;

	Ret=0;
	for(Index=0;Index<MapInfo->NoOfGroups;Index++)
		if(MapInfo->GroupSize[Index] > No)
			Ret++;

	return Ret;
}


void	MapRJRates(double* Rates, int *MappingVect, int MVLen, double *Prams)
{
	int	Index;
	int	Hit=FALSE;

	for(Index=0;Index<MVLen;Index++)
	{

		if(MappingVect[Index] == ZERORATENO)
		{
			Prams[Index] = 0.000000000000000000001;
			Hit = TRUE;

		}
		else
			Prams[Index] = Rates[MappingVect[Index]];
	}
/*
	if(Hit == TRUE)
	{

		for(Index=0;Index<MVLen;Index++)
			printf("%d\t", MappingVect[Index]);

		printf("\n");
		for(Index=0;Index<MVLen;Index++)
			printf("%f\t", Prams[Index]);

		exit(0);
	}
*/
}

int		NoOfMappings(RATES* Rates, int	RateNo)
{
	int	Ret=0;
	int	Index;

	for(Index=0;Index<Rates->NoOfFullRates;Index++)
		if(Rates->MappingVect[Index] == RateNo)
			Ret++;

	return Ret;
}

int		NoOfPramGroups(RATES* Rates, int *GroupID, int *GroupSize)
{
	int	Index=0;
	int	Ret=0;
	int	InPast;
	int	PIndex;

	for(Index=0;Index<Rates->NoOfFullRates;Index++)
	{
		InPast=FALSE;
		for(PIndex=0;PIndex<(Index);PIndex++)
		{
			if(Rates->MappingVect[Index] == Rates->MappingVect[PIndex])
				InPast = TRUE;
		}

		/* If we are in the zero bin don't count us */
		if(Rates->MappingVect[Index] == ZERORATENO)
			InPast = TRUE;

		if((NoOfMappings(Rates,Rates->MappingVect[Index]) >= 1) && (InPast == FALSE))
		{
			if(!((GroupID == NULL) || (GroupSize == NULL)))
			{
				GroupID[Ret] = Rates->MappingVect[Index];
				GroupSize[Ret] = NoOfMappings(Rates,Rates->MappingVect[Index]);
			}

			Ret++;
		}
	}

	return Ret;
}

void		PrintMapInfo(MAPINFO *MapInfo)
{
	int	Index=0;
	int	PIndex;

	printf("No Of Groups:\t%d\n", MapInfo->NoOfGroups);
	for(Index=0;Index<MapInfo->NoOfGroups;Index++)
	{
		printf("Group\t%d\t%d\t", MapInfo->GroupID[Index], MapInfo->GroupSize[Index]);

		for(PIndex=0;PIndex<MapInfo->GroupSize[Index];PIndex++)
			printf("%d\t", MapInfo->GroupPos[Index][PIndex]);

		printf("\n");
	}

	printf("Zero Group:\t%d\t", MapInfo->NoInZero);
	for(Index=0;Index<MapInfo->NoInZero;Index++)
		printf("%d\t", MapInfo->ZeroPos[Index]);
	printf("\n");

	printf("\n");
}

void		FreeMapInfo(MAPINFO *MapInfo)
{
	int	Index;

	for(Index=0;Index<MapInfo->NoOfGroups;Index++)
		free(MapInfo->GroupPos[Index]);
	free(MapInfo->GroupPos);
	
	free(MapInfo->GroupSize);
	free(MapInfo->GroupID);

	if(MapInfo->NoInZero > 0)
		free(MapInfo->ZeroPos);

	free(MapInfo);
}

void	FindZeroBin(RATES* Rates, MAPINFO* MapInfo)
{
	int		Index;
	int		ZVectIndex;

	MapInfo->NoInZero = 0;

	for(Index=0;Index<Rates->NoOfFullRates;Index++)
		if(Rates->MappingVect[Index] == ZERORATENO)
			MapInfo->NoInZero++;

	if(MapInfo->NoInZero > 0)
	{
		MapInfo->ZeroPos = (int*)malloc(sizeof(int) * MapInfo->NoInZero);
		if(MapInfo->ZeroPos == NULL)
			MallocErr();
	}
	else
		MapInfo->ZeroPos = NULL;

	ZVectIndex = 0;
	for(Index=0;Index<Rates->NoOfFullRates;Index++)
	{
		if(Rates->MappingVect[Index] == ZERORATENO)
		{
			MapInfo->ZeroPos[ZVectIndex] = Index;
			ZVectIndex++;
		}
	}
}

MAPINFO*	CreatMapping(RATES* Rates)
{
	MAPINFO*	Ret;
	int			Index;
	int			PIndex;
	int			Pos;

	Ret = (MAPINFO*)malloc(sizeof(MAPINFO));
	if(Ret == NULL)
		MallocErr();

	Ret->NoOfGroups = NoOfPramGroups(Rates, NULL, NULL);

	Ret->GroupSize	= (int*)malloc(sizeof(int)* Ret->NoOfGroups);
	Ret->GroupID	= (int*)malloc(sizeof(int)* Ret->NoOfGroups);
	if((Ret->GroupSize == NULL) || (Ret->GroupID == NULL))
		MallocErr();

	NoOfPramGroups(Rates, Ret->GroupID, Ret->GroupSize);

	Ret->GroupPos = (int**)malloc(sizeof(int*) * Ret->NoOfGroups);
	if(Ret->GroupPos == NULL)
		MallocErr();

	for(Index=0;Index<Ret->NoOfGroups;Index++)
	{
		Ret->GroupPos[Index] = (int*)malloc(sizeof(int) * Ret->GroupSize[Index]);
		if(Ret->GroupPos[Index] == NULL)
			MallocErr();

		Pos = 0;
		for(PIndex=0;PIndex<Rates->NoOfFullRates;PIndex++)
		{
			if(Rates->MappingVect[PIndex] == Ret->GroupID[Index])
			{
				Ret->GroupPos[Index][Pos] = PIndex;
				Pos++;
			}

		}
	}

	FindZeroBin(Rates, Ret);

	return Ret;
}

int		PickSplitGroup(RATES *Rates, MAPINFO *MapInfo)
{
	int	Ret;

	do
	{
		Ret = RandUSLong(Rates->RS) % MapInfo->NoOfGroups;
	}while(MapInfo->GroupSize[Ret] <= 1);

	return Ret;
}

int*	MakeSplitMask(RATES* Rates, MAPINFO *MapInfo, int GroupNo)
{
	int	*Ret=NULL;
	int	Index;
	int	G0,G1;
	int	GSize;

	GSize = MapInfo->GroupSize[GroupNo];

	Ret = (int*)malloc(sizeof(int) * GSize);
	if(Ret==NULL)
		MallocErr();

	G0 = G1 = 0;
	for(Index=0;Index<GSize-1;Index++)
	{
		if(RandDouble(Rates->RS) < 0.5)
		{
			Ret[Index] = 0;
			G0++;
		}
		else
		{
			Ret[Index] = 1;
			G1++;
		}
	}

	if((G0 > 0) && (G1 > 0))
	{
		if(RandDouble(Rates->RS) < 0.5)
			Ret[Index] = 0;
		else
			Ret[Index] = 1;
	}
	else
	{
		if(G0 == 0)
			Ret[Index] = 0;
		
		if(G1 == 0)
			Ret[Index] = 1;
	}

	return Ret;
}

double	GenDoubleIn(RANDSTATES*	RS, double Min, double Max)
{
	double	Diff;
	
	Diff = Max - Min;
	return (RandDouble(RS) * Diff) + Min;
}


void	RJSplit(RATES* Rates, OPTIONS* Opt)
{
	MAPINFO *MapInfo=NULL;
	int		GroupNo;
	int		SplitID;
	int		*SplitMask;
	int		NewID;
	double	*NewRates;
	int		G0, G1;
	int		Index;
	int		GSize;
	int		Pos;
	double	NewR;
	double	OldR;
	double	Mue;


	MapInfo = CreatMapping(Rates);

	if(FindGroupsGrater(1, MapInfo)  == 0)
	{
		FreeMapInfo(MapInfo);
		return;
	}

	GroupNo = PickSplitGroup(Rates, MapInfo);
	SplitID	= MapInfo->GroupID[GroupNo];

	GSize = MapInfo->GroupSize[GroupNo];
	SplitMask = MakeSplitMask(Rates, MapInfo, GroupNo);

	G0=G1=0;
	for(Index=0;Index<GSize;Index++)
	{
		if(SplitMask[Index] == 0)
			G0++;
		if(SplitMask[Index] == 1)
			G1++;
	}


	OldR = Rates->Rates[Rates->MappingVect[MapInfo->GroupPos[GroupNo][0]]];

	RJSplitPropRatio(Rates, Opt, MapInfo, G0, G1, OldR, SplitID);
/* Old
	Mue = GenDoubleIn(0, (G0 + G1) * OldR);
	NewR = (((G0 + G1) * OldR) - Mue)/(double)G1;
	OldR = Mue / (double)G0;
	*/

	Mue = GenDoubleIn(Rates->RS, -(G0*OldR), (G1*OldR));
	NewR = OldR + (Mue / (double)G0);
	OldR = OldR - (Mue / (double)G1);

	NewRates = (double*)malloc(sizeof(double) * (Rates->NoOfRates + 1));
	if(NewRates == NULL)
		MallocErr();

	memcpy(NewRates, Rates->Rates, sizeof(double)*Rates->NoOfRates);
	
	NewRates[Rates->MappingVect[MapInfo->GroupPos[GroupNo][0]]] = OldR;
	NewRates[Rates->NoOfRates] = NewR;

	NewID = MapInfo->NoOfGroups;
	for(Index=0;Index<GSize;Index++)
	{
		if(SplitMask[Index] == 1)
		{
			Pos = MapInfo->GroupPos[GroupNo][Index];
			Rates->MappingVect[Pos] = NewID;
		}
	}

	free(Rates->Rates);
	Rates->Rates = NewRates;

	Rates->NoOfRates = MapInfo->NoOfGroups + 1;

	free(SplitMask);
	FreeMapInfo(MapInfo);
}

/* Calcl the probablity of picking any particular split		*/
/* T1 the P of picking the split operator					*/
/* T2 the P of picking the "correct" group to split			*/ 
/* T3 the P of splitting the group corectly					*/
/* T4 the P of creating the correct rate after the split	*/
double	RJSplitRatio(RATES* Rates, OPTIONS* Opt, MAPINFO *MapInfo, int SizeG0, int SizeG1, double Rate, double Extra)
{
	double	NoSplitable;
	int		Index;
	double	T1,T2,T3,T4;
	int		NiNj;

	NiNj = SizeG0 + SizeG1;

	if(MapInfo->NoOfGroups + 1== Rates->NoOfFullRates)
		T1 = 1;
	else
		T1 = 0.5;

	NoSplitable = 0.0;
	for(Index=0;Index<MapInfo->NoOfGroups;Index++)
	{
		if(MapInfo->GroupSize[Index] > 1)
			NoSplitable++;
	}
	
	NoSplitable = NoSplitable + Extra;

	T2 = (1.0 / NoSplitable);


	T3 = (1.0 / ( pow(2, NiNj - 1) - 1.0));
 
	
	T4 = (1.0 / (Rate * (NiNj)));

	return T1 * T2 * T3 *T4;
}

double	RJMergeRatio(int NoOfGroups)
{
	double	Ret;
	int		KC2;

	if(NoOfGroups - 1 == 1)
		Ret = 1;
	else
		Ret = 0.5;

	KC2 = NoOfGroups * (NoOfGroups - 1);
	KC2 = KC2 / 2;
	Ret = Ret * (1.0 / (double)KC2);


	return Ret;
}

void	RJMergePropRatio(RATES* Rates, OPTIONS* Opt, MAPINFO *MapInfo, int G0ID, int G1ID, double Rate)
{
	double	Top;
	double	Bottom;
	double	Ret;
	int		SizeG0;
	int		SizeG1;
/*	int		Index; */

	SizeG0 = MapInfo->GroupSize[G0ID];
	SizeG1 = MapInfo->GroupSize[G1ID];
/*
	printf("Split: group: %d\t%d\tni %d\tnj %d\tRate %f\n", G0ID, G1ID, SizeG0, SizeG1, Rate);

	for(Index=0;Index<Rates->NoOfRates;Index++)
		printf("Rate %f\n", Rates->Rates[Index]);

	for(Index=0;Index<Rates->NoOfFullRates;Index++)
		printf("%d\t", Rates->MappingVect[Index]);
	printf("\n\n");
*/
	Rates->LogJacobion = (double)((double)(SizeG0 * SizeG1) / (double)(SizeG0 + SizeG1));
	Rates->LogJacobion = log(Rates->LogJacobion);

	Top		= RJSplitRatio(Rates, Opt, MapInfo, SizeG0, SizeG1, Rate, 1);
	Bottom	= RJMergeRatio(MapInfo->NoOfGroups);


	Ret = Top / Bottom;

	Rates->LnHastings = log(Ret) + Rates->LogJacobion;
}

void	RJSplitPropRatio(RATES* Rates, OPTIONS* Opt, MAPINFO *MapInfo, int G0Size, int G1Size, double Rate, int OGroup)
{
	double	Top;
	double	Bottom;
	double	Ret;
/*	int		Index;

	printf("Split: group: %d\t ni %d\tnj %d\tRate %f\n", OGroup, G0Size, G1Size, Rate);

	for(Index=0;Index<Rates->NoOfRates;Index++)
		printf("Rate %f\n", Rates->Rates[Index]);

	for(Index=0;Index<Rates->NoOfFullRates;Index++)
		printf("%d\t", Rates->MappingVect[Index]);
	printf("\n\n");

*/
	Rates->LogJacobion = (double)((double)(G0Size + G1Size) / (double)(G0Size * G1Size));
	Rates->LogJacobion = log(Rates->LogJacobion);


	Top		= RJMergeRatio(MapInfo->NoOfGroups+1);
	Bottom	= RJSplitRatio(Rates, Opt, MapInfo, G0Size, G1Size, Rate, 0);

	Ret = Top / Bottom;

	Rates->LnHastings = log(Ret) + Rates->LogJacobion;

}

void	RJMerge(RATES* Rates, OPTIONS* Opt)
{
	MAPINFO *MapInfo=NULL;
	int		G0, G1;
	int		G0ID, G1ID;
	int		Index;
	double	NewRate;
	int		DelPosInRates;
	double	*NewRates;
	int		ORPos;
	int		NRPos;

	MapInfo = CreatMapping(Rates);

	if(MapInfo->NoOfGroups == 1)
	{
		FreeMapInfo(MapInfo);
		return;
	}

	G0 = RandUSLong(Rates->RS) % MapInfo->NoOfGroups;
	do
		G1 = RandUSLong(Rates->RS) % MapInfo->NoOfGroups;
	while(G1 == G0);

	if(G1 < G0)
	{
		Index = G1;
		G1 = G0;
		G0 = Index;
	}


	G0ID = MapInfo->GroupID[G0];
	G1ID = MapInfo->GroupID[G1];

	DelPosInRates = Rates->MappingVect[MapInfo->GroupPos[G1][0]];


	NewRate = Rates->Rates[Rates->MappingVect[MapInfo->GroupPos[G0][0]]] * MapInfo->GroupSize[G0]; 
	NewRate += Rates->Rates[Rates->MappingVect[MapInfo->GroupPos[G1][0]]] * MapInfo->GroupSize[G1]; 
	NewRate = NewRate / (double)(MapInfo->GroupSize[G0] + MapInfo->GroupSize[G1]);
	Rates->Rates[Rates->MappingVect[MapInfo->GroupPos[G0][0]]] = NewRate;

	RJMergePropRatio(Rates, Opt, MapInfo, G0, G1, NewRate);


	for(Index=0;Index<Rates->NoOfFullRates;Index++)
	{
		if(Rates->MappingVect[Index] == G1ID)
			Rates->MappingVect[Index] = G0ID;
		
		if(Rates->MappingVect[Index] > G1ID)
			Rates->MappingVect[Index]--;
	}

	NewRates = (double*)malloc(sizeof(double) * (Rates->NoOfRates - 1));
	if(NewRates == NULL)
		MallocErr();

	NRPos = 0;
	for(ORPos=0;ORPos<Rates->NoOfRates;ORPos++)
	{
		if(ORPos != G1ID)
		{
			NewRates[NRPos] = Rates->Rates[ORPos];
			NRPos++;
		}
	}

	free(Rates->Rates);
	Rates->NoOfRates = MapInfo->NoOfGroups-1;
	Rates->Rates = NewRates;

	FreeMapInfo(MapInfo);
}


int		FindGropToReduce(RATES* Rates, MAPINFO* MapInfo)
{
	int		NoCandidates;
	int		No;

	NoCandidates = FindGroupsGrater(1, MapInfo);

	if(NoCandidates == 0)
		return -1;

	do
	{
		No = RandUSLong(Rates->RS) % MapInfo->NoOfGroups;
		if(MapInfo->GroupSize[No] > 1)
			return No;
	}while(1);

	return -1;
}


void	RJReduce(RATES* Rates, OPTIONS* Opt)
{
	MAPINFO *MapInfo=NULL;
	int		Group;
	int		TheOne;
	int		Pos;
	double	Top, Bottom;

	MapInfo = CreatMapping(Rates);

	Group = FindGropToReduce(Rates, MapInfo);

	if(Group == -1)
	{
		FreeMapInfo(MapInfo);
		return;
	}

	TheOne = RandUSLong(Rates->RS) % MapInfo->GroupSize[Group];

	Pos = MapInfo->GroupPos[Group][TheOne];

	Rates->MappingVect[Pos] = ZERORATENO;

	Top = 0.5;
	Top = Top * (1.0 / ((double)MapInfo->NoInZero + 1.0));
	Top = Top * (1.0 / (double)MapInfo->NoOfGroups);

	Bottom = 0.5;
	Bottom = Bottom * (1.0 / (double)FindGroupsGrater(1, MapInfo));
	Bottom = Bottom * (1.0 / (double)MapInfo->GroupSize[Group]);

	Rates->LnHastings = log(Top / Bottom);

	FreeMapInfo(MapInfo);
}

void	RJAugment(RATES* Rates, OPTIONS* Opt)
{
	double	Top, Bottom;
	MAPINFO *MapInfo=NULL;
	int		ZeroPos;
	int		Group;
	int		NoOfGroups;

	MapInfo = CreatMapping(Rates);		

	if(MapInfo->NoInZero == 0)
	{
		FreeMapInfo(MapInfo);
		return;
	}

	ZeroPos = RandUSLong(Rates->RS) % MapInfo->NoInZero;
	ZeroPos = MapInfo->ZeroPos[ZeroPos];
	Group = RandUSLong(Rates->RS) % MapInfo->NoOfGroups;

	Rates->MappingVect[ZeroPos] = Group;
	
	Top = 0.5;

	NoOfGroups = FindGroupsGrater(1, MapInfo);
	if(MapInfo->GroupSize[Group] == 1)
		NoOfGroups++;
	
	Top = Top * (1.0 / (double)NoOfGroups);
	Top = Top * (1.0 / (double)(MapInfo->GroupSize[Group] + 1));

	Bottom = 0.5;
	Bottom = Bottom * (1.0 / (double)MapInfo->NoInZero);
	Bottom = Bottom * (1.0 / (double)MapInfo->NoOfGroups);

	Rates->LnHastings = log(Top / Bottom);


	FreeMapInfo(MapInfo);
}
