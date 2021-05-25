/*
*  BayesTriats 3.0
*
*  copyright 2017
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

#include "SliceSampler.h"
#include "TypeDef.h"
#include "GenLib.h"

SLICESAMPLER*	CrateSliceSampler(int NoSteps)
{
	SLICESAMPLER* Ret;

	Ret = (SLICESAMPLER*)SMalloc(sizeof(SLICESAMPLER));
	
	Ret->NoSlices = 0;
	Ret->NoSteps = NoSteps;


	Ret->SliceX = (double*)SMalloc(sizeof(double) * NoSteps);
	Ret->SliceY = (double*)SMalloc(sizeof(double) * NoSteps);
	
	Ret->SliceMin = (double*)SMalloc(sizeof(double) * NoSteps);
	Ret->SliceMax = (double*)SMalloc(sizeof(double) * NoSteps);
	
	return Ret;
}

void FreeSliceSampler(SLICESAMPLER* SS)
{
	free(SS->SliceX);
	free(SS->SliceY);

	free(SS->SliceMax);
	free(SS->SliceMin);

	free(SS);
}

void	SSSetXPosVect(SLICESAMPLER* SS, double Min, double Max)
{
	double StepSize;
	int Index;

	StepSize = (Max - Min) / (SS->NoSteps - 1);

	for(Index=0;Index<SS->NoSteps;Index++)
		SS->SliceX[Index] = Min + (StepSize * Index);

//	printf("%f\t%f\n", Min, Max);

//	for(Index=0;Index<SS->NoSteps;Index++)
//		printf("%d\t%f\n", Index, SS->SliceX[Index]);

//	exit(0);
}



void	GetMinMaxY(double *List, int Size, double *Min, double *Max)
{
	int Index;

	*Min = List[0];
	*Max = List[0];

	for(Index=1;Index<Size;Index++)
	{
		if(List[Index] > *Max)
			*Max = List[Index];

		if(List[Index] < *Min)
			*Min = List[Index];
	}
}

double	GetXCrossPoint(double YPoint, int Index, double *YList, double *XList)
{
	double Ret;
	
	double Diff1, Diff2;

//	Slope = (YList[Index-1] - YList[Index]) / (XList[Index-1] - XList[Index]);
//	YDiff = YPoint - YList[Index-1];

	Diff1 = fabs(YPoint - YList[Index-1]);
	Diff2 = fabs(YPoint - YList[Index]);


	Ret = Diff1 / (Diff1 + Diff2);
	Ret = Ret * (XList[Index] - XList[Index-1]);
	Ret = Ret + XList[Index-1];

	return Ret;
}

double SSFindSliceStart(SLICESAMPLER *SS, double YPoint, int *Pos)
{
	if(*Pos == 0 && SS->SliceY[0] > YPoint)
		return SS->SliceX[0];

	while(*Pos < SS->NoSteps)
	{
		if(SS->SliceY[*Pos] > YPoint)
			return GetXCrossPoint(YPoint, *Pos, SS->SliceY, SS->SliceX);

		(*Pos)++;
	}

	return 0.0;
}

double SSFindSliceEnd(SLICESAMPLER *SS, double YPoint, int *Pos)
{
	while(*Pos < SS->NoSteps)
	{
		if(*Pos == SS->NoSteps - 1)
			return SS->SliceX[*Pos];
		
		if(SS->SliceY[*Pos] < YPoint)
			return GetXCrossPoint(YPoint, *Pos, SS->SliceY, SS->SliceX);

		(*Pos)++;
	}

	return 0.0;
}

double	FindNewPoint(SLICESAMPLER *SS, RANDSTATES *RS)
{
	double Sum, Point;
	int Index;

	Sum = 0.0;

	for(Index=0;Index<SS->NoSlices;Index++)
		Sum += SS->SliceMax[Index] - SS->SliceMin[Index];

	Point = Sum * RandDouble(RS);

	Sum = 0.0;

	for(Index=0;Index<SS->NoSlices;Index++)
	{
		if(Point < Sum + (SS->SliceMax[Index] - SS->SliceMin[Index]))
			return SS->SliceMin[Index] + (Point - Sum);

		Sum += SS->SliceMax[Index] - SS->SliceMin[Index];
	}

	for(Index=0;Index<SS->NoSlices;Index++)
		printf("%d\t%f\t%f\n", Index, SS->SliceX[Index], SS->SliceY[Index]);
	fflush(stdout);

	printf("Error in (%s::%d)\n", __FILE__, __LINE__);
	exit(0);

	return 1.0;
}


double	SSGetNewPoint(SLICESAMPLER *SS, RANDSTATES *RS, double POld)
{
	int CPos, Exit, NoSlices;
	double	SPos, EPos, Ret, YPoint, MinY, MaxY;


	Ret = 0;
	
	GetMinMaxY(SS->SliceY, SS->NoSteps, &MinY, &MaxY);

	do
	{
		YPoint = MinY + ((POld - MinY) * RandDouble(RS));
	} while(YPoint > MaxY);
	
	Exit = FALSE;

	CPos = 0;
	NoSlices = 0;

	do
	{
		SPos = SSFindSliceStart(SS, YPoint, &CPos);

		if(CPos >= SS->NoSteps)
			Exit = TRUE;
		else
		{
			EPos = SSFindSliceEnd(SS, YPoint, &CPos);

			SS->SliceMin[NoSlices] = SPos;
			SS->SliceMax[NoSlices] = EPos;
			NoSlices++;
			CPos+=1;
		}

	}while(Exit == FALSE);

	SS->NoSlices = NoSlices;

	Ret = FindNewPoint(SS, RS);
	return Ret;
}