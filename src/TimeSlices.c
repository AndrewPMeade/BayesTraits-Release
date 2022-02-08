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
#include <string.h>
#include <math.h>

#include "TypeDef.h"
#include "GenLib.h"
#include "TimeSlices.h"
#include "Trees.h"
#include "LocalTransform.h"
#include "Priors.h"

TIME_SLICE*		AllocTimeSlice(char *Name)
{
	TIME_SLICE* Ret;

	Ret = (TIME_SLICE*)SMalloc(sizeof(TIME_SLICE));
	Ret->Name = StrMake(Name);

	
	Ret->FixedScale = FALSE;
	Ret->FixedTime	= FALSE;

	Ret->Scale = -1.0;
	Ret->Time = -1.0;

	return Ret;
}

void	FreeTimeSlice(TIME_SLICE* TS)
{
	free(TS->Name);
	free(TS);
}

TIME_SLICE*		CloneTimeSlice(TIME_SLICE* TS)
{
	TIME_SLICE* Ret;

	Ret = AllocTimeSlice(TS->Name);

	Ret->FixedScale = TS->FixedScale;
	Ret->FixedTime = TS->FixedTime;

	Ret->Scale = TS->Scale;
	Ret->Time = TS->Time;

	return Ret;
}

TIME_SLICE*	AddTimeSlice(TIME_SLICES *TSlices, char *Name, double Time, double Scale)
{
	TIME_SLICE* TS;

	TS = AllocTimeSlice(Name);

	if(Time != -1.0)
	{
		TS->Time = Time;
		TS->FixedTime = TRUE;
	}

	if(Scale != -1.0)
	{
		TS->Scale = Scale;
		TS->FixedScale = TRUE;
	}

	TSlices->TimeSlices = (TIME_SLICE**)AddToList(&TSlices->NoTimeSlices, (void**)TSlices->TimeSlices, (void*)TS);

	return TS;
}

void	FreeTimeSlices(TIME_SLICES *TSlices)
{
	int Index;

	for(Index=0;Index<TSlices->NoTimeSlices;Index++)
		FreeTimeSlice(TSlices->TimeSlices[Index]);

	if(TSlices->TimeSlices != NULL)
		free(TSlices->TimeSlices);
	free(TSlices);
}

TIME_SLICES*	CreateTimeSlices(void)
{
	TIME_SLICES* Ret;

	Ret = (TIME_SLICES*)SMalloc(sizeof(TIME_SLICES));
	Ret->NoTimeSlices = 0;
	Ret->TimeSlices = NULL;

	return Ret;
}

TIME_SLICES*	CloneTimeSlices(TIME_SLICES* TSlices)
{
	TIME_SLICES* Ret;
	int Index;

	Ret = CreateTimeSlices();

	Ret->NoTimeSlices = TSlices->NoTimeSlices;

	if(Ret->NoTimeSlices > 0)
	{
		Ret->TimeSlices = (TIME_SLICE**)SMalloc(sizeof(TIME_SLICES*) * Ret->NoTimeSlices);
		for(Index=0;Index<Ret->NoTimeSlices;Index++)
			Ret->TimeSlices[Index] = CloneTimeSlice(TSlices->TimeSlices[Index]);
	}
	else
		Ret->TimeSlices = NULL;

	return Ret;
}

TIME_SLICE*		GetTimeSlice(TIME_SLICES *TSlices, char *Name)
{
	int Index;

	for(Index=0;Index<TSlices->NoTimeSlices;Index++)
	{
		if(strcmp(TSlices->TimeSlices[Index]->Name, Name) == 0)
			return TSlices->TimeSlices[Index];
	}

	return NULL;
}

void	PrintTimeSlices(FILE *Str, TIME_SLICES *TSlices)
{
	int Index;
	TIME_SLICE *TS;

	if(TSlices->NoTimeSlices == 0)
		return;

	fprintf(Str, "Number of time slices %d.\n", TSlices->NoTimeSlices);

	for(Index=0;Index<TSlices->NoTimeSlices;Index++)
	{
		TS = TSlices->TimeSlices[Index];
		fprintf(Str, "\t%s\t", TS->Name);

		if(TS->FixedScale == TRUE && TS->FixedTime == TRUE)
			fprintf(Str, "%f\t%f", TS->Time, TS->Scale);
		else
			if(TS->FixedTime == TRUE)
				fprintf(Str, "%f", TS->Time);

		fprintf(Str, "\n");
	}
}

double			RandScale(RANDSTATES *RS)
{
	double Ret;

	if(RandDouble(RS) < 0.5)
		return RandDouble(RS);

	do
	{
		Ret = 1.0 / RandDouble(RS);
	}while(Ret > 100);

	return Ret;
}

TIME_SLICES*	CreateRatesTimeSlices(RATES *Rates, TIME_SLICES *TSlices)
{
	TIME_SLICES	*Ret;
	TIME_SLICE	*T;
	int Index;
	

	Ret = CloneTimeSlices(TSlices);

	for(Index=0;Index<Ret->NoTimeSlices;Index++)
	{
		T = Ret->TimeSlices[Index];

		if(T->FixedTime == FALSE)
			T->Time = RandDouble(Rates->RS);

		if(T->FixedScale == FALSE)
			T->Scale = RandScale(Rates->RS);
	}

	return Ret;
}

void			CopyTimeSlices(TIME_SLICES *A, TIME_SLICES *B)
{
	int Index;
	TIME_SLICE	*TA, *TB;

	for(Index=0;Index<A->NoTimeSlices;Index++)
	{
		TB = B->TimeSlices[Index];
		TA = GetTimeSlice(A, TB->Name);

		TA->Scale = TB->Scale;
		TA->Time = TB->Time;
	}
}


int		CompTimeSlices(const void *CV1, const void *CV2)
{
	TIME_SLICE **TS1, **TS2;

	TS1 = (TIME_SLICE**)CV1;
	TS2 = (TIME_SLICE**)CV2;

	if((*TS1)->Time > (*TS2)->Time)
		return 1;

	if((*TS1)->Time < (*TS2)->Time)
		return -1;

	return 0;
}

void	PrintTimeSlice(FILE *Out, TIME_SLICE *TS)
{
	fprintf(Out, "\t%s\t%f\t%f\n", TS->Name, TS->Time, TS->Scale);
}

void	ReSetNodeScalars(TREE *Tree)
{
	int Index;

	for(Index=0;Index<Tree->NoNodes;Index++)
		Tree->NodeList[Index]->ScaleFactor = 0;
}

void	ApplyNodeScalars(TREE *Tree)
{
	int Index;

	for(Index=0;Index<Tree->NoNodes;Index++)
	{
		if(Tree->NodeList[Index]->Length != 0)
			Tree->NodeList[Index]->Length = Tree->NodeList[Index]->Length + Tree->NodeList[Index]->ScaleFactor;
	}
}

double	GetMaxRootToTip(TREE *Tree)
{
	double Ret;
	int Index;

	Ret = 0;
	for(Index=0;Index<Tree->NoNodes;Index++)
		if(Tree->NodeList[Index]->DistToRoot > Ret)
			Ret = Tree->NodeList[Index]->DistToRoot;

	return Ret;
}

double	FindPctBLOverLap(double BL_Start, double BL_End, double TS_Start, double TS_End)
{
	double A, Ret;
	

	if(TS_Start > BL_End)
		return 0.0;

	if(TS_End < BL_Start)
		return 0.0;
	
	if(TS_End > BL_End)
		TS_End = BL_End;

	if(TS_Start < BL_Start)
		TS_Start = BL_Start;

	A =  TS_End - TS_Start;

	Ret = A / (BL_End - BL_Start);



	return Ret;
}

double FindScaledBL(double Length, double P, double Scale)
{
	return ((1.0 - P) + (P * Scale)) * Length;
}

void	RunTimeSlice(TREE *Tree, double Start, double End, double Scale)
{
	int Index;
	double P, NLen;
	NODE N;

	for(Index=0;Index<Tree->NoNodes;Index++)
	{
		N = Tree->NodeList[Index];
		if(N->Ans != NULL)
		{
			P = FindPctBLOverLap(N->Ans->DistToRoot, N->DistToRoot, Start, End);
		//	N->ScaleFactor += (Scale * Pct);

			NLen = FindScaledBL(N->Length, P, Scale);
			

//			N->ScaleFactor = N->ScaleFactor * (NLen / N->Length);
//			N->ScaleFactor = (NLen / N->Length);
			N->ScaleFactor += (NLen - N->Length);

//			printf("%d\t%f\t%f\t%f\t%f\n", Index, N->Ans->DistToRoot, N->DistToRoot, N->ScaleFactor, P);
		}
	}

//	exit(0);
}

// Carnivores.trees
// 
//./Seq/Carnivores.trees ./Seq/Carnivores.txt < in.txt > sout.txt
//./Seq/3TaxaTest.trees ./Seq/Carnivores.txt < in.txt > sout.txt

// ./Seq/MammalBig.trees ./Seq/MammalBig.txt < in.txt > sout.txt

void	ApplyTimeSlices(RATES *Rates, TREES *Trees)
{
	TIME_SLICES *TS;
	TIME_SLICE	*Slice;
	TREE *Tree;
	double	Start, End, Scale, RootToTip;
	int Index;
	
	Tree = Trees->Tree[Rates->TreeNo];

	ReSetNodeScalars(Tree);

	RootToTip = GetMaxRootToTip(Tree);

	TS = Rates->TimeSlices;

//	for(Index=0;Index<TS->NoTimeSlices;Index++)
//		PrintTimeSlice(stdout, TS->TimeSlices[Index]);

	qsort(TS->TimeSlices, TS->NoTimeSlices, sizeof(TIME_SLICE*), CompTimeSlices);


	for(Index=0;Index<TS->NoTimeSlices;Index++)
	{
		Slice = TS->TimeSlices[Index];

		Scale = Slice->Scale;
		Start = Slice->Time;
		End = 1.0;
		if(Index != TS->NoTimeSlices - 1)
			End = TS->TimeSlices[Index+1]->Time;
		
		Start = Start * RootToTip;
		End = End * RootToTip;

//		printf("%s\t%f\t%f\t%f\n", Slice->Name, Start, End, Scale);

		RunTimeSlice(Tree, Start, End, Scale);

//		PrintTimeSlice(stdout, TS->TimeSlices[Index]);
	}

	ApplyNodeScalars(Tree);

//	SaveTrees("TimeSliced.trees", Trees);		exit(0);
}

void	PrintTimeSliceHeader(FILE *Str, TIME_SLICES *TS)
{
	int Index;
	TIME_SLICE *Slice;


	for(Index=0;Index<TS->NoTimeSlices;Index++)
	{
		Slice = TS->TimeSlices[Index];
		fprintf(Str, "%s - Time\t%s - Scale\t", Slice->Name, Slice->Name);
	}
}

void	PrintTimeSliceRates(FILE *Str, TIME_SLICES *TS_Opt, TIME_SLICES *TS_Rates)
{
	int Index;
	TIME_SLICE *Slice;

	for(Index=0;Index<TS_Opt->NoTimeSlices;Index++)
	{
		Slice = GetTimeSlice(TS_Rates, TS_Opt->TimeSlices[Index]->Name);
		fprintf(Str, "%f\t%f\t", Slice->Time, Slice->Scale);
	}

}

int				TimeSliceEstTime(TIME_SLICES *TS)
{
	int Index;
	TIME_SLICE *T;

	if(TS == NULL)
		return FALSE;

	for(Index=0;Index<TS->NoTimeSlices;Index++)
	{
		T = TS->TimeSlices[Index];
		if(T->FixedTime == FALSE)
			return TRUE;
	}

	return FALSE;
}

int				TimeSliceEstScale(TIME_SLICES *TS)
{
	int Index;
	TIME_SLICE *T;

	if(TS == NULL)
		return FALSE;

	for(Index=0;Index<TS->NoTimeSlices;Index++)
	{
		T = TS->TimeSlices[Index];
		if(T->FixedScale == FALSE)
			return TRUE;
	}

	return FALSE;
}

TIME_SLICE*	GetTimeSliceEstTime(RANDSTATES *RS, TIME_SLICES *TS)
{
	int Pos;

	while(TRUE)
	{
		Pos = RandUSInt(RS) % TS->NoTimeSlices;

		if(TS->TimeSlices[Pos]->FixedTime == FALSE)
			return TS->TimeSlices[Pos];
	}

	return NULL;
}

void	ChangeTimeSliceTime(RATES *Rates, SCHEDULE *Shed)
{
	double Dev, NTime;
	TIME_SLICE *TS;

	Shed->CurrentAT = Shed->TimeSliceTimeAT;
	Dev = Shed->CurrentAT->CDev;


	TS = GetTimeSliceEstTime(Rates->RS, Rates->TimeSlices);

	do
	{
		NTime = RandNormal(Rates->RS, TS->Time, Dev);

		if(NTime >= 0 && NTime <= 1.0)
		{
			TS->Time = NTime;
			return;
		}
	} while(TRUE);
}


TIME_SLICE*	GetTimeSliceEstScale(RANDSTATES *RS, TIME_SLICES *TS)
{
	int Pos;

	while(TRUE)
	{
		Pos = RandUSInt(RS) % TS->NoTimeSlices;

		if(TS->TimeSlices[Pos]->FixedScale == FALSE)
			return TS->TimeSlices[Pos];
	}

	return NULL;
}


void	ChangeTimeSliceScale(RATES *Rates, SCHEDULE *Shed)
{
	double NScale, Dev;
	TIME_SLICE *TS;

	Shed->CurrentAT = Shed->TimeSliceScaleAT;
	Dev = Shed->CurrentAT->CDev;

	TS = GetTimeSliceEstScale(Rates->RS, Rates->TimeSlices);

	NScale = ChangeLocalScale(Rates->RS, TS->Scale, Dev);

	Rates->LnHastings = CalcNormalHasting(TS->Scale, Dev);

	if(NScale > MAX_LOCAL_RATE || NScale < MIN_LOCAL_RATE)
		NScale = TS->Scale;

	TS->Scale = NScale;
}