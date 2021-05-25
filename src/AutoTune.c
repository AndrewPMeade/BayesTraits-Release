#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "typedef.h"
#include "genlib.h"
#include "AutoTune.h"

void		CalcRSqr(double *x, double *y, int Size, double *R2, double *Slope, double *Intercept);


int			GetAutoTuneSize(AUTOTUNE *AutoTune)
{
	if(AutoTune->No < AT_HSIZE)
		return AutoTune->No;

	return  AT_HSIZE;
}


AUTOTUNE*	CreatAutoTune(double Min, double Max)
{
	AUTOTUNE *Ret;
	int		Index;

	Ret = (AUTOTUNE*)malloc(sizeof(AUTOTUNE));
	if(Ret == NULL)
		return Ret;

	Ret->RateAcc = (double*)malloc(sizeof(double) * AT_HSIZE);
	Ret->RateDev = (double*)malloc(sizeof(double) * AT_HSIZE);
	
	if((Ret->RateAcc == NULL) || (Ret->RateDev == NULL))
		MallocErr();

	for(Index=0;Index<AT_HSIZE;Index++)
	{
		Ret->RateAcc[Index] = 0;
		Ret->RateDev[Index] = 0;
	}

	Ret->Min	= Min;
	Ret->Max	= Max;
	Ret->Target = ((Max - Min) * 0.5) + Min;
//	Ret->Target = 0.2;
	Ret->No		= 0;
	
	return Ret;
}

void		FreeAutoTune(AUTOTUNE *AutoTune)
{
	free(AutoTune->RateDev);
	free(AutoTune->RateAcc);
	free(AutoTune);
}

int			GetClosest(AUTOTUNE *AutoTune, double RD, double Acc)
{
	int Ret, Index, Size;
	
	Ret = -1;

	Size = GetAutoTuneSize(AutoTune);

	for(Index=0;Index<Size;Index++)
	{
		if(AutoTune->RateAcc[Index] < Acc)
		{

		}
	}

	return Ret;
}

double		BlindUpDate(AUTOTUNE *AutoTune, RANDSTATES *RS, double RD, double Acc)
{
	double Scale;

	if(Acc < AutoTune->Target)
		Scale = (RandDouble(RS) * 0.5) + 0.5;
	else
		Scale = RandDouble(RS) + 1;
		
	return RD * Scale;

/*	if(Acc < AutoTune->Target)
		return RD * (1.0 / AT_SCALE);
	else
		return RD * AT_SCALE;*/
}

int			InList(double *List, int Size, double RD)
{
	int Index;

	for(Index=0;Index<Size;Index++)
	{
		if(List[Index] == RD)
			return TRUE;
	}

	return FALSE;
}

void		AddAutoTune(AUTOTUNE *AutoTune, double RD, double Acc)
{
	int Pos;

	Pos = AutoTune->No % AT_HSIZE;

	if((AutoTune->No < AT_HSIZE)  || (AutoTune->Last != RD))
	{
		AutoTune->RateAcc[Pos] = Acc;
		AutoTune->RateDev[Pos] = RD;
	}

	AutoTune->Last = RD;
	AutoTune->No++;
}

int			AutoTuneValid(AUTOTUNE *AutoTune, double Acc)
{
	if((Acc >= AutoTune->Min) && (Acc <= AutoTune->Max))
		return TRUE;

	return FALSE;
}


void		PrintAutoTuneRates(AUTOTUNE *AutoTune)
{
	int Size, Index;

	Size = GetAutoTuneSize(AutoTune);

	for(Index=0;Index<Size;Index++)
		printf("%f\t%f\n", AutoTune->RateDev[Index], AutoTune->RateAcc[Index]);
	printf("\n");
}


double		AutoTuneNextRD(AUTOTUNE *AutoTune, RANDSTATES *RS, double RD, double Acc)
{
	double Ret;
	double R2, Slope, Int;
	double	OldRate;
	
//	return 0.548;
//	return 1;
	
//	return RandDouble(RS) * 2;
	

//	if(AutoTuneValid(AutoTune, Acc) == TRUE)
		AddAutoTune(AutoTune, RD, Acc);

	if(AutoTune->No > AT_HSIZE)
	{		
		if(AutoTuneValid(AutoTune, Acc) == FALSE)
			return  BlindUpDate(AutoTune, RS, RD, Acc);

	//	CalcRSqr(AutoTune->RateDev, AutoTune->RateAcc, AT_HSIZE, &R2, &Slope, &Int);
		CalcRSqr(AutoTune->RateAcc, AutoTune->RateDev, AT_HSIZE, &R2, &Slope, &Int);
		Ret = Int + (AutoTune->Target * Slope);
		
		/*	PrintAutoTuneRates(AutoTune);
		
			printf("%f\t%f\t%f\n", R2, Slope, Int);
			printf("New RD\t%f\n", Ret);

			printf("\n\n\n\n"); 
			fflush(stdout); */

		if(Ret < 0)
			Ret = BlindUpDate(AutoTune, RS, RD, Acc);
		
		return Ret;
	}
	
//	PrintAutoTuneRates(AutoTune);
	
	if(AutoTuneValid(AutoTune, Acc) == TRUE)
		return RD;

	Ret = BlindUpDate(AutoTune, RS, RD, Acc);
	
	return Ret;
}


