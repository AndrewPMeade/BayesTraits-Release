#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "typedef.h"
#include "genlib.h"
#include "AutoTune.h"

void		CalcRSqr(double *x, double *y, int Size, double *R2, double *Slope, double *Intercept);

void		ReSetAutoTune(AUTOTUNE *AutoTune)
{
	AutoTune->NoTried	= 0;
	AutoTune->NoAcc		= 0;
}

int			GetAutoTuneSize(AUTOTUNE *AutoTune)
{
	if(AutoTune->No < AT_HSIZE)
		return AutoTune->No;

	return  AT_HSIZE;
}

AUTOTUNE*	CreatAutoTune(char *Name, double InitDev, double Min, double Max)
{
	AUTOTUNE *Ret;
	int		Index;

	Ret = (AUTOTUNE*)malloc(sizeof(AUTOTUNE));
	if(Ret == NULL)
		return Ret;

	Ret->RateAcc = (double*)malloc(sizeof(double) * AT_HSIZE);
	Ret->RateDev = (double*)malloc(sizeof(double) * AT_HSIZE);
	
	if(Ret->RateAcc == NULL || Ret->RateDev == NULL)
		MallocErr();

	for(Index=0;Index<AT_HSIZE;Index++)
	{
		Ret->RateAcc[Index] = 0;
		Ret->RateDev[Index] = 0;
	}

	Ret->Last	= -1.0;

	Ret->Min	= Min;
	Ret->Max	= Max;
	Ret->Target = ((Max - Min) * 0.5) + Min;
	Ret->No		= 0;

	if(Name != NULL)
		Ret->Name = StrMake(Name);
	else
		Ret->Name = NULL;
	
	Ret->CDev	= InitDev;

	ReSetAutoTune(Ret);

	return Ret;
}

void		FreeAutoTune(AUTOTUNE *AutoTune)
{
	free(AutoTune->RateDev);
	free(AutoTune->RateAcc);

	if(AutoTune->Name != NULL)
		free(AutoTune->Name);

	free(AutoTune);
}


void		BlindUpDate(AUTOTUNE *AT, RANDSTATES *RS, double Acc)
{
	double Scale;

	if(Acc < AT->Target)
		Scale = (RandDouble(RS) * 0.5) + 0.5;
	else
		Scale = RandDouble(RS) + 1;
	
	AT->CDev = AT->CDev * Scale;
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

void		AddAutoTune(AUTOTUNE *AT, double Acc)
{
	int Pos;

	Pos = AT->No % AT_HSIZE;

	if(AT->Last != AT->CDev)
	{
		AT->RateAcc[Pos] = Acc;
		AT->RateDev[Pos] = AT->CDev;
		AT->No++;
	}

	AT->Last = AT->CDev;
}

int			AutoTuneValid(AUTOTUNE *AutoTune, double Acc)
{
	if(Acc >= AutoTune->Min && Acc <= AutoTune->Max)
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

double		AutoTuneCalcAcc(AUTOTUNE *AT)
{
	if (AT->NoTried == 0)
		return 0.0;

	return (double)AT->NoAcc / AT->NoTried;
}

void		AutoTuneUpDate(AUTOTUNE *AT, RANDSTATES *RS)
{
	double Ret;
	double R2, Slope, Int;
	double	Acc;

	if(AT->NoTried < AT_MIN_TRIED)
		return;

	Acc = AutoTuneCalcAcc(AT);
	ReSetAutoTune(AT);

	AddAutoTune(AT, Acc);

	if(AT->No > AT_HSIZE)
	{		
		if(AutoTuneValid(AT, Acc) == FALSE)
		{
			BlindUpDate(AT, RS, Acc);
			return;
		}

		CalcRSqr(AT->RateAcc, AT->RateDev, AT_HSIZE, &R2, &Slope, &Int);
		Ret = Int + (AT->Target * Slope);

		if(Ret < 0)
			BlindUpDate(AT, RS, Acc);
		
		return;
	}
	
	if(AutoTuneValid(AT, Acc) == TRUE)
		return;
	
	BlindUpDate(AT, RS, Acc);
}