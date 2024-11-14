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
#include <math.h>

#include "TypeDef.h"
#include "GenLib.h"
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

void		SetMaxDev(AUTOTUNE *AutoTune, double MaxDev)
{
	AutoTune->MaxDev = MaxDev;
}

AUTOTUNE*	CreatAutoTune(char *Name, double InitDev, double Min, double Max)
{
	AUTOTUNE *Ret;
	int		Index;

	Ret = (AUTOTUNE*)SMalloc(sizeof(AUTOTUNE));
	Ret->RateAcc = (double*)SMalloc(sizeof(double) * AT_HSIZE);
	Ret->RateDev = (double*)SMalloc(sizeof(double) * AT_HSIZE);
	
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
	Ret->MaxDev	= -1.0;

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


void		BlindUpDate(AUTOTUNE *AT, gsl_rng *RNG, double Acc)
{
	double Scale;

	if(Acc < AT->Target)
		Scale = (gsl_rng_uniform_pos(RNG) * 0.5) + 0.5;
	else
		Scale = gsl_rng_uniform_pos(RNG) + 1;
	
	AT->CDev = AT->CDev * Scale;

	if(AT->MaxDev != -1.0 && AT->CDev > AT->MaxDev)
		AT->CDev = AT->MaxDev;
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

void		AutoTuneUpDate(AUTOTUNE *AT, gsl_rng *RNG)
{
	double NDev;
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
			BlindUpDate(AT, RNG, Acc);
			return;
		}

		CalcRSqr(AT->RateAcc, AT->RateDev, AT_HSIZE, &R2, &Slope, &Int);
		NDev = Int + (AT->Target * Slope);

		if(NDev < 0)
			BlindUpDate(AT, RNG, Acc);
		else
			AT->CDev = NDev;

		if(AT->MaxDev != -1.0 && AT->CDev > AT->MaxDev)
			AT->CDev = AT->MaxDev;
		
		return;
	}
	
	if(AutoTuneValid(AT, Acc) == TRUE)
		return;
	
	BlindUpDate(AT, RNG, Acc);
}


