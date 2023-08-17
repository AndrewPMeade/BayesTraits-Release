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

#include "StableDist.h"
#include "GenLib.h"

STABLEDIST*	CreatStableDist(void)
{
	STABLEDIST* Ret;
	int Index;

	Ret = (STABLEDIST*)SMalloc(sizeof(STABLEDIST));

	Ret->Alpha = -1;
	Ret->Scale = -1;
	Ret->ScaledAlpha  = -1;

	Ret->P = (double**)SMalloc(sizeof(double*) * 4);
	for(Index=0;Index<4;Index++)
		Ret->P[Index] = (double*)SMalloc(sizeof(double) * 4);

	return Ret;
}

void		SetStableDist(STABLEDIST* SD, double Alpha, double Scale)
{
	SD->Alpha = Alpha;
	SD->Scale = Scale;

	SD->ScaledAlpha = Alpha * 20 - 1;
}

void		FreeStableDist(STABLEDIST* SD)
{
	int Index;
	
	for(Index=0;Index<4;Index++)
		free(SD->P[Index]);
	free(SD->P);
	free(SD);
}



double		StableDistPDF(STABLEDIST* SD, double x)
{
	return StableDistTPDF(SD, x, 1.0);
}

// Scale is Variance
// Lh is in log space
double		CaclNormalLogLh(double X, double Scale, double T)
{
	double Ret;
	double T1;

	Scale = Scale * T;

	Ret = log(1.0 / (sqrt(Scale) * 2.506628274631));

	T1 = -((X * X) / (2.0 * Scale));

	return Ret + T1;
}

double		StableDistTPDF(STABLEDIST* SD, double X, double t)
{
	return CaclNormalLogLh(X, SD->Scale, t);
}

