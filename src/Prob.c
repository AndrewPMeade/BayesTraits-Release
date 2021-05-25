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



#include <math.h>
#include <stdio.h>
#include <stdlib.h>


#include "TypeDef.h"
#include "Prob.h"

#ifndef M_PI
	#define M_PI       3.14159265358979323846
#endif


double	PDFNorm(double x, double Mean, double Var)
{
	double Ret, T1, T2;

	Ret = sqrt(Var) * sqrt(2*M_PI);
	Ret = 1.0 / Ret;

	T1 = (x - Mean) * (x - Mean);
	T2 = 2.0 * Var;

	T1 = T1 / T2;
	T2 = exp(-T1);

	Ret = Ret * T2;
	return Ret;
}

double		PDFExp(double X, double Mean)
{
	double Ret;

	Ret = Mean * exp(-Mean * X);

	return Ret;
}

double		PDFGamma(double X, double Shape, double Scale)
{
	double Ret;

	Ret = tgamma(Shape) * pow(Scale, Shape);

	Ret = 1.0 / Ret;

	Ret = Ret * pow(X, Shape-1) * exp(-(X/Scale));

	return Ret;
}

double	PDFSGamma(double x, double Alpha, double Beta)
{
	double Ret, s, T1, T2;

	s = 1 / ((Alpha - 1) * Beta);

	T1 = exp(-(x/s) / Beta);
	T2 = pow(x/s, -1 + Alpha);
	T2 = T1 * T2 * pow(Beta, -Alpha);
	Ret = T2 / tgamma(Alpha);

	Ret = Ret / s;

	Ret = Ret * VAR_RATES_PRIOR_SCALE;

	return Ret;
}

double		PDFInvGamma(double X, double Alpha, double Beta)
{
	double Ret;

	Ret = pow(Beta, Alpha) / tgamma(Alpha); 

	Ret = Ret * pow(X, -Alpha - 1.0);

	Ret = Ret * exp(-Beta / X);

	return Ret;
}

double		CDFNorm(double X, double Mean, double Var)
{
	double Ret;
	
	Ret = 0.5;

	Ret = (X - Mean) / (sqrt(Var) * sqrt(2.0));
	Ret = 0.5 * (1.0 + erf(Ret));
	
	return Ret;
}


