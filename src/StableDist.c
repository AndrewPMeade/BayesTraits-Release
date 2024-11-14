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

double		StableDistTPDF(double Scale, double X, double t)
{
	return CaclNormalLogLh(X, Scale, t);
}
