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

#include "typedef.h"
#include "genlib.h"
#include "Threaded.h"


#ifdef USE_MLK

#include "InvMKL.h"
#include "btlapack.h"

double		GetLnDet(double **Mat, int N)
{
	double Ret;
	int Index;
	
	Ret = 0;
	for(Index=0;Index<N;Index++)
	{
		if(Mat[Index][Index] < 0)
			Ret += log(-Mat[Index][Index]);
		else
			Ret += log(Mat[Index][Index]);
	}

	return Ret;
}

int			InvMLK(TREES *Trees, TREE *Tree)
{
	int Index, Info, WorkSize;
	TEMPCONVAR*	TempCon;
	double *Mat, *Temp, LnDet;

	TempCon = Trees->TempConVars;
	
	CopyMatrix(Tree->ConVars->InvV, Tree->ConVars->V);

	Mat = Tree->ConVars->InvV->me[0];
	Temp = TempCon->TMat->me[0];

	WorkSize = Trees->NoTaxa * Trees->NoTaxa;

	dgetrf(&Trees->NoTaxa, &Trees->NoTaxa, Mat, &Trees->NoTaxa, TempCon->T2, &Info);

	Tree->ConVars->LogDetOfV = GetLnDet(Tree->ConVars->InvV->me, Trees->NoTaxa);
	
	if(Info != 0)
		return Info;

	dgetri(&Trees->NoTaxa, Mat, &Trees->NoTaxa, TempCon->T2, Temp, &WorkSize, &Info);


	return Info;
}

#endif