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
#include <string.h>

#include "Threaded.h"
#include "TypeDef.h"
#include "GenLib.h"
#include "Continuous.h"
#include "CKappa.h"
#include "Part.h"

void	RecCalcKappaV(TREES* Trees, TREE *Tree, NODE N, PART *DiffPart, double Kappa, double SumLogPath)
{
	int x,y, XPos, YPos;
	double **Mat;
	double Dist;

	Mat = Tree->ConVars->V->me;
	Dist = SumLogPath;
		
	GetPartDiff(N->Ans->Part, N->Part, DiffPart);

	for(x=0;x<N->Part->NoTaxa;x++)
	{
		XPos = N->Part->Taxa[x];
		for(y=0;y<DiffPart->NoTaxa;y++)
		{
			YPos = DiffPart->Taxa[y];
			Mat[YPos][XPos] = Dist;
			Mat[XPos][YPos] = Dist;
		}
	}

	if(N->Tip == TRUE)
	{
		XPos = GetMapID(Trees, N->Taxa->No);
		Mat[XPos][XPos] = SumLogPath + pow(N->Length, Kappa);
		return;
	}

	Dist = Dist + pow(N->Length, Kappa);

	for(x=0;x<N->NoNodes;x++)
		RecCalcKappaV(Trees, Tree, N->NodeList[x], DiffPart, Kappa, Dist);
}

void	MakeKappaV(TREES* Trees, TREE *Tree, double Kappa)
{
	int		Index;
	NODE	N;
	PART	*CPart;

	CPart = CreatPart(Trees->NoTaxa);


	N = Tree->Root;
	for(Index=0;Index<N->NoNodes;Index++)
		RecCalcKappaV(Trees, Tree, N->NodeList[Index], CPart, Kappa, 0);

	FreePart(CPart);
}

