#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

#include "Threaded.h"
#include "typedef.h"
#include "genlib.h"
#include "continuous.h"
#include "ckappa.h"
#include "part.h"

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

