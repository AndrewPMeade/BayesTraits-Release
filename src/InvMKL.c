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