#ifdef BTOCL_DSC

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "typedef.h"
#include "trees.h"
#include "continuous.h"
#include "genlib.h"
#include "matrix.h"
// OpenCL headders

#include "btocl_discrete.h"

void ComputeQ(INVINFO *InvInfo, double **Q, int NOS);

double*	GenTVect(TREE *Tree, double RateMult, double Kappa)
{
	double *tvect;
	double Len;
	int Index;
	NODE N;

	tvect = (double*)malloc(sizeof(double) * Tree->NoNodes);
	if(tvect == NULL)
		MallocErr();
	
	// Trees->PList[N->ID]->me


	for(Index=0;Index<Tree->NoNodes;Index++)
	{
		N = Tree->NodeList[Index];
		
		Len = N->Length;
		if(Kappa != -1)
			Len = pow(Len, Kappa);

		Len = Len * RateMult;

		tvect[Index] = Len;
	}
	printf("\n");
	return tvect;
}

double	CreatPMatrix(double t, INVINFO *InvInfo, double **P, int NOS)
{
	int i,j,k;
	double E1, E2, TempD;
	double **TempM;
	double *Val;
	double	**Vec, **InvVec;

	// Allocate and populate entire Temp in GPU
	// if short of memory - do not use temporary matrix, may a temp vector with exp
	TempM = AllocMatMem(NOS, NOS);


	Val		= InvInfo->val;
	Vec		= InvInfo->vec->me;
	InvVec	= InvInfo->inv_vec->me;
	
	// If using temp, use a single call to update Temp with exp result.
	for(i=0;i<NOS;i++)
	{
		TempD = exp(t * Val[i]);
		for(j=0;j<NOS;j++)
			TempM[j][i] = Vec[j][i] * TempD;
	}

	E1 = 0.0;
	
	for(i=0;i<NOS;i++)
	{
		E2 = 0;
		for(j=0;j<NOS;j++)
		{
			P[i][j] = 0;
			for(k=0;k<NOS;k++)
				P[i][j] += TempM[i][k] * InvVec[k][j];
			// P[i][j] calculated
			// assert P[i][j] = exp(t * InvInfo->A->me[i][j])
			//printf("(%lf eqt %lf) ",P[i][j], exp(t * InvInfo->Q->me[i][j]) );
			E2 += P[i][j];

			if(P[i][j] < 0)
				return 1.0;
		}

		E2 = E2 - 1.0;
		E1 += E2 * E2;
	}
	
	FreeMatMem(TempM);

	return E1;
}

int		SetAllPMatrixOpenCL(RATES *Rates, TREES *Trees, OPTIONS *Opt, double RateMult, double Kappa)
{
	double *TVect;
	TREE *Tree;
	int Index, NOS;
	NODE N;
	double Err;
	
	//btdebug_enter("pmatrix");
	
	Tree = Trees->Tree[Rates->TreeNo];
	TVect = GenTVect(Tree, RateMult, Kappa);
	
	//ComputeQ(Trees->InvInfo, Trees->InvInfo->Q->me, NOS);
	
	NOS = Trees->NoOfStates;
	
	ComputeQ(Trees->InvInfo, Trees->InvInfo->Q->me, NOS);
	
	// Index 0 is root (what about NID??) - 0 it seems
	//printf("Addresses\n");
	for(Index=1;Index<Tree->NoNodes;Index++)
	{
		N = Tree->NodeList[Index];
		//printf("(Index %d diff %d) ",Index, Trees->PList[N->ID]->me-rootaddr);
		Err = CreatPMatrix(TVect[Index], Trees->InvInfo, Trees->PList[N->ID]->me , NOS);
		if(Err > 0.001)
		{
			free(TVect);
			return TRUE;
		}

//		PrintMatrix(Trees->PList[N->ID], "p=", stdout);
	}
	//printf("\n");
//	exit(0);

	free(TVect);
	
	//btdebug_exit("pmatrix");
	
	return FALSE;
}

void ComputeQ(INVINFO *InvInfo, double **Q, int NOS)
{
	int i,j,k;
	double TempD;
	double **TempM;
	double *Val;
	double	**Vec, **InvVec;

	// Allocate and populate entire Temp in GPU
	// if short of memory - do not use temporary matrix, may a temp vector with exp
	TempM = AllocMatMem(NOS, NOS);

	Val		= InvInfo->val;
	Vec		= InvInfo->vec->me;
	InvVec	= InvInfo->inv_vec->me;
	
	// If using temp, use a single call to update Temp with exp result.
	for(i=0;i<NOS;i++)
	{
		TempD = Val[i];
		for(j=0;j<NOS;j++)
			TempM[j][i] = Vec[j][i] * TempD;
	}

	for(i=0;i<NOS;i++)
	{

		for(j=0;j<NOS;j++)
		{
			Q[i][j] = 0;
			for(k=0;k<NOS;k++)
				Q[i][j] += TempM[i][k] * InvVec[k][j];
		}


	}

	FreeMatMem(TempM);

	return;

}

#endif



