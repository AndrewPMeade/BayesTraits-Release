#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "genlib.h"
#include "matrix.h"
#include "ContrastRegCalcReg.h"
#include "contrasts.h"
#include "typedef.h"
#include "linalg.h"

void	SetContrastMultiRegData(MATRIX *M, int N,  int NoX, REG_BETA_SPACE*	RegSpace)
{
	int x, y;

	for(x=0;x<N;x++)
		RegSpace->Uy->me[0][x] = M->me[x][0];
	
	for(x=0;x<N;x++)
	{
		for(y=0;y<NoX;y++)
			RegSpace->Ux->me[x][y] = M->me[x][y+1];
	}
}

//double*	ContrastMultiReg(double *Y, double **X, int N,  int NoX, int TestCorrel)
double*	ContrastMultiReg(MATRIX *M, int TestCorrel)
{
	int NoX, N;
	double *Ret;
	REG_BETA_SPACE*	RegSpace;
	int		Index, Err;

	N = M->NoOfRows;
	NoX = M->NoOfCols-1;

	RegSpace =	InitRegBetaSpace(NoX, N);
	
	Ret = (double*)malloc(sizeof(double) * NoX);
	if(Ret == NULL)
		MallocErr();

	SetContrastMultiRegData(M, N,  NoX, RegSpace);
	
	Transpose(RegSpace->Ux, RegSpace->TUx);

	MatrixMult(RegSpace->TUx, RegSpace->Ux, RegSpace->Prod1);

	if(NoX == 1)
	{
		if(TestCorrel == FALSE)
			RegSpace->InvUx->me[0][0] = 1;
		else
			RegSpace->InvUx->me[0][0] = 1.0 / RegSpace->Prod1->me[0][0];
	}
	else
	{
		if(TestCorrel == FALSE)
			SetIdentityMatrix(RegSpace->InvUx);
		else
		{
			Err = InvertMatrix(RegSpace->Prod1->me, NoX, RegSpace->TempDVect, RegSpace->TempIVect, RegSpace->InvUx->me);

			if(Err != NO_ERROR)
			{
				printf("Matrix singular: %s::%d\n", __FILE__, __LINE__);
				exit(0);
			}
		}
	}

	MatrixMult(RegSpace->TUx, RegSpace->Uy, RegSpace->Prod2);

	MatrixMult(RegSpace->InvUx, RegSpace->Prod2, RegSpace->Prod3);

	for(Index=0;Index<NoX;Index++)
		Ret[Index] = RegSpace->Prod3->me[0][Index];

	FreeRegBetaSpace(RegSpace);
	
	return Ret;
}