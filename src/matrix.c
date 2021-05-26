#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "typedef.h"
#include "genlib.h"
#include "matrix.h"


/* Some matrix manipulasion rotueses */
/* The representatin of me is Row By Coloum*/


void	FreeMatMem(double** Mat)
{
	free(Mat[0]);
	free(Mat);
}

double** AllocMatMem(int NoOfRows, int NoOfCols)
{
	double**	Ret=NULL;
	int		RIndex=0;
	int		CIndex=0;

	Ret = (double**)malloc(NoOfRows * sizeof(double*));
	if(Ret==NULL)
		MallocErr();

	Ret[0] = (double*)malloc(NoOfRows * NoOfCols * sizeof(double));

	if(Ret[0]==NULL)
		MallocErr();

	for(RIndex=1;RIndex<NoOfRows;RIndex++)
		Ret[RIndex] = Ret[0] + RIndex * NoOfCols;

	return Ret;
}

MATRIX*	AllocMatrix(int NoOfRows, int NoOfCols)
{
	MATRIX* Ret=NULL;
	int		RIndex=0;
	int		CIndex=0;
	
	Ret = (MATRIX*)malloc(sizeof(MATRIX));
	if(Ret == NULL)
		MallocErr();

	Ret->NoOfCols = NoOfCols;
	Ret->NoOfRows = NoOfRows;

	Ret->me = AllocMatMem(NoOfRows, NoOfCols);

	for(RIndex=0;RIndex<NoOfRows;RIndex++)
		for(CIndex=0;CIndex<NoOfCols;CIndex++)
			Ret->me[RIndex][CIndex] = 0;
	
	return Ret;
}

void	FreeMatrix(MATRIX *Matrix)
{
	free(Matrix->me[0]);
	free(Matrix->me);
	free(Matrix);
}

void	CopyMatrix(MATRIX *A, MATRIX *B)
{
	int	Total;

	Total = B->NoOfCols * B->NoOfRows;

	memcpy(A->me[0], B->me[0], Total * sizeof(double));

	A->NoOfCols = B->NoOfCols;
	A->NoOfRows = B->NoOfRows;
}

void	PrintMatrix(MATRIX *Matrix, char* Headder, FILE*	Str)
{
	int	RIndex, CIndex;

	fprintf(Str, "Matrix: %s\n", Headder);

	fprintf(Str, "%d\t%d\n", Matrix->NoOfRows, Matrix->NoOfCols);

	for(RIndex=0;RIndex<Matrix->NoOfRows;RIndex++)
	{
		for(CIndex=0;CIndex<Matrix->NoOfCols;CIndex++)
			fprintf(Str, "%10.10f\t", Matrix->me[RIndex][CIndex]);
		fprintf(Str, "\n");
	}


}


void	PrintMathematicaMatrix(MATRIX *Matrix, char* Headder, FILE*	Str)
{
	int	RIndex, CIndex;

	fprintf(Str, "%s", Headder);

	fprintf(Str, "{");
    
	for(RIndex=0;RIndex<Matrix->NoOfRows;RIndex++)
	{
		fprintf(Str, "{");
		for(CIndex=0;CIndex<Matrix->NoOfCols;CIndex++)
		{
			fprintf(Str, "%10.10f", Matrix->me[RIndex][CIndex]);
			if(CIndex!=Matrix->NoOfCols-1)
				fprintf(Str, ",");

		}
		fprintf(Str, "}");

		if(RIndex!=Matrix->NoOfRows-1)
			fprintf(Str, ",");
		fprintf(Str, "\n");
	}
	
	fprintf(Str, "};\n");

}

void	KroneckerProduct(MATRIX *A, MATRIX *B, MATRIX *C)
{
	int X,Y;
	int	XOut, YOut;
	int	XOffSet, YOffSet;

	for(XOut=0;XOut<A->NoOfCols;XOut++)
	{
		for(YOut=0;YOut<A->NoOfRows;YOut++)
		{
			XOffSet = (XOut * B->NoOfCols);
			YOffSet = (YOut * B->NoOfRows);
			
			for(X=0;X<B->NoOfCols;X++)
			{
				for(Y=0;Y<B->NoOfRows;Y++)
				{
					C->me[YOffSet + Y][XOffSet + X] = A->me[YOut][XOut] * B->me[Y][X];
				}
			}
		}
	}
}

void	VectByMatrixMult(double *Vect, MATRIX *Mat, double *Ret)
{
	int		x,y;
	double	Temp;
	
	for(y=0;y<Mat->NoOfCols;y++)
	{
		Temp = 0;
		for(x=0;x<Mat->NoOfRows;x++)	
			Temp = Temp + (Mat->me[x][y] * Vect[x]);

	
		Ret[y] = Temp;
	}
}

void	MatrixByVectMult(MATRIX* Mat, double *Vect, double *Ret)
{
	int		x,y;
	double	Temp;

	for(x=0;x<Mat->NoOfRows;x++)	
	{
		Temp = 0;
		for(y=0;y<Mat->NoOfCols;y++)
			Temp = Temp + (Mat->me[x][y] * Vect[y]);

		Ret[x] = Temp;
	}
}
/*
void	MatrixMult(MATRIX *A, MATRIX *B, MATRIX *Prod)
{
	int		Row, Col;
	int		Index;
	double	Temp;

	if(A->NoOfRows != Prod->NoOfRows)
	{
		printf("1. Matrix is not of the correct order to multiple\n");
		exit(0);
	}

	if(B->NoOfCols != Prod->NoOfCols)
	{
		printf("2. Matrix is not of the correct order to multiple\n");
		exit(0);
	}

	if(A->NoOfCols!= B->NoOfRows)
	{
		printf("3. Matrix is not of the correct order to multiple\n");
		exit(0);
	}

	for(Row=0;Row<A->NoOfRows;Row++)
	{
		for(Col=0;Col<A->NoOfCols;Col++)
		{
			Temp = 0;
			for(Index=0;Index<B->NoOfCols;Index++)
			{
				Temp += A->me[Row][Col] * B->me[Row][Index];
			}
			Prod->me[Row][Col] = Temp;
		}
	}	
}
*/


void	MatrixMult(MATRIX *A, MATRIX *B, MATRIX *Prod)
{
	int		k,i,j;

	if(A->NoOfRows != Prod->NoOfRows)
	{
		printf("1. Matrix is not of the correct order to multiple\n");
		exit(0);
	}

	if(B->NoOfCols != Prod->NoOfCols)
	{
		printf("2. Matrix is not of the correct order to multiple\n");
		exit(0);
	}
	
	for (k=0; k<B->NoOfCols; k++)
	{
		for (i=0; i<A->NoOfRows; i++) 
		{
			Prod->me[i][k] = 0.0;
			for (j=0; j<A->NoOfCols; j++)
				Prod->me[i][k] = Prod->me[i][k] + A->me[i][j] * B->me[j][k];
		}
	}
/*

	for(Row=0;Row<B->NoOfRows;Row++)
	{
		for(Col=0;Col<B->NoOfCols;Col++)
		{
			Temp = 0;
			for(Index=0;Index<A->NoOfCols;Index++)
			{
			//	Temp += B->me[Row][Col] * A->me[Row][Index];
				Temp += B->me[Index][Col] * A->me[Row][Index];
			}
			Prod->me[Row][Col] = Temp;
		}
	}	*/
}



void	PrintMathematicaVect(double* Vect, int N, char* Banna, FILE *Str)
{
	int	Index;
	fprintf(Str, "%s", Banna);
	fprintf(Str, "{");
	for(Index=0;Index<N;Index++)
	{
		fprintf(Str, "%10.10f", Vect[Index]);
		if(Index!=N-1)
			fprintf(Str, ",");
	}

	fprintf(Str, "};\n");
}

void	SetIdentityMatrix(MATRIX *M)
{
	int x,y;

	for(x=0;x<M->NoOfRows;x++)
		for(y=0;y<M->NoOfCols;y++)
			M->me[x][y] = 0;

	for(x=0;x<M->NoOfRows;x++)
		M->me[x][x] = 1;

}

void	Transpose(MATRIX *T, MATRIX *TransT)
{
	int x,y;

	if((TransT->NoOfCols != T->NoOfRows) || (TransT->NoOfRows != T->NoOfCols))
	{
		printf("Cannot Transpose matrix, as the destinastion matrix is of the order\n");
		exit(0);
	}

	for(x=0;x<T->NoOfRows;x++)
	{
		for(y=0;y<T->NoOfCols;y++)
			TransT->me[y][x] = T->me[x][y];
	}
}

void	ScaleMatrix(MATRIX *M, double Scalar)
{
	int r,c;

	for(r=0;r<M->NoOfRows;r++)
		for(c=0;c<M->NoOfCols;c++)
			M->me[r][c] *= Scalar;
}

double	VectByVectMult(double *Vect1, double *Vect2, int Size)
{
	double Ret;
	int		Index;

	Ret = 0;
	for(Index=0;Index<Size;Index++)
		Ret += Vect1[Index] * Vect2[Index];

	return Ret;
}
