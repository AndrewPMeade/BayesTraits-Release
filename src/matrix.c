#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "genlib.h"
#include "matrix.h"
#include "linalg.h"

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

MATRIX**	AllocMultiMatrixLinMem(int NoMat, int NoRows, int NoCols)
{
	MATRIX **Ret;
	MATRIX *Mat;
	double *Mem;
	int		Index, RIndex;
	
	Ret = (MATRIX**)malloc(sizeof(MATRIX*) * NoMat);
	if(Ret == NULL)
		MallocErr();

	Mem = (double*)malloc(sizeof(double) * NoRows * NoCols * NoMat);
	if(Mem == NULL)
		MallocErr();	
	
	for(Index=0;Index<NoMat;Index++)
	{
		Ret[Index] = (MATRIX*)malloc(sizeof(MATRIX));
		if(Ret[Index] == NULL)
			MallocErr();

		Mat = Ret[Index];
		Mat->NoOfCols = NoCols;
		Mat->NoOfRows = NoRows;

		Mat->me = (double**)malloc(sizeof(double*) * NoRows);
		if(Mat->me == NULL)
			MallocErr();

		Mat->me[0] = &Mem[Index * NoRows * NoCols];
		for(RIndex=1;RIndex<NoRows;RIndex++)
			Mat->me[RIndex] = Mat->me[0] + (RIndex * NoCols);
	}

	return Ret;
}

void		FreeMultiMatrixLinMem(MATRIX **MList, int NoMat)
{
	MATRIX	*Mat;
	double *Mem;
	int		MIndex;

	Mem = MList[0]->me[0];

	for(MIndex=0;MIndex<NoMat;MIndex++)
	{
		Mat = MList[MIndex];

		free(Mat->me);
		free(Mat);
	}

	free(Mem);
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
			fprintf(Str, "%30.30f\t", Matrix->me[RIndex][CIndex]);
		fprintf(Str, "\n");
	}
}

void	PrintMatrixBinary(MATRIX *Matrix, char* Headder, FILE*	Str)
{
	int	RIndex, CIndex;

	fprintf(Str, "Matrix: %s\n", Headder);
	fprintf(Str, "%d\t%d\n", Matrix->NoOfRows, Matrix->NoOfCols);

	for(RIndex=0;RIndex<Matrix->NoOfRows;RIndex++)
	{
		for(CIndex=0;CIndex<Matrix->NoOfCols;CIndex++)
		{
			PrintDoubleHex(Str, Matrix->me[RIndex][CIndex]);
			fprintf(Str, "\t");
		}
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

void	PrintMathematicaTFMatrix(MATRIX *Matrix, char* Headder, FILE*	Str)
{
	int	RIndex, CIndex;

	fprintf(Str, "%s", Headder);

	fprintf(Str, "{");
    
	for(RIndex=0;RIndex<Matrix->NoOfRows;RIndex++)
	{
		fprintf(Str, "{");
		for(CIndex=0;CIndex<Matrix->NoOfCols;CIndex++)
		{
			fprintf(Str, "%10.10f", log(1+Matrix->me[RIndex][CIndex]));
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

// Helper function used by VectByKroneckerMult below.
// Assumptions: sigma and mat are square matrices, column-major.
// However, since sigma and mat are both symmetric, this is not a problem - both (linear)buffers are the same.
// vres: pointer to result vector buffer, size = sigma_dim*mat_dim
// v: pointer to vector buffer, size = sigma_dim*mat_dim
// sigma: pointer to buffer storing a sigma_dim by sigma_dim matrix
// 
// vres = v x (sigma kx mat)
int kronecker_vectmultColumn(double* vres, double* v, double* sigma, int sigma_dim,double*  mat, int mat_dim) {
	
	double *pmat;	
	int mat_col, mat_row;
	double *psigma;
	double *pvres_start, *pvres;
	double *pv, *pv_start;
	int sumacc_idx, extra_idx;
	double* sumacc, a, sum;
//	double *vstart;
	
	sumacc = (double*)malloc(sizeof(double)*sigma_dim);
		
	pmat = mat;
	pvres_start = vres;
	
	for(mat_col=0; mat_col < mat_dim; mat_col++) {	
		for(sumacc_idx=0; sumacc_idx < sigma_dim; sumacc_idx++) 
			sumacc[sumacc_idx] = 0.0;
		pv_start = v;
		for(mat_row=0; mat_row < mat_dim; mat_row++) {
			a = *pmat;
			pv = pv_start;
			for(sumacc_idx=0; sumacc_idx < sigma_dim; sumacc_idx++) {	
				sumacc[sumacc_idx] += a*(*pv);
				pv += mat_dim;
			}	
			pv_start++;
			pmat++;
		}
				
		// update result
		pvres = pvres_start;
		psigma = sigma; 
		for(extra_idx=0; extra_idx < sigma_dim; extra_idx++) {
			sum = 0.0;
			for(sumacc_idx=0; sumacc_idx < sigma_dim; sumacc_idx++) {
				sum+= *psigma*sumacc[sumacc_idx];
				psigma++;
			}
			*pvres = sum;
			pvres += mat_dim;
		}	
		pvres_start++;
	}
	
	
	free(sumacc);
	return 0;
}

// Simpler version (column major) for sigma_dim = 1
int kronecker_vectmultColumnOne(double* vres, double* v, double* sigma, int sigma_dim,double*  mat, int mat_dim) {
	int mat_col, mat_row;     // included for clarity - no need
	double *pvres, *pv;
//	double sum;
	double *pmat;
	
	double sigma_constant;
	
	if (sigma_dim > 1) {
		return kronecker_vectmultColumn(vres, v, sigma, sigma_dim, mat, mat_dim);
	}
	
	sigma_constant = *sigma;	
	pvres = vres;
	pmat = mat;		
	for(mat_col=0; mat_col < mat_dim; mat_col++) {	
		pv = v;
		*pvres = 0.0;	
		for(mat_row = 0; mat_row < mat_dim; mat_row++) {
			*pvres += (*pmat)*(*pv);
			pv++;
			pmat++;
		}
		(*pvres) *= sigma_constant;
		pvres++;
	}
	return 0;
}

void	VectByKroneckerMult(double *Vect, MATRIX *A, MATRIX* B, double *Ret) 
{
	if (A->NoOfRows == 1)
		kronecker_vectmultColumnOne(Ret, Vect, A->me[0], A->NoOfRows, B->me[0], B->NoOfRows);
	else
		kronecker_vectmultColumn(Ret, Vect, A->me[0], A->NoOfRows, B->me[0], B->NoOfRows);
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

void	MatrixSelfMult(MATRIX *A, MATRIX *Prod)
{
	int		k,i,j;

	for (k=0; k<A->NoOfRows; k++)
	{
		for (i=0; i<A->NoOfRows; i++) 
		{
			Prod->me[i][k] = 0.0;
			for (j=0; j<A->NoOfCols; j++)
				Prod->me[i][k] += A->me[i][j] * A->me[k][j];
		}
	}	
}

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
				Temp += B->me[Index][Col] * A->me[Row][Index];
			
			Prod->me[Row][Col] = Temp;
		}
	}	
*/
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

double*	FlattenMatix(MATRIX *M)
{
	double *Ret;
	int		x,y,i;

	i = 0;

	Ret = (double*)malloc(sizeof(double) * (M->NoOfRows * M->NoOfCols));
	if(Ret == NULL)
		MallocErr();

	for(x=0;x<M->NoOfRows;x++)
	{
		for(y=0;y<M->NoOfCols;y++)
			Ret[i++] = M->me[x][y];
	}

	return Ret;
}

void	FillMatrix(MATRIX *M, double Value)
{
	int x,y;

	for(x=0;x<M->NoOfRows;x++)
		for(y=0;y<M->NoOfCols;y++)
			M->me[x][y] = Value;
	
}

MATRIX_INVERT*	CreatMatInvertInfo(int N)
{
	MATRIX_INVERT*	Ret;

	Ret = (MATRIX_INVERT*)malloc(sizeof(MATRIX_INVERT));
	if(Ret == NULL)
		MallocErr();

	Ret->Inv = AllocMatrix(N, N);
	Ret->TempD = (double*)malloc(sizeof(double) * N);
	Ret->TempI = (int*)malloc(sizeof(double) * N);

	return Ret;
}

void			FreeMatInvertInfo(MATRIX_INVERT *M)
{
	FreeMatrix(M->Inv);
	free(M->TempD);
	free(M->TempI);
}

int				Matrix_Invert(MATRIX *Mat, MATRIX_INVERT *Info)
{
	int Ret;

	Ret = InvertMatrixAndDet(Mat->me, Mat->NoOfCols, Info->TempD, Info->TempI, Info->Inv->me, &Info->Det);
	Info->Det = exp(Info->Det);
	
	return Ret;
}