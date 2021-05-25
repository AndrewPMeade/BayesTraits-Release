#if !defined (MATRIXHEADDER)
#define MATRIXHEADDER

#include <stdio.h>

typedef	struct
{
	double**	me;
	int			NoOfCols;
	int			NoOfRows;
} MATRIX;

MATRIX*	AllocMatrix(int NoOfRows, int NoOfCols);
void	FreeMatrix(MATRIX *Matrix);
void	PrintMatrix(MATRIX *Matrix, char* Headder, FILE*	Str);
void	PrintMatrixBinary(MATRIX *Matrix, char* Headder, FILE*	Str);
void	CopyMatrix(MATRIX *A, MATRIX *B);
double** AllocMatMem(int NoOfRows, int NoOfCols);
void	FreeMatMem(double** Mat);
void	KroneckerProduct(MATRIX *A, MATRIX *B, MATRIX *C);
void	VectByMatrixMult(double *Vect, MATRIX *Mat, double *Ret);
void	MatrixByVectMult(MATRIX* Mat, double *Vect, double *Ret);
double	VectByVectMult(double *Vect1, double *Vect2, int Size);
void	MatrixMult(MATRIX *A, MATRIX *B, MATRIX *Prod);
void	SetIdentityMatrix(MATRIX *M);
void	Transpose(MATRIX *T, MATRIX *TransT);
void	ScaleMatrix(MATRIX *M, double Scalar);

void	FillMatrix(MATRIX *M, double Value);

void	PrintMathematicaMatrix(MATRIX *Matrix, char* Headder, FILE*	Str);
void	PrintMathematicaVect(double* Vect, int N, char* Banna, FILE *Str);

void	PrintMathematicaTFMatrix(MATRIX *Matrix, char* Headder, FILE* Str);

double*	FlattenMatix(MATRIX *M);

#endif
