#ifndef BTLIN_DMATRIX_H
#define BTLIN_DMATRIX_H

/* Double Matrix - Column-major */
typedef struct {
  double* m;
  int nrows;
  int ncols;
} DMATRIX;


// ************************************************************************


DMATRIX* btlin_allocDMATRIX(int r, int c);
DMATRIX* btlin_randomDMATRIX(int r, int c,int max);
DMATRIX* btlin_copyallocDMATRIX(DMATRIX* dm);
void btlin_freeDMATRIX(DMATRIX* m);
void btlin_printDMATRIX(DMATRIX* m);

double btlin_getDMATRIX(DMATRIX* m, int i, int j);
void btlin_setDMATRIX(DMATRIX* m, int i, int j, double v);

void btlin_copyMATRIXtoDMATRIX(double** m,DMATRIX* dm);
void btlin_copyDMATRIXtoMATRIX(DMATRIX* dm,double** m);
void btlin_copyDMATRIX(DMATRIX* left, DMATRIX* right);
void btlin_copylowDMATRIX(DMATRIX* dm);

void btlin_makeSymmetricDMATRIX(char uplo, DMATRIX* dm);


/* ******** Algorithms ********* */
int btlin_choleskyDMATRIX(DMATRIX* m, double* det);
int btlin_invcholeskyDMATRIX(DMATRIX* m, double *det);

#endif