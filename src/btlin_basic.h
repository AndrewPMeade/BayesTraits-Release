#ifndef BTLIN_H
#define BTLIN_H


void btlin_copylow(double* m, int n);
void btlin_print(double* m, int nrows, int ncols);
void btlin_printR(double* m, int nrows, int ncols);
void btlin_makeSymmetric(char uplo, double* a, int n);



#endif