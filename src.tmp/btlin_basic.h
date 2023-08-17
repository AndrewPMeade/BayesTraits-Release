#ifndef BTLIN_H
#define BTLIN_H

double* btlin_matrix_newcopy(double* a, int m, int n);
void btlin_copylow(double* m, int n);
void btlin_print(double* m, int nrows, int ncols);
void btlin_printR(double* m, int nrows, int ncols);
void btlin_makeSymmetric(char uplo, double* a, int n);

void btlin_transpose(double* r, double* s, int n);

#endif