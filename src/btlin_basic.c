#include "btlin_basic.h"

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <ctype.h>

// Copy
double* btlin_matrix_newcopy(double* a, int m, int n) {
	double *ret,*p;
	int i,size;
	size = m*n;
	ret = (double*)malloc(sizeof(double)*size);
	// Alternatives: Copy in chunks (2^n blocks)
	p = ret;
	for(i=0; i < size; i++) {
		(*p) = (*a);
		p++; a++;
	}
		//ret[i] = a[i];
	return ret;
}

// Copy Lower triangle to upper triangle of square matrix
void btlin_copylow(double* m, int n) 
{
  int i,j;
  for(i=0; i < n; i++) {
    for(j=i+1; j < n; j++) {
      m[j*n+i] = m[i*n+j];
    }
  }
}

// Print matrix (column major)
void btlin_print(double* m, int nr, int nc) {
  int i,j, index;
  
  printf("MATRIX %d x %d\n",nr,nc);
  for(i=0; i < nr; i++) {
    index = i;
    printf("[%d] ",i);
    for(j=0; j < nc; j++) {
      printf("%lf ",m[index]);
      index += nr;
    }
    printf("\n");
  }
}

// Print matrix (row major)
void btlin_printR(double* m, int nr, int nc) {
  int i,j, index;
  
  printf("MATRIX %d x %d\n",nr,nc);
  index=0;
  for(i=0; i < nr; i++) {
    printf("[%d] ",i);
    for(j=0; j < nc; j++) {
      printf("%lf ",m[index]);
      index++;
    }
    printf("\n");
  }
}

void btlin_makeSymmetric(char uplo, double* a, int n) {
  int i,l,u,col,dim,colxnr;

  dim = n*n;
  uplo = tolower(uplo);
  colxnr = n;
  if (uplo == 'l') {  // copies lower triangular to upper triangulr
    for (i=1; i < dim; i+=(n+1),colxnr+=n) {
      for(l=i,u=i+n-1; l < colxnr; l++, u+=n) {
	a[u]= a[l];
      } 
    }
  } else { // copy upper to lower
    for (i=1; i < dim; i+=(n+1),colxnr+=n) {
      for(l=i,u=i+n-1; l < colxnr; l++, u+=n) {
	a[l]= a[u];
      } 
    }
  }
}

// r = transpose(s)
void btlin_transpose(double* r, double* s, int n) {
	int i,j;
	double *pr, *ps, *ps_start, temp;
	pr = r;
	ps_start = s;
	for(i=0; i < n; i++) {
		ps = ps_start;
		for(j=0; j < n; j++) {
			pr = ps;
			pr++;
			ps += n;
		}
		ps_start++;
	}
}




