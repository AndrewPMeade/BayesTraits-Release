/*
*  BayesTriats 3.0
*
*  copyright 2017
*
*  Andrew Meade
*  School of Biological Sciences
*  University of Reading
*  Reading
*  Berkshire
*  RG6 6BX
*
* BayesTriats is free software: you can redistribute it and/or modify
* it under the terms of the GNU General Public License as published by
* the Free Software Foundation, either version 3 of the License, or
* (at your option) any later version.
* 
* This program is distributed in the hope that it will be useful,
* but WITHOUT ANY WARRANTY; without even the implied warranty of
* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
* GNU General Public License for more details.
* 
* You should have received a copy of the GNU General Public License
* along with this program.  If not, see <http://www.gnu.org/licenses/>
*
*/



#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <time.h>
#include <ctype.h>

#include "btlin_dmatrix.h"
#include "btlin_alg.h"

DMATRIX* btlin_allocDMATRIX(int r, int c) {
  DMATRIX* ret = (DMATRIX*)malloc(sizeof(DMATRIX));
  ret->nrows = r;
  ret->ncols = c;
  ret->m = (double*)malloc(sizeof(double)*r*c);
  return ret;
}


DMATRIX* btlin_randomDMATRIX(int r, int c,int max) {
  int i;
  DMATRIX* ret = btlin_allocDMATRIX(r,c);
//  srand(time(0));
  for(i = 0; i < r*c; i++) {
    ret->m[i] = (rand()%(2*max+1) - max)*1.0;

  }
  return ret;
}

void btlin_freeDMATRIX(DMATRIX* dm) {
  free(dm->m);
  free(dm);
}

void btlin_printDMATRIX(DMATRIX* dm) {
  int i,j, index;
  int nr, nc;
  nr = dm->nrows;
  nc = dm->ncols;
  printf("DMATRIX %d x %d\n",nr,nc);
  for(i=0; i < nr; i++) {
    index = i;
    printf("[%d] ",i);
    for(j=0; j < nc; j++) {
      printf("%lf ",dm->m[index]);
      index += nr;
    }
    printf("\n");
  }
}

double btlin_getDMATRIX(DMATRIX* dm, int i, int j) {
  return dm->m[j*(dm->nrows)+i];
}

void btlin_setDMATRIX(DMATRIX* dm, int i, int j, double v) {
  dm->m[j*(dm->nrows)+i] = v;
}

void btlin_copyMATRIXtoDMATRIX(double** m, DMATRIX* dm) 
{
  int i,j,nr,nc;
  nr = dm->nrows;
  nc = dm->ncols;
  for(i=0; i < nr; i++) {
    for(j=0; j < nc; j++) {
      dm->m[j*nr+i] = m[i][j];
    }
  }
}


void btlin_copyDMATRIXtoMATRIX(DMATRIX* dm, double** m) 
{
  int i,j,nr,nc;
  nr = dm->nrows;
  nc = dm->ncols;
  for(i=0; i < nr; i++) {
    for(j=0; j < nc; j++) {
      m[i][j] = dm->m[j*nr+i];
    }
  }
}

DMATRIX* btlin_copyallocDMATRIX(DMATRIX* dm) 
{
  DMATRIX* ret;
 
  ret = btlin_allocDMATRIX(dm->nrows,dm->ncols);
  // perhaps it would be better to perform this copying in chucks of e.g 1024
  memcpy(ret->m,dm->m,sizeof(double)*dm->nrows*dm->ncols);

  return ret;
}

// left := right (buffer values)
void btlin_copyDMATRIX(DMATRIX* left, DMATRIX* right) {
	int lsize, rsize, size;
	lsize = left->nrows*left->ncols;
	rsize = right->nrows*right->ncols;
	size = (lsize < rsize)?lsize:rsize;
	memcpy(left->m,right->m,sizeof(double)*size);
	return;
}

void btlin_copylowDMATRIX(DMATRIX* dm) 
{
  int i,j,nr,nc;
  nr = dm->nrows;
  nc = dm->ncols;
  for(i=0; i < nr; i++) {
    for(j=i+1; j < nc; j++) {
      dm->m[j*nr+i] = dm->m[i*nr+j];
    }
  }
}

void btlin_makeSymmetricDMATRIX(char uplo, DMATRIX* dm)
{
//	int col;
  int i,l,u,n,dim,colxnr;
  double *a = dm->m;
  n = dm->nrows;
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


/*  ****************** Algorithms *****************   */

// Cholesky unblocked
int btlin_choleskyDMATRIX(DMATRIX* dm, double* det) {
	return btlin_cholesky(dm->m, dm->nrows,det);
}

int btlin_invcholeskyDMATRIX(DMATRIX* dm, double *det) {
	return btlin_invcholesky(dm->m,dm->nrows,det);
}









