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



#include "btlapack.h"
#include "typedef.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <time.h>

#ifdef BTLAPACK

#include "btlapack.h"
#include "btlin_alg.h"

#include "btdebug.h"


// square matrix
int btlapack_cholesky(double* a, int n, double* det) {
  int info;
  char uplo;

  uplo = 'L';
  dpotrf_(&uplo, &n, a, &n, &info);

  if (det != 0) {
	*det = btlapack_ltri_det(a,n,n);
	*det *= 2.0;
  }

  return info;
}


int btlapack_invcholesky(double* a, int n, double *det) {
	int info;
	char uplo;

  //btdebug_enter("btlapackchfact");
  
  uplo = 'L';
  dpotrf_(&uplo, &n, a, &n, &info);
  
  //btdebug_exit("btlapackchfact");
  
  //btlapack_printDMATRIX(dm);

  if (info!=0) return info;

  if (det != 0) {
	*det = btlapack_ltri_det(a,n,n);
	*det *= 2.0;
  }

  
  // proceed with dpotri
  //btdebug_enter("btlapackcholeskysub");
  dpotri_(&uplo, &n, a, &n, &info);
  //btdebug_exit("btlapackcholeskysub");
  
  //btlapack_copylowDMATRIX(dm);
  
  btlin_makeSymmetric('L',a,n);

  return info;


}




int btlapack_choleskyU(double* a, int n, double* det) {
  int info;
  char uplo;

  uplo = 'U';
  dpotrf_(&uplo, &n, a, &n, &info);
  

  if (det != 0) {
	*det = btlapack_ltri_det(a,n,n);
	*det *= 2.0;
  }
  
  return info;
}


int btlapack_cholesky2(double* a, int n, double* det) {
  int info;
  char uplo;

  uplo = 'L';
  dpotf2_(&uplo, &n, a, &n, &info);

  if (det != 0) {
	*det = btlapack_ltri_det(a,n,n);
	*det *= 2.0;
  }
  
  return info;
}



int btlapack_invlu(double* a, int n, double* det) {
  int *IPIV = (int*)malloc(sizeof(int)*(n+1));
  int LWORK = n*n;
  double *WORK = (double*)malloc(sizeof(double)*LWORK);
  int info,i;

  dgetrf_(&n,&n,a,&n,IPIV,&info);

  //printf("LU decomposition:\n");
  
  // btlapack_printDMATRIX(dm);

  if (det != NULL) {
	*det = 0.0;
	//printf("Diag: ");
	for(i=0;i<LWORK;i+=(n+1)) {
		printf("%10.10lf ",a[i]);
		if (a[i] < 0)
			*det += log(-a[i]);
		else
			*det += log(a[i]);
	}
  //printf("\n");
  }
  
  dgetri_(&n,a,&n,IPIV,WORK,&LWORK,&info);

  free(IPIV);
  free(WORK);

  return info;
}




int btlapack_ldl(double* a, int n, double* det) {
  int *IPIV = (int*)malloc(sizeof(int)*(n+1));
  int LWORK = n*n;
  double *WORK = (double*)malloc(sizeof(double)*LWORK);
  int info,i,rc;
  double detblock,c;
  char uplo;

  uplo = 'L';

  //extern void dsytrf_(char* uplo, int* N, double* A, int* lda, int* IPIV, double* WORK, int* LWORK, int* INFO);

  dsytrf_(&uplo, &n, a, &n,IPIV,WORK,&LWORK, &info);

  //if (info==0) {
  //  printf("Factored matrix\n");
  //  btlapack_printDMATRIX(dm);
  //  printf("P = [");
  //  for(i=0; i < n; i++)
  //    printf(" %d",IPIV[i]);
  //  printf(" ]\n");
  //} else {
  //  printf("ERROR with LDL factorization\n");
  //}
  
  // Calculating the determinant
    if (det != NULL) {
	*det = 0.0;
	rc = 0; // row/column
	for(i=0;i<n*n;) {
		//    printf("i=%d n=%d ",i,n);
		if (IPIV[rc] < 0) { // 2-by-2 pivot block
			detblock = a[i];
			c = a[i+1];  // symmetric
			i += (n+1); rc++;
			detblock = detblock*a[i] - c*c; // a*d - b*c 
		} else {
			detblock = a[i];
    }
    i+= (n+1); rc++;
    if (detblock < 0)
      *det += log(-detblock);
    else
      *det += log(detblock);
	}
  }
  

  free(IPIV);
  free(WORK);

  return info;
}


int btlapack_invldl(double* a, int n, double* det) {
  int *IPIV = (int*)malloc(sizeof(int)*(n+1));
  int LWORK = n*n;
  double *WORK = (double*)malloc(sizeof(double)*LWORK);
  int info,i,rc;
  double detblock,c;
  char uplo;

  uplo = 'L';

  dsytrf_(&uplo, &n, a, &n,IPIV,WORK,&LWORK, &info);

    //printf("Factored matrix\n");
    //btlapack_printDMATRIX(dm);
    //printf("P = [");
    //for(i=0; i < n; i++)
    //printf(" %d",IPIV[i]);
    //printf(" ]\n");

  if (info!=0) {
    printf("ERROR with LDL factorization\n");
  }

  if (det != NULL) {
	*det = 0.0;
	rc = 0; // row/column
	for(i=0;i<n*n;) {
		//    printf("i=%d n=%d ",i,n);
		if (IPIV[rc] < 0) { // 2-by-2 pivot block
			detblock = a[i];
			c = a[i+1];  // symmetric
			i += (n+1); rc++;
			detblock = detblock*a[i] - c*c; // a*d - b*c 
		} else {
			detblock = a[i];
    }
    i+= (n+1); rc++;
    if (detblock < 0)
      *det += log(-detblock);
    else
      *det += log(detblock);
	}
  }

  // now do the inverse.
  dsytri_(&uplo, &n, a, &n, IPIV, WORK, &info);

  btlin_makeSymmetric('L',a,n);

  free(IPIV);
  free(WORK);

  return info;
}



int btlapack_invldlW(double* a, double* w, int n, double* det) {
  int *IPIV = (int*)malloc(sizeof(int)*(n+1));
  int LWORK = n*n;
  //double *WORK = (double*)malloc(sizeof(double)*LWORK);
  int info,i,rc;
  double detblock,c;
  char uplo;

  uplo = 'L';

  btdebug_enter("lapackfactor");
  dsytrf_(&uplo, &n, a, &n,IPIV,w,&LWORK, &info);
  btdebug_exit("lapackfactor");

    //printf("Factored matrix\n");
    //btlapack_printDMATRIX(dm);
    //printf("P = [");
    //for(i=0; i < n; i++)
    //printf(" %d",IPIV[i]);
    //printf(" ]\n");

  if (info!=0) {
    printf("ERROR with LDL factorization\n");
  }

  if (det != NULL) {
	*det = 0.0;
	rc = 0; // row/column
	for(i=0;i<n*n;) {
		//    printf("i=%d n=%d ",i,n);
		if (IPIV[rc] < 0) { // 2-by-2 pivot block
			detblock = a[i];
			c = a[i+1];  // symmetric
			i += (n+1); rc++;
			detblock = detblock*a[i] - c*c; // a*d - b*c 
		} else {
			detblock = a[i];
    }
    i+= (n+1); rc++;
    if (detblock < 0)
      *det += log(-detblock);
    else
      *det += log(detblock);
	}
  }

  // now do the inverse.
  btdebug_enter("lapackfinish");
  dsytri_(&uplo, &n, a, &n, IPIV, w, &info);
  btdebug_exit("lapackfinish");

  btlin_makeSymmetric('L',a,n);

  free(IPIV);
  //free(WORK);

  return info;
}






// new naming/argument style
int btlapack_ducholesky(char uplo,double* a,int n,int lda, double* det) {
	int info; 
	
	dpotf2_(&uplo, &n, a, &lda, &info);
	return info;
}



int btlapack_ducholeskydet(char uplo,double* a,int n,int lda, double* det) {
	int info; 
	//double* pdiag;
	
	//printf("btlapack_ducholesky corner %lf n %d lda %d  \n",*a, n,lda);
	
	dpotf2_(&uplo, &n, a, &lda, &info);


	

	//*det = 0.0;
	// warning: use lda instead of n
	//pdiag = a;
	//for(i=0; i<n; i++) {
	//	if (*pdiag < 0)
	//		*det += log(-*pdiag);
	//	else
	//		*det += log(*pdiag);
	//	pdiag += (lda+1);
	//}

	if (det != NULL) {
		*det = btlapack_ltri_det(a,n,lda);
		*det *= 2.0;
	}
	
	return info;
}


// Inverse lower triangular
int btlapack_ltri_inv2(double* a, int n, int lda, double* det) {
	int info;
	char uplo,diag;

	diag = 'N';
	uplo = 'L';

	// dtrti2_(char* uplo, char* diag, int *n, double* A, int* lda, int* INFO);
	dtrti2_(&uplo,&diag,&n,a,&lda,&info);

	if (det != NULL) {
		*det = btlapack_ltri_det(a,n,lda);
	}

	return info;
}

int btlapack_ltri_inv(double* a, int n, int lda, double* det) {
	int info;
	char uplo,diag;

	diag = 'N';
	uplo = 'L';

	// dtrti2_(char* uplo, char* diag, int *n, double* A, int* lda, int* INFO);
	dtrtri_(&uplo,&diag,&n,a,&lda,&info);

	if (det != NULL) {
		*det = btlapack_ltri_det(a,n,lda);
	}

	return info;

	
}

#endif  // if BTLAPACK defined

// computes the log of the determinant of a lower triangular matrix
double btlapack_ltri_det(double* a, int n, int lda) 
{
	double  det, *pdiag;
	int i;

	det = 0.0;
	pdiag = a;
	for(i=0; i<n; i++) 
	{
		if (*pdiag < 0)
			det += log(-*pdiag);
		else
			det += log(*pdiag);
		pdiag += (lda+1);
	}

	return det;
}


