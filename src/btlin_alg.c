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



#include "btlin_alg.h"

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "btdebug.h"

/* **** Local Auxiliary **** */
int btlin_dbcholupdcol(double* pcol,int next_n,double* pdiag,int nb,int lda);
int btlin_dbcholupdmat(double* pmat,int next_n,double* pcol,int nb, int lda);

// A <- alpha*L*M
int btlin_ltri_LowerByMatrix(double alpha, double* pL,double* pA,int m,int n,int lda);
// A <- alpha*M*L
int btlin_ltri_MatrixByLower(double alpha, double* pA,double* pL,int m,int n,int lda);

// uses column major traversal
// size of v = size of vres = sigma_dim * mat_dim
// Special case: for sigma_dim = 1
int btlin_kronecker_vectmultOne(double* vres, double* v, double* sigma, int sigma_dim,double*  mat, int mat_dim) {
	int mat_col, mat_row;     // included for clarity - no need
	double *pvres, *pv;
	double sum;
	double *pmat;
	
	double sigma_constant;
	
	if (sigma_dim > 1) {
		return btlin_kronecker_vectmult(vres, v, sigma, sigma_dim, mat, mat_dim);
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



// uses column major traversal
// size of v = size of vres = sigma_dim * mat_dim
int btlin_kronecker_vectmult(double* vres, double* v, double* sigma, int sigma_dim,double*  mat, int mat_dim) {
	int sigma_col, sigma_row; // included for clarity - no need
	int mat_col, mat_row;     // included for clarity - no need
	int vres_idx, v_idx;      // can be replaced with pointers
	double sum, subsum;
	double *pmat,*pmat_start;
	double *psigma, *psigma_start;
	
	vres_idx = 0;
	psigma_start = sigma;
	for(sigma_col=0; sigma_col < sigma_dim; sigma_col++) {
		pmat_start = mat;
		for(mat_col=0; mat_col < mat_dim; mat_col++) {	
			v_idx = 0;
			sum = 0.0;
			psigma = psigma_start;
			for(sigma_row=0; sigma_row < sigma_dim; sigma_row++) {
				subsum = 0.0;
				// compute subsum
				pmat = pmat_start;
				for(mat_row = 0; mat_row < mat_dim; mat_row++) {
					// subsum += mat[mat_row][mat_col]*v[v_idx];
					//printf("sub += %lf * %lf ",*pmat,v[v_idx]);
					subsum += (*pmat)*(v[v_idx]);
					v_idx++;
					pmat++;
				}
				// sum += subsum* sigma[sigma_row][sigma_col];
				//printf("subsum*%lf\n",*psigma);
				sum += subsum*(*psigma);
				psigma++;
			}
			vres[vres_idx] = sum;
			vres_idx++;
			pmat_start += mat_dim;
		}
			
		psigma_start += sigma_dim;
	}	

	return 0;
}

// Assumptions: sigma and mat are square matrices, column-major.
// This version assumes that the matrices use a column-major representation. 
// However, since sigma and mat
// are both symmetric, this is not a problem - both (linear)buffers are the same.
// vres: pointer to result vector buffer, size = sigma_dim*mat_dim
// v: pointer to vector buffer, size = sigma_dim*mat_dim
// sigma: pointer to buffer storing a sigma_dim by sigma_dim matrix
// 
// vres = v x (sigma kx mat)
int btlin_kronecker_vectmult2(double* vres, double* v, double* sigma, int sigma_dim,double*  mat, int mat_dim) {
	
	double *pmat;	
	int mat_col, mat_row;
	double *psigma;
	double *pvres_start, *pvres;
	double *pv, *pv_start;
	int sumacc_idx, extra_idx;
	double* sumacc, a, sum;
	double *vstart;
	
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

// This version uses row-major layout.
// used to compare against column-major version.
int btlin_kronecker_vectmult2Row(double* vres, double* v, double* sigma, int sigma_dim,double*  mat, int mat_dim) {
	
	double *pmat, *pmat_start;	
	int mat_col, mat_row;
	double *psigma, *psigma_start;
	double *pvres_start, *pvres;
	double *pv, *pv_start;
	int sumacc_idx, extra_idx;
	double* sumacc, a, sum;
	double *vstart;
	
	sumacc = (double*)malloc(sizeof(double)*sigma_dim);
		
	pmat_start = mat;
	pvres_start = vres;
	
	for(mat_col=0; mat_col < mat_dim; mat_col++) {	
		for(sumacc_idx=0; sumacc_idx < sigma_dim; sumacc_idx++) 
			sumacc[sumacc_idx] = 0.0;
		pv_start = v;
		pmat = pmat_start;
		for(mat_row=0; mat_row < mat_dim; mat_row++) {
			a = *pmat;
			pv = pv_start;
			for(sumacc_idx=0; sumacc_idx < sigma_dim; sumacc_idx++) {	
				sumacc[sumacc_idx] += a*(*pv);
				pv += mat_dim;
			}	
			pv_start++;
			pmat += mat_dim;  // add lda
		}
		pmat_start++;
				
		// update result
		pvres = pvres_start;
		psigma_start = sigma; 
		for(extra_idx=0; extra_idx < sigma_dim; extra_idx++) {
			sum = 0.0;
			psigma = psigma_start;
			for(sumacc_idx=0; sumacc_idx < sigma_dim; sumacc_idx++) {
				sum+= *psigma*sumacc[sumacc_idx];
				psigma += sigma_dim;
			}
			psigma_start++;
			*pvres = sum;
			pvres += mat_dim;
			
		}	
		pvres_start++;
	}
	
	free(sumacc);
	return 0;
}



// Cholesky decomposition
// Standard algorithm - assuming column-major and lower triangular
// obscure due to one-dimensional indexes
int btlin_cholesky(double* m, int n, double* det) {
	int i,j,k,diagidx, idx1, idx2, idx, subdiagidx;
	double diag, a, sum;
	*det = 0;
	// may want to check for squareness
	
	diagidx = 0;
	for(i=0; i < n; i++) {	
		diag = m[diagidx];
		if (diag < 0) {
			printf("Error index %d negative diagonal %lf\n",diagidx,diag);		
			return 1;
		}
		diag = m[diagidx] = sqrt(diag);	
		if (diag < 0)
			*det += log(-diag);
		else
			*det += log(diag);
		// update column		
		idx = diagidx;
		for(j=i+1; j < n; j++) {
			idx++;			
			m[idx] = m[idx]/diag;
		}	
		// update the remaining submatrix
		
		// jump to next diagonal element
		idx1 = diagidx;
		diagidx += n+1; // update diagonal diagonal
		subdiagidx = diagidx; // diagonal indexes of submatrix (next loop)
		for(j=i+1;j < n; j++) {
			idx = subdiagidx;
			idx1++;
			idx2 = idx1;
			a = m[idx1];
			for(k=j; k < n; k++) { 				
				m[idx] -= a*m[idx2];				
				idx2++;
				idx++;
			}
			subdiagidx += n+1;
		}
	}
	*det *= 2.0;
	return 0;
}



int btlin_invcholesky(double* m, int n, double* det) {
	int col,row, diagidx, rowidx, colidx, k, kidx;
	double diag, sum, *pdiag;
	// assert m->nrows=m->ncols

	
	btdebug_enter("btlincholeskyfact");
	btlin_cholesky(m,n,det);
	// build L\L**T (symmetric)
	btdebug_exit("btlincholeskyfact");
	btdebug_enter("btlincholeskycopy");
	btlin_copylow(m,n);  // make it symmetric
	btdebug_exit("btlincholeskycopy");
	//printf("After Cholesky:\n");
	//btlapack_printDMATRIX(dm);
	
	btdebug_enter("btlincholeskysub");
	
	// first substitution - compute X into lower triangle
	// LX = I
	// traverse each column from top to bottom (start from leftmost in seq version).
	diagidx = 0;
	for(col=0; col < n; col++) {
		//printf("COLUMN %d diagidx %d\n",col,diagidx);
		diag = 1/m[diagidx];  // compute but don't store
		sum = 0;
		for(row=col+1; row < n; row++) {
			//printf("  row %d\n", row);
			// go to column col
			rowidx = n*row+col; // row in L. Traversed in its transpose (column down)
			// sum=diag*l(row,k) = diag*(k,row)
			sum = m[rowidx] / m[diagidx];
			//printf("sum = %lf = %lf / %lf\n",sum, m[rowidx],m[diagidx]);
			colidx = diagidx;
			// compute sum of product: l(row,k)*x(k,col),...,l(row,row)*x(row,col)
			rowidx++;
			colidx++;
			for(k = 0; k < row-col-1; k++) {				
				sum += m[colidx]*m[rowidx];
				//printf("sum = %lf += %lf * %lf\n",sum, m[colidx],m[rowidx]);
				rowidx++;
				colidx++;
			}
			// for k = row, conpute x(row,col)
			if (colidx != diagidx) {  // do not update the diagonal
				//printf("    updating %d rowidx %d\n",colidx,rowidx);
				m[colidx] = -sum / m[rowidx];
			}
			// note: - we don't need to divide by m[rowidx], just store sum.
			// we know that in the final matrix, we can get the value by
			// dividing  m[row,col] by m[row,row]
			// plus the real value in the diag is 1/m[x,x]
			// Check if this can be useful
		}
		diagidx += n+1;
	}
	
	//printf("After first substitution:\n");
	//btlapack_printDMATRIX(dm);
	
	// second substitution - compute inverse
	// column: left to right. rows: bottom-up
	// this approach may access L**T inefficiently but it may be the right one
	// for the GPU version
	for(col=0; col < n; col++) {	
		for(row=n-1; row >= col; row--) {
			// compute sum
			sum = 0;
			rowidx = n*(n-1)+row;
			colidx = (col+1)*n - 1;
			for(k=n-1; k >row; k--) {
				sum += m[rowidx]*m[colidx];
				//  row: go from right to left
				rowidx -= n;
				colidx--;
			}
			// colidx should point the 'unkown'
			// if colidx diagonal
			if (row==col) { // or rowidx = colidx, then X(i,i) = 1/L(i,i)
				m[colidx] = ((1.0/m[colidx])-sum) / m[rowidx];
			} else {
				m[colidx] = (m[colidx]-sum) / m[rowidx];
			}
		}
	}
	
	
 	btdebug_exit("btlincholeskysub");

	return 0;

}

// ********* Cholesky Blocked version **********

#ifdef BTLAPACK
// Cholesky decomposition, blocked version, calling btlapack subroutines
int btlin_bcholesky(double* m, int n, double* det) {
	int nb, block_offset; // block dimension
	int next_n, info, lda;
	char uplo = 'L';
	double block_det;
	double* pdiag, *pcol, *pmat;
	
	*det = 0.0;
	nb = 2; // IMPORTANT
	
	// may want to check for squareness
	lda = n;
	pdiag = m;	
	block_offset = nb*(lda+1);
		
	//printf("block offset %d\n",block_offset);	
	//printf("diag %lf\n",*pdiag);
	
	pcol = pdiag + nb;
	pmat = pdiag+block_offset;
	
	
	while (n > nb) {
		//btlapack_printDMATRIX(dm);
		//printf("remaining: %d  nb %d\n",n,nb);
		next_n = n-nb;
		//btlapack_printDMATRIX(dm);
		info = btlapack_ducholeskydet('L',pdiag,nb,lda,&block_det);
		if (info != 0) {
			printf("Error ducholeskydet %d\n",info);
			return info;
		}
		*det += block_det;
		
		//printf("det=%lf\n",block_det);
		// update column
		info = btlin_dbcholupdcol(pcol,next_n,pdiag,nb,lda);
		
		//printf("after column update\n");
		//btlapack_printDMATRIX(dm);
				
		// update mat
		//btlapack_printDMATRIX(dm);
		//printf("----> DOING UPDATE\n");
		info = btlin_dbcholupdmat(pmat,next_n,pcol,nb,lda);
		//btlapack_printDMATRIX(dm);
		// update values for next iteration
		n = next_n;
		pdiag = pmat;
		pcol += block_offset;
		pmat += block_offset;
	}
	//btlapack_printDMATRIX(dm);
	// du =  double-unblocked
	info = btlapack_ducholeskydet('L',pdiag,n,lda,&block_det);
	*det += block_det;
	//printf("det=%lf\n",block_det);
	return info;
	

}

// Calculates the inverse of a lower triangular matrix
// Blocked version. uses Fortran dtrti2 (unblocked version)
int btlin_ltri_invhybrid(double* a, int n, double* det) {

	// The structure of this code is very similar to Fortran dtrtri
	// The main difference is that the block inverse is taken before the second product
	// The idea is that they can be done in parallel
	
	int diag, old_diag, nb, info;
	int current_nb, nblockcolumn;
	double *pdiag,*old_pdiag, *pblockcolumn;
	

	nb = 2;
	if (n <= nb) {
		current_nb = n;
	} else {
		current_nb = n % nb;
		if (current_nb == 0)
			current_nb = nb;
	}
	diag = n - current_nb;
	pblockcolumn = pdiag = a + (n*diag +diag);
	nblockcolumn = current_nb;
	// Inverse of first block
	info = btlapack_ltri_inv2(pdiag,current_nb,n,NULL);
	if (info != 0) {
		printf("Error btlapack_ltri_inv2 %d\n",info);
		exit(0);
		return info;
	}
	old_diag = diag;
	diag -= nb; // always go up a full block
	
	// outer loop
	while(diag >= 0) {
		pblockcolumn = pdiag - n*nb;
		old_pdiag = pdiag;
		pdiag = pblockcolumn-nb;
		// block update 1: X  <- L1^(-1) * X
		// in-place - No extra workspace required (limited parallelism)
		
		//printf("--- After first inverse\n"); 
		//btlin_print(a,n,n);
		
		//printf("Before lowerbypmatrix %lf %lf\n",*old_pdiag,*pblockcolumn);
		
		btlin_ltri_LowerByMatrix(-1.0,old_pdiag,pblockcolumn,nblockcolumn,nb,n);
		
		//btlin_print(a,n,n);
		
		// Inverse of new block
		info = btlapack_ltri_inv2(pdiag,nb,n,NULL);
		if (info != 0) {
			printf("Error btlapack_ltri_inv2 %d\n",info);
			exit(0);
			return info;
		}
		
		//printf("After local inverse\n");
		//btlin_print(a,n,n);
		
		// update column block 2
		//printf("Before matrixbylower %lf %lf\n",*pblockcolumn,*pdiag);
		btlin_ltri_MatrixByLower(1.0,pblockcolumn,pdiag,nblockcolumn,nb,n);
		
		//btlin_print(a,n,n);
		
		// update variables
		nblockcolumn += nb;
		old_diag = diag;
		diag -= nb;
		
		//exit(0);
		
	}
	if (det != NULL) {
		*det = btlin_ltri_det(a,n,n);		
	}

	return 0;
}



#endif // BTLAPACK special case

// ******** Auxiliary procedures ***********

int btlin_dbcholupdcol(double* pcol,int nrows,double* pl,int nb,int lda) {
	int row, col, i;
	double sum, *local, *pl_start, *pcol_iter;
	double *pcol_start;
	
	pl_start = pl;

	local = (double*)malloc(sizeof(double)*nb);
	
	for(col=0; col < nb; col++) {  // col:[0,nb-1]
		//pcol_start = pcol+lda*col;
		// obtained after exiting the internal loop (ends in correct pointer)
		
		// load localL from pL 
		pl = pl_start;
		for(i=0; i <= col; i++) {  // col is column number
			local[i] = *pl;
			pl +=lda;
		}
		
		pl_start++; // next row
		sum = 0;
		pcol_start = pcol;
		for(row=0; row < nrows; row++) { //row: [0,nrows-1]
			pcol_iter = pcol_start;
			//printf("row %d\n",row);
			//compute sum
			sum = 0;
			for(i=0; i < col;i++) {
				sum+= local[i]*(*pcol_iter);
				pcol_iter += lda;
			}
			// compute unknown (row,col)
			
			*pcol_iter =(*pcol_iter - sum)/local[i];
			pcol_start++;
		}
	}

	free(local);
	
	return 0;
}


int btlin_dbcholupdmat(double* pmat,int dim_mat,double* pcol,int nb, int lda) {
	int row,col,i;
	double *local_idx, *local_start;
	double *pmat_start, *pmat_index, *pcol_idx, *pcol_start;
	double* local, sum;

	
	local = (double*)malloc(sizeof(double)*nb);
	local_idx = pcol; // top-left corner of L
	local_start = pcol;
	pmat_start = pmat;
	pcol_start = pcol;
	for(col=0; col < dim_mat; col++) {
		// load col into local array
		pcol_start = pcol+col;
		local_idx = local_start;
		for(i=0; i < nb; i++) {
			local[i] = *local_idx;
			local_idx += lda;
		}
		
		pmat_index = pmat_start;		
		for(row=col; row < dim_mat; row++) {
			sum = 0.0;
			// compute sum
			pcol_idx = pcol_start;
			//printf("sum :");
			for(i=0;i<nb;i++) {
				//printf(" %lf*%lf ",local[i],(*pcol_idx));
				sum += local[i]*(*pcol_idx);
				pcol_idx += lda;
			}
			//printf("--> old %lf sum %lf", *pmat_index, sum);
			*pmat_index -= sum;
			//printf(" new %lf\n",*pmat_index); 
			pmat_index++;
			pcol_start++;
		}
		local_start++;  
		pmat_start += lda+1; // follow diagonal
		pcol_start++;
	}
	
	free(local);

	return 0;
}


// computes the log of the determinant of a lower triangular matrix
double btlin_ltri_det(double* a, int n, int lda) {
	double  det, *pdiag;
	int i;

	det = 0.0;
	pdiag = a;
	//printf("Diag elements: ");
	for(i=0; i<n; i++) {
		//printf("%lf ",*pdiag);
		if (*pdiag < 0)
			det += log(-*pdiag);
		else
			det += log(*pdiag);
		pdiag += (lda+1);
	}
	//printf("-- total %lf\n",det);

	return det;
}

// **************************************************
// Auxiliary functions related to btlin_ltri_inv
// **************************************************

// A <- alpha*L*M
int btlin_ltri_LowerByMatrix(double alpha, double* pL,double* pA,int m,int n,int lda) {
	int rowA, colA, k;
	double *pL_start, *pL_worker;
	double *pA_start, *pA_worker;
	double *pA_result, pA_result_start;
	double temp;
	
	//printf("------------ m %d n %d\n",m,n);

	// DIM(a) = m x n
	pA_start = pA; // start from bottom left corner
	for(colA=0; colA < n; colA++) {
		pA_result = pA_start+m-1;
		pL_start = pL+m-1;
		for(rowA=m-1; rowA >= 0; rowA--) {
			temp = 0.0;
			pA_worker = pA_start;
			pL_worker = pL_start;
			for(k = 0; k <= rowA; k++) {
				//printf("product L %lf A %lf\n",*pL_worker,*pA_worker);
				temp += (*pL_worker) * (*pA_worker);
				pA_worker++;
				pL_worker += lda;
			}
			//printf("!!!updating %lf with %lf\n",*pA_result,alpha*temp);
			*pA_result = alpha * temp;
			pA_result--;
			pL_start--;
		}
		pA_start += lda;
	}
	
	return 0;
}


// A <- alpha*M*L
int btlin_ltri_MatrixByLower(double alpha, double* pA,double* pL,int m,int n,int lda) {
	int rowA, colA, k;
	double *pL_start, *pL_worker;
	double *pA_start, *pA_worker, *pA_result;
	double temp;

	//printf("m %d n %d\n",m,n);
	
	pA_start = pA+(n-1)*lda; // go to top right corner
	for(rowA=0; rowA < m; rowA++) {
		pA_result = pA;
		pL_start = pL+n-1;  // go to bottom left corner
		for(colA=0; colA < n; colA++) {
			temp = 0.0;
			pA_worker = pA_start;
			pL_worker = pL_start;
			for(k = colA; k < n; k++) {
				//printf("product L %lf A %lf\n",*pL_worker,*pA_worker);
				temp += (*pL_worker) * (*pA_worker);
				pA_worker -= lda;
				pL_worker--;
			}
			//printf("!!!updating %lf with %lf\n",*pA_result,alpha*temp);
			*pA_result = alpha*temp;
			pA_result += lda;
			pL_start += lda;
		}
		// row finished
		pA_start++; // go to next row
		pA++; // re-using PA as a marker for top left;
	}
	
	return 0;
}





