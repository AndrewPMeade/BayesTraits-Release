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



#ifdef BTOCL

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "btocl_dmatrix.h"
#include "btlapack.h"


int btocl_dbcholupdcol(cl_command_queue queue, cl_kernel kernel, cl_mem buffer, int col_idx, int nrows, int diag_idx, int cbs, int lda);

// ------------------------ BLOCKED version Lower Triangular-----------------------------

int btocl_choleskyDMATRIX(cl_mem buffer, DMATRIX* dm, double* det, int cbs, int mbs) {
	double* m;
	int n,lda;

	// may want to check for squareness
	m = dm->m;
	lda = n = dm->nrows;


	return btocl_cholesky(buffer,m,n,det,cbs,mbs);
}

int btocl_invcholeskyDMATRIX(cl_mem buffer, DMATRIX* dm, double* det, int cbs, int mbs) {
	double* m;
	int n,lda;

	// may want to check for squareness
	m = dm->m;
	lda = n = dm->nrows;

	return btocl_invcholesky(buffer,m,n,det,cbs,mbs);

}

// ---------------------------------------------------------------------------------------




/* *********** OLD code below ************ */
/* *********** OLD code below ************ */
/* *********** OLD code below ************ */
/* *********** OLD code below ************ */


/*  ********************* OLD code below ******************* */

void btocl_cholupdmat_kernel(int id, double* matrix, double* diagval, int diagidx, int lda, int n, int dimBlocks);


int btocl_invdetcholeskyDMATRIX1(cl_mem buffer, DMATRIX* dm, double* det) {
	// need extra buffer to store new diagonal
	int i,j, k,diagidx, subdiagidx, idx, idx1, idx2,n;
	double diag, *m;
	double* diagArray;
	//OCL variables
	cl_kernel kernel1, kernel2;
	cl_int err;
	cl_context context;
	cl_command_queue queue;
	cl_mem buffer_diag;
	size_t global_work_size;

	*det = 0;
	n = dm->nrows;
	m = dm->m;
	diagArray = (double*)malloc(sizeof(double)*n);

	context = btocl_getContext();
	queue =  btocl_getCommandQueue();
	kernel1 = btocl_getKernel(BTOCL_CHOLUPDCOL);

	if (kernel1 == NULL) {
		printf("Error: Couldn't load kernel %s\n",btocl_getKernelName(BTOCL_CHOLUPDCOL));
		exit(1);
	} else {
		printf("loaded %s\n",btocl_getKernelName(BTOCL_CHOLUPDCOL));
	}
	kernel2 = btocl_getKernel(BTOCL_CHOLUPDMAT);
	if (kernel1 == NULL) {
		printf("Error: Couldn't load kernel %s\n",btocl_getKernelName(BTOCL_CHOLUPDMAT));
		exit(1);
	} else {
		printf("loaded %s\n",btocl_getKernelName(BTOCL_CHOLUPDMAT));
	}

	// buffer to store diagonal
	buffer_diag = clCreateBuffer(context, CL_MEM_WRITE_ONLY, sizeof(double), NULL, &err);
	if (err < 0) {
		printf("Couldn't create diagonal buffer\n");
		exit(1);
	}

	//btdebug_enter("btoclcopy1");
	// copy matrix to buffer - may not need this - probably copied during creation
	clEnqueueWriteBuffer(queue,buffer,CL_TRUE,0,n*n*sizeof(double),&m[0],0,0,NULL);
	//btdebug_exit("btoclcopy1");

	diagidx = 0;
	diag = m[diagidx];	// first diagonal
	global_work_size = n; // size of column to be updated
	for(i=0; i < n; i++) {
		//printf("Column %d\n",i);
		//printf("new diag %lf\n",diag);
		// diag = m[diagidx]; -- computed in previous iteration
		if (diag < 0) {
			printf("Error index %d negative diagonal %lf\n",diagidx,diag);
			return 1;
		}
		global_work_size--;
		diag  = sqrt(diag);
		//printf("sq root diag %lf\n",diag);
		diagArray[i] = diag; // square root of diagonal
		if (diag < 0)
			*det += log(-diag);
		else
			*det += log(diag);

		// call update column kernel
		//perhaps it could be extended to have three phases:
		// update, copy to local and do multiplication

		//idx = diagidx;
		//for(j=i+1; j < n; j++) {
		//	idx++;
		//	m[idx] = m[idx]/diag;
		//}

		//printf("before update - diag %lf\n",diag);
		//btlapack_printDMATRIX(dm);

		// copy parameters
		// Pass Arguments
		if (global_work_size > 0) {
			if ((err = clSetKernelArg(kernel1,0,sizeof(cl_mem), &buffer)) < 0) {
				printf("Couldnt set first argument\n");
				exit(1);
			}
			if ((err = clSetKernelArg(kernel1,1,sizeof(diagidx), &diagidx)) < 0) {
				printf("Couldnt set second argument\n");
				exit(1);
			}
			if ((err = clSetKernelArg(kernel1,2,sizeof(diag), &diag)) < 0) {
				printf("Couldnt set third argument\n");
				exit(1);
			}
			// copy colum to buffer
			//btdebug_enter("btoclkernel1");
			// schedule kernel
			clEnqueueNDRangeKernel(queue,kernel1,1,NULL,&global_work_size,NULL,0,NULL,NULL);
			// copy back
			//btdebug_exit("btoclkernel1");

			if ((err = clSetKernelArg(kernel2,0,sizeof(cl_mem), &buffer)) < 0) {
				printf("Couldnt set first argument\n");
				exit(1);
			}
			if ((err = clSetKernelArg(kernel2,1,sizeof(cl_mem), &buffer_diag)) < 0) {
				printf("Couldnt set second argument\n");
				exit(1);
			}
			if ((err = clSetKernelArg(kernel2,2,sizeof(diagidx), &diagidx)) < 0) {
				printf("Couldnt set third argument\n");
				exit(1);
			}
			if ((err = clSetKernelArg(kernel2,3,sizeof(n), &n)) < 0) {  // lda
				printf("Couldnt set fifth argument\n");
				exit(1);
			}
			if ((err = clSetKernelArg(kernel2,4,sizeof(global_work_size), &global_work_size)) < 0) {
				printf("Couldnt set forth argument\n");
				exit(1);
			}
			//btdebug_enter("btoclkernel2");
			clEnqueueNDRangeKernel(queue,kernel2,1,NULL,&global_work_size,NULL,0,NULL,NULL);
			//btdebug_exit("btoclkernel2");
			// copy diagonal_buffer to diag
			//btdebug_enter("btocldiag");
			clEnqueueReadBuffer(queue,buffer_diag,CL_TRUE,0,sizeof(double),&diag,0,0,NULL);
			//btdebug_exit("btocldiag");
		}
		//printf("after update, global_work size %d\n",global_work_size);
		//btlapack_printDMATRIX(dm);

		diagidx += n+1; // update diagonal index


	}
	// copy matrix from device to host
	clEnqueueReadBuffer(queue,buffer,CL_TRUE,0,n*n*sizeof(double),&m[0],0,0,NULL);

	// copy diagonals
	diagidx=0;
	for(j = 0; j < n; j++) {
		m[diagidx]=diagArray[j];
		diagidx += (n+1);
	}

	*det *= 2.0;
	return 0;
}

// This version does not return the diagonal after each iteration
// it computes it after the matrix update.


int btocl_invdetcholeskyDMATRIX2(cl_mem buffer, DMATRIX* dm, double* det) {
	// need extra buffer to store new diagonal
	int i,j, k,diagidx, subdiagidx, idx, idx1, idx2,n;
	double diag, a, *m;
	//OCL variables
	cl_kernel kernel1, kernel2;
	cl_int err;
	cl_context context;
	cl_command_queue queue;

	size_t global_work_size;

	*det = 0;
	n = dm->nrows;
	m = dm->m;

	context = btocl_getContext();
	queue =  btocl_getCommandQueue();
	kernel1 = btocl_getKernel(BTOCL_CHOLUPDCOL);

	if (kernel1 == NULL) {
		printf("Error: Couldn't load kernel %s\n",btocl_getKernelName(BTOCL_CHOLUPDCOL));
		exit(1);
	} else {
		printf("loaded %s\n",btocl_getKernelName(BTOCL_CHOLUPDCOL));
	}
	kernel2 = btocl_getKernel(BTOCL_CHOLUPDMAT);
	if (kernel1 == NULL) {
		printf("Error: Couldn't load kernel %s\n",btocl_getKernelName(BTOCL_CHOLUPDMAT));
		exit(1);
	} else {
		printf("loaded %s\n",btocl_getKernelName(BTOCL_CHOLUPDMAT));
	}
	diagidx = 0;
	diag = m[diagidx];
	if (diag < 0) {
			printf("Error index %d negative diagonal %lf\n",diagidx,diag);
			return 1;
	} else {
		m[diagidx] = diag;
	}



	//btdebug_enter("btoclcopy1");
	// copy matrix to buffer - may not need this - probably copied during creation
	clEnqueueWriteBuffer(queue,buffer,CL_TRUE,0,n*n*sizeof(double),&m[0],0,0,NULL);
	//btdebug_exit("btoclcopy1");


	global_work_size = n; // size of column to be updated
	for(i=0; i < n; i++) {
		//printf("Column %d\n",i);
		//printf("new diag %lf\n",diag);
		// diag = m[diagidx]; -- computed in previous iteration

		global_work_size--;

		//printf("sq root diag %lf\n",diag);


		// call update column kernel
		//perhaps it could be extended to have three phases:
		// update, copy to local and do multiplication

		//idx = diagidx;
		//for(j=i+1; j < n; j++) {
		//	idx++;
		//	m[idx] = m[idx]/diag;
		//}

		//printf("before update - diag %lf\n",diag);
		//btlapack_printDMATRIX(dm);

		// copy parameters
		// Pass Arguments
		if (global_work_size > 0) {
			if ((err = clSetKernelArg(kernel1,0,sizeof(cl_mem), &buffer)) < 0) {
				printf("Couldnt set first argument\n");
				exit(1);
			}
			if ((err = clSetKernelArg(kernel1,1,sizeof(diagidx), &diagidx)) < 0) {
				printf("Couldnt set second argument\n");
				exit(1);
			}

			// copy colum to buffer
			//btdebug_enter("btoclkernel1");
			// schedule kernel
			clEnqueueNDRangeKernel(queue,kernel1,1,NULL,&global_work_size,NULL,0,NULL,NULL);
			// copy back
			//btdebug_exit("btoclkernel1");

			if ((err = clSetKernelArg(kernel2,0,sizeof(cl_mem), &buffer)) < 0) {
				printf("Couldnt set first argument\n");
				exit(1);
			}

			if ((err = clSetKernelArg(kernel2,1,sizeof(diagidx), &diagidx)) < 0) {
				printf("Couldnt set second argument\n");
				exit(1);
			}
			if ((err = clSetKernelArg(kernel2,2,sizeof(n), &n)) < 0) {  // lda
				printf("Couldnt set third argument\n");
				exit(1);
			}
			if ((err = clSetKernelArg(kernel2,3,sizeof(global_work_size), &global_work_size)) < 0) {
				printf("Couldnt set forth argument\n");
				exit(1);
			}
			//btdebug_enter("btoclkernel2");
			clEnqueueNDRangeKernel(queue,kernel2,1,NULL,&global_work_size,NULL,0,NULL,NULL);
			//btdebug_exit("btoclkernel2");
			// copy diagonal_buffer to diag
			//btdebug_enter("btocldiag");

			//btdebug_exit("btocldiag");
		}
		//printf("after update, global_work size %d\n",global_work_size);
		//btlapack_printDMATRIX(dm);

		diagidx += n+1; // update diagonal index


	}
	// copy matrix from device to host
	clEnqueueReadBuffer(queue,buffer,CL_TRUE,0,n*n*sizeof(double),&m[0],0,0,NULL);

	// copy diagonals
	diagidx=0; *det = 0.0;
	for(j = 0; j < n; j++) {
		diag = m[diagidx];
		if (diag < 0)
			*det += log(-diag);
		else
			*det += log(diag);
		diagidx += (n+1);
	}

	*det *= 2.0;
	return 0;
}

// This version divides the matrix into blocks NBxNB
// the kernel that updates the matrix does it in blocks NBxNB
int btocl_invdetcholeskyDMATRIXunblocked(cl_mem buffer, DMATRIX* dm, double* det) {
	// need extra buffer to store new diagonal
	int i,j, k,diagidx, subdiagidx, idx, idx1, idx2, dimBlocks;
	int n,nrem;
	double diag, a, *m;
	double* diagArray;
	const int bs = 32;
	//OCL variables
	cl_kernel kernel1, kernel2;
	cl_int err;
	cl_context context;
	cl_command_queue queue;
	cl_mem buffer_diag;
	size_t kernel1_work_size,kernel2_work_size;

	*det = 0;
	nrem = n = dm->nrows;
	m = dm->m;
	diagArray = malloc(sizeof(double)*n);

	context = btocl_getContext();
	queue =  btocl_getCommandQueue();
	kernel1 = btocl_getKernel(BTOCL_CHOLUPDCOL);

	if (kernel1 == NULL) {
		printf("Error: Couldn't load kernel %s\n",btocl_getKernelName(BTOCL_CHOLUPDCOL));
		exit(1);
	} else {
		printf("loaded %s\n",btocl_getKernelName(BTOCL_CHOLUPDCOL));
	}
	kernel2 = btocl_getKernel(BTOCL_CHOLUPDMAT);
	if (kernel1 == NULL) {
		printf("Error: Couldn't load kernel %s\n",btocl_getKernelName(BTOCL_CHOLUPDMAT));
		exit(1);
	} else {
		printf("loaded %s\n",btocl_getKernelName(BTOCL_CHOLUPDMAT));
	}

	// buffer to store diagonal
	buffer_diag = clCreateBuffer(context, CL_MEM_WRITE_ONLY, sizeof(double), NULL, &err);
	if (err < 0) {
		printf("Couldn't create diagonal buffer\n");
		exit(1);
	}

	//btdebug_enter("btoclcopy1");
	// copy matrix to buffer - may not need this - probably copied during creation
	clEnqueueWriteBuffer(queue,buffer,CL_TRUE,0,n*n*sizeof(double),&m[0],0,0,NULL);
	//btdebug_exit("btoclcopy1");

	diagidx = 0;
	diag = m[diagidx];	// first diagonal
	kernel1_work_size = n; // size of column to be updated
	for(i=0; i < n; i++) {
		//printf("Column %d\n",i);
		//printf("new diag %lf\n",diag);
		// diag = m[diagidx]; -- computed in previous iteration
		if (diag < 0) {
			printf("Error index %d negative diagonal %lf\n",diagidx,diag);
			return 1;
		}
		nrem--;
		dimBlocks = nrem/bs;
		if (nrem%bs != 0) dimBlocks++;
		kernel1_work_size = nrem;
		kernel2_work_size = dimBlocks*(dimBlocks+1)/2;
		//printf("Column %d blockwidth %d wgsize %d\n",i,dimBlocks,kernel2_work_size);
		diag  = sqrt(diag);
		//printf("sq root diag %lf\n",diag);
		diagArray[i] = diag; // square root of diagonal
		if (diag < 0)
			*det += log(-diag);
		else
			*det += log(diag);

		// call update column kernel
		//perhaps it could be extended to have three phases:
		// update, copy to local and do multiplication

		//idx = diagidx;
		//for(j=i+1; j < n; j++) {
		//	idx++;
		//	m[idx] = m[idx]/diag;
		//}

		//printf("before update - diag %lf\n",diag);
		//btlapack_printDMATRIX(dm);

		// copy parameters
		// Pass Arguments
		if (kernel1_work_size > 0) {
			if ((err = clSetKernelArg(kernel1,0,sizeof(cl_mem), &buffer)) < 0) {
				printf("Couldnt set first argument\n");
				exit(1);
			}
			if ((err = clSetKernelArg(kernel1,1,sizeof(diagidx), &diagidx)) < 0) {
				printf("Couldnt set second argument\n");
				exit(1);
			}
			if ((err = clSetKernelArg(kernel1,2,sizeof(diag), &diag)) < 0) {
				printf("Couldnt set third argument\n");
				exit(1);
			}
			// copy colum to buffer
			//btdebug_enter("btoclkernel1");
			// schedule kernel
			clEnqueueNDRangeKernel(queue,kernel1,1,NULL,&kernel1_work_size,NULL,0,NULL,NULL);
			// copy back
			//btdebug_exit("btoclkernel1");

			if ((err = clSetKernelArg(kernel2,0,sizeof(cl_mem), &buffer)) < 0) {
				printf("Couldnt set first argument\n");
				exit(1);
			}
			if ((err = clSetKernelArg(kernel2,1,sizeof(cl_mem), &buffer_diag)) < 0) {
				printf("Couldnt set second argument\n");
				exit(1);
			}
			if ((err = clSetKernelArg(kernel2,2,sizeof(diagidx), &diagidx)) < 0) {
				printf("Couldnt set third argument\n");
				exit(1);
			}
			if ((err = clSetKernelArg(kernel2,3,sizeof(n), &n)) < 0) {  // lda
				printf("Couldnt set forth argument\n");
				exit(1);
			}
			if ((err = clSetKernelArg(kernel2,4,sizeof(nrem), &nrem)) < 0) {
				printf("Couldnt set fifth argument\n");
				exit(1);
			}
			if ((err = clSetKernelArg(kernel2,5,sizeof(dimBlocks), &dimBlocks)) < 0) {  // n
				printf("Couldnt set sixth argument\n");
				exit(1);
			}
			//btdebug_enter("btoclkernel2");
			clEnqueueNDRangeKernel(queue,kernel2,1,NULL,&kernel2_work_size,NULL,0,NULL,NULL);
			//btdebug_exit("btoclkernel2");
			// copy diagonal_buffer to diag
			//btdebug_enter("btocldiag");
			clEnqueueReadBuffer(queue,buffer_diag,CL_TRUE,0,sizeof(double),&diag,0,0,NULL);
			//btdebug_exit("btocldiag");
		}
		//printf("after update, global_work size %d\n",global_work_size);
		//btlapack_printDMATRIX(dm);

		diagidx += n+1; // update diagonal index


	}
	// copy matrix from device to host
	clEnqueueReadBuffer(queue,buffer,CL_TRUE,0,n*n*sizeof(double),&m[0],0,0,NULL);

	// copy diagonals
	diagidx=0;
	for(j = 0; j < n; j++) {
		m[diagidx]=diagArray[j];
		diagidx += (n+1);
	}

	*det *= 2.0;
	return 0;
}


// ------------ Blocked version - use Upper Triangular 
// is not used, 

int btocl_invdetcholeskyDMATRIXUpper(cl_mem buffer, DMATRIX* dm, double* det) {
	int block_offset; // block dimension
	int n, next_n, info, lda;
	char uplo = 'U';
	double block_det;
	double* pdiag, *pcol, *pmat;
	double diag_idx, col_idx, mat_idx;
	// for buffer copying
	int copy_idx;
	double* pcopy;
	int i;
	double* m;
	size_t copy_offset;
	int copy_size;
	int iter;
	int dimBlocks;
	const int bs = 32;  // block size updmat
	const int nb = 32; // block size cholesky parition

	cl_kernel kernel1, kernel2;
	cl_int err;
	cl_context context;
	cl_command_queue queue;


	context = btocl_getContext();
	queue =  btocl_getCommandQueue();
	kernel1 = btocl_getKernel(BTOCL_CHOLUPDCOL);

	if (kernel1 == NULL) {
		printf("Error: Couldn't load kernel %s\n",btocl_getKernelName(BTOCL_CHOLUPDCOL));
		exit(1);
	} else {
		printf("loaded %s\n",btocl_getKernelName(BTOCL_CHOLUPDCOL));
	}
	kernel2 = btocl_getKernel(BTOCL_CHOLUPDMAT);
	if (kernel1 == NULL) {
		printf("Error: Couldn't load kernel %s\n",btocl_getKernelName(BTOCL_CHOLUPDMAT));
		exit(1);
	} else {
		printf("loaded %s\n",btocl_getKernelName(BTOCL_CHOLUPDMAT));
	}

	*det = 0.0;


	// may want to check for squareness
	m = dm->m;
	lda = n = dm->nrows;
	pdiag = m;
	diag_idx = 0;
	block_offset = nb*(lda+1);

	pcol    = pdiag + nb;
	col_idx = diag_idx + nb;
	pmat    = pdiag + block_offset;
	mat_idx = diag_idx + block_offset;

	copy_idx = diag_idx;
	pcopy = pdiag;

	// copy matrix to buffer
	clEnqueueWriteBuffer(queue,buffer,CL_TRUE,0,n*n*sizeof(double),m,0,0,NULL);

	//printf("START\n");
	//btlapack_printDMATRIX(dm);
	copy_size = nb*lda*sizeof(double);
	iter = 0;
	while (n > nb) {
		iter++;

		next_n = n-nb;

		//btlapack_printDMATRIX(dm);
		info = btlapack_ducholeskydet(uplo,pdiag,nb,lda,&block_det);
		if (info != 0) {
			printf("Error ducholeskydet %d\n",info);
			return info;
		}
		*det += block_det;

		//printf("After local cholesky\n");
		//btlapack_printDMATRIX(dm);

		copy_offset = copy_idx*sizeof(double);
		//printf("write buffer copy_idx %d copy_offset %d\n",copy_idx,copy_offset);
		err = clEnqueueWriteBuffer(queue,buffer,CL_TRUE,copy_offset,copy_size,&m[copy_idx],0,0,NULL);
		if(err != 0) {
			printf("Error writing buffer\n");
			exit(0);
		}

		//clEnqueueWriteBuffer(queue,buffer,CL_TRUE,0,lda*lda*sizeof(double),m,0,0,NULL);

		// -- update column
		// info = btlin_dbcholupdcol(pcol,next_n,pdiag,nb,lda);
		info = btocl_dbcholupdcol(queue,kernel1,buffer,col_idx,next_n,diag_idx,nb,lda);

		// debug - copy eveything to have a look
		clEnqueueReadBuffer(queue,buffer,CL_TRUE,0,lda*lda*sizeof(double),m,0,0,NULL);
		//btlapack_printDMATRIX(dm);

		// -- update mat
		dimBlocks = next_n/bs;
		if (next_n%bs != 0) dimBlocks++;

		// Calling this function "btocl_invdetcholeskyDMATRIXUpper" is not called
		// the funciton below has an incorrect number of paramters, its been commneted out and a exit
//		info = btocl_dbcholupdmat(queue,kernel2,buffer,mat_idx,next_n,col_idx,nb,lda,dimBlocks);

		exit(1);

		//btlapack_printDMATRIX(dm);

		// read diag+column from GPU to buffer

		// update values for next iteration
		n = next_n;
		pdiag = pmat;
		diag_idx = mat_idx;
		pcol += block_offset;
		col_idx += block_offset;
		pmat += block_offset;
		mat_idx += block_offset;

		// buffer copying
		pcopy += nb*lda;
		copy_idx += nb*lda;
		if (n < nb)
			copy_size = n*lda*sizeof(double);
		copy_offset = copy_idx*sizeof(double);

		//printf("writing from buffer to matrix. copy_idx %d size %d\n",copy_idx,copy_size/sizeof(double));
		clEnqueueReadBuffer(queue,buffer,CL_TRUE,copy_offset,copy_size,&m[copy_idx],0,0,NULL);
		//clEnqueueReadBuffer(queue,buffer,CL_TRUE,0,lda*lda*sizeof(double),m,0,0,NULL);


	}
	//btlapack_printDMATRIX(dm);

	// du =  double-unblocked
	info = btlapack_ducholeskydet(uplo,pdiag,n,lda,&block_det);
	*det += block_det;
	//printf("det=%lf\n",block_det);
	return info;

}


// ------- Blocked version auxiliary functions -----------------------------



// ---------- end Blocked version auxoliary functions --------------------------

// PRACTICE version
// call non-GPU "kernel" - unblocked
int btocl_invdetcholeskyDMATRIXdebug(cl_mem buffer, DMATRIX* dm, double* det) {
	// need extra buffer to store new diagonal
	int i,j, k,diagidx, subdiagidx, idx, idx1, idx2, dimBlocks;
	int n,nrem;
	double diag, a, *m;
	double* diagArray;
	//OCL variables
	cl_kernel kernel1, kernel2;
	cl_int err;
	cl_context context;
	cl_command_queue queue;
	cl_mem buffer_diag;
	size_t kernel1_work_size,kernel2_work_size;

	*det = 0;
	nrem = n = dm->nrows;
	m = dm->m;
	diagArray = malloc(sizeof(double)*n);

	context = btocl_getContext();
	queue =  btocl_getCommandQueue();
	kernel1 = btocl_getKernel(BTOCL_CHOLUPDCOL);

	if (kernel1 == NULL) {
		printf("Error: Couldn't load kernel %s\n",btocl_getKernelName(BTOCL_CHOLUPDCOL));
		exit(1);
	} else {
		printf("loaded %s\n",btocl_getKernelName(BTOCL_CHOLUPDCOL));
	}
	kernel2 = btocl_getKernel(BTOCL_CHOLUPDMAT);
	if (kernel1 == NULL) {
		printf("Error: Couldn't load kernel %s\n",btocl_getKernelName(BTOCL_CHOLUPDMAT));
		exit(1);
	} else {
		printf("loaded %s\n",btocl_getKernelName(BTOCL_CHOLUPDMAT));
	}

	// buffer to store diagonal
	buffer_diag = clCreateBuffer(context, CL_MEM_WRITE_ONLY, sizeof(double), NULL, &err);
	if (err < 0) {
		printf("Couldn't create diagonal buffer\n");
		exit(1);
	}

	//btdebug_enter("btoclcopy1");
	// copy matrix to buffer - may not need this - probably copied during creation
	clEnqueueWriteBuffer(queue,buffer,CL_TRUE,0,n*n*sizeof(double),&m[0],0,0,NULL);
	//btdebug_exit("btoclcopy1");

	diagidx = 0;
	diag = m[diagidx];	// first diagonal
	kernel1_work_size = n; // size of column to be updated
	for(i=0; i < n; i++) {
		//printf("Column %d\n",i);
		//printf("new diag %lf\n",diag);
		// diag = m[diagidx]; -- computed in previous iteration
		if (diag < 0) {
			printf("Error index %d negative diagonal %lf\n",diagidx,diag);
			return 1;
		}
		nrem--;
		dimBlocks = nrem/16;
		if (nrem%16 != 0) dimBlocks++;
		kernel1_work_size = nrem;
		kernel2_work_size = dimBlocks*(dimBlocks+1)/2;
		//printf("Column %d blockwidth %d wgsize %d\n",i,dimBlocks,kernel2_work_size);
		diag  = sqrt(diag);
		//printf("sq root diag %lf\n",diag);
		diagArray[i] = diag; // square root of diagonal
		if (diag < 0)
			*det += log(-diag);
		else
			*det += log(diag);

		// call update column kernel
		//perhaps it could be extended to have three phases:
		// update, copy to local and do multiplication

		idx = diagidx;
		for(j=i+1; j < n; j++) {
			idx++;
			m[idx] = m[idx]/diag;
		}

		//printf("before update - diag %lf\n",diag);
		//btlapack_printDMATRIX(dm);

		// copy parameters
		// Pass Arguments
		if (kernel1_work_size > 0) {




			for(j=0; j < kernel2_work_size;j++)
				btocl_cholupdmat_kernel(j,m,diagArray,diagidx,n,nrem,dimBlocks);

			//btdebug_exit("btoclkernel2");
			// copy diagonal_buffer to diag
			//btdebug_enter("btocldiag");
			clEnqueueReadBuffer(queue,buffer_diag,CL_TRUE,0,sizeof(double),&diag,0,0,NULL);
			//btdebug_exit("btocldiag");
		}
		//printf("after update, global_work size %d\n",global_work_size);
		//btlapack_printDMATRIX(dm);

		diagidx += n+1; // update diagonal index
		diag = m[diagidx];

	}
	// copy matrix from device to host
	clEnqueueReadBuffer(queue,buffer,CL_TRUE,0,n*n*sizeof(double),&m[0],0,0,NULL);

	// copy diagonals
	diagidx=0;
	for(j = 0; j < n; j++) {
		m[diagidx]=diagArray[j];
		diagidx += (n+1);
	}

	*det *= 2.0;
	return 0;
}


//int btocl_choleskyDMATRIX(DMATRIX* m, double* det) {
//
//	return 0;
//}

// C version of kernel, adapted for testing e.g. global_id is passed as argument
void btocl_cholupdmat_kernel(int id, double* matrix, double* diagval, int diagidx, int lda, int n, int dimBlocks) {
	double a;
	int i,j, idx,idx2, idxstart, idx2start, idxa;
	int row, col, idtmp;

	int lastrow, lastcol,last;

	// number of workers = dimBlocs(dimBlocks+1)/2

	// each kernel takes care of a 16x16 square of the matrix
	// Numbering(id): top-botton/left-right top left corner=0;
	// get row,col
	col=0; row=0; idtmp = id;
	while(idtmp >= dimBlocks) {
		idtmp -= dimBlocks;
		dimBlocks--;
		col++;
	}
	row = col+idtmp;
	//printf("kernel %d n %drow %d col %d - ",id,n,row,col);
	col = 16*col;
	row = 16*row;
	last = n-1;
	lastrow = (row+16-1);
	if (lastrow > last) lastrow = last;
	lastcol = (col+16-1);
	if (lastcol > last) lastcol = last;

	idx2start = diagidx+row+1;
	idxa = diagidx+col+1;
	idxstart = diagidx+lda*(1+col)+row+1;  // recal, row >= col
	// if diagonal block
	//printf(" r %d c %d lastr %d lastc %d idx2s %d idxs %d idxa %d\n",row,col,lastrow,lastcol,idx2start, idxstart,idxa);

	for(i=col; i <= lastcol;i++){
		idx = idxstart;
		idx2 = idx2start;
		a = matrix[idxa];
		//printf("   new col idx %d idx2 %d idxa %d a %lf\n",idx,idx2,idxa,a);
		for(j=row;j<=lastrow;j++) {
			//printf("idx %d old value %lf idx2 %d v2 %lf- ",idx,matrix[idx],idx2,matrix[idx2]);
			matrix[idx] -= a*matrix[idx2];
			//printf("new value  v %lf\n",matrix[idx]);
			idx++; idx2++;
		}
		idxstart += lda;
		idxa++;
	}

	if (id==0) {
		*diagval = matrix[diagidx+lda+1];
	}


}


#endif // if BTOCL defined


