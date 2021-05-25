

#ifdef BTOCL

#include "btocl_lin.h"

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

// Local
// ****** Auxiliary Cholesky ******
int btocl_dbcholupdcol(cl_command_queue q, cl_kernel k, cl_mem m,int pcol,int nrows, int pdiag, int nb, int lda);
int btocl_dbcholupdmat(cl_command_queue q, cl_ushort kernel_type, cl_kernel k, cl_mem m,int pmat, int dim_mat, int pcol,int nb, int lda, int mbs);
// *** Auxiliary Lower Triangular Inverse
int btocl_ltri_LowerByMatrix(cl_ushort kernel_type,cl_kernel kernel,cl_mem buffer,double alpha,int L_idx,int M_idx,int MT_idx, int m, int n,int lda);
int btocl_ltri_MatrixByLower(cl_ushort kernel_type,cl_kernel kernel,cl_mem buffer,double alpha,int MT_idx,int M_idx,int L_idx, int m, int n, int lda);
int btocl_ltri_LTbyL(cl_mem buffer,int n);


int btocl_invcholesky(cl_mem buffer, double* m, int n,  double* det,int cbs, int mbs) {
	char uplo='L';
	int info;
	cl_command_queue queue;	

	queue =  btocl_getCommandQueue();

	btocl_cholesky(buffer,m,n,det,cbs,mbs);
	// buffer and m are identical due to memory swaps in algorithm
	
	dpotri_(&uplo, &n, m, &n, &info);
	
	// copy lower triangle 
	//btlin_makeSymmetric('L',m,n);
	// at this point buffer is not symmetric
	
	// copy lower triangle 
	//printf("transposing....\n");
	// this is only a half/one-sided transposition: must change name
	clEnqueueWriteBuffer(queue,buffer,CL_TRUE,0,n*n*sizeof(double),m,0,0,NULL);
	btocl_writeLtoU(buffer,0,n,n,8);	// mbs
	clEnqueueReadBuffer(queue,buffer,CL_TRUE,0,n*n*sizeof(double),m,0,0,NULL);

	return 0;

}

// This version uses opencl kernels only
int btocl_invcholesky_pure(cl_mem buffer, double* m, int n,  double* det,int cbs, int mbs) {
	char uplo='L';
	int info;
	cl_command_queue queue;	
	double det2;


	queue =  btocl_getCommandQueue();

	info = btocl_cholesky(buffer,m,n,det,cbs,mbs);
	// buffer and m are identical due to memory swaps in algorithm
	if (info != 0) {
		printf("Error Cholesky Factorization\n");
		return info;
	}

	
	//dpotri_(&uplo, &n, m, &n, &info);
	
	// invert L
	btocl_ltri_inv(buffer, m, n, &det2, 128);
	printf("End ltri_inv\n");
	// L**-T * L^(-1)
	btocl_ltri_LTbyL(buffer,n);
	printf("End ltri_lybyl\n"); 

	clEnqueueReadBuffer(queue,buffer,CL_TRUE,0,n*n*sizeof(double),m,0,0,NULL);

	// Buffer and m should both have the same values: a symmetric matrix = m <- m^(-1)
	return 0;

}

// ------------------------ BLOCKED version Lower Triangular-----------------------------

//int btocl_choleskyDMATRIX_blocked(cl_mem buffer, BTLAPACK_DMATRIX* dm, double* det, int cbs, int mbs) {

// Cholesky Factorisation of matrix mat (mat_dim x mat_dim)
int btocl_cholesky(cl_mem buffer, double* mat, int mat_dim, double* det,int cbs, int mbs) {
	int block_offset; // block dimension
	int next_n, info;
	char uplo = 'L';
	double block_det;
	double* pdiag, *pcol, *pmat;
	int diag_idx, col_idx, mat_idx;
	// for buffer copying
	int copy_idx;
	double* pcopy;
	size_t copy_offset;
	int copy_size;
	int iter;
	//int dimBlocks;
	int n;
	
	int lda;
	double* m;
	
	//mbs = 8;  // block size updmat
	//cbs = 8; // block size cholesky parition
	
	cl_kernel kernel1, kernel2;
	cl_int err;
	cl_context context;
	cl_command_queue queue;
	cl_ushort kernel_type1, kernel_type2;

	context = btocl_getContext();
	queue =  btocl_getCommandQueue();	
		
	kernel_type1 = BTOCL_CHOLUPDCOL_PURE;
	kernel1 = btocl_getKernel(kernel_type1);
		
	if (kernel1 == NULL) {
		printf("Error: Couldn't load kernel %s\n",btocl_getKernelName(kernel_type1));
		exit(1);
	} else {
		printf("loaded %s\n",btocl_getKernelName(kernel_type1));
	}
	
	// The matrix update kernel needs to know the block size of the Cholesky part
	// TODO: Use one kernel with compilation directive if possible
	switch (cbs) {
		case 4: kernel_type2 = BTOCL_CHOLUPDMAT_P4; break;
		case 8: kernel_type2 = BTOCL_CHOLUPDMAT_P8; break;
		case 16: kernel_type2 = BTOCL_CHOLUPDMAT_P16; break;
		case 32: kernel_type2 = BTOCL_CHOLUPDMAT_P32; break;
		case 64: kernel_type2 = BTOCL_CHOLUPDMAT_P64; break;
		case 128: kernel_type2 = BTOCL_CHOLUPDMAT_P128; break;
		default: kernel_type2 = BTOCL_CHOLUPDMAT_P8; // default to 8
	}
	// New kernel that does not take advantage of symmetry
	// Used instead of others because (row,col) calculations were causing kernel timeout problems
	kernel_type2 = BTOCL_CHOLUPDMAT_PURE;
	
	printf("cbs %d  using kernel %s\n",cbs,btocl_getKernelName(kernel_type2));
	kernel2 = btocl_getKernel(kernel_type2);
	if (kernel1 == NULL) {
		printf("Error: Couldn't load kernel %s\n",btocl_getKernelName(kernel_type2));
		exit(1);
	} else {
		printf("loaded %s\n",btocl_getKernelName(kernel_type2));
	}
	
	*det = 0.0;
	m = mat;
	lda = n = mat_dim;
	pdiag = mat;	
	diag_idx = 0;
	block_offset = cbs*(lda+1);
		
	pcol    = pdiag + cbs;
	col_idx = diag_idx + cbs;
	pmat    = pdiag + block_offset;
	mat_idx = diag_idx + block_offset;
	
	copy_idx = diag_idx;
	pcopy = pdiag;
	
	// copy matrix to buffer
	clEnqueueWriteBuffer(queue,buffer,CL_TRUE,0,n*n*sizeof(double),m,0,0,NULL);
	//printf("Size of buffer: %d\n",n*n*sizeof(double));
	
	//printf("START\n");
	//btlapack_printDMATRIX(dm);
	copy_size = cbs*lda*sizeof(double);
	iter = 0;
	while (n > cbs) {
		//printf("Iteration %d n %d\n",iter,n);
		iter++;				
		next_n = n-cbs;		
		
		//btlapack_printDMATRIX(dm);
		info = btlapack_ducholeskydet(uplo,pdiag,cbs,lda,&block_det);
		if (info != 0) {
			printf("Error ducholeskydet %d\n",info);
			return info;
		}
		*det += block_det;
		//printf("det %d\n",*det);
		
		copy_offset = copy_idx*sizeof(double);
		//printf("write buffer copy_idx %d copy_offset %d copy size %d\n",copy_idx,copy_offset,copy_size);
		//printf("Valid addres?? %lf last %lf \n",m[copy_idx],m[copy_idx+cbs*lda-1]);
		err=clFinish(queue);
		if (err != CL_SUCCESS) {
			printf("error clFinish\n");
			exit(0);
		}
		err = clEnqueueWriteBuffer(queue,buffer,CL_TRUE,copy_offset,copy_size,&m[copy_idx],0,0,NULL);
		if(err != CL_SUCCESS) {
			printf("BTOCL Cholesky: Error writing buffer\n");
			btocl_printRuntimeError(err);
			exit(0);
		}
						
		//clEnqueueWriteBuffer(queue,buffer,CL_TRUE,0,lda*lda*sizeof(double),m,0,0,NULL);
		
		// -- update column
		// info = btlin_dbcholupdcol(pcol,next_n,pdiag,nb,lda);
		//printf("UpdChol-----");
		info = btocl_dbcholupdcol(queue,kernel1,buffer,col_idx,next_n,diag_idx,cbs,lda);
		if(err != CL_SUCCESS) {
			printf("BTOCL Cholesky: Error updating block column\n");
			btocl_printRuntimeError(err);
			exit(0);
		}
		//printf("finished\n");
		
		// debug - copy everything to have a look
		// do we need this?
		//err = clEnqueueReadBuffer(queue,buffer,CL_TRUE,0,lda*lda*sizeof(double),m,0,0,NULL);
		//if (err != CL_SUCCESS) {
		//	printf("Error Read buffer ALL for debug\n");
		//	exit(0);
		//}
		//btlapack_printDMATRIX(dm);
		
		// -- update mat		
		//info = btlin_dbcholupdmat(pmat,next_n,pcol,nb,lda);
		//printf("Params mat_idx %d next_n %d col_idx %d cbs %d\n",mat_idx,next_n,col_idx,cbs);
		//printf("UpdMat -----");
		info = btocl_dbcholupdmat(queue,kernel_type2,kernel2,buffer,mat_idx,next_n,col_idx,cbs,lda,mbs);
		if(err != CL_SUCCESS) {
			printf("BTOCL Cholesky: Error updating remaining matrix\n");
			btocl_printRuntimeError(err);
			exit(0);
		}
		
		//printf("finished\n");
		
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
		pcopy += cbs*lda;
		copy_idx += cbs*lda;
		if (n < cbs)
			copy_size = n*lda*sizeof(double);
		copy_offset = copy_idx*sizeof(double);
		
		//printf("writing from buffer to matrix. copy_idx %d size %d\n",copy_idx,copy_size/sizeof(double));
		err = clEnqueueReadBuffer(queue,buffer,CL_TRUE,copy_offset,copy_size,&m[copy_idx],0,0,NULL);
		//clEnqueueReadBuffer(queue,buffer,CL_TRUE,0,lda*lda*sizeof(double),m,0,0,NULL);
		if (err != CL_SUCCESS) {
			printf("BTOCL Cholesky: Error reading buffer\n");
			btocl_printRuntimeError(err);
			btocl_printBufferInfo(buffer);
			exit(0);
		}
			
	}
	//btlapack_printDMATRIX(dm);
	
	// du =  double-unblocked
	info = btlapack_ducholeskydet(uplo,pdiag,n,lda,&block_det);
	*det += block_det;
	copy_offset = copy_idx*sizeof(double);
	//printf("write buffer copy_idx %d copy_offset %d\n",copy_idx,copy_offset);
	err = clEnqueueWriteBuffer(queue,buffer,CL_TRUE,copy_offset,copy_size,&m[copy_idx],0,0,NULL);
	if(err != 0) {
		printf("BTOCL Cholesky: Error writing buffer\n");
		btocl_printRuntimeError(err);
		exit(0);
	}
	
	// write everything back to array
	clEnqueueReadBuffer(queue,buffer,CL_TRUE,0,lda*lda*sizeof(double),m,0,0,NULL);
	
	return info;
	
}

// ---------------------------------------------------------------------------------------




// ------------------------- Auxiliary Cholesky ----------------------------


int btocl_dbcholupdcol(cl_command_queue queue, cl_kernel kernel, cl_mem buffer, int col_idx, int nrows, int diag_idx, int cbs, int lda) {
	size_t kernel_work_size;
	cl_int err;
	if ((err = clSetKernelArg(kernel,0,sizeof(cl_mem), &buffer)) < 0) {
			printf("Couldnt set first argument\n");
			exit(1);
	}
	if ((err = clSetKernelArg(kernel,1,sizeof(col_idx), &col_idx)) < 0) {
		printf("Couldnt set second argument\n");
	exit(1);
	}
	if ((err = clSetKernelArg(kernel,2,sizeof(nrows), &nrows)) < 0) {
		printf("Couldnt set third argument\n");
		exit(1);
	}
	if ((err = clSetKernelArg(kernel,3,sizeof(diag_idx), &diag_idx)) < 0) {
		printf("Couldnt set fourth argument\n");
		exit(1);
	}
	if ((err = clSetKernelArg(kernel,4,sizeof(cbs), &cbs)) < 0) {
		printf("Couldnt set fifth argument\n");
		exit(1);
	}
	if ((err = clSetKernelArg(kernel,5,sizeof(lda), &lda)) < 0) {
		printf("Couldnt set sixth argument\n");
		exit(1);
	}
	
	kernel_work_size = nrows;
	//printf("Launching kernel with %d workintems\n",kernel_work_size);
	//printf("args col_idx %d diag_idx %d nb %d lda %d\n",col_idx,diag_idx,nb,lda);
	err = clEnqueueNDRangeKernel(queue,kernel,1,NULL,&kernel_work_size,NULL,0,NULL,NULL);
	if (err != CL_SUCCESS) {
		printf("Error cholesky update column\n");
		btocl_printRuntimeError(err);
	}
	err = clFinish(queue);
	if (err != CL_SUCCESS) {
		printf("Error cholesky update column - finish\n");
		btocl_printRuntimeError(err);
	}
	return 0;
}


// mbs: block size to be used to divide "mat"
int btocl_dbcholupdmat(cl_command_queue queue, cl_ushort kernel_type, cl_kernel kernel, cl_mem buffer, int mat_idx, int mat_dim, int col_idx, int nb, int lda, int mbs) {
	size_t kernel_work_size;
	cl_int err;
	int dimBlocks; // number of blocks per side
	
	if (kernel_type == BTOCL_CHOLUPDMAT_PURE) {
		//kernel_work_size = mat_dim*(mat_dim+1)/2;
		kernel_work_size = mat_dim*mat_dim;
	} else {
		dimBlocks = mat_dim/mbs;
		if (mat_dim%mbs != 0) dimBlocks++;
		kernel_work_size = dimBlocks*(dimBlocks+1)/2;
	}
	
	if ((err = clSetKernelArg(kernel,0,sizeof(cl_mem), &buffer)) < 0) {
			printf("Couldnt set first argument\n");
			exit(1);
	}
	if ((err = clSetKernelArg(kernel,1,sizeof(mat_idx), &mat_idx)) < 0) {
		printf("Couldnt set second argument\n");
	exit(1);
	}
	if ((err = clSetKernelArg(kernel,2,sizeof(mat_dim), &mat_dim)) < 0) {
		printf("Couldnt set third argument\n");
		exit(1);
	}
	if ((err = clSetKernelArg(kernel,3,sizeof(col_idx), &col_idx)) < 0) {
		printf("Couldnt set fourth argument\n");
		exit(1);
	}
	if ((err = clSetKernelArg(kernel,4,sizeof(nb), &nb)) < 0) {
		printf("Couldnt set fifth argument\n");
		exit(1);
	}
	if ((err = clSetKernelArg(kernel,5,sizeof(lda), &lda)) < 0) {
		printf("Couldnt set sixth argument\n");
		exit(1);
	}
	if (kernel_type != BTOCL_CHOLUPDMAT_PURE) {
		if ((err = clSetKernelArg(kernel,6,sizeof(mbs), &mbs)) < 0) {
			printf("Couldnt set seventh argument\n");
			exit(1);
		}
		if ((err = clSetKernelArg(kernel,7,sizeof(dimBlocks), &dimBlocks)) < 0) {
			printf("Couldnt set eight argument\n");
			exit(1);
		}
	} 
	err = clEnqueueNDRangeKernel(queue,kernel,1,NULL,&kernel_work_size,NULL,0,NULL,NULL);
	//__kernel void btocl_cholupdmat_kernel(__global double* matrix, int mat_idx, int dim_mat, int col_idx, int nb, int lda, int dimBlocks)
	if (err != CL_SUCCESS) {
		printf("Error cholesky update Mat\n"); 
		btocl_printRuntimeError(err);
		exit(0);
	}
	err = clFinish(queue);
	if (err != CL_SUCCESS) {
		printf("Error cholesky update Mat - clFinish\n");
		btocl_printRuntimeError(err);
		exit(0);
	}
	
	return 0;
}


// ------------------------- Transpose Buffer ------------------------------
int btocl_writeLtoU(cl_mem buffer,int mat_idx, int mat_dim, int lda, int mbs) {
	size_t kernel_work_size;
	cl_int err;
	int dimBlocks; // number of blocks per side
	cl_command_queue queue;
	cl_kernel kernel;
	cl_context context;
	cl_ushort kernel_type;
	
	
	context = btocl_getContext();
	queue =  btocl_getCommandQueue();	
	
	
	kernel_type = BTOCL_TRANSPOSE;
	kernel = btocl_getKernel(kernel_type);
		
	if (kernel == NULL) {
		printf("Error: Couldn't load kernel %s\n",btocl_getKernelName(kernel_type));
		exit(1);
	} else {
		printf("loaded %s\n",btocl_getKernelName(kernel_type));
	}
	
	dimBlocks = mat_dim/mbs;
	if (mat_dim%mbs != 0) dimBlocks++;
	kernel_work_size = dimBlocks*(dimBlocks+1)/2;
	
	if ((err = clSetKernelArg(kernel,0,sizeof(cl_mem), &buffer)) < 0) {
			printf("Couldnt set first argument\n");
			exit(1);
	}
	if ((err = clSetKernelArg(kernel,1,sizeof(mat_idx), &mat_idx)) < 0) {
		printf("Couldnt set second argument\n");
	exit(1);
	}
	if ((err = clSetKernelArg(kernel,2,sizeof(mat_dim), &mat_dim)) < 0) {
		printf("Couldnt set third argument\n");
		exit(1);
	}
	
	if ((err = clSetKernelArg(kernel,3,sizeof(lda), &lda)) < 0) {
		printf("Couldnt set fourth argument\n");
		exit(1);
	}
	if ((err = clSetKernelArg(kernel,4,sizeof(mbs), &mbs)) < 0) {
		printf("Couldnt set fifth argument\n");
		exit(1);
	}
	if ((err = clSetKernelArg(kernel,5,sizeof(dimBlocks), &dimBlocks)) < 0) {
		printf("Couldnt set sixth argument\n");
		exit(1);
	}
	
	clEnqueueNDRangeKernel(queue,kernel,1,NULL,&kernel_work_size,NULL,0,NULL,NULL);
	//__kernel void btocl_cholupdmat_kernel(__global double* matrix, int mat_idx, int dim_mat, int col_idx, int nb, int lda, int dimBlocks)
	
	return 0;

}

// Kronecker Product with vector multiplication
// vres = v * [ sigma *K* mat ]
// dim(vres) = dim(v) = dim(sigma)*dim(mat)
// dim(sigma) small
int btocl_kronecker_vectmult(cl_mem vres_buffer, cl_mem v_buffer, cl_mem sigma_buffer, int sigma_dim, cl_mem mat_buffer, int mat_dim) {
	cl_command_queue queue;
	cl_kernel kernel;
	cl_context context;
	cl_ushort kernel_type;
	int argnum, err;
	size_t globalws[3];
	size_t localws[3];  // number of wokitems per dimension
	int dim, usewg, local_memSize;  
	int LS=1;
	context = btocl_getContext();
	queue =  btocl_getCommandQueue();
	
	// set kernel type
	//kernel_type = BTOCL_KRON_BASIC;
	//kernel_type = BTOCL_KRON_ACCUM;
	if (sigma_dim == 1)
		kernel_type = BTOCL_KRON_LOG1;
	else
		kernel_type = BTOCL_KRON_ACCUM;
	
	kernel = btocl_getKernel(kernel_type);		
	if (kernel == NULL) {
		printf("Error: Couldn't load kernel %s\n",btocl_getKernelName(kernel_type));
		exit(1);
	} else {
		printf("loaded %s\n",btocl_getKernelName(kernel_type));
	}

	// set arguments
	argnum = 0;
	if ((err = clSetKernelArg(kernel,argnum++,sizeof(cl_mem), &vres_buffer)) < 0)
		return 0;
	if ((err = clSetKernelArg(kernel,argnum++,sizeof(cl_mem), &v_buffer)) < 0)
		return 0;
	if ((err = clSetKernelArg(kernel,argnum++,sizeof(cl_mem), &sigma_buffer)) < 0)
		return 0;
	if ((err = clSetKernelArg(kernel,argnum++,sizeof(int), &sigma_dim)) < 0)
		return 0;
	if ((err = clSetKernelArg(kernel,argnum++,sizeof(cl_mem), &mat_buffer)) < 0)
		return 0;
	if ((err = clSetKernelArg(kernel,argnum++,sizeof(int), &mat_dim)) < 0)
		return 0;
	if (kernel_type == BTOCL_KRON_ACCUM) {
		LS = 128;   // 2^k
		local_memSize = sizeof(double)*LS;
		if ((err = clSetKernelArg(kernel,argnum++,local_memSize, NULL)) < 0) // accumulator
			return 0;
	} else if (kernel_type == BTOCL_KRON_LOG1) {
		local_memSize = sizeof(double)*sigma_dim;
		if ((err = clSetKernelArg(kernel,argnum++,local_memSize, NULL)) < 0) // log accumulator
			return 0;
		if ((err = clSetKernelArg(kernel,argnum++,sizeof(int), &LS)) < 0)
		return 0;
	}
		
	dim = 1;
	usewg = 0;
	globalws[0] = mat_dim*sigma_dim;
	if (kernel_type == BTOCL_KRON_ACCUM) {
		usewg = 1;
		localws[0] = sigma_dim;
	} else if (kernel_type == BTOCL_KRON_LOG1) {
		//printf("LOG version!");
		globalws[0] = mat_dim*LS;
		usewg = 1;
		localws[0] = LS;
	} 
	
	
	// set variable arguments
	// none yet
	
	// call kernel
	if (usewg) {  // yes
	  //printf("mat %d sigma %d global %lu local %lu\n",mat_dim,sigma_dim,globalws[0],localws[0]);
	  clEnqueueNDRangeKernel(queue,kernel,dim,NULL,globalws,localws,0,NULL,NULL);
	} else {
	  clEnqueueNDRangeKernel(queue,kernel,dim,NULL,globalws,NULL,0,NULL,NULL);
	}

	// don't need to copy result buffer - stored in vres_buffer

	return 1;
} 

int  btocl_ltri_inv(cl_mem buffer, double* pMat, int mat_dim, double* det,int nb) {
	// The structure of this code is very similar to Fortran dtrtri
	// The main difference is that the block inverse is taken before the second product
	// The idea is that they can be done in parallel
	
	int diag, old_diag, info, i, err;
	int current_nb, nblockcolumn;
	double *pdiag,*old_pdiag, *pblockcolumn;
	int diag_idx, old_diag_idx, blockcolumn_idx;
	double *pMat_test;
	size_t offset_copy;

	// Select kernels
	cl_command_queue queue;	
	cl_context context;
	cl_kernel kernel1, kernel2;
	cl_ushort kernel_type1, kernel_type2;
	
	pMat_test = (double*)malloc(sizeof(double)*mat_dim*mat_dim);


	
	kernel_type1= BTOCL_LTRI_LBYM;
	kernel_type2 = BTOCL_LTRI_MBYL;
	
	kernel1 = btocl_getKernel(kernel_type1);	
	if (kernel1 == NULL) {
		printf("Error: Couldn't load kernel %s\n",btocl_getKernelName(kernel_type1));
		exit(1);
	}
	
	kernel2 = btocl_getKernel(kernel_type2);	
	if (kernel2 == NULL) {
		printf("Error: Couldn't load kernel %s\n",btocl_getKernelName(kernel_type2));
		exit(1);
	}
	
	queue =  btocl_getCommandQueue();
	context = btocl_getContext();
	
	//  ------------------ blocked algorithm ------------------
	//nb = 2;  // 128/256
	if (mat_dim <= nb) {
		current_nb = mat_dim;
	} else {
		current_nb = mat_dim % nb;
		if (current_nb == 0)
			current_nb = nb;
	}
	diag = mat_dim - current_nb;
	pblockcolumn = pdiag = pMat + (mat_dim*diag +diag);
	blockcolumn_idx = diag_idx = mat_dim*diag +diag;
	nblockcolumn = current_nb;
	// Inverse of first block
	info = btlapack_ltri_inv2(pdiag,current_nb,mat_dim,NULL);
	if (info != 0) {
		printf("Error btlapack_ltri_inv2 %d\n",info);
		exit(0);
		return info;
	}
	// Copy result of new inverse to buffer
	// if buffer is empty, copy whole matrix
	offset_copy = sizeof(double)*diag*mat_dim;
	clEnqueueWriteBuffer(queue,buffer,CL_TRUE, offset_copy, current_nb*mat_dim*sizeof(double),
	pMat+diag*mat_dim,0,0,NULL);
	
	old_diag = diag;
	diag -= nb; // always go up a full block
	
	//printf("--- After first inverse\n"); 
	//btlin_print(pMat,mat_dim,mat_dim);
	
	// outer loop
	while(diag >= 0) {
		pblockcolumn = pdiag - mat_dim*nb;
		blockcolumn_idx = diag_idx - mat_dim*nb;
		old_pdiag = pdiag;
		old_diag_idx = diag_idx;
		pdiag = pblockcolumn - nb;
		diag_idx = blockcolumn_idx - nb;
		// block update 1: X  <- L1^(-1) * X
		// in-place - No extra workspace required (limited parallelism)
		
		//printf("--- After first inverse\n"); 
		//btlin_print(pMat,mat_dim,mat_dim);
		
		//printf("Before lowerbypmatrix %lf %lf\n",*old_pdiag,*pblockcolumn);
		
		// test: copy matrix to buffer
		//clEnqueueWriteBuffer(queue,buffer,CL_TRUE,0,mat_dim*mat_dim*sizeof(double),pMat,0,0,NULL);
		
		//btlin_ltri_LowerByMatrix(-1.0,old_pdiag,pblockcolumn,nblockcolumn,nb,mat_dim);
		
		// GPU version
		btocl_ltri_LowerByMatrix(kernel_type1,kernel1,buffer,-1.0,old_diag_idx,blockcolumn_idx,old_diag_idx-1,nblockcolumn,nb,mat_dim);
				
		// Inverse of new block
		info = btlapack_ltri_inv2(pdiag,nb,mat_dim,NULL);
		if (info != 0) {
			printf("Error btlapack_ltri_inv2 %d\n",info);
			exit(0);
			return info;
		}
		
		// Copy result of new inverse to buffer
		//printf("Copying new inverse from %d elements %d\n",diag*mat_dim,nb*mat_dim);
		offset_copy = sizeof(double)*diag*mat_dim;
		// blocking
		clEnqueueWriteBuffer(queue,buffer,CL_TRUE, offset_copy, nb*mat_dim*sizeof(double),
		pMat+diag*mat_dim,0,0,NULL);
		
		
		// update column block 2
		//printf("Before matrixbylower %lf %lf\n",*pblockcolumn,*pdiag);
		//btlin_ltri_MatrixByLower(1.0,pblockcolumn,pdiag,nblockcolumn,nb,mat_dim);		
		//printf("---> Calling GPU MbyL M %d\n", blockcolumn_idx);
		err = btocl_ltri_MatrixByLower(kernel_type2,kernel2,buffer,1.0,old_diag_idx-1,diag_idx,blockcolumn_idx,nblockcolumn,nb,mat_dim);
		if (err != 0) {
			printf("Error duting MbyL, argnum %d\n",err);
			exit(0);
		}
		
		// test: copy buffer to test matrix
		// commented out - it should be ok
		//clEnqueueReadBuffer(queue,buffer,CL_TRUE,0,mat_dim*mat_dim*sizeof(double),
		// pMat_test,0,0,NULL);
		//printf("Matrix after MbyL\n");
		//btlin_print(pMat,mat_dim,mat_dim);
		//printf("------------Buffer after GPU MbyL\n");
		//btlin_print(pMat_test,mat_dim,mat_dim);
		
		
		//btlin_print(pMat,n,n);
		
		// update variables
		nblockcolumn += nb;
		old_diag = diag;
		diag -= nb;
		
		//exit(0);
		
	}
	
	// Copy results 
	clEnqueueReadBuffer(queue,buffer,CL_TRUE,0,mat_dim*mat_dim*sizeof(double),pMat,0,0,NULL);
	
	// debug
	free(pMat_test);
	
	// Code that computes the determinant of the inverse of L
	if (det != NULL) {
		*det = 0.0;
		pdiag = pMat;
	//printf("Diag elements: ");
		for(i=0; i< mat_dim; i++) {
		//printf("%lf ",*pdiag);
			if (*pdiag < 0)
				*det += log(-*pdiag);
			else
				*det += log(*pdiag);
			pdiag += (mat_dim+1);
		}
		//*det = btlin_ltri_det(pMat, mat_dim, mat_dim);
	}
	
	
	return 0;
}


// Local
int btocl_ltri_LowerByMatrix(cl_ushort kernel_type,cl_kernel kernel,cl_mem buffer,double alpha,
		int L_idx,int M_idx,int MT_idx, int m, int n,int lda) {
	cl_command_queue queue;	
	int argnum,err,dim;
	size_t globalws[3];
	size_t localws[3];

	queue =  btocl_getCommandQueue();
		
	// set arguments
	argnum = 0;
	if ((err = clSetKernelArg(kernel,argnum++,sizeof(cl_mem), &buffer)) < 0)
		return 0;
	if ((err = clSetKernelArg(kernel,argnum++,sizeof(double), &alpha)) < 0)
		return 0;
	if ((err = clSetKernelArg(kernel,argnum++,sizeof(int), &L_idx)) < 0)
		return 0;
	if ((err = clSetKernelArg(kernel,argnum++,sizeof(int), &M_idx)) < 0)
		return 0;
	if ((err = clSetKernelArg(kernel,argnum++,sizeof(int), &MT_idx)) < 0)
		return 0;
	if ((err = clSetKernelArg(kernel,argnum++,sizeof(int), &m)) < 0)
		return 0;
	if ((err = clSetKernelArg(kernel,argnum++,sizeof(int), &n)) < 0)
		return 0;
	if ((err = clSetKernelArg(kernel,argnum++,sizeof(int), &lda)) < 0)
		return 0;
		
		
	dim = 1;
	globalws[0] = m*n;
	localws[0] = n;
	
	//printf("GPU m %d n %d\n",m,n);
	clEnqueueNDRangeKernel(queue,kernel,dim,NULL,globalws,localws,0,NULL,NULL);
					
	return 0;
}

int btocl_ltri_MatrixByLower(cl_ushort kernel_type,cl_kernel kernel,cl_mem buffer,double alpha,int MT_idx,int L_idx,int M_idx, int m, int n, int lda){
	cl_command_queue queue;	
	int argnum,err,dim;
	size_t globalws[3];
	size_t localws[3];
	
	queue =  btocl_getCommandQueue();
	
	// set arguments
	argnum = 0;
	if ((err = clSetKernelArg(kernel,argnum++,sizeof(cl_mem), &buffer)) < 0)
		return argnum;
	if ((err = clSetKernelArg(kernel,argnum++,sizeof(double), &alpha)) < 0)
		return argnum;
	if ((err = clSetKernelArg(kernel,argnum++,sizeof(int), &MT_idx)) < 0)
		return argnum;
	if ((err = clSetKernelArg(kernel,argnum++,sizeof(int), &L_idx)) < 0)
		return argnum;
	if ((err = clSetKernelArg(kernel,argnum++,sizeof(int), &M_idx)) < 0)
		return argnum;
	if ((err = clSetKernelArg(kernel,argnum++,sizeof(int), &m)) < 0)
		return argnum;
	if ((err = clSetKernelArg(kernel,argnum++,sizeof(int), &n)) < 0)
		return argnum;
	if ((err = clSetKernelArg(kernel,argnum++,sizeof(int), &lda)) < 0)
		return argnum;
		
		
	dim = 1;
	globalws[0] = m*n;
	localws[0] = n;
	
	//printf("MbyL ---- GPU m %d n %d\n",m,n);
	clEnqueueNDRangeKernel(queue,kernel,dim,NULL,globalws,localws,0,NULL,NULL);
	
	return 0;
}

/* ***************************************
Description: Computes L^T * L, where L is a lower triangular matrix. 
The algorithm works in-place and generates a complete symmetric matrix (L^T by L) - as opposed to generating a lower or triangular versions.

Input:
- buffer: contains the input matrix
- n: matrix dimension
Output:
- buffer: symmetric matrix containing the result of L^t by L.

  ****************************************** */
int btocl_ltri_LTbyL(cl_mem buffer,int n) {
	cl_command_queue queue;
	cl_context context;
	cl_mem buffer_diag;
	cl_int err;
	cl_kernel kernel1, kernel2;
	cl_ushort kernel_type1, kernel_type2;
	int argnum,dim, mat_idx, mat_dim, lda;
	size_t globalws[3];
	//size_t localws[3];

	
	context = btocl_getContext();
	buffer_diag = clCreateBuffer(context, CL_MEM_READ_WRITE, sizeof(double)*n, NULL, &err);
	queue =  btocl_getCommandQueue();
	
	kernel_type1 = BTOCL_LTRI_LTBYL;
	kernel_type2 = BTOCL_WRITEUTOL_DIAG;
	
	kernel1 = btocl_getKernel(kernel_type1);	
	if (kernel1 == NULL) {
		printf("Error: Couldn't load kernel %s\n",btocl_getKernelName(kernel_type1));
		return 1;
	}
	
	kernel2 = btocl_getKernel(kernel_type2);	
	if (kernel2 == NULL) {
		printf("Error: Couldn't load kernel %s\n",btocl_getKernelName(kernel_type2));
		return 1;
	}
	
	// ************* kernel1 : perform multiplication *******************
	
	dim = 1;
	globalws[0] = n*n;
	//localws[0] =...
	mat_idx = 0;
	mat_dim = lda = n;
	
	// Arguments to kernel1
	argnum = 0;
	if ((err = clSetKernelArg(kernel1,argnum++,sizeof(cl_mem), &buffer)) < 0)
		return argnum;
	if ((err = clSetKernelArg(kernel1,argnum++,sizeof(cl_mem), &buffer_diag)) < 0)
		return argnum;
	if ((err = clSetKernelArg(kernel1,argnum++,sizeof(int), &mat_idx)) < 0)
		return argnum;
	if ((err = clSetKernelArg(kernel1,argnum++,sizeof(int), &mat_dim)) < 0)
		return argnum;
	if ((err = clSetKernelArg(kernel1,argnum++,sizeof(int), &lda)) < 0)
		return argnum;

	err = clEnqueueNDRangeKernel(queue,kernel1,dim,NULL,globalws,NULL,0,NULL,NULL);

	//double* ones = (double*)malloc(sizeof(double)*n);
	//int i;
	//for(i=0;i<n;i++) ones[i]=1;
	//err = clEnqueueWriteBuffer(queue,buffer_diag,CL_TRUE,0,sizeof(double)*n,ones,0,0,NULL);
	//free(ones);
	
	// ************** kernel2 : Copy Diagonal back and do half-transpose (U -> L)
	argnum = 0;
	if ((err = clSetKernelArg(kernel2,argnum++,sizeof(cl_mem), &buffer)) < 0)
		return argnum;
	if ((err = clSetKernelArg(kernel2,argnum++,sizeof(cl_mem), &buffer_diag)) < 0)
		return argnum;
	if ((err = clSetKernelArg(kernel2,argnum++,sizeof(int), &mat_idx)) < 0)
		return argnum;
	if ((err = clSetKernelArg(kernel2,argnum++,sizeof(int), &mat_dim)) < 0)
		return argnum;
	if ((err = clSetKernelArg(kernel2,argnum++,sizeof(int), &lda)) < 0)
		return argnum;
		// Arguments to kernel1

	
	clEnqueueNDRangeKernel(queue,kernel2,dim,NULL,globalws,NULL,0,NULL,NULL);
	
	printf("Finished ltbyl %d\n",n);
	
	
	clReleaseMemObject(buffer_diag);

	return 0;
}



#endif  // if BTOCL is defined
