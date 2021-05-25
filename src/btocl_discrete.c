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
#include <string.h>

#include "typedef.h"
#include "trees.h"
#include "Continuous.h"
#include "genlib.h"
#include "matrix.h"
// OpenCL headers



#include "btocl_discrete.h"
#include "btlin_basic.h"

#include "btdebug.h"

typedef struct {
	// Maximum dimension is 3 - we mainly use 1
	size_t globalws[3];
	size_t localws[3];  // number of wokitems per dimension
	int dim;  //determines actual length of globalws and localws
	int usewg;
	int num_nodes_wg;  // number of nodes per workgroup
	int num_rows_wg;   // if workgroups must group by rows
	size_t num_wg;
	int check_error;  // indicates of kernel checks pmatrix error
	int transpose_inv;
	cl_ushort kernel_type;
	// Params info for PartialLh computation
	int action; // 0=Collect, 1=Collect/Compute 2=Compute (tips)
	int passedCut;
	int lastGroup;
	int num_nodes;
} KernelInfo;

/* Set PMatrix functions */
void btocl_AllocInvInfo(cl_context context, INVINFO* InvInfo, int NOS, int MaxNodes);
int	btocl_SetAllPMatrixKernel(RATES *Rates, TREES *Trees, OPTIONS *Opt, double RateMult, cl_ushort kernel_type);
void btocl_loadInvInfoBuffers(cl_command_queue queue, INVINFO *InvInfo, int NOS, int NoNodes, KernelInfo k);
void printTreeTraversal(TREE* Tree);
int btocl_computeExpComponent(TREES *Trees, TREE* Tree);

/* New */
void btocl_AllocPMatrixInfo(TREES* Trees) {
	cl_context context;
	int NOS, max_nnodes;
	cl_int err;

	context = btocl_getContext();
	NOS = Trees->NoStates;
	max_nnodes = Trees->MaxNodes;

	//printf("**********************************\n");
	//printf("NOS %d max_nnodes %d\n",NOS,max_nnodes);

	// cl_mem buffer_pmatrix  --  NoNodes*NOS*NOS
	Trees->buffer_pmatrix  = btocl_clCreateBuffer("Pmatrix buffer",context, CL_MEM_READ_WRITE,
		sizeof(double)*max_nnodes*NOS*NOS, NULL, &err);

	Trees->buffer_exp_eigen  = btocl_clCreateBuffer("Eigenvalue matrix buffer",context, CL_MEM_READ_WRITE,
		sizeof(double)*max_nnodes*NOS, NULL, &err);


	Trees->buffer_error  = btocl_clCreateBuffer("Pmatrix error buffer",context, CL_MEM_WRITE_ONLY,
		sizeof(int), NULL, &err);

	btocl_AllocInvInfo(context, Trees->InvInfo[0], NOS, max_nnodes);

	// host
	Trees->check_pmatrix = (double*)malloc(sizeof(double)*max_nnodes*NOS*NOS);  // used for debugging and error checking
	Trees->perror = (int*)malloc(sizeof(int));
	return;
}



void btocl_FreePMatrixInfo(TREES* Trees) {
	INVINFO* InvInfo = Trees->InvInfo[0];

	// host
	free(Trees->perror);
	free(Trees->check_pmatrix);  // used for debugging and error checking

	free(InvInfo->vect_t);
	free(InvInfo->vect_id);
	free(InvInfo->vect_test);

	// opencl
	clReleaseMemObject(InvInfo->buffer_vec);
	clReleaseMemObject(InvInfo->buffer_inv_vec);
	clReleaseMemObject(InvInfo->buffer_val);
	clReleaseMemObject(InvInfo->buffer_t);
	clReleaseMemObject(InvInfo->buffer_id);
	clReleaseMemObject(InvInfo->buffer_temp);

	clReleaseMemObject(Trees->buffer_error);
	clReleaseMemObject(Trees->buffer_exp_eigen);
	clReleaseMemObject(Trees->buffer_pmatrix);

}



void btocl_AllocInvInfo(cl_context context, INVINFO* InvInfo, int NOS, int MaxNodes) {
	int err;

	// Host arrays for intermediate information
		// these two could be part of InvInfo, to avoid malloc and frees
	InvInfo->vect_t = (double*)malloc(sizeof(double) * MaxNodes);
	if(InvInfo->vect_t == NULL)
		MallocErr();
	InvInfo->vect_id = (int*)malloc(sizeof(int) * MaxNodes);
	if(InvInfo->vect_id == NULL)
		MallocErr();

	// ONLY for TESTING
	InvInfo->vect_test = (double*)malloc(sizeof(double) * MaxNodes* NOS * NOS);
	if(InvInfo->vect_test == NULL)
		MallocErr();

	// The OpenCL buffers

	// cl_mem vec_buffer  -  NOS*NOS  read-only
	InvInfo->buffer_vec  = btocl_clCreateBuffer("EigenVector Matrix buffer",context, CL_MEM_READ_ONLY,
		sizeof(double)*NOS*NOS, NULL, &err);

	// cl_mem buffer_inv_vec  -  NOS*NOS  read-only
	InvInfo->buffer_inv_vec  = btocl_clCreateBuffer("Eigenvector inverse matrix buffer",context, CL_MEM_READ_ONLY,
		sizeof(double)*NOS*NOS, NULL, &err);

	// 	cl_mem buffer_val  -  NOS vector read-only
	InvInfo->buffer_val  = btocl_clCreateBuffer("Eigenvalues vector buffer",context, CL_MEM_READ_ONLY,
		sizeof(double)*NOS, NULL, &err);

	// 	cl_mem buffer_t  --  NoNodes vector  read-only  -- current tree
	// allocate upper bound MaxNodes
	InvInfo->buffer_t  = btocl_clCreateBuffer("(Discrete) t-len vector buffer",context, CL_MEM_READ_ONLY,
		sizeof(double)*MaxNodes, NULL, &err);

	// cl_mem buffer_i  -- NoNodes read only
	// allocate upper bound MaxNodes
	InvInfo->buffer_id  = btocl_clCreateBuffer("Node id vector matrix",context, CL_MEM_READ_ONLY,
		sizeof(int)*MaxNodes, NULL, &err);

	// cl_mem buffer_temp -- NoNodes*NOS*NOS read-write
	InvInfo->buffer_temp  = btocl_clCreateBuffer("Temporary matrix buffer", context, CL_MEM_READ_WRITE,
		sizeof(double)*MaxNodes*NOS*NOS, NULL, &err);

}


// Loads error and eigendecomposition buffers
void btocl_loadInvInfoBuffers(cl_command_queue queue, INVINFO *InvInfo, int NOS, int NoNodes, KernelInfo k) {
	double transpose[16]; // holds NOS=2 or 4 transpose
	double* invV;
	// vec --> buffer_vec;
	clEnqueueWriteBuffer(queue,InvInfo->buffer_vec,CL_TRUE,0,NOS*NOS*sizeof(double),InvInfo->vec->me[0],0,0,NULL);
	invV = InvInfo->inv_vec->me[0];
	if (k.transpose_inv) {
		if (NOS==2) { // transpose inv_vec --> buffer_inv_vec
			transpose[0] = invV[0]; transpose[3] = invV[3];  // diagonal
			transpose[1] = invV[2]; transpose[2] = invV[1];
		} else if (NOS==4) {
			transpose[0] = invV[0];  transpose[1] = invV[4];  transpose[2] = invV[8];   transpose[3] = invV[12];
			transpose[4] = invV[1];  transpose[5] = invV[5];  transpose[6] = invV[9];   transpose[7] = invV[13];
			transpose[8] = invV[2];  transpose[9] = invV[6];  transpose[10] = invV[10]; transpose[11] = invV[14];
			transpose[12] = invV[3]; transpose[13] = invV[7]; transpose[14] = invV[11]; transpose[15] = invV[15];
		} else {
			btlin_transpose(transpose,invV,NOS);
		}
		// here we should have code that transposes any square matrix
		clEnqueueWriteBuffer(queue,InvInfo->buffer_inv_vec,CL_TRUE,0,NOS*NOS*sizeof(double),&transpose[0],0,0,NULL);
	} else { // inv_vec --> buffer_inv_vec
		clEnqueueWriteBuffer(queue,InvInfo->buffer_inv_vec,CL_TRUE,0,NOS*NOS*sizeof(double),InvInfo->inv_vec->me[0],0,0,NULL);
	}
	// val   --> buffer_val
	clEnqueueWriteBuffer(queue,InvInfo->buffer_val,CL_TRUE,0,NOS*sizeof(double),InvInfo->val,0,0,NULL);
	// vect_t -> buffer_t;    // NoNodes  read-only - copy NoNodes only
	clEnqueueWriteBuffer(queue,InvInfo->buffer_t,CL_TRUE,0,NoNodes*sizeof(double),InvInfo->vect_t,0,0,NULL);
	// vect_id -> buffer_id;    // NoNodes  read-only - copy NoNodes only
	clEnqueueWriteBuffer(queue,InvInfo->buffer_id,CL_TRUE,0,NoNodes*sizeof(int),InvInfo->vect_id,0,0,NULL);

}


// Assumption: result to test is in I->vect_Test
// unless we only need to check pmatrix
void printSetPMatrixResults(RATES *Rates, TREES *Trees, OPTIONS *Opt, double RateMult) {
	double* pmatrix;
	INVINFO* I;
	TREE* Tree;
	int NOS;

	I = Trees->InvInfo[0];
	Tree = Trees->Tree[Rates->TreeNo];
	pmatrix = Trees->PList[0]->me[0];
	NOS = Trees->NoStates;

	//printf("Result=\n");
	//btlin_printR(I->vect_test+1*NOS*NOS,NOS,NOS);

	printf("Real Value**=\n");
	btlin_printR(pmatrix+1*NOS*NOS,NOS,NOS);


}

// Prints Pmatrix Node
void printPMatrixNode(RATES *Rates, TREES *Trees, OPTIONS *Opt, double RateMult,int node) {
	double* pmatrix;
	INVINFO* I;
	TREE* Tree;
	int NOS;

	I = Trees->InvInfo[0];
	Tree = Trees->Tree[Rates->TreeNo];
	pmatrix = Trees->PList[0]->me[0];
	NOS = Trees->NoStates;

	//printf("Result=\n");
	//btlin_printR(I->vect_test+1*NOS*NOS,NOS,NOS);

	printf("PMatrix**=\n");
	btlin_printR(pmatrix+node*NOS*NOS,NOS,NOS);

	printf("V=\n");
	btlin_printR(I->vec->me[0],NOS,NOS);
	printf("invV=\n");
	btlin_printR(I->inv_vec->me[0],NOS,NOS);
	printf("Lambda=\n");
	btlin_printR(I->val,1,NOS);
	printf("t = %lf\n",I->vect_t[node]);

}

void printPMatrixCell(TREES *Trees, TREE* Tree,int node,int row,int col) {
	double* pmatrix;
	INVINFO* I;

	int NOS,i;
	double *p,*pv,*pinv,*plambda,t;

	double accum = 0.0;
	double temp, temp2;
	double *ptemp;

	I = Trees->InvInfo[0];
	pmatrix = Trees->PList[0]->me[0];
	NOS = Trees->NoStates;

	//printf("Result=\n");
	//btlin_printR(I->vect_test+1*NOS*NOS,NOS,NOS);
	printf("****** printing info ****\n");
	printf("PMatrix value = %lf\n", *(pmatrix+node*NOS*NOS +NOS*row+col));
	printf("Check PMatrix value = %lf\n", *(Trees->check_pmatrix  +node*NOS*NOS +NOS*row+col));
	//btlin_printR(pmatrix+node*NOS*NOS,NOS,NOS);

	printf("V=\n");
	btlin_printR(I->vec->me[0]+NOS*row,1,NOS);

	printf("invV=\n");
	//btlin_printR(I->inv_vec->me[0],NOS,NOS);
	p = I->inv_vec->me[0]+col;
	for(i=0; i<NOS;i++) {
		printf("%lf ",*p);
		p += NOS;
	}
	printf("\n");
	printf("Lambda=\n");
	btlin_printR(I->val,1,NOS);
	printf("t = %lf\n",I->vect_t[node]);

	printf("Computing on-the-fly...\n");

	accum = 0.0;
	pv = I->vec->me[0]+NOS*row;
	pinv = I->inv_vec->me[0]+col;
	plambda = I->val;
	t = I->vect_t[node];
	for(i = 0; i < NOS; i++) {
		temp = exp(t*(*plambda));
		temp2 = (*pv)*temp;
		accum +=  temp2*(*pinv);
		pv++;
		pinv += NOS;
		plambda++;
	}
	printf("Value = %lf\n",accum);

	ptemp = (double*)malloc(sizeof(double)*NOS);
	pv = I->vec->me[0]+NOS*row;
	pinv = I->inv_vec->me[0]+col;
	plambda = I->val;
	t = I->vect_t[node];
	for(i = 0; i < NOS; i++) {
		temp = exp(t*plambda[i]);
		ptemp[i] = (*pv) * temp;
		pv++;
	}
	accum = 0.0;
	for(i = 0; i < NOS; i++) {
		accum +=  ptemp[i]*(*pinv);
		pinv += NOS;
	}
	printf("Value2 = %lf\n",accum);
	free(ptemp);

}


void	GenTIDVectors(TREE *Tree, INVINFO* InvInfo, double RateMult)
{

	double Len;
	int Index;
	NODE N;
	double* t;
	int* id;

	t = InvInfo->vect_t;
	id = InvInfo->vect_id;

	for(Index=0;Index<Tree->NoNodes;Index++)
	{
		N = Tree->NodeList[Index];

		Len = N->Length * RateMult;

		t[Index] = Len;
		id[Index] = N->ID;
	}
}

double	CreatPMatrix(double t, INVINFO *InvInfo, double **P, int NOS, double** TempM)
{
	int i,j,k;
	double E1, E2, TempD;

	double *Val;
	double	**Vec, **InvVec;

	Val		= InvInfo->val;
	Vec		= InvInfo->vec->me;
	InvVec	= InvInfo->inv_vec->me;

	//printf("t=<%lf>\n",t);
	// If using temp, use a single call to update Temp with exp result.
	for(i=0;i<NOS;i++)
	{
		TempD = exp(t * Val[i]);
		for(j=0;j<NOS;j++)
			TempM[j][i] = Vec[j][i] * TempD;
	}

	E1 = 0.0;

	for(i=0;i<NOS;i++)
	{
		E2 = 0;
		for(j=0;j<NOS;j++)
		{
			P[i][j] = 0;
			for(k=0;k<NOS;k++)
				P[i][j] += TempM[i][k] * InvVec[k][j];
			//printf("[%lf]",P[i][j]);
			// P[i][j] calculated
			// assert P[i][j] = exp(t * InvInfo->A->me[i][j])
			//printf("(%lf eqt %lf) ",P[i][j], exp(t * InvInfo->Q->me[i][j]) );
			E2 += P[i][j];

			if(P[i][j] < 0)
				return 1.0;
		}

		E2 = E2 - 1.0;
		E1 += E2 * E2;
	}

	return E1;
}

// Sets the actual parameters of the expqt kernels
// The function returns the total number of parameters processed (even if the last one failed) and
// sets 'err' to the error code returned by opencl when setting the last parameter.
// The function stops if OpenCl fails to set a parameter (err < 0).

// The functions sets the 4th argument to the number of the next argument to be set (last argument plus one).
int setExpqtArgs(cl_kernel kernel, TREES* Trees, TREE* Tree, KernelInfo k, cl_ushort kernel_type, int* err) {
	INVINFO* I;
	int NOS, NoNodes;
	int argnum, local_memSize, num_cells;

	I = Trees->InvInfo[0];
	NOS = Trees->NoStates;
	NoNodes =  Tree->NoNodes;

	argnum = 0;
	if ((*err = clSetKernelArg(kernel,argnum++,sizeof(cl_mem), &(I->buffer_vec))) < 0)
		return argnum;
	if ((*err = clSetKernelArg(kernel,argnum++,sizeof(cl_mem), &(I->buffer_val))) < 0)
		return argnum;
	if ((*err = clSetKernelArg(kernel,argnum++,sizeof(cl_mem), &(I->buffer_inv_vec))) < 0)
		return argnum;
	if ((*err = clSetKernelArg(kernel,argnum++,sizeof(cl_mem), &(I->buffer_t))) < 0)
		return argnum;
	if ((*err = clSetKernelArg(kernel,argnum++,sizeof(cl_mem), &(I->buffer_id))) < 0)
		return argnum;
	if ((*err = clSetKernelArg(kernel,argnum++,sizeof(cl_mem), &(Trees->buffer_pmatrix))) < 0)
		return argnum;

	if (kernel_type == BTOCL_EXPQT_NOTEMP) {
		if ((*err = clSetKernelArg(kernel,argnum++,sizeof(NOS), &NOS)) < 0)
			return argnum;
		num_cells = NOS*NOS*NoNodes;
		if ((*err = clSetKernelArg(kernel,argnum++,sizeof(num_cells), &num_cells)) < 0)
			return argnum;
	} else if ((kernel_type == BTOCL_EXPQT3) || (kernel_type == BTOCL_EXPQT3_NOS4)) {
		if ((*err = clSetKernelArg(kernel,argnum++,sizeof(cl_mem), &(Trees->buffer_error))) < 0)
			return argnum;
		num_cells = NOS*NOS*(NoNodes-1); // skipping root
		//printf("num cells = %d\n", num_cells);
		if ((*err = clSetKernelArg(kernel,argnum++,sizeof(num_cells), &num_cells)) < 0)
			return argnum;
	} else if ((kernel_type == BTOCL_EXPQT3_LOCAL) || (kernel_type == BTOCL_EXPQT3_LOCALERR) ) {
		if ((*err = clSetKernelArg(kernel,argnum++,sizeof(cl_mem), &(Trees->buffer_error))) < 0)
			return argnum;
		num_cells = NOS*NOS*(NoNodes-1); // skipping root
		//printf("num cells = %d\n", num_cells);
		if ((*err = clSetKernelArg(kernel,argnum++,sizeof(num_cells), &num_cells)) < 0)
			return argnum;
		local_memSize = sizeof(double)*NOS*(k.num_rows_wg);
		if ((*err = clSetKernelArg(kernel,argnum++,local_memSize, NULL)) < 0)  // local_invV
			return argnum;
	} else if (kernel_type == BTOCL_EXPQT) {
		if ((*err = clSetKernelArg(kernel,argnum++,sizeof(NOS), &NOS)) < 0)
			return argnum;
		if ((*err = clSetKernelArg(kernel,argnum++,sizeof(cl_mem), &(I->buffer_temp))) < 0)
			return argnum;
	} else if ((kernel_type == BTOCL_EXPQT_ROWNOS2) || (kernel_type == BTOCL_EXPQT_ROWNOS4)) {
		if ((*err = clSetKernelArg(kernel,argnum++,sizeof(NoNodes), &NoNodes)) < 0) {
			return argnum;
		}
		if ((*err = clSetKernelArg(kernel,argnum++,sizeof(k.num_nodes_wg), &k.num_nodes_wg)) < 0) {
			return argnum;
		}
		if ((*err = clSetKernelArg(kernel,argnum++,sizeof(cl_mem), &(Trees->buffer_error))) < 0)
			return argnum;
	} else if (kernel_type == BTOCL_EXPQT_LOCAL) {
		if ((*err = clSetKernelArg(kernel,argnum++,sizeof(NOS), &NOS)) < 0)
			return argnum;
		// int nonodes, int num_nodes_wg,  __local double* local_invV, __local double* local_V
		if ((*err = clSetKernelArg(kernel,argnum++,sizeof(cl_mem), &(I->buffer_temp))) < 0)
			return argnum;
		if ((*err = clSetKernelArg(kernel,argnum++,sizeof(NoNodes), &NoNodes)) < 0)
			return argnum;
		if ((*err = clSetKernelArg(kernel,argnum++,sizeof(k.num_nodes_wg), &k.num_nodes_wg)) < 0)
			return argnum;

		local_memSize = sizeof(double)*NOS*NOS;
		if ((*err = clSetKernelArg(kernel,argnum++,local_memSize, NULL)) < 0)  // local_invV
			return argnum;

		if ((*err = clSetKernelArg(kernel,argnum++,local_memSize, NULL)) < 0)  // local_V
			return argnum;
	} else if (kernel_type == BTOCL_EXPQT_LOCALERR) {
		if ((*err = clSetKernelArg(kernel,argnum++,sizeof(NOS), &NOS)) < 0)
			return argnum;
		// int nonodes, int num_nodes_wg,  __local double* local_invV, __local double* local_V
		if ((*err = clSetKernelArg(kernel,argnum++,sizeof(cl_mem), &(I->buffer_temp))) < 0)
			return argnum;
		printf("Args - NoNodes = %d\n",NoNodes);
		if ((*err = clSetKernelArg(kernel,argnum++,sizeof(NoNodes), &NoNodes)) < 0)
			return argnum;
		if ((*err = clSetKernelArg(kernel,argnum++,sizeof(k.num_nodes_wg), &k.num_nodes_wg)) < 0)
			return argnum;

		local_memSize = sizeof(double)*NOS*NOS;
		if ((*err = clSetKernelArg(kernel,argnum++,local_memSize, NULL)) < 0)  // local_invV
			return argnum;

		if ((*err = clSetKernelArg(kernel,argnum++,local_memSize, NULL)) < 0)  // local_V
			return argnum;

		if ((*err = clSetKernelArg(kernel,argnum++,sizeof(cl_mem), &(Trees->buffer_error))) < 0)
			return argnum;
	}

	return argnum;
}

void setExpqtKernelInfo(TREES* Trees,TREE* Tree,cl_ushort kernel_type,KernelInfo* k) {
	int NOS,NoNodes;
	int wg_size;
	NOS = Trees->NoStates;
	NoNodes =  Tree->NoNodes;
	// IMPORTANT - check before compiling

	k->dim = 1;
	k->check_error = 0; // default. Kernel does not check for pmatrix error
	k->transpose_inv = 0;  // default: do not use transpose

	k->kernel_type = kernel_type;

	if (kernel_type == BTOCL_EXPQT) {
		//case BTOCL_EXPQT_NOTEMP:  -- no workgroup version
		k->usewg = FALSE;
		k->globalws[0] = NOS*NOS*NoNodes; // one worker per cell
	} else if (kernel_type == BTOCL_EXPQT_NOTEMP) {
		k->usewg = TRUE;
		k->localws[0] = 512; // 'n' workitems per workgroup
		// NoNodes*NOS*NOS total number of cells
		k->num_wg = NoNodes*NOS*NOS / k->localws[0];
		if (NoNodes*NOS*NOS %  k->localws[0] != 0 )
			k->num_wg++;
		k->globalws[0] = k->localws[0]*k->num_wg;
	} else if ((kernel_type == BTOCL_EXPQT3) || (kernel_type == BTOCL_EXPQT3_LOCAL) ||
				(kernel_type == BTOCL_EXPQT3_LOCALERR)) {
		k->usewg = TRUE;
		k->check_error = 1;
		wg_size = 256;  // ideal size
		k->num_rows_wg = wg_size / NOS;
		k->localws[0] = k->num_rows_wg * NOS;  // real size
		k->num_wg = (NoNodes-1)*NOS*NOS / k->localws[0];
		// we are skipping the root node -> (NoNodes-1)
		if ((NoNodes-1)*NOS*NOS %  k->localws[0] != 0 )
			k->num_wg++;
		k->globalws[0] = k->localws[0]*k->num_wg;
		//printf("NOS %d NoNOdes %d local %d global %d\n",NOS,NoNodes,k->localws[0],k->globalws[0]);
	} else if (kernel_type == BTOCL_EXPQT3_NOS4) {  // uses private memory instead of local memory
		k->usewg = TRUE;
		k->check_error = 1;
		k->transpose_inv = 1;
		wg_size = 64;  // ideal size
		k->num_rows_wg = wg_size / NOS;
		k->localws[0] = k->num_rows_wg * NOS;  // real size
		k->num_wg = (NoNodes-1)*NOS*NOS / k->localws[0];
		// we are skipping the root node -> (NoNodes-1)
		if ((NoNodes-1)*NOS*NOS %  k->localws[0] != 0 )
			k->num_wg++;
		k->globalws[0] = k->localws[0]*k->num_wg;
		//printf("NOS %d NoNOdes %d local %d global %d\n",NOS,NoNodes,k->localws[0],k->globalws[0]);
	} else if ((kernel_type == BTOCL_EXPQT_ROWNOS2) || (kernel_type == BTOCL_EXPQT_ROWNOS4)) {
		k->usewg = TRUE;
		k->check_error = 1;
		k->transpose_inv = 1;
		wg_size = 32;  // multiple of 2 or 4
		k->num_nodes_wg = wg_size / NOS;  // should be an exact integer divivsion
		// k->num_nodes_wg =40;  // change this to optimise "warp"
		k->num_wg = NoNodes / k->num_nodes_wg;
		if (NoNodes % k->num_nodes_wg != 0)
			k->num_wg++;
		k->localws[0] = NOS*k->num_nodes_wg;        // row/column based
		k->globalws[0] = k->localws[0]*k->num_wg;
		//printf("NOS %d NoNOdes %d local %d global %d\n",NOS,NoNodes,k->localws[0],k->globalws[0]);
	} else if ((kernel_type == BTOCL_EXPQT_LOCALERR) || (kernel_type == BTOCL_EXPQT_LOCAL)) {
		if (kernel_type == BTOCL_EXPQT_LOCALERR) {
			k->check_error = 1;
		}
		k->usewg = TRUE;
		// start from desired wg_size
		wg_size = 1024;  // target. Final size will depend on nos
		k->num_nodes_wg = wg_size / (NOS*NOS);
		//k->num_nodes_wg = 8;  // number of nodes per workgroup
		if (k->num_nodes_wg ==0) { // This may blowup the workitem allocation per wg
			printf("Insuffient workgroup size\n");
			// check against max wg size to avoid blowup
			k->num_nodes_wg = 1;  // wg_size increases
		}
		k->num_wg = (NoNodes-1) / k->num_nodes_wg;
		if ((NoNodes-1) % k->num_nodes_wg != 0)  // not exact, then go one up
			k->num_wg++;
		k->localws[0] = NOS*NOS*k->num_nodes_wg; // workitems per workgroup
		k->globalws[0] = k->localws[0]*k->num_wg;
		//printf("Num nodes wg %d Num workgroups %d\n", k->num_nodes_wg, k->num_wg );
		//printf("local %d global %d\n",k->localws[0],k->globalws[0]);
	} else {
		printf("EXPQT: Unknown Kernel\n");
		exit(0);
	}


}

// Used for debugging
// Copies PList -> check_pmatrix
void copy_pmatrix(TREES* Trees) {
	int n,row,col;
	int nos, max_nnodes;
	double *pmatrix,*p;

	nos = Trees->NoStates;
	pmatrix = Trees->PList[0]->me[0];
	p = Trees->check_pmatrix;
	max_nnodes = Trees->MaxNodes;

	//printf("Copying Pmatrix to check_pmatrix\n");

	// All nodes: max_nnodes
	for(n=0; n < max_nnodes; n++) {
		//printf("n=%d\n",n);
		//printf("Copy node %d ",n);
		for(row=0; row<nos; row++) {
			for(col=0; col < nos; col++) {
				//printf("{%lf}",*pmatrix);
				*p = *pmatrix;
				p++; pmatrix++;
			}
		}

	}
	printf("\n");
}

// Used for debugging
// Compares PList with check_pmatrix
void compare_pmatrix(TREES* Trees, TREE* Tree) {
	int n,row,col;
	int nos, max_nnodes;
	double* pmatrix,*p;
	int* IVect;

	IVect = Trees->InvInfo[0]->vect_id;

	nos = Trees->NoStates;
	pmatrix = Trees->PList[0]->me[0];
	p = Trees->check_pmatrix;
	max_nnodes = Trees->MaxNodes;

	// Skip root
	pmatrix += nos*nos;
	p += nos*nos;
	for(n=1; n < max_nnodes; n++)
		for(row=0; row < nos; row++)
			for(col=0; col < nos; col++) {
				if (fabs(*p - *pmatrix) > 0.0001) {
					printPMatrixCell(Trees, Tree,n,row,col);
					printf("-------- ERROR!\n");
					printf("NOS %d\n", nos);
					printf("Error node %d row %d col %d\n",n,row,col);
					printf("pmatrix %lf new %lf\n",*p,*pmatrix);
					printf("Number of Nodes = %d\n",Tree->NoNodes);
					printf("Number of MAx Nodes = %d\n",max_nnodes);
					pmatrix -= col; p -= col;
					for(col=0; col < nos; col++) {
						printf("%lf (%lf) ",*pmatrix++, *p++);
					}

					exit(0);
				}
				p++; pmatrix++;
			}
}

int checkPMatrixError(double* pmatrix, INVINFO* InvInfo, int NoNodes, int NOS) {
	int i,nid, row,col;
	double* pmatrix_ptr;
	double E1, E2;
	int* IVect;
	int pmatrix_idx;

	//printf("here!\n");

	IVect = InvInfo->vect_id;

	//printf("NoNodes=%d\n",NoNodes);
	for(i = 1; i < NoNodes; i++) {
		//printf("i=%d - ",i);
		nid = IVect[i];
		//printf("nid=%d\n",nid);
		pmatrix_ptr = pmatrix + NOS*NOS*nid;
		pmatrix_idx =  NOS*NOS*nid;
		E1 = 0.00;  // pmatrix[nid] accumulated error
		for(row=0; row<NOS;row++) {
			for(col=0; col<NOS; col++) {
				if (*pmatrix_ptr < 0.0) {  // negative probability
					printf("i=%d - ",i);
					//printf("nid=%d\n",nid);
					printf("row %d col %d\n",row,col);
					printf("negative %lf\n",*pmatrix_ptr);
					// print row
					//pmatrix_ptr  -= col;
					//for(col=0; col < NOS; col++) printf("%lf ",*pmatrix_ptr++);
					//printf("\n");
					//printf("NEGATIVE!!!\n"); exit(0);
					return TRUE;
					//return pmatrix_idx;
				}
				E2 += *pmatrix_ptr;
				pmatrix_ptr++;
				pmatrix_idx++;
			}
			// check row
			E2 = E2 - 1.0;
			E1 += E2 * E2;
		}
		// Check pmatrix error
		if (E1 > 0.001) {
			printf("greater than 0.001\n");
				//printf("Real Value=\n");
				//btlin_printR(pmatrix+nid*NOS*NOS,NOS,NOS);
			return TRUE;
		}
	}
	return FALSE;
}

// Sequential version
int	btocl_SetAllPMatriCPU(RATES *Rates, TREES *Trees, OPTIONS *Opt, double RateMult)
{
	double *TVect;
	int *IVect;
	TREE *Tree;
	int Index, NOS;
	double Err;
	double **TempM;

	//btdebug_enter("pmatrix");
	NOS = Trees->NoStates;

	TempM = AllocMatMem(NOS, NOS);

	Tree = Trees->Tree[Rates->TreeNo];
	GenTIDVectors(Tree, Trees->InvInfo[0], RateMult);
	TVect = Trees->InvInfo[0]->vect_t;
	IVect = Trees->InvInfo[0]->vect_id;

	for(Index=1;Index<Tree->NoNodes;Index++)
	{
		//N = Tree->NodeList[Index];
		Err = CreatPMatrix(TVect[Index], Trees->InvInfo[0], Trees->PList[IVect[Index]]->me , NOS, TempM);
		if(Err > 0.001)
		{
			printf("quitting pmatrix calculation due to error\n");
			return TRUE;
		}
		//	PrintMatrix(Trees->PList[N->ID], "p=", stdout);
	}

	FreeMatMem(TempM);
	//btdebug_exit("pmatrix");
	//printf("NoNodes %d MaxNoNodes %d NOS %d\n",Tree->NoNodes, Trees->MaxNodes,NOS);
	//exit(0);



	return FALSE;
}

// Entry Point for GPU version
int	btocl_SetAllPMatrix(RATES *Rates, TREES *Trees, OPTIONS *Opt, double RateMult) {
	int NOS;
	cl_ushort kernel_type;
	int err;
	double* pmatrix;
	TREE* Tree;
	//int standard_err;

	NOS = Trees->NoStates;

	switch(NOS) {
		case 2:
			kernel_type = BTOCL_EXPQT_ROWNOS2;
			//kernel_type = BTOCL_EXPQT_NOTEMP;
			//kernel_type = BTOCL_EXPQT3_LOCALERR;
			//kernel_type = BTOCL_EXPQT3_LOCALERR;
			//kernel_type = BTOCL_EXPQT3;
			//kernel_type = BTOCL_EXPQT_LOCALERR;
			break;
		case 4:
			//kernel_type = BTOCL_EXPQT_ROWNOS4;
			//kernel_type = BTOCL_EXPQT_NOTEMP;
			//kernel_type = BTOCL_EXPQT3_NOS4;
			kernel_type = BTOCL_EXPQT3_LOCALERR;
			//kernel_type = BTOCL_EXPQT_LOCALERR;
			//kernel_type = BTOCL_EXPQT3_LOCALERR;
			break;
		default:
			//printf("Using EXPQT3_LOCAL...\n");
			kernel_type = BTOCL_EXPQT3_LOCALERR;
			//kernel_type = BTOCL_EXPQT3;
			//kernel_type = BTOCL_EXPQT_LOCALERR;  -- used as default?????
			//kernel_type = BTOCL_EXPQT_LOCAL;
			//kernel_type = BTOCL_EXPQT_LOCALERR;  // modify to check error
	}


	// removed feb 2014
	//standard_err = btocl_SetAllPMatriCPU(Rates, Trees, Opt, RateMult);

	pmatrix = Trees->PList[0]->me[0];
	Tree = Trees->Tree[Rates->TreeNo];
	//standard_err = checkPMatrixError(pmatrix, Trees->InvInfo, Tree->NoNodes, NOS);
	//if (standard_err == TRUE) { printf("-->Standard version: error\n"); }

	// copy result to check_pmatrix
	//copy_pmatrix(Trees);

	//printf("start setpmatrix %d --\n",Trees->count++);
	//printf("start setpmatrix  --");
	//printf("Kernel %s\n",btocl_getKernelName(kernel_type));
	err = btocl_SetAllPMatrixKernel(Rates,Trees,Opt,RateMult,kernel_type);

	//if (((standard_err != FALSE) && (err == FALSE)) ||
	//	((standard_err == FALSE) && (err != FALSE))) {
	//	printf("STOP!!\n");
	//	exit(0);
	//}

	// removed feb 2014
	//printf("Comparing......\n");
	//if ((standard_err == FALSE) && (err == FALSE)) {
	//	printf("Comparing......\n");
	//	compare_pmatrix(Trees,Tree);
	//}
	//printf("Comparison OK!!\n");

	//printf(" end setpmatrix  --\n");

	//exit(0);


	return err;

}


// General function that sets up kernel arguments and properties for OpenCL version of SetAllPMatrix.
int	btocl_SetAllPMatrixKernel(RATES *Rates, TREES *Trees, OPTIONS *Opt, double RateMult, cl_ushort kernel_type) {
	TREE *Tree;
	int NOS, MaxNodes, NoNodes;
	INVINFO* I;
	double* pmatrix;

	cl_command_queue queue;
	cl_kernel kernel;
	cl_context context;

	KernelInfo k;

	int err;
	// int errCorrect;
	int argnum;

	//btdebug_enter("pmatrix");
	I = Trees->InvInfo[0];
	NOS = Trees->NoStates;
	MaxNodes = Trees->MaxNodes;

	Tree = Trees->Tree[Rates->TreeNo];
	NoNodes =  Tree->NoNodes;
	GenTIDVectors(Tree, Trees->InvInfo[0], RateMult);

	queue =  btocl_getCommandQueue();
	context = btocl_getContext();

	kernel = btocl_getKernel(kernel_type);
	if (kernel == NULL) {
		printf("Error: Couldn't load kernel %s\n",btocl_getKernelName(kernel_type));
		exit(1);
	}

	//printf("Loaded %s\n",btocl_getKernelName(kernel_type));
	//print("Rate TreeNo %d\n",Rates->TreeNo);

	//printf("NOS %d MaxNodes %d\n",NOS,MaxNodes);

	// decide on workgroups and workitems
	// Update k
	setExpqtKernelInfo(Trees,Tree,kernel_type,&k);

	// the rest
	btocl_loadInvInfoBuffers(queue,Trees->InvInfo[0], NOS, Tree->NoNodes,k);

	// testing
	//btdebug_enter("expeigen");
	//btocl_computeExpComponent(Trees,Tree);
	//btdebug_exit("expeigen");
	//exit(0);

	// error is always the last parameter (if used)
	*(Trees->perror) = 0;
	if (k.check_error) {
		clEnqueueWriteBuffer(queue,Trees->buffer_error,CL_TRUE,0,sizeof(int),(Trees->perror),0,0,NULL);
	}

	argnum = setExpqtArgs(kernel,Trees,Tree,k,kernel_type,&err);
	if (err < 0) {
		printf("Error. Couldn't set argument: %d\n",argnum-1);
	}

	//printf("NOS %d Nodes %d\n",NOS,NoNodes);
	//printf("usewg %d dim %d numnodeswg %d numwg %d  local %d global %d\n",k.usewg,k.dim,
	//k.num_nodes_wg,k.num_wg,k.localws[0], k.globalws[0]);


	//printf("launched kernel...");
	if (k.usewg == FALSE) {
		err = clEnqueueNDRangeKernel(queue,kernel,k.dim,NULL,k.globalws,NULL,0,NULL,NULL);
	} else {
		err = clEnqueueNDRangeKernel(queue,kernel,k.dim,NULL,k.globalws,k.localws,0,NULL,NULL);
	}
	//printf("returned from kernel\n");
	if (err != CL_SUCCESS) {
		printf("Error kernel execution ExppQt\n");
		exit(0);
	}

	pmatrix = Trees->PList[0]->me[0];

	// Read pMatrix buffer
	// Trees->PList
	//clEnqueueReadBuffer(queue,I->buffer_pmatrix,CL_TRUE,0,
	//			sizeof(double)*MaxNodes*NOS*NOS, I->vect_test,0,0,NULL);
	clEnqueueReadBuffer(queue,Trees->buffer_pmatrix,CL_TRUE,0,
				sizeof(double)*MaxNodes*NOS*NOS, pmatrix,0,0,NULL);


	if (k.check_error) {
		// read the error back
		clEnqueueReadBuffer(queue,Trees->buffer_error,CL_TRUE,0,
				sizeof(int), Trees->perror,0,0,NULL);
		err = *(Trees->perror);
	} else {
		err = checkPMatrixError(pmatrix, I,NoNodes,NOS);
	}

	if (err != FALSE) {
		//printf("-->btocl version: error\n");


		//printf("Checking with CPU error function\n");
		//err = checkPMatrixError(pmatrix, I,NoNodes,NOS);
		//if (err == FALSE) {
		//	printf("NO ERROR!!...quitting\n");
		//	exit(0);
		//}
	}

	//printf("ERROR=%d\n", err);
	// or call check error function

	//printSetPMatrixResults(Rates, Trees, Opt, RateMult);
	// exit(0);



	//return

	//printf("Node 1\n");
	//printPMatrixNode(Rates, Trees, Opt, RateMult, 1);

	//errCorrect = checkPMatrixError(pmatrix, I,NoNodes,NOS);
	//printf("Correct error = %d\n",errCorrect);

	//printf("finished error checking %d\n",scounter++);

	//int i, pidx,local_idx;
	//if (err != FALSE) {
	//	//printf("extra ERROR in pmatrix\n");
	//	printf("err with type = %d\n",err);
	//	err = err/10;
	//	pidx = err;
	//	printf("err cell = %d\n",err);
	//	printf("Pmatrix value %lf ( %lf )\n",pmatrix[err], Trees->check_pmatrix[err]);
	//	err = err - (err%NOS);
	//	printf("Row :\n");
	//	for(i=0; i < NOS;i++) {
	//		printf("%lf ( %lf )",pmatrix[err],  Trees->check_pmatrix[err]);
	//		err++;
	//	}
	//	printf("\n");
	//	local_idx = pidx % (NOS*NOS);
	//	printPMatrixCell(Trees, Tree,pidx/(NOS*NOS), local_idx/NOS,local_idx%NOS);

	//	exit(0);
	//}

	//printf("Rate TreeNo %d\n",Rates->TreeNo);

	if (err != FALSE) {
		err = TRUE;
	}


	return err;
	//return FALSE;
}

// Assumption: loadInvInfo has been called - OpencCL buffers loaded
// Trees: Number of states, Eigenvalues
// Tree: Number of Nodes
int btocl_computeExpComponent(TREES *Trees, TREE *Tree) {

	int NOS, MaxNodes, NoNodes;
	INVINFO* I;
	cl_ushort kernel_type;

	cl_command_queue queue;
	cl_kernel kernel;
	cl_context context;


	int argnum,err, dim,i;
	size_t globalws;

	double *expeigen,*p;

	// ** end variable declaration

	NOS = Trees->NoStates;
	NoNodes =  Tree->NoNodes;
	I = Trees->InvInfo[0];
	MaxNodes = Trees->MaxNodes;  // real size of pmatrix


	queue =  btocl_getCommandQueue();
	context = btocl_getContext();

	kernel_type = BTOCL_EXP_EIGEN;
	kernel = btocl_getKernel(kernel_type);
	if (kernel == NULL) {
		printf("Error: Couldn't load kernel %s\n",btocl_getKernelName(kernel_type));
		exit(1);
	}

	// Kernel Arguments
	argnum = 0;
	if ((err = clSetKernelArg(kernel,argnum++,sizeof(cl_mem), &(I->buffer_val))) < 0)
		return err;
	if ((err = clSetKernelArg(kernel,argnum++,sizeof(cl_mem), &(I->buffer_t))) < 0)
		return err;
	if ((err = clSetKernelArg(kernel,argnum++,sizeof(NOS), &NOS)) < 0)
		return err;
	if ((err = clSetKernelArg(kernel,argnum++,sizeof(cl_mem), &(Trees->buffer_exp_eigen))) < 0)
		return err;

	dim = 1;
	globalws = NOS*NoNodes;
	//printf("NOS %d NoNodes %d MaxNodes %d Number worktimes %d\n",NOS,NoNodes, MaxNodes,globalws);
	err = clEnqueueNDRangeKernel(queue,kernel,dim,NULL,&globalws,NULL,0,NULL,NULL);

	if (err != CL_SUCCESS) {
		printf("Error: Enqueue kernel\n");
		return err;
	}

	//expeigen = (double*)malloc(NOS*MaxNodes);

	//err = clEnqueueReadBuffer(queue,Trees->buffer_exp_eigen,CL_TRUE,0,
	//			sizeof(double)*MaxNodes*NOS,expeigen ,0,0,NULL);
	//if (err != CL_SUCCESS) {
	//	printf("Error: Couldn't read exp_eigen buffer\n");
	//	btocl_printRuntimeError(err);
	//	return err;
	//}


	//p = expeigen;
	//for(i=0; i < 20; i++,p++) {
	//	printf("%lf ",*p);
	//}
	//printf("\n");




	//btocl_free_runtime();

	//printf("freeing\n");
	//free(expeigen);
	//printf("done\n");



	return 0; // success
}


#endif



