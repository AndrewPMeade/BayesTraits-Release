#ifdef BTOCL_CON

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "typedef.h"
#include "trees.h"
#include "continuous.h"

// OpenCL headers
#include "btocl_continuous.h"
#include "btocl_lin.h"

// Open CL to find the invers of V and log det of inv V
// V is in Tree->ConVars->V 
// InvV should be stroed in Tree->ConVars->InvV
// Log Det of V should be stroed in Tree->ConVars->LogDetOfV
void	btocl_FindInvV(TREES *Trees, TREE* Tree) 
{ 

	int err;
	int cbs;  // Cholesky Block size
	//printf("btocl findinvv\n");
	CopyMatrix(Tree->ConVars->InvV, Tree->ConVars->V);
		
	cbs = 256;
	if (Tree->ConVars->InvV->NoOfRows < 200) {
		cbs = 64;
	}
	//cbs = 32;
	
	//btdebug_enter("btoclcholesky");
	err = btocl_invcholesky(Tree->ConVars->buffer_invV, Tree->ConVars->InvV->me[0],Tree->ConVars->InvV->NoOfRows, &Tree->ConVars->LogDetOfV,cbs,32);
	//err = btocl_invcholesky_pure(Tree->ConVars->buffer_invV, Tree->ConVars->InvV->me[0],Tree->ConVars->InvV->NoOfRows, &Tree->ConVars->LogDetOfV,cbs,32);
	//btdebug_exit("btoclcholesky");	
	
	//printf("LogDetOfV=%f;\n", Tree->ConVars->LogDetOfV);
	//btlin_print(Tree->ConVars->InvV->me[0],Tree->ConVars->InvV->NoOfRows,Tree->ConVars->InvV->NoOfRows);
	
	if(err != 0)
	{
		printf("V Matrix inversion error in %s %d\n", __FILE__, __LINE__);
		PrintMathematicaMatrix(Tree->ConVars->V, "V=", stdout);
		exit(0);
	}	
	
}

void    btocl_VectByKroneckerMult(TREE* Tree) {

	// VectByKroneckerMult(Tree->ConVars->ZA, Tree->ConVars->InvSigma,
	//                     Tree->ConVars->InvV,Tree->ConVars->ZATemp);
	int mat_dim, sigma_dim;
	CONVAR * convar;
	cl_command_queue queue;

	convar = Tree->ConVars;
	mat_dim = convar->InvV->NoOfRows;
	sigma_dim = convar->InvSigma->NoOfRows;
	
	queue = btocl_getCommandQueue();
	
	// Load invSigma
	clEnqueueWriteBuffer(queue,convar->buffer_invSigma,CL_TRUE,0,sigma_dim*sigma_dim*sizeof(double),convar->InvSigma->me[0],0,0,NULL);
	// Load ZA
	clEnqueueWriteBuffer(queue,convar->buffer_ZA,CL_TRUE,0,mat_dim*sigma_dim*sizeof(double),convar->ZA,0,0,NULL);
	
	btocl_kronecker_vectmult(convar->buffer_ZATemp, convar->buffer_ZA, convar->buffer_invSigma, 
	sigma_dim, convar->buffer_invV, mat_dim);
	
	// Read ZATemp
	clEnqueueReadBuffer(queue,convar->buffer_ZATemp,CL_TRUE,0,mat_dim*sigma_dim*sizeof(double),
	convar->ZATemp,0,0,NULL);
}


void	btocl_AllocConVar(CONVAR* ConVar, TREES *Trees)
{
	cl_context context;
	int err;
	
	context = btocl_getContext();
	
	ConVar->buffer_invV = clCreateBuffer(context, CL_MEM_READ_WRITE, 
		sizeof(double)*(Trees->NoOfTaxa)*(Trees->NoOfTaxa), NULL, &err);		
	if (err != 0) {
		printf("Error allocating OpenCL buffer for InvV\n");
		exit(0);
	}
	// Assuming that ConVar has been allocated
	ConVar->buffer_invSigma = clCreateBuffer(context, CL_MEM_READ_ONLY, 
		sizeof(double)*(ConVar->InvSigma->NoOfRows)*(ConVar->InvSigma->NoOfRows), NULL, &err);		
	if (err != 0) {
		printf("Error allocating OpenCL buffer for InvV\n");
		exit(0);
	}
	ConVar->buffer_ZA = clCreateBuffer(context, CL_MEM_READ_ONLY, 
		sizeof(double)*(ConVar->InvV->NoOfRows)*(ConVar->InvSigma->NoOfRows), NULL, &err);		
	if (err != 0) {
		printf("Error allocating OpenCL buffer for InvV\n");
		exit(0);
	}
	ConVar->buffer_ZATemp = clCreateBuffer(context, CL_MEM_WRITE_ONLY, 
		sizeof(double)*(ConVar->InvV->NoOfRows)*(ConVar->InvSigma->NoOfRows), NULL, &err);		
	if (err != 0) {
		printf("Error allocating OpenCL buffer for InvV\n");
		exit(0);
	}
	

}

void	btocl_FreeConVar(CONVAR* ConVar)
{
	clReleaseMemObject(ConVar->buffer_invV);
	clReleaseMemObject(ConVar->buffer_invSigma);
	clReleaseMemObject(ConVar->buffer_ZA);
	clReleaseMemObject(ConVar->buffer_ZATemp);
}
#endif



