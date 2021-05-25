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

#include "TypeDef.h"
#include "GenLib.h"
#include "Trees.h"
#include "Continuous.h"

// OpenCL headers
#include "btocl_continuous.h"
#include "btocl_lin.h"

// Open CL to find the invers of V and log det of inv V
// V is in Tree->ConVars->V
// InvV should be stroed in Tree->ConVars->InvV
// Log Det of V should be stroed in Tree->ConVars->LogDetOfV
int	btocl_FindInvV(TREES *Trees, TREE* Tree)
{

	int err;
	int cbs;  // Cholesky Block size
	//printf("btocl findinvv\n");
	CopyMatrix(Tree->ConVars->InvV, Tree->ConVars->V);

	cbs = 256;
	while (cbs > Tree->ConVars->InvV->NoOfRows) {
		cbs /= 2;
	}
	//printf("CBS = %d\n",cbs);

	//btdebug_enter("btoclcholesky");
	//err = btocl_invcholesky(Tree->ConVars->buffer_invV, Tree->ConVars->InvV->me[0],Tree->ConVars->InvV->NoOfRows, &Tree->ConVars->LogDetOfV,cbs,32);
	err = btocl_invcholesky_pure(Tree->ConVars->buffer_invV, Tree->ConVars->InvV->me[0],Tree->ConVars->InvV->NoOfRows, &Tree->ConVars->LogDetOfV,cbs,32);
	//btdebug_exit("btoclcholesky");

	//printf("LogDetOfV=%f;\n", Tree->ConVars->LogDetOfV);
	//btlin_print(Tree->ConVars->InvV->me[0],Tree->ConVars->InvV->NoOfRows,Tree->ConVars->InvV->NoOfRows);

	if(err != 0)
	{
		return FALSE;
		printf("V Matrix inversion error in %s %d\n", __FILE__, __LINE__);
		PrintMathematicaMatrix(Tree->ConVars->V, "V=", stdout);
		exit(0);
	}
	return TRUE;
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



	if (sigma_dim == 1) {
		//btocl_kronecker_vectmult_one(cl_mem vres_buffer, cl_mem v_buffer, double sigma, cl_mem mat_buffer, int mat_dim)
		btocl_kronecker_vectmult_one(convar->buffer_ZATemp, convar->buffer_ZA, convar->InvSigma->me[0][0],
		convar->buffer_invV, mat_dim);
	} else {
		btocl_kronecker_vectmult(convar->buffer_ZATemp, convar->buffer_ZA, convar->buffer_invSigma,
		sigma_dim, convar->buffer_invV, mat_dim);
	}


	// Read ZATemp
	clEnqueueReadBuffer(queue,convar->buffer_ZATemp,CL_TRUE,0,mat_dim*sigma_dim*sizeof(double),
	convar->ZATemp,0,0,NULL);
}


void	btocl_AllocConVar(CONVAR* ConVar, TREES *Trees)
{
	cl_context context;
	int err;

	context = btocl_getContext();

	ConVar->buffer_invV = clCreateBuffer(context, CL_MEM_READ_WRITE, sizeof(double)*(Trees->NoTaxa)*(Trees->NoTaxa), NULL, &err);
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



