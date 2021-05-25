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
	//printf("btocl findinvv\n");
	CopyMatrix(Tree->ConVars->InvV, Tree->ConVars->V);
		
	if (err < 0) {
		printf("Couldn't create matrix buffer\n");
		exit(1);
	}
	
	//btdebug_enter("btoclcholesky");
	err = btocl_invcholesky(Tree->ConVars->buffer_invV, Tree->ConVars->InvV->me[0],Tree->ConVars->InvV->NoOfRows, &Tree->ConVars->LogDetOfV,4,2);
	//btdebug_exit("btoclcholesky");	
	
	//printf("LogDetOfV=%f;\n", Tree->ConVars->LogDetOfV);
	btlin_print(Tree->ConVars->InvV->me[0],Tree->ConVars->InvV->NoOfRows,Tree->ConVars->InvV->NoOfRows);
	
	if(err != 0)
	{
		printf("V Matrix inverstion error in %s %d\n", __FILE__, __LINE__);
		PrintMathematicaMatrix(Tree->ConVars->V, "V=", stdout);
		exit(0);
	}	
	
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

}

void	btocl_FreeConVar(CONVAR* ConVar)
{
	clReleaseMemObject(ConVar->buffer_invV);
}
#endif



