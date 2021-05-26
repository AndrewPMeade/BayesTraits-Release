#ifdef BTLAPACK

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "typedef.h"
#include "trees.h"
#include "continuous.h"


// OpenCL headers
#include "btlapack_interface.h"




// Open CL to find the invers of V and log det of inv V
// V is in Tree->ConVars->V 
// InvV should be stroed in Tree->ConVars->InvV
// Log Det of V should be stroed in Tree->ConVars->LogDetOfV
void	btlapack_FindInvV(TREES *Trees, TREE* Tree) 
{ 
	int		err;
	
	//btdebug_enter("inverse");
	
	CopyMatrix(Tree->ConVars->InvV, Tree->ConVars->V);

	//printf("size %d\n",Tree->ConVars->InvV->NoOfRows);
	
	//err = btlapack_invldlW(Tree->ConVars->InvV->me[0], Trees->TempConVars->TMat->me[0]  , Tree->ConVars->InvV->NoOfRows, &Tree->ConVars->LogDetOfV);
	err = btlapack_invcholesky(Tree->ConVars->InvV->me[0], Tree->ConVars->InvV->NoOfRows, &Tree->ConVars->LogDetOfV);

	//printf("LogDetOfV=%f;\n", Tree->ConVars->LogDetOfV);
	//btlin_print(Tree->ConVars->InvV->me[0],Tree->ConVars->InvV->NoOfRows,Tree->ConVars->InvV->NoOfRows);

	if(err != 0)
	{
		printf("V Matrix inverstion error in %s %d\n", __FILE__, __LINE__);
		PrintMathematicaMatrix(Tree->ConVars->V, "V=", stdout);
		exit(0);
	}	
	
	//btdebug_exit("inverse");
}

void	btlapack_InitConTree(TREES *Trees, TREE* Tree)
{

}

void	btlapack_FreeConTree(TREES* Trees, TREE* Tree)
{
	
}

#endif  // if BTLAPACK defined



