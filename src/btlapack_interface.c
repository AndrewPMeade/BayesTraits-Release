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

#include "TypeDef.h"

#ifdef BTLAPACK

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "Trees.h"
#include "Continuous.h"
#include "GenLib.h"


// OpenCL headers
#include "btlapack_interface.h"




// Open CL to find the invers of V and log det of inv V
// V is in Tree->ConVars->V 
// InvV should be stroed in Tree->ConVars->InvV
// Log Det of V should be stroed in Tree->ConVars->LogDetOfV
int	btlapack_FindInvV(TREES *Trees, TREE* Tree) 
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
		return FALSE;
		printf("V Matrix inverstion error in %s %d\n", __FILE__, __LINE__);
		PrintMathematicaMatrix(Tree->ConVars->V, "V=", stdout);
		exit(0);
	}	

	return TRUE;
	
	//btdebug_exit("inverse");
}

void	btlapack_InitConTree(TREES *Trees, TREE* Tree)
{

}

void	btlapack_FreeConTree(TREES* Trees, TREE* Tree)
{
	
}

#endif  // if BTLAPACK defined



