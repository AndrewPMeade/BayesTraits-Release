/*
*  BayesTriats 4.0
*
*  copyright 2022
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



#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

#include "TypeDef.h"
#include "GenLib.h"
#include "RJLocalScalar.h"
#include "Priors.h"

TRANSFORM_TYPE	NameToRJLocalType(char *Name, int *Err)
{
	int Index;

	MakeLower(Name);

	*Err = FALSE;

	for(Index=0;Index<NO_RJ_LOCAL_SCALAR;Index++)
	{
		if(StrICmp(Name, RJ_LOCAL_SCALAR_NAMES[Index]) == 0)
			return (TRANSFORM_TYPE)Index;
	}

	*Err = TRUE;
	return (TRANSFORM_TYPE)0;
}


int	UseRJLocalScalars(OPTIONS *Opt)
{
	int Index;

	for(Index=0;Index<NO_RJ_LOCAL_SCALAR;Index++)
	{
		if(Opt->UseRJLocalScalar[Index] == TRUE)
			return TRUE;
	}

	return FALSE;
}


PRIOR*	GetPriorFromRJRatesScalar(OPTIONS *Opt, TRANSFORM_TYPE Type)
{
	if(Type == VR_KAPPA)
		return GetPriorFromName("Kappa", Opt->AllPriors, Opt->NoAllPriors);

	if(Type == VR_LAMBDA)
		return GetPriorFromName("Lambda", Opt->AllPriors, Opt->NoAllPriors);

	if(Type == VR_DELTA)
		return GetPriorFromName("Delta", Opt->AllPriors, Opt->NoAllPriors);

	if(Type == VR_OU)
		return GetPriorFromName("OU", Opt->AllPriors, Opt->NoAllPriors);

	if(Type == VR_NODE)
		return GetPriorFromName("VRNode", Opt->AllPriors, Opt->NoAllPriors);

	if(Type == VR_BL)
		return GetPriorFromName("VRBranch", Opt->AllPriors, Opt->NoAllPriors);

	if(Type == VR_FABRIC_BETA)
		return GetPriorFromName("FabricBeta", Opt->AllPriors, Opt->NoAllPriors);

	printf("Unknown transform type");
	exit(1);
	return NULL;
}