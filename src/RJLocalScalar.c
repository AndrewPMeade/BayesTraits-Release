
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

#include "typedef.h"
#include "genlib.h"
#include "RJLocalScalar.h"
#include "priors.h"

TRANSFORM_TYPE	NameToRJLocalType (char *Name)
{
	int Index;

	MakeLower(Name);

	for(Index=0;Index<NO_RJ_LOCAL_SCALAR;Index++)
	{
		if(strcmp(Name, RJ_LOCAL_SCALAR_NAMES[Index]) == 0)
			return (TRANSFORM_TYPE)Index;
	}

	return (TRANSFORM_TYPE)-1;
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

	printf("Unknown transform type");
	exit(1);
	return NULL;
}