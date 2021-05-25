
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

#include "typedef.h"
#include "genlib.h"
#include "RJLocalScalar.h"

RJ_VARRATE_TYPE	NameToRJLocalType (char *Name)
{
	int Index;

	MakeLower(Name);

	for(Index=0;Index<NO_RJ_LOCAL_SCALAR;Index++)
	{
		if(strcmp(Name, RJ_LOCAL_SCALAR_NAMES[Index]) == 0)
			return (RJ_VARRATE_TYPE)Index;
	}

	return (RJ_VARRATE_TYPE)-1;
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

PRIORS*	GetPriorFromRJRatesScalar(OPTIONS *Opt, RJ_VARRATE_TYPE Type)
{
	return Opt->RJLocalScalarPriors[Type];
}