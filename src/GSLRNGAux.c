#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "GSLRNGAux.h"
#include "TypeDef.h"


int RandPosInProportion(gsl_rng *RNG, double *List, int No)
{
	double	Sum, SSF;
	int		Index;
	double	Point;

	Sum = 0;
	for(Index=0;Index<No;Index++)
		Sum += List[Index];

	Point = gsl_rng_uniform_pos(RNG);
	SSF = 0;
	for(Index=0;Index<No;Index++)
	{
		if(Point <= (SSF + List[Index]) / Sum)
			return Index;

		SSF += List[Index];
	}

	return 0;
}

