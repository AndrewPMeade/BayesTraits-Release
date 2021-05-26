#ifndef SLICESAMPLER_H
#define SLICESAMPLER_H

#include "TypeDef.h"

SLICESAMPLER*	CrateSliceSampler(int NoSlices);
void			FreeSliceSampler(SLICESAMPLER* SS);

void			SSSetXPosVect(SLICESAMPLER* SS, double Min, double Max);

double			SSGetNewPoint(SLICESAMPLER *SS, RANDSTATES *RS, double POld);

#endif