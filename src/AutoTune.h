#ifndef AUTOTUNE_H
#define AUTOTUNE_H

#include "typedef.h"

#define		AT_HSIZE	128
#define		AT_SCALE	2.0

AUTOTUNE*	CreatAutoTune(double Min, double Max);
void		FreeAutoTune(AUTOTUNE *AutoTune);

double		AutoTuneNextRD(AUTOTUNE *AutoTune, RANDSTATES *RS, double RD, double Acc);
int			AutoTuneValid(AUTOTUNE *AutoTune, double Acc);
#endif