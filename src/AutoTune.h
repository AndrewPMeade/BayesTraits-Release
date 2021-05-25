#ifndef AUTOTUNE_H
#define AUTOTUNE_H

#include "typedef.h"

#define		AT_HSIZE		128
#define		AT_SCALE		2.0
#define		AT_MIN_TRIED	10					

typedef struct 
{
	int No;
	
	double	CDev;

	int		NoAcc, NoTried;

	double	Min, Max, Target;
	double	Last;	
	
	double	*RateAcc;
	double	*RateDev;
} AUTOTUNE;

AUTOTUNE*	CreatAutoTune(double InitDev, double Min, double Max);
void		FreeAutoTune(AUTOTUNE *AutoTune);

void		AutoTuneUpDate(AUTOTUNE *AutoTune, RANDSTATES *RS);
void		ReSetAutoTune(AUTOTUNE *AutoTune);
#endif