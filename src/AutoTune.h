#ifndef AUTOTUNE_H
#define AUTOTUNE_H

#include "TypeDef.h"

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
	double	MaxDev;	
	double	*RateAcc;
	double	*RateDev;

	char	*Name;
		
} AUTOTUNE;

AUTOTUNE*	CreatAutoTune(char *Name, double InitDev, double Min, double Max);
void		FreeAutoTune(AUTOTUNE *AutoTune);

void		SetMaxDev(AUTOTUNE *AutoTune, double MaxDev);

void		AutoTuneUpDate(AUTOTUNE *AutoTune, RANDSTATES *RS);

double		AutoTuneCalcAcc(AUTOTUNE *AT);

void		ReSetAutoTune(AUTOTUNE *AutoTune);

#endif