#ifndef SCHEDULE_HEADDER
#define SCHEDULE_HEADDER

#include "TypeDef.h"

SCHEDULE*	CreatSchedule(OPTIONS *Opt, gsl_rng *RNG, TREES *Trees);
void		FreeeSchedule(SCHEDULE* Shed);

void		UpDateShedAcc(int Acc, SCHEDULE* Shed);

void		UpDateSchedule(OPTIONS *Opt, SCHEDULE* Shed, gsl_rng *RNG);

double		GetAccRate(int Op, SCHEDULE* Shed);

int		FindNoOfAutoCalibRates(OPTIONS *Opt, TREES *Trees);
void		PrintShed(OPTIONS* Opt, SCHEDULE* Shed, FILE* Str);

void		SetCustomShed(SCHEDULE* Shed);

void		SetShedOpFreq(SCHEDULE*	Shed, int No, double Val);

CUSTOM_SCHEDULE*	AllocCustomSchedule(void);
void				FreeCustomSchedule(CUSTOM_SCHEDULE*	CShed);
void				PrintCustomSchedule(FILE *Str, int NoCShed, CUSTOM_SCHEDULE **ShedList);

void				SetCustomSchedule(OPTIONS *Opt, FILE *ShedFile, size_t Itters, SCHEDULE* Shed);

void				SetRJLockedModel(SCHEDULE* Shed);

#endif
