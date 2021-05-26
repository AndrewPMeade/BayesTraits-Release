#ifndef SCHEDULE_HEADDER
#define SCHEDULE_HEADDER

#include "typedef.h"

SCHEDULE*	CreatSchedule(OPTIONS *Opt, RANDSTATES *RS);
void		FreeeSchedule(SCHEDULE* Shed);

void		UpDateShedAcc(int Acc, SCHEDULE* Shed);

void		UpDateSchedule(OPTIONS *Opt, SCHEDULE* Shed, RANDSTATES *RS);

int			UseAutoTune(OPTIONS *Opt);

double		GetAccRate(int Op, SCHEDULE* Shed);

int			FindNoOfAutoCalibRates(OPTIONS *Opt);

#endif
