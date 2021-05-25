#ifndef TIME_SLICES_H
#define TIME_SLICES_H

#include "typedef.h"
#include "genlib.h"

TIME_SLICE*		AddTimeSlice(TIME_SLICES *TSlices, char *Name, double Time, double Scale);
TIME_SLICE*		AllocTimeSlice(char *Name);
TIME_SLICE*		GetTimeSlice(TIME_SLICES *TSlices, char *Name);

TIME_SLICES*	CreateTimeSlices(void);
void			FreeTimeSlices(TIME_SLICES *TSlices);

void			CopyTimeSlices(TIME_SLICES *A, TIME_SLICES *B);

void			PrintTimeSlices(FILE *Str, TIME_SLICES *TSlices);

TIME_SLICES*	CreateRatesTimeSlices(RATES *Rates, TIME_SLICES *TSlices);

void			ApplyTimeSlices(RATES *Rates, TREES *Trees);

void			PrintTimeSliceHeader(FILE *Str, TIME_SLICES *TS);
void			PrintTimeSliceRates(FILE *Str, TIME_SLICES *TS_Opt, TIME_SLICES *TS_Rates);

int				TimeSliceEstTime(TIME_SLICES *TS);
int				TimeSliceEstScale(TIME_SLICES *TS);

void			ChangeTimeSliceTime(RATES *Rates, SCHEDULE *Shed);
void			ChangeTimeSliceScale(RATES *Rates, SCHEDULE *Shed);


#endif
