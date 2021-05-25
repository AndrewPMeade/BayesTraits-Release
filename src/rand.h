#ifndef __RANDOM
#define __RANDOM

#include "typedef.h"


unsigned long	RandomLong(void);
double			RandomReal(void);
long			RandLong(long Low, long High);
unsigned		Poisson(double ExpectedValue);
double			NegExp(double ExpectedValue);
int				SetSeed(void);
void			IntSetSeed(int S);
float			Other(void);
float			box_muller(float m, float s);
double			nrand(void);
unsigned long	GetSeed(void);

/* A good random number genetrater */
//double			GenRand(void);
//int				SmallRand(unsigned int Mag);
//int				IntRand();

void			DirichletRandomVariable(double *alp, double *z, int n);

RANDSTATES*		CreateRandStates(void);
void			FreeRandStates(RANDSTATES* RandS);
double			GenRandState(RANDSTATES* RandS);

#endif
