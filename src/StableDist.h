#ifndef STABLEDIST_H
#define STABLEDIST_H

#define LINEAR_INTERPOLATION


typedef struct
{
	double	Alpha;
	double	Scale;
	double	ScaledAlpha;

	double	**P;

} STABLEDIST;

STABLEDIST*	CreatStableDist(void);
void		FreeStableDist(STABLEDIST* SD);

double		StableDistPDF(STABLEDIST* SD, double x);
double		StableDistTPDF(double Scale, double x, double t);

double		CaclNormalLogLh(double X, double Scale, double T);

#endif