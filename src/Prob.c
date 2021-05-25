#include <math.h>
#include <stdio.h>
#include <stdlib.h>


#include "typedef.h"
#include "Prob.h"

#ifndef M_PI
	#define M_PI       3.14159265358979323846
#endif

extern double gamma(double a);

//	extern double erf(double a);

extern double beta(double a, double b);
extern double igamc ( double, double );
extern double igam (double a, double x);
extern double incbet(double aa, double bb, double xx );

/*
double erf(double x)
{
	double t, y;
	double a1 =  0.254829592;
	double a2 = -0.284496736;
	double a3 =  1.421413741;
	double a4 = -1.453152027;
	double a5 =  1.061405429;
	double p  =  0.3275911;

    // Save the sign of x
    int sign = 1;
    if (x < 0)
        sign = -1;
    x = fabs(x);

    // A&S formula 7.1.26
    t = 1.0/(1.0 + p*x);
    y = 1.0 - (((((a5*t + a4)*t) + a3)*t + a2)*t + a1)*t*exp(-x*x);

    return sign*y;
}
*/

double	PDFNorm(double x, double Mean, double Var)
{
	double Ret, T1, T2;

	Ret = sqrt(Var) * sqrt(2*M_PI);
	Ret = 1.0 / Ret;

	T1 = (x - Mean) * (x - Mean);
	T2 = 2.0 * Var;

	T1 = T1 / T2;
	T2 = exp(-T1);

	Ret = Ret * T2;
	return Ret;
}

double		PDFExp(double X, double Mean)
{
	double Ret;

	Ret = Mean * exp(-Mean * X);

	return Ret;
}
/*
double		PDFGamma(double X, double Alpha, double Beta)
{
	double Ret;

	Ret = pow(Beta, Alpha) / gamma(Alpha);
	Ret = Ret * pow(X, Alpha-1.0) * exp(-Beta * X);

	return Ret;
}
*/

double		PDFGamma(double X, double Shape, double Scale)
{
	double Ret;

	Ret = gamma(Shape) * pow(Scale, Shape);

	Ret = 1.0 / Ret;

	Ret = Ret * pow(X, Shape-1) * exp(-(X/Scale));

	return Ret;
}

double		PDFBeta(double X, double Alpha, double Beta)
{
	double Ret;

	Ret = pow(X, Alpha - 1.0);

	Ret = Ret * pow(1.0 - X, Beta - 1.0);

	Ret = Ret / beta(Alpha, Beta);

	return Ret;
}

double		PDFChi(double X, double Alpha, double Beta)
{
//	double Ret;

	printf("hs not been implmented.\n");
	exit(0);

	return 1.0;
}

double		PDFInvGamma(double X, double Alpha, double Beta)
{
	double Ret;

	Ret = pow(Beta, Alpha) / gamma(Alpha); 

	Ret = Ret * pow(X, -Alpha - 1.0);

	Ret = Ret * exp(-Beta / X);

	return Ret;
}

double		CDFNorm(double X, double Mean, double Var)
{
	double Ret;
	
	Ret = 0.5;

	Ret = (X - Mean) / (sqrt(Var) * sqrt(2.0));
	Ret = 0.5 * (1.0 + erf(Ret));
	
	return Ret;
}

double		CDFExp(double X, double Alpha)
{
	double Ret;

	Ret = 1.0 - exp(-Alpha * X);
	
	return Ret;
}

double		CDFGamma(double X, double Shape, double Scale)
{
	return igam(Shape, X/Scale);
}

double		CDFBeta(double X, double Alpha, double Beta)
{
	return incbet(Alpha, Beta, X);
}

double		CDFChi(double X, double Alpha, double Beta);

double CDFInvGamma(double X, double Alpha, double Beta)
{
	return igamc(Alpha, Beta / X);
}

