#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include "typedef.h"
#include "genlib.h"
#include "1dopt.h"
#include "likelihood.h"
#include "rand.h"

#define R 0.61803399 
#define C (1.0-R)
#define SHFT2(a,b,c) (a)=(b);(b)=(c);
#define SHFT3(a,b,c,d) (a)=(b);(b)=(c);(c)=(d);

#define	GOLD 1.6618034
#define GLIMIT 100.00
#define TINY 1.0e-20

static double maxarg1,maxarg2;

#define FMAX(a,b) (maxarg1=(a),maxarg2=(b),(maxarg1) > (maxarg2) ?\
        (maxarg1) : (maxarg2))



#define SIGN(a,b) ((b) >= 0.0 ? fabs(a) : -fabs(a))

RATES*		FRates;
TREES*		FTrees;
OPTIONS*	FOpt;


double golden (int maxiter, int count, double ax, double bx, double cx, double (*f)(double), double tol, double *xmin)
{
	double f1,f2,x0,x1,x2,x3;

	if (maxiter<=0)
	{
		*xmin = bx;
		count = 0;
		return 0;
	}
	else if (maxiter==1)
	{
		*xmin = bx;
		count = 1;
		return (*f)(bx);
	}

	x0 = ax; 
	x3 = cx;
	if (fabs(cx-bx) > fabs(bx-ax)) 
	{ 
		x1 = bx;
		x2 = bx+C*(cx-bx); 
	} 
	else 
	{
		x2 = bx;
		x1 = bx-C*(bx-ax);
	}
	f1 = (*f)(x1); 
				
	f2 = (*f)(x2);
	count = 2;
	while (count<maxiter && fabs(x3-x0)>tol*(fabs(x1)+fabs(x2))) 
	{
		if (f2 < f1) 
		{ 
			SHFT3(x0,x1,x2,R*x1+C*x3) 
			SHFT2(f1,f2,(*f)(x2)) 
		} 
		else 
		{
			SHFT3(x3,x2,x1,R*x2+C*x0)
			SHFT2(f2,f1,(*f)(x1)) 
		}
		count++;
	} 
	if (f1 < f2) 
	{ 
		*xmin=x1;
		return f1;
	} 
	else 
	{
		*xmin=x2;
		return f2;
	}
}

double	Opt1DFunc(double Pram)
{
	double	Ret;

	if(FOpt->EstDelta == TRUE)
		FRates->Delta = Pram;

	if(FOpt->EstKappa == TRUE)
		FRates->Kappa = Pram;

	if(FOpt->EstLambda == TRUE)
		FRates->Lambda = Pram;

	if(FOpt->DataType == DISCRETE)
		FRates->Rates[0] = Pram;

	Ret = Likelihood(FRates, FTrees, FOpt);

	return -Ret;
}

void	mnbrak(double *ax, double *bx, double *cx, double *fa, double *fb, double *fc, double (*func)(double))
{
	double	ulim, u, r, q, fu, dum;

	*fa = (*func)(*ax);
	*fb = (*func)(*bx);

	if(*fb > *fa)
	{
		SHFT3(dum, *ax, *bx, dum);
		SHFT3(dum, *fb, *fa, dum);
	}

	*cx = (*bx) + GOLD * (*bx - *ax);

	*fc = (*func)(*cx);

	while(*fb > *fc)
	{
		r = (*bx - *ax) * (*fb - *fc);
		q = (*bx - *cx) * (*fb - *fa);
		u = (*bx) - ((*bx - *cx)*q-(*bx-*ax)*r) / (2.0 * SIGN(FMAX(fabs(q-r), TINY), q-r));
		ulim = (*bx) + GLIMIT * (*cx - *bx);

		if((*bx-u) * (u-*cx) > 0.0)
		{
			fu = (*func)(u);
			if(fu < *fc)
			{
				*ax = (*bx);
				*bx = u;
				*fa = (*fb);
				*fb = fu;
				return;
			}
			else if (fu > *fb)
			{
				*cx = u;
				*fc = fu;
				return;
			}
			u = (*cx) + GOLD * (*cx - *bx);
			fu = (*func)(u);
		}
		else if ((*cx - u) * (u-ulim) > 0.0)
		{
			fu = (*func)(u);
			if(fu < *fc)
			{
				SHFT3(*bx, *cx, u, *cx + GOLD*(*cx - *bx))
				SHFT3(*fb, *fc, fu, (*func)(u))
			}
		}
		else if((u-ulim) * (ulim-*cx) >= 0.0)
		{
			u = ulim;
			fu = (*func)(u);
		} 
		else
		{
			u = (*cx) + GOLD * (*cx-*bx);
			fu = (*func)(u);
		}
		SHFT3(*ax, *bx, *cx, u)
		SHFT3(*fa, *fb, *fc, fu)
	}
}
double	Opt1DDesc(OPTIONS *Opt, RATES *Rates, TREES *Trees)
{
	double	BestLh;
	double	BestP;
	double	Current;
	double	xmin;


	BestLh = golden(500, 1, 0, 0.05, 0.1, Opt1DFunc, 0.00001, &xmin);
	BestP = xmin;

	Current = golden(500, 1, 0, 0.5, 1, Opt1DFunc, 0.00001, &xmin);
	if((Current != -ERRLH) && (BestLh > Current))
	{
		BestP = xmin;
		BestLh= Current;
	}

	Current = golden(500, 1, 0, 5, 10, Opt1DFunc, 0.00001, &xmin);
	if((Current != -ERRLH) && (BestLh > Current))
	{
		BestP = xmin;
		BestLh= Current;
	}

	Current = golden(500, 1, 0, 50, 100, Opt1DFunc, 0.00001, &xmin);
	if((Current != -ERRLH) && (BestLh > Current))
	{
		BestP = xmin;
		BestLh= Current;
	}

	Current = golden(500, 1, 0, 500, 1000, Opt1DFunc, 0.00001, &xmin);
	if((Current != -ERRLH) && (BestLh > Current))
	{
		BestP = xmin;
		BestLh= Current;
	}

	return BestP;
}

void	Opt1D(OPTIONS *Opt, RATES *Rates, TREES *Trees)
{
	double	xmin;
	FTrees	= Trees;
	FRates	= Rates;
	FOpt	= Opt;

	if(Opt->EstDelta == TRUE)
	{
		golden(500, 1, 0, 2, 3, Opt1DFunc, 0.000001, &xmin);
		Rates->Delta = xmin;
	}

	if(Opt->EstKappa == TRUE)
	{
		golden(500, 1, 0, 1.5, 3, Opt1DFunc, 0.00001, &xmin);
		Rates->Kappa = xmin;
	}

	if(Opt->EstLambda == TRUE)
	{
		golden(500, 1, 0, 0.5, 1 , Opt1DFunc, 0.00001, &xmin);
		Rates->Lambda = xmin;
	}

	if(Opt->DataType == DISCRETE)
	{
	/*	
		golden(500, 1, 0, 500, 1000, Opt1DFunc, 0.00001, &xmin);
		Rates->Rates[0] = xmin; 
	*/
		Rates->Rates[0] = Opt1DDesc(Opt, Rates, Trees);
	}
}


/* Function with brackiting */
/*
void	Opt1D(OPTIONS *Opt, RATES *Rates, TREES *Trees)
{
	double	xmin;
	double	ax,bx,cx;
	double	fa,fb,fc;

	FTrees	= Trees;
	FRates	= Rates;
	FOpt	= Opt;

	ax = 0;

	if(Opt->EstDelta == TRUE)
	{
		ax = TINY;
		bx = 8;

		mnbrak(&ax, &bx, &cx, &fa, &fb, &fc, Opt1DFunc);
		golden(500, 1, ax, bx, cx, Opt1DFunc, 0.001, &xmin);
		Rates->Delta = xmin;
	}

	if(Opt->EstKappa == TRUE)
	{
		ax = TINY;
		bx = 3;
		mnbrak(&ax, &bx, &cx, &fa, &fb, &fc, Opt1DFunc);
		golden(500, 1, ax, bx, cx, Opt1DFunc, 0.001, &xmin);
		Rates->Kappa = xmin;
	}

	if(Opt->EstLambda == TRUE)
	{
		ax = TINY;
		bx = 1;
		mnbrak(&ax, &bx, &cx, &fa, &fb, &fc, Opt1DFunc);

		golden(500, 1, ax, bx, cx, Opt1DFunc, 0.001, &xmin);
		Rates->Lambda = xmin;
	}

	if(Opt->DataType == DISCRETE)
	{
		ax = TINY;
		bx = 50;
		cx =0;

		mnbrak(&ax, &bx, &cx, &fa, &fb, &fc, Opt1DFunc);

		golden(500, 1, ax, bx, cx, Opt1DFunc, 0.001, &xmin);

		Rates->Rates[0] = xmin;
	}
}
*/