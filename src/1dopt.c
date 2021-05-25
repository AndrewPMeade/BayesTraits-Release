#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include "typedef.h"
#include "genlib.h"
#include "1dopt.h"
#include "likelihood.h"
#include "RandLib.h"

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

	if(FOpt->EstOU == TRUE)
		FRates->OU = Pram;

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

double	GetConP(OPTIONS *Opt, RATES *Rates)
{
	if(Opt->EstDelta == TRUE)
		return Rates->Delta;

	if(Opt->EstKappa == TRUE)
		return Rates->Kappa;

	if(Opt->EstLambda == TRUE)
		return Rates->Lambda;

	if(Opt->EstOU == TRUE)
		return Rates->OU;

	exit(0);
}

void	SetConP(OPTIONS *Opt, RATES *Rates, double P)
{
	if(Opt->EstDelta == TRUE)
	{
		Rates->Delta = P;
		return;
	}

	if(Opt->EstKappa == TRUE)
	{
		Rates->Kappa = P;
		return;
	}

	if(Opt->EstLambda == TRUE)
	{
		Rates->Lambda = P;
		return;
	}

	if(Opt->EstOU == TRUE)
	{
		Rates->OU = P;
		return;
	}

	exit(0);
}

void	UpDateBestLh(OPTIONS *Opt, RATES *Rates, double *BestLh, double *BestP)
{
	if(Rates->Lh > *BestLh)
	{
		*BestLh = Rates->Lh;
		*BestP = GetConP(Opt, Rates);
	}
}

double	SetDefConP(OPTIONS *Opt)
{
	if(Opt->EstDelta == TRUE)
		return 1.0;
	
	if(Opt->EstKappa == TRUE)
		return 1.0;

	if(Opt->EstLambda == TRUE)
		return 1.0;

	if(Opt->EstOU == TRUE)
		return MIN_OU + MIN_OU;

	return -1;
}

double	SetMinConP(OPTIONS *Opt)
{
	if(Opt->EstDelta == TRUE)
		return MIN_DELTA;
	
	if(Opt->EstKappa == TRUE)
		return MIN_KAPPA;

	if(Opt->EstLambda == TRUE)
		return MIN_LAMBDA;

	if(Opt->EstOU == TRUE)
		return MIN_OU;

	return -1;
}


double	SetMaxConP(OPTIONS *Opt)
{
	if(Opt->EstDelta == TRUE)
		return MAX_DELTA;
	
	if(Opt->EstKappa == TRUE)
		return MAX_KAPPA;

	if(Opt->EstLambda == TRUE)
		return MAX_LAMBDA;

	if(Opt->EstOU == TRUE)
		return MAX_OU;

	return -1;
}

double	MidPoint(double Min, double Max)
{
	double Ret;

	Ret = (Max - Min) * 0.5;
	return Ret + Min;
}

double	SetAveConP(OPTIONS *Opt)
{
	if(Opt->EstDelta == TRUE)
		return MidPoint(MIN_DELTA, MAX_DELTA);
	
	if(Opt->EstKappa == TRUE)
		return MidPoint(MIN_KAPPA, MAX_KAPPA);
		
	if(Opt->EstLambda == TRUE)
		return MidPoint(MIN_LAMBDA, MAX_LAMBDA);
		
	if(Opt->EstOU == TRUE)
		return MidPoint(MIN_OU, MAX_OU);

	return -1;
}

double	RandPoint(RANDSTATES *RS, double Min, double Max)
{
	double Ret;
	
	Ret = (Max - Min) * RandDouble(RS);
	return Ret + Min;
}

double	GetRandRates(OPTIONS *Opt, RATES *Rates)
{

	if(Opt->EstDelta == TRUE)
		return RandPoint(Rates->RS, MIN_DELTA, MAX_DELTA);
	
	if(Opt->EstKappa == TRUE)
		return RandPoint(Rates->RS, MIN_KAPPA, MAX_KAPPA);
		
	if(Opt->EstLambda == TRUE)
		return RandPoint(Rates->RS, MIN_LAMBDA, MAX_LAMBDA);
		
	if(Opt->EstOU == TRUE)
		return RandPoint(Rates->RS, MIN_OU, MAX_OU);

	return -1;
}

void	SetRandConRates(OPTIONS *Opt, RATES *Rates, TREES *Trees, double *BestLh, double *BestP)
{
	*BestP = SetDefConP(Opt);
	SetConP(Opt, Rates, *BestP);
	*BestLh = Likelihood(Rates, Trees, Opt);
	return;
	
	SetConP(Opt, Rates, GetRandRates(Opt, Rates));
	Rates->Lh = Likelihood(Rates, Trees, Opt);
	UpDateBestLh(Opt, Rates, BestLh, BestP);
	
	return;

	SetConP(Opt, Rates, SetMinConP(Opt));
	Rates->Lh = Likelihood(Rates, Trees, Opt);
	UpDateBestLh(Opt, Rates, BestLh, BestP);

	SetConP(Opt, Rates, SetMaxConP(Opt));
	Rates->Lh = Likelihood(Rates, Trees, Opt);
	UpDateBestLh(Opt, Rates, BestLh, BestP);
	
	SetConP(Opt, Rates, SetAveConP(Opt));
	Rates->Lh = Likelihood(Rates, Trees, Opt);
	UpDateBestLh(Opt, Rates, BestLh, BestP); 
	
	
	
/*	*BestP = GetRandRates(Opt, Rates);
	SetConP(Opt, Rates, *BestP);
	*BestLh = Likelihood(Rates, Trees, Opt);
	return;

	SetConP(Opt, Rates, SetMinConP(Opt));
	Rates->Lh = Likelihood(Rates, Trees, Opt);
	UpDateBestLh(Opt, Rates, BestLh, BestP);

	SetConP(Opt, Rates, SetMaxConP(Opt));
	Rates->Lh = Likelihood(Rates, Trees, Opt);
	UpDateBestLh(Opt, Rates, BestLh, BestP);
	
	SetConP(Opt, Rates, SetAveConP(Opt));
	Rates->Lh = Likelihood(Rates, Trees, Opt);
	UpDateBestLh(Opt, Rates, BestLh, BestP); */
}

void	Opt1D(OPTIONS *Opt, RATES *Rates, TREES *Trees)
{
	double	Tol, xmin, Start;
	double	BestLh, BestP;
	int Index;

	Tol = 0.000001;

	Tol = 0.00001;
	FTrees	= Trees;
	FRates	= Rates;
	FOpt	= Opt;


	if(Opt->DataType == DISCRETE)
	{
		Rates->Rates[0] = Opt1DDesc(Opt, Rates, Trees);
		return;
	}
	/*
	BestP = 0;
	for(Index=0;Index<1000;Index++)
	{
		Rates->OU = BestP;
		BestLh = Likelihood(Rates, Trees, Opt);

		printf("%f\t%f\n", BestP, BestLh);fflush(stdout);
		BestP += 0.01;
	}

	exit(0);
	*/

	SetRandConRates(Opt, Rates, Trees, &BestLh, &BestP);
	
	for(Index=0;Index<Opt->MLTries;Index++)
	{
		if(Index > 0)
			Start = GetRandRates(Opt, Rates);
		else
			Start = BestP; 

		SetConP(Opt, Rates, Start);

		if(Opt->EstDelta == TRUE)
		{
			golden(500, 1, MIN_DELTA, Start, MAX_DELTA, Opt1DFunc, Tol, &xmin);
		}

		if(Opt->EstKappa == TRUE)
		{
			golden(500, 1, MIN_KAPPA, Start, MAX_KAPPA, Opt1DFunc, Tol, &xmin);
		}

		if(Opt->EstLambda == TRUE)
		{
			golden(500, 1, MIN_LAMBDA,Start, MAX_LAMBDA , Opt1DFunc, Tol, &xmin);
		}

		if(Opt->EstOU == TRUE)
		{
			golden(500, 1, MIN_OU, Start, MAX_OU, Opt1DFunc, Tol, &xmin);
		}

		SetConP(Opt, Rates, xmin);
		Rates->Lh = Likelihood(Rates, Trees, Opt);

//		printf("ETest\t%d\t%f\t%f\t%f\n", Index, Start, Rates->Lh, xmin);fflush(stdout);

		if(Rates->Lh > BestLh)
		{
			BestLh = Rates->Lh;
			BestP = xmin;
		}
	}

	SetConP(Opt, Rates, BestP);

	Rates->Lh = Likelihood(Rates, Trees, Opt);

//	printf("Best\t%d\t%f\t%f\n", Index, Rates->Lh, BestP);
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
