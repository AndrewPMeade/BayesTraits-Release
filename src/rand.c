#include <math.h>
#include <time.h>
#include <stdio.h>
#include <stdlib.h>

#include "genlib.h"
#include "rand.h"
#include "typedef.h"

#define		N		25
#define		MAGIC	7
#define		MAXRAND	0x7FFF


#ifndef TRUE
	#define		TRUE	0
#endif

#ifndef FALSE
	#define		FALSE	1
#endif


static const long A = 48271L;
static const long M = 2147483647L;
static const long Q = 2147483647L / 48271L;
static const long R = 2147483647L % 48271L;


unsigned long	Seed;
unsigned long	KeepSeed;

double			GenRand(void);
unsigned long	RandomLong(void);
double			RandomReal(void);
long			RandLong(long Low, long High);
unsigned		Poisson(double ExpectedValue);
double			NegExp(double ExpectedValue);
int				SetSeed(void);
float			Other(void);
float			box_muller(float m, float s);

/* initial 25 seeds, change as you wish */
static unsigned long x[N] =
{ 
	0x95f24dab, 0x0b685215, 0xe76ccae7, 0xaf3ec239, 0x715fad23,
	0x24a590ad, 0x69e4b5ef, 0xbf456141, 0x96bc1b7b, 0xa7bdf825,
	0xc1de75b7, 0x8858a9c9, 0x2da87693, 0xb657f9dd, 0xffdc8a9f,
	0x8121da71, 0x8b823ecb, 0x885d05f5, 0x4e20cd47, 0x5a9ad5d9,
	0x512c0c03, 0xea857ccd, 0x4cc1d30f, 0x8891a8a1, 0xa6b7aadb
};

/* this is magic vector `a', don't change */
static unsigned long mag01[2]=
{ 
	0x0, 0x8ebfd028 
};


/* Set both the random seends and print it out for futer use */
int SetSeed(void)
{
	int	S;
	int	i;
	
	S =  (unsigned)time(NULL);

	#ifndef JNIRUN
		printf("Rand Seed\t%d\n", S);
	#endif
	
	srand(S);
	KeepSeed = Seed = S;

	for(i=0;i<N;i++)
		x[i] = RandomLong();

	return S;
}

void	IntSetSeed(int S)
{
	int i;

	srand(S);
	KeepSeed = Seed = S;

	for(i=0;i<N;i++)
		x[i] = RandomLong();
}

unsigned long GetSeed(void)
{
	return KeepSeed;
}

unsigned long RandomLong(void)
{
	long TmpSeed = A * ( Seed % Q ) - R * ( Seed / Q );
	if( TmpSeed >= 0 )
		Seed = TmpSeed;
	else
		Seed = TmpSeed + M;

	return Seed;
}

double RandomReal(void)
{
    return RandomLong( ) / 2147483647.0;
}

/* Not best algorithm for linear congruential generators */
long RandLong( long Low, long High )
{
    return RandomLong( ) % ( High - Low + 1 ) + Low;
}

unsigned Poisson( double ExpectedValue )
{
    double Limit = exp( -ExpectedValue );
    double Product = RandomReal( );

    int Count;
    for( Count = 0; Product > Limit; Count++ )
        Product *= RandomReal( );

    return Count;
}

double NegExp( double ExpectedValue )
{
    return - ExpectedValue * log(RandomReal());
}

float	Other(void)
{
         float x1, x2, w, y1, y2;
 
         do {
                 x1 = (float)2.0 * (float)RandomReal() - (float)1.0;
                 x2 = (float)2.0 * (float)RandomReal() - (float)1.0;
                 w = x1 * x1 + x2 * x2;
         } while ( w >= 1.0 );

         w = (float)sqrt( (-2.0 * log( w ) ) / w );
         y1 = x1 * w;
         y2 = x2 * w;
		 return y1;
}

float box_muller(float m, float s)	/* normal random variate generator */
{				        /* mean m, standard deviation s */
	float x1, x2, w, y1;
	static float y2;
	static int use_last = 0;

	if (use_last)		        /* use value from previous call */
	{
		y1 = y2;
		use_last = 0;
	}
	else
	{
		do {
			x1 = (float)2.0 * (float)RandomReal() - (float)1.0;
			x2 = (float)2.0 * (float)RandomReal() - (float)1.0;
			w = x1 * x1 + x2 * x2;
		} while ( w >= 1.0 );

		w = (float)sqrt( (-2.0 * log( w ) ) / w );
		y1 = x1 * w;
		y2 = x2 * w;
		use_last = 1;
	}

	return( m + y1 * s );
}


/* A C-program for TT800 : July 8th 1996 Version */
/* by M. Matsumoto, email: matumoto@math.keio.ac.jp */
/* genrand() generate one pseudorandom number with double precision */
/* which is uniformly distributed on [0,1]-interval */
/* for each call.  One may choose any initial 25 seeds */
/* except all zeros. */

/* See: ACM Transactions on Modelling and Computer Simulation, */
/* Vol. 4, No. 3, 1994, pages 254-266. */
#ifdef adkasdaslj

double GenRand(void)
{
	unsigned long y;
	static int k = 0;
	int kk;

//	return (double)(rand() % 10000) / 10000;
	/* generate N words at one time */
	if (k==N) 
	{ 
		for (kk=0;kk<N-MAGIC;kk++) 
			x[kk] = x[kk+MAGIC] ^ (x[kk] >> 1) ^ mag01[x[kk] % 2];

		for (; kk<N;kk++) 
			x[kk] = x[kk+(MAGIC-N)] ^ (x[kk] >> 1) ^ mag01[x[kk] % 2];

		k=0;
    }

	y = x[k];
	y ^= (y << 7) & 0x2b5b2500; /* s and b, magic vectors */
	y ^= (y << 15) & 0xdb8b0000; /* t and c, magic vectors */
	y &= 0xffffffff; /* you may delete this line if word size = 32 */

	y ^= (y >> 16); /* added to the 1994 version */
	k++;
	return( (double) y / (unsigned long) 0xffffffff);
}
#endif
/*
int		OneInTen(void)
{
	unsigned int	r;

	r = (int)(GenRand() * 10000) % 10000;

	if(r>=9000)	return 9;
	if(r>=8000)	return 8;
	if(r>=7000)	return 7;
	if(r>=6000)	return 6;
	if(r>=5000)	return 5;
	if(r>=4000)	return 4;
	if(r>=3000)	return 3;
	if(r>=2000)	return 2;
	if(r>=1000)	return 1;
	return 0;
}
*/
/*
int		SmallRand(unsigned int Mag)
{
	while(Mag > 0)
	{
		Mag--;
		if(OneInTen()!=0)
			return FALSE;
		if(Mag == 1)
			if(OneInTen() == 0) return TRUE;
			else return FALSE;
	}

	printf("Eror\n");
	exit(1);
	return -1;
}
*/

/*
int	IntRand()
{
	double Ret;

	Ret = GenRand() * MAXRAND;
	
	return (int)Ret;
}
*/

/*
double	LnRGamma (double alpha)
{

	double x = alpha, f = 0.0, z;

	if (x < 7)
		{
		f = 1.0;
		z = x-1.0;
		while (++z < 7.0)
			f *= z;
		x = z;
		f = -log(f);
		}
	z = 1.0/(x*x);
	return  f + (x-0.5)*log(x) - x + 0.918938533204673 +
		(((-0.000595238095238*z+0.000793650793651)*z-0.002777777777778)*z +
		0.083333333333333)/x;

}
*/


double RndGamma1 (double  s)
{

	double	r, x=0.0, small=1e-37, w;
	static double a, p, uf, ss=10.0, d;

	if (s != ss)
		{
		a  = 1.0 - s;
		p  = a / (a+s*exp(-a));
		uf = p * pow(small/a,s);
		d  = a * log(a);
		ss = s;
		}
	for (;;)
		{
		r = RandomReal();
		if (r > p)
			x = a - log((1.0 - r) / (1.0 - p)), w = a * log(x) - d;
		else if (r>uf)
			x = a * pow(r/p,1/s), w = x;
		else
			return (0.0);
		r = RandomReal();
		if (1.0-r <= w && r > 0.0)
		if (r*(w+1.0) >= 1.0 || -log(r) <= w)
			continue;
		break;
		}
	return (x);

}


double	RndGamma2 (double s)
{
	double	r ,d, f, g, x;
	static double b, h, ss=0;

	if (s!=ss)
		{
		b  = s-1.0;
		h  = sqrt(3.0*s-0.75);
		ss = s;
		}
	for (;;)
		{
		r = RandomReal();
		g = r-r*r;
		f = (r-0.5)*h/sqrt(g);
		x = b+f;
		if (x <= 0.0)
			continue;
		r = RandomReal();
		d = 64*r*r*g*g*g;
		if (d*x < x-2.0*f*f || log(d) < 2*(b*log(x/b)-f))
			break;
		}
	return (x);
}

double RndGamma (double s)
{

	double r=0.0;

	if (s <= 0.0)
		puts ("jgl gamma..");
	else if (s < 1.0)
		r = RndGamma1 (s);
	else if (s > 1.0)
		r = RndGamma2 (s);
	else
		r =- log(RandomReal());
	return (r);

}

void DirichletRandomVariable (double *alp, double *z, int n)
{

	int		i;
	double	sum;

	sum = 0.0;
	for(i=0; i<n; i++)
		{
		z[i] = RndGamma (alp[i]) / 1.0;
		sum += z[i];
		}
	for(i=0; i<n; i++)
		z[i] /= sum;

}			 
							 
RANDSTATES*	CreateRandStates(void)
{
	RANDSTATES*		Ret;
	int				Index;

	Ret = (RANDSTATES*)malloc(sizeof(RANDSTATES));
	if(Ret == NULL)
		MallocErr();
	Ret->States = (unsigned long*) malloc(sizeof(long) * N);
	if(Ret->States == NULL)
		MallocErr();

	for(Index=0;Index<N;Index++)
		Ret->States[Index] = RandomLong();

	Ret->K = 0;
	return Ret;
}

/* A C-program for TT800 : July 8th 1996 Version */
/* by M. Matsumoto, email: matumoto@math.keio.ac.jp */
/* genrand() generate one pseudorandom number with double precision */
/* which is uniformly distributed on [0,1]-interval */
/* for each call.  One may choose any initial 25 seeds */
/* except all zeros. */

/* See: ACM Transactions on Modelling and Computer Simulation, */
/* Vol. 4, No. 3, 1994, pages 254-266. */

double GenRandState(RANDSTATES*	RandS)
{
	unsigned long y;
	int kk;

	
	/* generate N words at one time */
	if (RandS->K == N) 
	{ 
		for (kk=0;kk<N-MAGIC;kk++) 
			RandS->States[kk] = RandS->States[kk+MAGIC] ^ (RandS->States[kk] >> 1) ^ mag01[RandS->States[kk] % 2];

		for (; kk<N;kk++) 
			RandS->States[kk] = RandS->States[kk+(MAGIC-N)] ^ (RandS->States[kk] >> 1) ^ mag01[RandS->States[kk] % 2];

		RandS->K = 0;
    }

	y = RandS->States[RandS->K];
	y ^= (y << 7) & 0x2b5b2500;		/* s and b, magic vectors */
	y ^= (y << 15) & 0xdb8b0000;	/* t and c, magic vectors */
	y &= 0xffffffff;				/* you may delete this line if word size = 32 */

	y ^= (y >> 16);					/* added to the 1994 version */
	RandS->K++;
	return( (double) y / (unsigned long) 0xffffffff);
}

void			FreeRandStates(RANDSTATES* RandS)
{
	free(RandS->States);
	free(RandS);
}

double nrand(RANDSTATES* RandS)
{
 /* gives a distribution with mean 0 and std 1. 

  To change the mean to M, simply add M to whatever
  is returned

  To change the std to S, simply multiply whatever is returned 
  by S. Do the mult first.

  Eg: this returns Z (the thing in the return line)

  for mean M and std S, instead return M + S*Z
*/
  double a, b;
  double pi = 3.14159265358979323846;
 
  a = GenRandState(RandS);
  b = GenRandState(RandS);

 return(  sqrt(-2.0*log(a)) * cos(2*pi*b));
}