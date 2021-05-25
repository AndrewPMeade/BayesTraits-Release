
/* The following routines were copied from Yang (1994) */
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

#include "gamma.h"

#define MAINCAT 4

#ifndef max
	#define min(x,y)							((x) < (y) ? (x) : (y))	
#endif

#ifndef max
	#define max(x,y)							((x) > (y) ? (x) : (y))	
#endif

#define square(a)							((a)*(a))
#define POINTGAMMA(prob,alpha,beta) 		PointChi2(prob,2.0*(alpha))/(2.0*(beta))
#define PAI2								6.283185307
#define max2(a,b)	 						((a)>(b)?(a):(b))

static double   PointNormal (double prob);
static double   CdfNormal (double x);
double   LnGamma (double alpha);
double   IncompleteGamma (double x, double alpha, double LnGamma_alpha);
static double   LBinormal (double h1, double h2, double r);
/* static double   LBinormal (double h, double k, double r); */
static double   CdfBinormal (double h, double k, double r);
static double   PointChi2 (double prob, double v);
static double   Vha (double h, double k);
static double   Tha (double h1, double h2, double a1, double a2);

/* ------------------------------------------------------------------------------
|                                                                               |
|  Returns z so That Prob{x<z} = prob where x ~ N(0,1) and                      |
|  (1e-12) < prob < 1-(1e-12).  Returns (-9999) if in error.                    |
|                                                                               |
|  Odeh, R. E. and J. O. Evans.  1974.  The percentage points of the normal     |
|     distribution.  Applied Statistics, 22:96-97 (AS70)                        |
|                                                                               |
|  Newer methods:                                                               |
|                                                                               |
|  Wichura, M. J.  1988.  Algorithm AS 241: The percentage points of the        |
|     normal distribution.  37:477-484.                                         |
|  Beasley, JD & S. G. Springer.  1977.  Algorithm AS 111: The percentage       |
|     points of the normal distribution.  26:118-121.                           |
|                                                                               |
-------------------------------------------------------------------------------*/   
double PointNormal (double prob)
{

	double 		a0 = -0.322232431088, a1 = -1.0, a2 = -0.342242088547, a3 = -0.0204231210245,
 					a4 = -0.453642210148e-4, b0 = 0.0993484626060, b1 = 0.588581570495,
 					b2 = 0.531103462366, b3 = 0.103537752850, b4 = 0.0038560700634,
 					y, z = 0, p = prob, p1;

	p1 = (p<0.5 ? p : 1-p);
	if (p1<1e-20) 
	   return (-9999);
	y = sqrt (log(1/(p1*p1)));   
	z = y + ((((y*a4+a3)*y+a2)*y+a1)*y+a0) / ((((y*b4+b3)*y+b2)*y+b1)*y+b0);
	return (p<0.5 ? -z : z);

}





/* ------------------------------------------------------------------------------
|                                                                               |
|  Hill, I. D.  1973.  The normal integral.  Applied Statistics, 22:424-427.    |
|     (AS66)                                                                    |
|                                                                               |
|  Adapted by Z. Yang, March 1994.  Hill's routine is quite bad, and I          |
|  have not consulted:                                                          |
|                                                                               |
|  Adams, A. G.  1969.  Algorithm 39.  Areas under the normal curve.            |
|     Computer J. 12:197-198.                                                   |
|                                                                               |
-------------------------------------------------------------------------------*/   
double CdfNormal (double x)

{

	int 			invers = 0;
	double 		p, limit = 10.0, t = 1.28, y = x*x/2.0;

	if (x < 0) 
		{  
		invers = 1;  
		x  *= -1.0; 
		}
	if (x > limit)  
		return (invers ? 0 : 1);
	if (x < t)  
		p = 0.5 - x * (0.398942280444 - 0.399903438504 * y /
			(y + 5.75885480458 - 29.8213557808 /
			(y + 2.62433121679 + 48.6959930692 /
			(y + 5.92885724438))));
	else 
		p = 0.398942280385 * exp(-y) /
			(x - 3.8052e-8 + 1.00000615302 /
			(x + 3.98064794e-4 + 1.98615381364 /
			(x - 0.151679116635 + 5.29330324926 /
			(x + 4.8385912808 - 15.1508972451 /
			(x + 0.742380924027 + 30.789933034 /
			(x + 3.99019417011))))));
	return  invers ? p : 1-p;

}





/*-------------------------------------------------------------------------------
|                                                                               |
|  Returns ln(gamma(alpha)) for alpha > 0, accurate to 10 decimal places.       |
|  Stirling's formula is used for the central polynomial part of the procedure. |
|                                                                               |
|  Pike, M. C. and I. D. Hill.  1966.  Algorithm 291: Logarithm of the gamma    |
|     function.  Communications of the Association for Computing                |
|     Machinery, 9:684.                                                         |
|                                                                               |
-------------------------------------------------------------------------------*/   
double LnGamma (double alpha)

{

	double 		x = alpha, f = 0.0, z;

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





/*-------------------------------------------------------------------------------
|                                                                               |
|  Returns the incomplete gamma ratio I(x,alpha) where x is the upper           |
|  limit of the integration and alpha is the shape parameter.  Returns (-1)     |
|  if in error.                                                                 |
|  LnGamma_alpha = ln(Gamma(alpha)), is almost redundant.                      |
|  (1) series expansion     if (alpha>x || x<=1)                                |
|  (2) continued fraction   otherwise                                           |
|                                                                               |
|  RATNEST FORTRAN by                                                           |
|  Bhattacharjee, G. P.  1970.  The incomplete gamma integral.  Applied         |
|     Statistics, 19:285-287 (AS32)                                             |
|                                                                               |
-------------------------------------------------------------------------------*/   
double IncompleteGamma (double x, double alpha, double LnGamma_alpha)

{

	int 			i;
	double 		p = alpha, g = LnGamma_alpha,
					accurate = 1e-8, overflow = 1e30,
					factor, gin = 0.0, rn = 0.0, a = 0.0, b = 0.0, an = 0.0, 
					dif = 0.0, term = 0.0, pn[6];

	if (x == 0.0) 
		return (0.0);
	if (x < 0 || p <= 0) 
		return (-1.0);

	factor = exp(p*log(x)-x-g);   
	if (x>1 && x>=p) 
		goto l30;
	gin = 1.0;  
	term = 1.0;  
	rn = p;
	l20:
		rn++;
		term *= x/rn;   
		gin += term;
		if (term > accurate) 
			goto l20;
		gin *= factor/p;
		goto l50;
	l30:
		a = 1.0-p;   
		b = a+x+1.0;  
		term = 0.0;
		pn[0] = 1.0;  
		pn[1] = x;  
		pn[2] = x+1;  
		pn[3] = x*b;
		gin = pn[2]/pn[3];
	l32:
		a++;  
		b += 2.0;  
		term++;   
		an = a*term;
		for (i=0; i<2; i++) 
			pn[i+4] = b*pn[i+2]-an*pn[i];
		if (pn[5] == 0) 
			goto l35;
		rn = pn[4]/pn[5];   
		dif = fabs(gin-rn);
		if (dif>accurate) 
			goto l34;
		if (dif<=accurate*rn) 
			goto l42;
	l34:
		gin = rn;
	l35:
		for (i=0; i<4; i++) 
			pn[i] = pn[i+2];
		if (fabs(pn[4]) < overflow) 
			goto l32;
		for (i=0; i<4; i++) 
			pn[i] /= overflow;
		goto l32;
	l42:
		gin = 1.0-factor*gin;
	l50:
		return (gin);

}





/*-------------------------------------------------------------------------------
|                                                                               |
|  Returns z so That Prob{x<z} = prob where x is Chi2 distributed with df=v.    |
|  Returns -1 if in error.   0.000002 < prob < 0.999998.                        |
|                                                                               |
|  RATNEST FORTRAN by                                                           |
|  Best, D. J. and D. E. Roberts.  1975.  The percentage points of the          |
|     Chi2 distribution.  Applied Statistics 24:385-388.  (AS91)                |
|                                                                               |
|  Converted into C by Ziheng Yang, Oct. 1993.                                  |
|                                                                               |
-------------------------------------------------------------------------------*/   
double PointChi2 (double prob, double v)

{

	double 		e = 0.5e-6, aa = 0.6931471805, p = prob, g,
					xx, c, ch, a = 0.0, q = 0.0, p1 = 0.0, p2 = 0.0, t = 0.0, 
					x = 0.0, b = 0.0, s1, s2, s3, s4, s5, s6;

	if (p < 0.000002 || p > 0.999998 || v <= 0.0) 
		return (-1.0);
	g = LnGamma (v/2.0);
	xx = v/2.0;   
	c = xx - 1.0;
	if (v >= -1.24*log(p)) 
		goto l1;
	ch = pow((p*xx*exp(g+xx*aa)), 1.0/xx);
	if (ch-e<0) 
		return (ch);
	goto l4;
	l1:
		if (v > 0.32) 
			goto l3;
		ch = 0.4;   
		a = log(1.0-p);
	l2:
		q = ch;  
		p1 = 1.0+ch*(4.67+ch);  
		p2 = ch*(6.73+ch*(6.66+ch));
		t = -0.5+(4.67+2.0*ch)/p1 - (6.73+ch*(13.32+3.0*ch))/p2;
		ch -= (1.0-exp(a+g+0.5*ch+c*aa)*p2/p1)/t;
		if (fabs(q/ch-1.0)-0.01 <= 0.0) 
			goto l4;
		else                       
			goto l2;
	l3: 
		x = PointNormal (p);
		p1 = 0.222222/v;   
		ch = v*pow((x*sqrt(p1)+1.0-p1), 3.0);
		if (ch > 2.2*v+6.0)  
			ch = -2.0*(log(1.0-p)-c*log(0.5*ch)+g);
	l4:
		q = ch;   
		p1 = 0.5*ch;
		if ((t = IncompleteGamma (p1, xx, g)) < 0.0) 
			{
			printf ("\nerr IncompleteGamma");
			return (-1.0);
			}
		p2 = p-t;
		t = p2*exp(xx*aa+g+p1-c*log(ch));   
		b = t/ch;  
		a = 0.5*t-b*c;
		s1 = (210.0+a*(140.0+a*(105.0+a*(84.0+a*(70.0+60.0*a))))) / 420.0;
		s2 = (420.0+a*(735.0+a*(966.0+a*(1141.0+1278.0*a))))/2520.0;
		s3 = (210.0+a*(462.0+a*(707.0+932.0*a)))/2520.0;
		s4 = (252.0+a*(672.0+1182.0*a)+c*(294.0+a*(889.0+1740.0*a)))/5040.0;
		s5 = (84.0+264.0*a+c*(175.0+606.0*a))/2520.0;
		s6 = (120.0+c*(346.0+127.0*c))/5040.0;
		ch += t*(1+0.5*t*s1-b*c*(s1-b*(s2-b*(s3-b*(s4-b*(s5-b*s6))))));
		if (fabs(q/ch-1.0) > e) 
			goto l4;
		return (ch);

}





/*-------------------------------------------------------------------------------
|                                                                               |
|  Discretization of gamma distribution with equal proportions in each          |
|  category.                                                                    |
|                                                                               |
-------------------------------------------------------------------------------*/   
int DiscreteGamma (double *freqK, double *rK, double alfa, double beta, int K, int median)

{

	int 			i;
	double 		gap05 = 1.0/(2.0*K), t, factor = alfa/beta*K, lnga1;

	if (median) 
		{
		
		for (i=0; i<K; i++) 
			rK[i] = POINTGAMMA((i*2.0+1)*gap05, alfa, beta);
		for (i=0,t=0; i<K; i++) 
			t += rK[i];
		for (i=0; i<K; i++)     
			rK[i] *= factor/t;
		}
	else 
		{
		/* The one we use. */
		lnga1 = LnGamma(alfa+1);
		for (i=0; i<K-1; i++) 
			freqK[i] = POINTGAMMA((i+1.0)/K, alfa, beta);
		for (i=0; i<K-1; i++) 
			freqK[i] = IncompleteGamma(freqK[i]*beta, alfa+1, lnga1);
		rK[0] = freqK[0]*factor;
		rK[K-1] = (1-freqK[K-2])*factor;
		for (i=1; i<K-1; i++)  
			rK[i] = (freqK[i]-freqK[i-1])*factor;
		}
	for (i=0; i<K; i++) 
		freqK[i]=1.0/K;
	return (0);

}





/*-------------------------------------------------------------------------------
|                                                                               |
|  Auto-discrete-gamma distribution of rates over sites, K equal-probable       |
|  categories, with the mean for each category used.                            |
|  This routine calculates M[], freqK[] and rK[], using alfa, rho and K.        |
|                                                                               |
-------------------------------------------------------------------------------*/   
int AutodGamma (double **M, double *freqK, double *rK, double *rho1, double alfa, double rho, int K)

{

	int			i, j, i1, i2;
	double		*point=freqK, x, y, large=20, v1;
	
	/* if (fabs(rho)>1-1e-4) error ("rho out of range"); */
	for (i=0; i<K-1; i++) 
		point[i] = PointNormal ((i + 1.0) / K);
	for (i=0; i<K; i++) 
		{
		for (j=0; j<K; j++) 
			{
			x = (i < K-1 ? point[i]:large);
			y = (j < K-1 ? point[j]:large);
			M[i][j] = CdfBinormal (x, y, rho);
			}
		}
	for (i1=0; i1<2*K-1; i1++) 
		{
		for (i2=0; i2<K*K; i2++) 
			{
			i = i2 / K; 
			j = i2 % K;
			if (i+j != 2*(K-1)-i1) 
				continue;
			y = 0;
			if (i > 0) 
				y -= M[i-1][j];
			if (j > 0) 
				y -= M[i][j-1];
			if (i > 0 && j > 0) 
				y += M[i-1][j-1];
			M[i][j] = (M[i][j]+y) * K;
			if (M[i][j] < 0) 
				printf("M[%d][%d] = %12.8f < 0\n", i+1, j+1, M[i][j]);
			}
		}
	DiscreteGamma (freqK, rK, alfa, alfa, K, 0);

	for (i=0, v1=*rho1=0; i<K; i++) 
		{
		v1 += rK[i] * rK[i] * freqK[i];
		for (j=0; j<K; j++)
			*rho1 += freqK[i] * M[i][j] * rK[i] * rK[j];
		}
	v1 -= 1;
	*rho1 = (*rho1-1) / v1;
	
	return (0);
	
}





#if 0
/*-------------------------------------------------------------------------------
|                                                                               |
|  L(h,k,r) = prob(x>h, y>k)                                                    |
|           = V(h,a1) + V(k,a2) + 1 - 0.5*{F(h) + F(k)} - acos(r)/(2*PAI)       |
|                                                                               |
|  where x and y are standard binormal, with r=corr(x,y)                        |
|                                                                               |
|  The accuracy is poor if ((h or k is large) and (r is near to -1)), when      |
|  L(h,k,r) is small.                                                           |
|                                                                               |
-------------------------------------------------------------------------------*/   
double LBinormal(double h, double k, double r)

{

	double 		limit = 5.0, rsmall = 1e-10, t;

	if (h < -limit && k < -limit) 
		return (1.0);
	if (h > limit || k > limit) 
		return (0.0);
	if (h < -limit) 
		return (1.0-CdfNormal(k));
	if (k < -limit) 
		return (1.0-CdfNormal(h));
	if (fabs(r) < rsmall) 
		return ((1.0-CdfNormal(h))*(1.0-CdfNormal(k)));
	if (r > 1.0-rsmall) 
		return (1.0-CdfNormal(max(h,k)));
	if (r < -(1.0-rsmall)) 
		{
		if (h+k >= 0.0) 
			return(0.0);
		else 
			return (1.0-CdfNormal(h)-CdfNormal(k));
		}
	t = Vha(h,(k-r*h)/sqrt(1-r*r)) + Vha(k,(h-r*k)/sqrt(1.0-r*r)) + 
		1.0 - 0.5*(CdfNormal(h)+CdfNormal(k)) - acos(r)/PAI2;
	return max(t,0);

}




/*-------------------------------------------------------------------------------
|                                                                               |
|  Calculate the CDF of a standard bivariate normal distribution                |
|  r = corr(x,y)                                                                |
|                                                                               |
|  Johnson, N. L. and S. Kotz.  1972.  Distributions in statitics: continuous   |
|     multivariate distributions.  John Wiley & sons, New York.  pp. 93-100     |
|                                                                               |
|  F(h,k;r) =  prob(x<h, y<k)                                                   |
|                                                                               |
|           =   INT      INT     f(x,y; r) dy dx                                |
|            (-inf,h)  (-inf,k)                                                 |
|                                                                               |
|           = 0.5*{F(h)+F(k)-d(h,k)} - T(h,a1) - T(k,a2)                        |
|                                                                               |
-------------------------------------------------------------------------------*/   
double CdfBinormal (double h, double k, double r)

{

	double 		limit = 5.0, rsmall = 1e-10, a, dhk, t;

	if (h > limit && k > limit) 
		return (1.0);
	if (h < -limit || k < -limit) 
		return (0.0);
	if (h > limit) 
		return (CdfNormal(k));
	if (k > limit) 
		return (CdfNormal(h));
	if (fabs(r) < rsmall) 
		return (CdfNormal(h)*CdfNormal(k));
	if (r > 1.0-rsmall) 
		return (CdfNormal(min(h,k)));
	if (r < -(1-rsmall)) 
		{
		if (h+k>=0) 
			return (CdfNormal(h)+CdfNormal(k)-1);
		else       
			return 0.0;
		}
	a =  sqrt(1.0-r*r);
	dhk = !(h*k > 0.0 || (h*k == 0.0 && h+k >= 0.0));
	t = 0.5*(CdfNormal(h)+CdfNormal(k)-dhk) - 
		Tha(h,1,k-h*r, h*a) - Tha(k,1,h-k*r, k*a);
	return max(t,0);

}
#else
double CdfBinormal (double h1, double h2, double r)
{
/* F(h1,h2,r) = prob(x<h1, y<h2), where x and y are standard binormal, 
*/
   return (LBinormal(h1,h2,r)+CdfNormal(h1)+CdfNormal(h2)-1);
}    

double LBinormal (double h1, double h2, double r)
{
/* L(h1,h2,r) = prob(x>h1, y>h2), where x and y are standard binormal, 
   with r=corr(x,y),  error < 2e-7.
      Drezner Z., and G.O. Wesolowsky (1990) On the computation of the
      bivariate normal integral.  J. Statist. Comput. Simul. 35:101-107.
*/
   int i;
   double x[]={0.04691008, 0.23076534, 0.5, 0.76923466, 0.95308992};
   double w[]={0.018854042, 0.038088059, 0.0452707394,0.038088059,0.018854042};
   double Lh=0, r1, r2, r3, rr, aa, ab, h3, h5, h6, h7, h12;

   h12=(h1*h1+h2*h2)/2;
   if (fabs(r)>=0.7) {
      r2=1-r*r;   r3=sqrt(r2);
      if (r<0) h2*=-1;
      h3=h1*h2;   h7=exp(-h3/2);
      if (fabs(r)!=1) {
         h6=fabs(h1-h2);   h5=h6*h6/2; h6/=r3; aa=.5-h3/8;  ab=3-2*aa*h5;
         Lh = .13298076*h6*ab*(1-CdfNormal(h6))
            - exp(-h5/r2)*(ab+aa*r2)*0.053051647;
         for (i=0; i<5; i++) {
            r1=r3*x[i];  rr=r1*r1;   r2=sqrt(1-rr);
            Lh-=w[i]*exp(-h5/rr)*(exp(-h3/(1+r2))/r2/h7-1-aa*rr);
         }
      }
      if (r>0) Lh = Lh*r3*h7+(1-CdfNormal(max2(h1,h2)));
      else if (r<0) Lh = (h1<h2?CdfNormal(h2)-CdfNormal(h1):0) - Lh*r3*h7;
   }
   else {
      h3=h1*h2;
      if (r!=0) 
         for (i=0; i<5; i++) {
            r1=r*x[i]; r2=1-r1*r1;
           Lh+=w[i]*exp((r1*h3-h12)/r2)/sqrt(r2);
         }
      Lh=(1-CdfNormal(h1))*(1-CdfNormal(h2))+r*Lh;
   }
   return (Lh);
}    

#endif




double Vha (double h, double k)

{

	double 		tv1 = 1e-35, t, sign;

	if (fabs(h) < tv1) 
		{
		if (h == 0.0) 
			sign = (k >= 0.0 ? 1.0 : -1.0);
		else      
			sign = ((k >= 0.0 && h > 0.0 || k < 0.0 && h < 0.0) ? 1.0 : -1.0);
		t = 0.25*sign;
		}
	else 
		t = atan(k/h)/PAI2;
	return  t - Tha(h, 1, k, h);

}




/*-------------------------------------------------------------------------------
|                                                                               |
|  Calculate Owen's (1956) T(h,a) function, -inf <= h, a <= inf,                |
|  where h = h1/h2, a = a1/a2, from the program of:                             |
|                                                                               |
|  Young, J. C. and C. E. Minder.  1974.  Algorithm AS 76.  An integral         |
|     useful in calculating non-central t and bivariate normal                  |
|     probabilities.  Appl. Statist., 23:455-457.  [Correction: Appl.           |
|     Statist., 28:113 (1979).  Remarks: Appl. Statist. 27:379 (1978),          |
|     28: 113 (1979), 34:100-101 (1985), 38:580-582 (1988)]                     |
|                                                                               |
|  See also:                                                                    |
|                                                                               |
|  Johnson, N. L.  and S. Kotz.  1972.  Distributions in statistics:            |
|     multivariate distributions.  Wiley and Sons.  New York.  pp. 93-100.      |
|                                                                               |
-------------------------------------------------------------------------------*/   
double Tha (double h1, double h2, double a1, double a2)

{

	int 			ng = 5, i;
	double 		U[] = {0.0744372, 0.2166977, 0.3397048, 0.4325317, 0.4869533},
					R[] = {0.1477621, 0.1346334, 0.1095432, 0.0747257, 0.0333357},
					pai2 = 6.283185307, tv1 = 1e-35, tv2 = 15.0, tv3 = 15.0, tv4 = 1e-5,
					a, h, rt, t, x1, x2, r1, r2, s, k, sign = 1.0;

	if (fabs(h2) < tv1) 
		return (0.0);
	h = h1 / h2;
	if (fabs(a2) < tv1) 
		{
		t = CdfNormal(h);
		if (h >= 0.0) 
			t = (1.0 - t) / 2.0;
		else      
			t /= 2.0;
		return (t*(a1 >= 0.0 ? 1.0 : -1.0));
		}
	a = a1 / a2;
	if (a < 0.0) 
		sign = -1.0;  
	a = fabs(a);  
	h = fabs(h);   
	k = h*a;
	if (h > tv2 || a < tv1) 
		return (0.0);
	if (h < tv1) 
		return (atan(a)/pai2*sign);
	if (h < 0.3 && a > 7.0) /* (Boys RJ, 1989) */
		{             
		x1 = exp(-k*k/2.0)/k;
		x2 = (CdfNormal(k)-0.5)*sqrt(pai2);
		t = 0.25 - (x1+x2)/pai2*h + ((1.0+2.0/(k*k))*x1+x2)/(6.0*pai2)*h*h*h;
		return (max(t,0)*sign);
		}
	t = -h*h / 2.0;  
	x2 = a;  
	s = a*a;
	if (log(1.0+s)-t*s >= tv3) 
		{
		x1 = a/2;  
		s /= 4.0;
	for (;;) /* truncation point by Newton iteration */
		{        
		x2 = x1 + (t*s+tv3-log(s+1.0)) / (2.0*x1*(1.0/(s+1.0)-t));
		s = x2*x2;
		if (fabs(x2-x1) < tv4) 
			break;
		x1 = x2;
		}
	}
	for (i=0,rt=0; i<ng; i++) /* Gauss quadrature */
		{          
		r1 = 1.0+s*square(0.5+U[i]);
		r2 = 1.0+s*square(0.5-U[i]);
		rt+= R[i]*(exp(t*r1)/r1 + exp(t*r2)/r2);
		}
	return (max(rt*x2/pai2,0)*sign);

}
/*
int	main(void)
{
	double f[MAINCAT], rK[MAINCAT];
	double alpha;
	int	i;
	alpha = 10;
	

	DiscreteGamma (&f[0], &rK[0], alpha, alpha, MAINCAT, 0);
	
	for(i=0;i<MAINCAT;i++)
		printf("%f\t%f\n", f[i], rK[i]);

	return 0;
}
*/
