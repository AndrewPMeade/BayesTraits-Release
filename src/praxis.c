#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include "praxis.h"
#include "minfit.h"
#include "genlib.h"
#include "RandLib.h"

void	BlankPraxisState(PRAXSTATE*	S)
{
	S->i = S->j = S->k = S->k2 = S->nl = S->nf = S->kl = S->kt = 0;

	S->s= S->sl= S->dn= S->dmin= S->fx= S->f1= S->lds= S->ldt= S->sf= S->df=
	S->qf1= S->qd0= S->qd1= S->qa= S->qb= S->qc= S->m2= S->m4= S->small= S->vsmall= 
	S->large= S->vlarge= S->ldfac= S->t2= 0.0;


	for(S->i=0;S->i<PNSIZE;S->i++)
		S->d[S->i] = S->y[S->i] = S->z[S->i] = S->q0[S->i] = S->q1[S->i] = 0.0;
	
	for(S->i=0;S->i<PNSIZE;S->i++)
		for(S->j=0;S->j<PNSIZE;S->j++)
			S->v[S->i][S->j]  = 0.0;


	S->tol	= SQREPSILON,
	S->scbd	= 1.0;
	S->step	= 1.0;
	S->ktm	= 1;
	S->prin	= 2;

	S->NoLhCalls = 0;
}

PRAXSTATE*	AllocPraxState(void)
{
	PRAXSTATE*	Ret;
	
	Ret = (PRAXSTATE*)malloc(sizeof(PRAXSTATE));
	if(Ret == NULL)
		MallocErr();

	BlankPraxisState(Ret);

	return Ret;
}

void	FreePracxStates(PRAXSTATE* PS)
{
	free(PS);
}

/* --------------------------------------------------------------------------- */

/*
void	BlankPraxisGlobal(void)
{
	i= j= k= k2= nl= nf= kl= kt= 0;
	s= sl= dn= dmin= fx= f1= lds= ldt= sf= df= qf1= qd0= qd1= qa= qb= qc= m2= m4= small= vsmall= large= vlarge= ldfac= t2= 0.0;


	for(i=0;i<N;i++)
		d[i] = y[i] = z[i] = q0[i] = q1[i] = 0.0;
	
	for(i=0;i<N;i++)
		for(j=0;j<N;j++)
			v[i][j]  = 0.0;
}
*/

/* --------------------------------------------------------------------------- */
void sort(PRAXSTATE* S)
/* d and v in descending order */
{
 int k, i, j;
 double s;

 for (i=0; i<S->n-1; i++) {
     k = i; s = S->d[i];
     for (j=i+1; j<S->n; j++) {
         if (S->d[j] > s) {
            k = j;
            s = S->d[j];
         }
     }
     if (k > i) {
        S->d[k] = S->d[i];
        S->d[i] = s;
        for (j=0; j<S->n; j++) {
            s = S->v[j][i];
            S->v[j][i] = S->v[j][k];
            S->v[j][k] = s;
        }
     }
 }
}

/* --------------------------------------------------------------------------- */
void print(PRAXSTATE* S)
/* print a line of traces */
{
 printf("\n");
 printf("... chi square reduced to ... %20.10e\n", S->fx);
 printf("... after %u function calls ...\n", S->nf);
 printf("... including %u linear searches ...\n", S->nl);
 vecprint("... current values of x ...", S->x, S->n);
}

/* --------------------------------------------------------------------------- */
void matprint(char *s,double v[PNSIZE][PNSIZE], int n)
{
 int k, i;

 printf("%s\n", s);
 for (k=0; k<n; k++) {
     for (i=0; i<n; i++) {
         printf("%20.10e ", v[k][i]);
     }
     printf("\n");
 }
}

/* --------------------------------------------------------------------------- */
void vecprint(char *s, double x[PNSIZE], int n)
{
 int i;

 printf("%s\n", s);
 for (i=0; i<n; i++)
     printf("%20.10e ", x[i]);
 printf("\n");
}

double flin(PRAXSTATE* S, double l, int j)
{
 int i;

#ifndef MSDOS
 double tflin[PNSIZE];
#endif

 if (j != -1) {		/* linear search */
    for (i=0; i<S->n; i++)
        tflin[i] = S->x[i] + l *S->v[i][j];
 }
 else {			/* search along parabolic space curve */
      S->qa = l*(l-S->qd1)/(S->qd0*(S->qd0+S->qd1));
      S->qb = (l+S->qd0)*(S->qd1-l)/(S->qd0*S->qd1);
      S->qc = l*(l+S->qd0)/(S->qd1*(S->qd0+S->qd1));
      for (i=0; i<S->n; i++)
          tflin[i] = S->qa*S->q0[i]+S->qb*S->x[i]+S->qc*S->q1[i];
 }
 S->nf++;

	/* Possable Err woth a check */

//	return (*fun)(tflin);
	return S->fun(S, tflin);
}

/* --------------------------------------------------------------------------- */
void min1(PRAXSTATE* S, int j, int nits, double *d2, double *x1, double f1, int fk)
{
 int k, i, dz;
 double x2, xm, f0, f2, fm, d1, t2, s, sf1, sx1;

 sf1 = f1; sx1 = *x1;
 k = 0; xm = 0.0; fm = f0 = S->fx; dz = *d2 < S->macheps;
 /* find step size*/
 s = 0;
 for (i=0; i<S->n; i++)
     s += S->x[i]*S->x[i];
 s = sqrt(s);
 if (dz)
    t2 = S->m4*sqrt(fabs(S->fx)/S->dmin + s*S->ldt) + S->m2*S->ldt;
 else
    t2 = S->m4*sqrt(fabs(S->fx)/(*d2) + s*S->ldt) + S->m2*S->ldt;
 s = s*S->m4 + S->t;
 if (dz && t2 > s) t2 = s;
 if (t2 < S->small) t2 = S->small;
 if (t2 > 0.01*S->h) t2 = 0.01 * S->h;
 if (fk && f1 <= fm) {
    xm = *x1;
    fm = f1;
 }
 if (!fk || fabs(*x1) < t2) {
    *x1 = (*x1 > 0 ? t2 : -t2);
     f1 = flin(S, *x1, j);
 }
 if (f1 <= fm) {
    xm = *x1;
    fm = f1;
 }
next:
 if (dz) {
    x2 = (f0 < f1 ? -(*x1) : 2*(*x1));
    f2 = flin(S, x2, j);
    if (f2 <= fm) {
       xm = x2;
       fm = f2;
    }
    *d2 = (x2*(f1-f0) - (*x1)*(f2-f0))/((*x1)*x2*((*x1)-x2));
 }
 d1 = (f1-f0)/(*x1) - *x1**d2; dz = 1;
 if (*d2 <= S->small) {
    x2 = (d1 < 0 ? S->h : -S->h);
 }
 else {
    x2 = - 0.5*d1/(*d2);
 }
 if (fabs(x2) > S->h)
    x2 = (x2 > 0 ? S->h : -S->h);
test:
 f2 = flin(S, x2, j);
 if ((k < nits) && (f2 > f0)) {
    k++;
    if ((f0 < f1) && (*x1*x2 > 0.0))
       goto next;
    x2 *= 0.5;
    goto test;
 }
 S->nl++;
 if (f2 > fm) x2 = xm; else fm = f2;
 if (fabs(x2*(x2-*x1)) > S->small) {
    *d2 = (x2*(f1-f0) - *x1*(fm-f0))/(*x1*x2*(*x1-x2));
 }
 else {
    if (k > 0) *d2 = 0;
 }
 if (*d2 <= S->small) *d2 = S->small;
 *x1 = x2; S->fx = fm;
 if (sf1 < S->fx) {
    S->fx = sf1;
    *x1 = sx1;
 }
 if (j != -1)
    for (i=0; i<S->n; i++)
        S->x[i] += (*x1)*S->v[i][j];
}

/* --------------------------------------------------------------------------- */
void quadprax(PRAXSTATE* S)
/* look for a minimum along the curve q0, q1, q2 */
{
 int i;
 double l, s;

 s = S->fx; S->fx = S->qf1; S->qf1 = s; S->qd1 = 0.0;
 for (i=0; i<S->n; i++) {
     s = S->x[i]; l = S->q1[i]; S->x[i] = l; S->q1[i] = s;
     S->qd1 = S->qd1 + (s-l)*(s-l);
 }
 s = 0.0; S->qd1 = sqrt(S->qd1); l = S->qd1;
 if (S->qd0>0.0 && S->qd1>0.0 &&S->nl>=3*S->n*S->n) {
    min1(S, -1, 2, &s, &l, S->qf1, 1);
    S->qa = l*(l-S->qd1)/(S->qd0*(S->qd0+S->qd1));
    S->qb = (l+S->qd0)*(S->qd1-l)/(S->qd0*S->qd1);
    S->qc = l*(l+S->qd0)/(S->qd1*(S->qd0+S->qd1));
 }
 else {
    S->fx = S->qf1; S->qa = S->qb = 0.0; S->qc = 1.0;
 }
 S->qd0 = S->qd1;
 for (i=0; i<S->n; i++) {
     s = S->q0[i]; S->q0[i] = S->x[i];
     S->x[i] = S->qa*s + S->qb*S->x[i] + S->qc*S->q1[i];
 }
}

PRAXSTATE*	IntiPraxis(double(*_fun)(void*, double*), double *_x, int _n, int pr, int il, int ktm2, int LMaxFun)
{
	PRAXSTATE*	PState;

	PState = AllocPraxState();

	 /* init global extern variables and parameters */
//	PState->macheps = EPSILON; 
	PState->macheps = EPSILON * EPSILON * EPSILON * EPSILON; 
	PState->h		= PState->step;
	PState->t		= PState->tol;

	PState->n		= _n;
	PState->x		= _x; 
	PState->fun		= _fun;

	PState->prin	= pr;
	PState->illc	= il;
	PState->ktm		= ktm2;
	PState->maxfun	= LMaxFun;

	PState->small	= PState->macheps*PState->macheps;
	PState->vsmall	= PState->small*PState->small;
	PState->large	= 1.0/PState->small; 
	PState->vlarge	= 1.0/PState->vsmall;

	PState->m2		= sqrt(PState->macheps);
	PState->m4		= sqrt(PState->m2);
	PState->ldfac	= (PState->illc ? 0.1 : 0.01);

	PState->nl		= PState->kt = 0;
	PState->nf		= 1;

	PState->t2		= PState->small + fabs(PState->t);
	PState->t		= PState->t2; 
	PState->dmin	= PState->small;

	if (PState->h < 100.0*PState->t)
		PState->h = 100.0*PState->t;
 
	PState->ldt = PState->h;

	return PState;
}

/* --------------------------------------------------------------------------- */
double		praxis(PRAXSTATE* PState)
{
	PState->fx		= (*PState->fun)(PState, PState->x);
	PState->qf1		= PState->fx;


	for (PState->i=0; PState->i<PState->n; PState->i++)
		for (PState->j=0; PState->j<PState->n; PState->j++)
			PState->v[PState->i][PState->j] = (PState->i == PState->j ? 1.0 : 0.0);

	PState->d[0] = 0.0;
	PState->qd0 = 0.0;
 
	for (PState->i=0; PState->i<PState->n; PState->i++) 
		PState->q1[PState->i] = PState->x[PState->i];
 
	if(PState->prin > 1)
	{
		printf("\n------------- enter function praxis -----------\n");
		printf("... current parameter settings ...\n");
		printf("... scaling ... %20.10e\n", PState->scbd);
		printf("...   tol   ... %20.10e\n", PState->t);
		printf("... maxstep ... %20.10e\n", PState->h);
		printf("...   illc  ... %20u\n", PState->illc);
		printf("...   ktm   ... %20u\n", PState->ktm);
		printf("... maxfun  ... %20u\n", PState->maxfun);
	}

	if (PState->prin) 
		print(PState);

	mloop:
		PState->sf = PState->d[0];
		PState->s = PState->d[0] = 0.0;

	/* minimize along first direction */
	min1(PState, 0, 2, &PState->d[0], &PState->s, PState->fx, 0);
 if (PState->s <= 0.0)
    for (PState->i=0; PState->i < PState->n; PState->i++)
        PState->v[PState->i][0] = -PState->v[PState->i][0];
 if ((PState->sf <= (0.9 * PState->d[0])) || ((0.9 * PState->sf) >= PState->d[0]))
    for (PState->i=1; PState->i<PState->n; PState->i++)
        PState->d[PState->i] = 0.0;
 for (PState->k=1; PState->k<PState->n; PState->k++) {
     for (PState->i=0; PState->i<PState->n; PState->i++)
         PState->y[PState->i] = PState->x[PState->i];
     PState->sf = PState->fx;
     PState->illc = PState->illc || (PState->kt > 0);
next:
     PState->kl = PState->k;
     PState->df = 0.0;
     if (PState->illc) {        /* random step to get off resolution valley */
        for (PState->i=0; PState->i<PState->n; PState->i++)
		{
            PState->z[PState->i] = (0.1 * PState->ldt + PState->t2 * pow(10.0,(double)PState->kt)) * (RandDouble(PState->Rates->RS) - 0.5);
            PState->s = PState->z[PState->i];
            for (PState->j=0; PState->j < PState->n; PState->j++)
                PState->x[PState->j] += PState->s * PState->v[PState->j][PState->i];
        }
 
		/* Old  
		fx = (*fun)(x);
        */
		PState->fx = PState->fun(PState, PState->x);
		
		PState->nf++;
     }
     /* minimize along non-conjugate directions */
     for (PState->k2=PState->k; PState->k2<PState->n; PState->k2++) {
         PState->sl = PState->fx;
         PState->s = 0.0;
         min1(PState, PState->k2, 2, &PState->d[PState->k2], &PState->s, PState->fx, 0);
         if (PState->illc) {
            double szk = PState->s + PState->z[PState->k2];
            PState->s = PState->d[PState->k2] * szk*szk;
         }
         else
            PState->s = PState->sl - PState->fx;
         if (PState->df < PState->s) {
            PState->df = PState->s;
            PState->kl = PState->k2;
         }
     }
     if (!PState->illc && (PState->df < fabs(100.0 * PState->macheps * PState->fx))) {
        PState->illc = 1;
        goto next;
     }
     if ((PState->k == 1) && (PState->prin > 1))
        vecprint("\n... New Direction ...",PState->d,PState->n);
     /* minimize along conjugate directions */
     for (PState->k2=0; PState->k2<=PState->k-1; PState->k2++) {
         PState->s = 0.0;
         min1(PState, PState->k2, 2, &PState->d[PState->k2], &PState->s, PState->fx, 0);
     }
     PState->f1 = PState->fx;
     PState->fx = PState->sf;
     PState->lds = 0.0;
     for (PState->i=0; PState->i<PState->n; PState->i++)
	 {
         PState->sl = PState->x[PState->i];
         PState->x[PState->i] = PState->y[PState->i];
         PState->y[PState->i] = PState->sl - PState->y[PState->i];
         PState->sl = PState->y[PState->i];
         PState->lds = PState->lds + PState->sl*PState->sl;
     }
     PState->lds = sqrt(PState->lds);
     if (PState->lds > PState->small) {
        for (PState->i=PState->kl-1; PState->i>=PState->k; PState->i--)
		{
            for (PState->j=0; PState->j < PState->n; PState->j++)
                PState->v[PState->j][PState->i+1] = PState->v[PState->j][PState->i];
                PState->d[PState->i+1] = PState->d[PState->i];
            }
            PState->d[PState->k] = 0.0;
            for (PState->i=0; PState->i < PState->n; PState->i++)
                PState->v[PState->i][PState->k] = PState->y[PState->i] / PState->lds;
            min1(PState, PState->k, 4, &PState->d[PState->k], &PState->lds, PState->f1, 1);
            if (PState->lds <= 0.0) {
               PState->lds = -PState->lds;
               for (PState->i=0; PState->i<PState->n; PState->i++)
                   PState->v[PState->i][PState->k] = -PState->v[PState->i][PState->k];
            }
     }
     PState->ldt = PState->ldfac * PState->ldt;
     if (PState->ldt < PState->lds)
        PState->ldt = PState->lds;
     if (PState->prin > 1)
        print(PState);
     PState->t2 = 0.0;
     for (PState->i=0; PState->i<PState->n; PState->i++)
         PState->t2 += PState->x[PState->i]*PState->x[PState->i];
     PState->t2 = PState->m2 * sqrt(PState->t2) + PState->t;
     if (PState->ldt > (0.5 * PState->t2))
        PState->kt = 0;
     else
        PState->kt++;
     if (PState->kt > PState->ktm)
        goto fret;
 }
 /*  try quadratic extrapolation in case */
 /*  we are stuck in a curved valley */
 quadprax(PState);
 PState->dn = 0.0;
 for (PState->i=0; PState->i<PState->n; PState->i++)
 {
     PState->d[PState->i] = 1.0 / sqrt(PState->d[PState->i]);
     if (PState->dn < PState->d[PState->i])
        PState->dn = PState->d[PState->i];
 }
 if (PState->prin > 2)
    matprint("\n... New Matrix of Directions ...",PState->v,PState->n);
 for (PState->j=0; PState->j<PState->n; PState->j++) {
     PState->s = PState->d[PState->j] / PState->dn;
     for (PState->i=0; PState->i < PState->n; PState->i++)
         PState->v[PState->i][PState->j] *= PState->s;
 }
 if (PState->scbd > 1.0) {       /* scale axis to reduce condition number */
    PState->s = PState->vlarge;
    for (PState->i=0; PState->i<PState->n; PState->i++) {
        PState->sl = 0.0;
        for (PState->j=0; PState->j < PState->n; PState->j++)
            PState->sl += PState->v[PState->i][PState->j]*PState->v[PState->i][PState->j];
        PState->z[PState->i] = sqrt(PState->sl);
        if (PState->z[PState->i] < PState->m4)
           PState->z[PState->i] = PState->m4;
        if (PState->s > PState->z[PState->i])
           PState->s = PState->z[PState->i];
    }
    for (PState->i=0; PState->i<PState->n; PState->i++) {
        PState->sl = PState->s / PState->z[PState->i];
        PState->z[PState->i] = 1.0 / PState->sl;
        if (PState->z[PState->i] > PState->scbd) {
           PState->sl = 1.0 / PState->scbd;
           PState->z[PState->i] = PState->scbd;
        }
    }
 }
 for (PState->i=1; PState->i<PState->n; PState->i++)
     for (PState->j=0; PState->j<=PState->i-1; PState->j++) {
         PState->s = PState->v[PState->i][PState->j];
         PState->v[PState->i][PState->j] = PState->v[PState->j][PState->i];
         PState->v[PState->j][PState->i] = PState->s;
     }
 minfit(PState->n, PState->macheps, PState->vsmall, PState->v, PState->d);
 if (PState->scbd > 1.0) {
    for (PState->i=0; PState->i<PState->n; PState->i++) {
        PState->s = PState->z[PState->i];
        for (PState->j=0; PState->j<PState->n; PState->j++)
            PState->v[PState->i][PState->j] *= PState->s;
    }
    for (PState->i=0; PState->i<PState->n; PState->i++) {
        PState->s = 0.0;
        for (PState->j=0; PState->j<PState->n; PState->j++)
            PState->s += PState->v[PState->j][PState->i]*PState->v[PState->j][PState->i];
        PState->s = sqrt(PState->s);
        PState->d[PState->i] *= PState->s;
        PState->s = 1.0 / PState->s;
        for (PState->j=0; PState->j<PState->n; PState->j++)
            PState->v[PState->j][PState->i] *= PState->s;
    }
 }
 for (PState->i=0; PState->i<PState->n; PState->i++) {
     if ((PState->dn * PState->d[PState->i]) > PState->large)
        PState->d[PState->i] = PState->vsmall;
     else if ((PState->dn * PState->d[PState->i]) < PState->small)
        PState->d[PState->i] = PState->vlarge;
     else
        PState->d[PState->i] = pow(PState->dn * PState->d[PState->i],-2.0);
 }
 sort(PState);               /* the new eigenvalues and eigenvectors */
 PState->dmin = PState->d[PState->n-1];
 if (PState->dmin < PState->small)
    PState->dmin = PState->small;
 PState->illc = (PState->m2 * PState->d[0]) > PState->dmin;
 if ((PState->prin > 2) && (PState->scbd > 1.0))
    vecprint("\n... Scale Factors ...",PState->z,PState->n);
 if (PState->prin > 2)
    vecprint("\n... Eigenvalues of A ...",PState->d,PState->n);
 if (PState->prin > 2)
    matprint("\n... Eigenvectors of A ...",PState->v,PState->n);
 if ((PState->maxfun > 0) && (PState->nl > PState->maxfun)) {
    if (PState->prin)
       printf("\n... maximum number of function calls reached ...\n");
    goto fret;
 }
 goto mloop; 	 /* back to main loop */

fret:
 if (PState->prin > 0) {
    vecprint("\n... Final solution is ...", PState->x, PState->n);
    printf("\n... ChiSq reduced to %20.10e ...\n", PState->fx);
    printf("... after %20u function calls.\n", PState->nf);
 }

	return(PState->fx);
}

/* --------------------------------------------------------------------------- */


