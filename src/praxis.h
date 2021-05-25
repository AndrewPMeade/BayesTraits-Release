/*
	f u n c t i o n     p r a x i s

 praxis is a general purpose routine for the minimization of a
 function in several variables. the algorithm used is a modifi-
 cation of conjugate gradient search method by powell. the changes
 are due to r.p. brent, who gives an algol-w program, which served
 as a basis for this function.

 references:
     - powell, m.j.d., 1964. an efficient method for finding
       the minimum of a function in several variables without
       calculating derivatives, computer journal, 7, 155-162
     - brent, r.p., 1973. algorithms for minimization without
       derivatives, prentice hall, englewood cliffs.

     problems, suggestions or improvements are always wellcome
                       karl gegenfurtner   07/08/87
														 c - version

 usage: min = praxis(fun, x, n);

  fun        the function to be minimized. fun is called from
				 praxis with x and n as arguments
  x          a double array containing the initial guesses for
				 the minimum, which will contain the solution on
				 return
  n          an integer specifying the number of unknown
				 parameters
  min        praxis returns the least calculated value of fun

 some additional global variables control some more aspects of
 the inner workings of praxis. setting them is optional, they
 are all set to some reasonable default values given below.

	prin      controls the printed output from the routine.
				 0 -> no output
				 1 -> print only starting and final values
				 2 -> detailed map of the minimization process
				 3 -> print also eigenvalues and vectors of the
						search directions
				 the default value is 1

  tol        is the tolerance allowed for the precision of the
				 solution. praxis returns if the criterion
				 2 * ||x[k]-x[k-1]|| <= sqrt(macheps) * ||x[k]|| + tol
				 is fulfilled more than ktm times.
				 the default value depends on the machine precision

  ktm        see just above. default is 1, and a value of 4 leads
				 to a very(!) cautious stopping criterion.

  step       is a steplength parameter and should be set equal
				 to the expected distance from the solution.
				 exceptionally small or large values of step lead to
				 slower convergence on the first few iterations
				 the default value for step is 1.0

  scbd       is a scaling parameter. 1.0 is the default and
				 indicates no scaling. if the scales for the different
				 parameters are very different, scbd should be set to
				 a value of about 10.0.

  illc       should be set to true (1) if the problem is known to
				 be ill-conditioned. the default is false (0). this
				 variable is automatically set, when praxis finds
				 the problem to be ill-conditioned during iterations.

  maxfun     is the maximum number of calls to fun allowed. praxis
				 will return after maxfun calls to fun even when the
				 minimum is not yet found. the default value of 0
				 indicates no limit on the number of calls.
				 this return condition is only checked every n
				 iterations.  */


#ifndef PRAXIS
#define PRAXIS

#define EPSILON 1.0e-8
#define SQREPSILON 1.0e-16


#define N 90

double rndom(void);
void sort(void);
void print(void);
void matprint(char *s,double v[N][N], int n);
void vecprint(char *s, double x[N], int n);
double flin(double l, int j);

#ifdef min1
       #undef min1
       void min1(int j, int nits, double *d2, double *x1, double f1, int fk);

#endif

void quadprax(void);
double praxis(double(*_fun)(double*), double *_x, int _n, int pr, int il, int ktm2, int LMaxFun);

void	BlankPraxisGlobal(void);
#endif