#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <time.h>

#ifdef CLIK_P
	#include <cilk/cilk.h>
	#include <cilk/cilk_api.h>
#endif

#include "typedef.h"
#include "Threaded.h"

#ifdef USE_MLK
	#include "InvMKL.h"
#endif

int		GetThreadNo(void)
{
#ifdef OPENMP_THR
	return omp_get_thread_num();
#endif	
	return 0;
}

int		GetMaxThreads(void)
{
#ifdef OPENMP_THR
	return omp_get_num_procs();
#endif
	
	return 1;
}

void	SetNoOfThreads(int No)
{
#ifdef USE_MLK
	mkl_set_num_threads(No);
#endif

#ifdef OPENMP_THR
	omp_set_num_threads(No);
	return; 
#endif	

#ifdef CLIK_P
	char *TStr;

	TStr = (char*)SMalloc(sizeof(char) * 64);
	sprintf(TStr, "%d", No);

	if (0 != __cilkrts_set_param("nworkers",TStr))
	{
		printf("Failed to set worker count\n");
		exit(1);
	}

	free(TStr);
	 __cilkrts_init();
	return;
#endif

	return ;
}

double	GetSeconds(void)
{
#ifndef OPENMP_THR
	return  (double)time(NULL);
#else
	return omp_get_wtime();
#endif
}
