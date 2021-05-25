#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

#include "typedef.h"
#include "threaded.h"


int		GetThreadNo(void)
{
#ifndef THREADED
	return 0;
#else
	return omp_get_thread_num();
#endif	
}

int		GetMaxThreads(void)
{
#ifndef THREADED
	return 1;
#else
	return omp_get_num_procs();
#endif	
}

void	SetNoOfThreads(int No)
{
#ifndef THREADED
	return ;
#else
	omp_set_num_threads(No);
	return; 
#endif	
}


double	GetSeconds(void)
{
#ifndef THREADED
	return  (double)time(NULL);
#else
	return omp_get_wtime();
#endif
}