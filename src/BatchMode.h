#ifndef BATCHMODE_H
#define BATCHMODE_H


// Run the batch command using openMP
//#define PBATCH

#ifdef PBATCH
#include <omp.h>
#endif

void	BatchRun(char *BatchFN);

#endif