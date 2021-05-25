#ifndef BTOCL_RUNTIME_KERNELS_H
#define BTOCL_RUNTIME_KERNELS_H

#ifdef BTOCL

#include "btocl_runtime.h"

// added to get access to TREES and OPT
//#include "typedef.h"


// Don't forget "/" at the end
#define BTOCL_KERNELS "C:/MinGW/msys/1.0/home/User/code/btocl_kernels/"
// #define BTOCL_KERNELS "/Users/igor/z/academia/gpuproject/btJan26/btocl_kernels/"

// ********** Kernel indexes  ***********
// Important: Do not repeat, index < NUM_KERNELS
// it's ok to skip indexes
#define BTOCL_CHOLUPDCOL_PURE 0
#define BTOCL_CHOLUPDMAT_P4 1
#define BTOCL_CHOLUPDMAT_P8 2
#define BTOCL_CHOLUPDMAT_P16 3
#define BTOCL_CHOLUPDMAT_P32 4
#define BTOCL_CHOLUPDMAT_P64 5
#define BTOCL_CHOLUPDMAT_P128 6
#define BTOCL_CHOLUPDCOL 7
#define BTOCL_CHOLUPDMAT 8
#define BTOCL_TRANSPOSE 9
#define BTOCL_LTRI_MBYL 10
#define BTOCL_LTRI_LBYM 11
#define BTOCL_EXPQT 12
#define BTOCL_EXPQT_ROWNOS2 13
#define BTOCL_EXPQT_LOCAL 14
#define BTOCL_EXPQT_ROWNOS4 15
#define BTOCL_PLH 16
#define BTOCL_PLHNODE 17
#define BTOCL_PLHROW 18
#define BTOCL_PLHROWG 19
#define BTOCL_PLHROWALLG 20
#define BTOCL_PLHROWGL 21
#define BTOCL_PLHNODE_NOS2 22
#define BTOCL_PLHROWFULLG 23
#define BTOCL_PLHREDUCEGL 24
#define BTOCL_EXPQT_NOTEMP 25
#define BTOCL_EXPQT3 26
#define BTOCL_EXPQT3_LOCAL 27
#define BTOCL_EXPQT3_LOCALERR 28
#define BTOCL_EXPQT3_NOS4 29
#define BTOCL_EXPQT_LOCALERR 30
#define BTOCL_KRON_BASIC 31
#define BTOCL_KRON_ACCUM 32
#define BTOCL_KRON_LOG1 33
#define BTOCL_CHOLUPDMAT_PURE 34
#define BTOCL_LTRI_LTBYL 35
#define BTOCL_WRITEUTOL_DIAG 36

// ********* end kernel indexes ************

void initialise_kernel_info();
cl_int btocl_load_all(int cont,int disc, int nos, int nsites);
// internal
int btocl_load_continuousKernels(BTOCL_RUNTIME* rt);
int btocl_load_discreteKernels(BTOCL_RUNTIME* rt);

int load_kernel_file(unsigned short kernel_idx, char* kernel_name, char* kernel_file, BTOCL_RUNTIME* rt);
cl_int createProgramFromBuffers(cl_context context, cl_program* program, BTOCL_RUNTIME* rt);
cl_int freeProgramBuffers(BTOCL_RUNTIME* rt);
cl_int loadKernels(BTOCL_RUNTIME* rt);
cl_int createProgramVerbatim(cl_context context,cl_program *program);


#endif // if BTOCL defined

#endif

