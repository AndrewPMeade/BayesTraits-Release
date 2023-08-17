#ifndef BTOCL_KERNELS_BAYESTRAITS_H
#define BTOCL_KERNELS_BAYESTRAITS_H

#ifdef BTOCL

#include "btocl_runtime.h"

// Use environment variable!!!

// Don't forget "/" at the end
//#define BTOCL_KERNELS_DIR "C:/MinGW/msys/1.0/home/User/code/btocl_kernels/"
#define BTOCL_KERNELS_DIR "C:/GPUProject/code/btocl_kernels/xx"
// #define BTOCL_KERNELS_DIR "/Users/igor/z/academia/gpuproject/btJan26/btocl_kernels/"

// ********** Kernel indexes  ***********
// Important: Do not repeat, index < NUM_KERNELS
// it's ok to skip indexes

/*
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
#define BTOCL_KRON_ACCUM_SIGMA1 34
#define BTOCL_KRON_TILES1 35
#define BTOCL_CHOLUPDMAT_PURE 37
#define BTOCL_LTRI_LTBYL 38
#define BTOCL_WRITEUTOL_DIAG 39
#define BTOCL_LTRI_LTBYL_B 40
#define BTOCL_LTRI_LTBYL_B0 41
#define BTOCL_EXP_EIGEN 45
*/


extern unsigned short BTOCL_CHOLUPDCOL_PURE;
extern unsigned short BTOCL_CHOLUPDMAT_P4;
extern unsigned short BTOCL_CHOLUPDMAT_P8;
extern unsigned short BTOCL_CHOLUPDMAT_P16;
extern unsigned short BTOCL_CHOLUPDMAT_P32;
extern unsigned short BTOCL_CHOLUPDMAT_P64;
extern unsigned short BTOCL_CHOLUPDMAT_P128;
extern unsigned short BTOCL_CHOLUPDCOL;
extern unsigned short BTOCL_CHOLUPDMAT;
extern unsigned short BTOCL_TRANSPOSE;
extern unsigned short BTOCL_LTRI_MBYL;
extern unsigned short BTOCL_LTRI_LBYM;
extern unsigned short BTOCL_EXPQT;
extern unsigned short BTOCL_EXPQT_ROWNOS2;
extern unsigned short BTOCL_EXPQT_LOCAL;
extern unsigned short BTOCL_EXPQT_ROWNOS4;
extern unsigned short BTOCL_PLH;
extern unsigned short BTOCL_PLHNODE;
extern unsigned short BTOCL_PLHROW;
extern unsigned short BTOCL_PLHROWG;
extern unsigned short BTOCL_PLHROWALLG;
extern unsigned short BTOCL_PLHROWGL;
extern unsigned short BTOCL_PLHNODE_NOS2;
extern unsigned short BTOCL_PLHROWFULLG;
extern unsigned short BTOCL_PLHREDUCEGL;
extern unsigned short BTOCL_EXPQT_NOTEMP;
extern unsigned short BTOCL_EXPQT3;
extern unsigned short BTOCL_EXPQT3_LOCAL;
extern unsigned short BTOCL_EXPQT3_LOCALERR;
extern unsigned short BTOCL_EXPQT3_NOS4;
extern unsigned short BTOCL_EXPQT_LOCALERR;
extern unsigned short BTOCL_KRON_BASIC;
extern unsigned short BTOCL_KRON_ACCUM;
extern unsigned short BTOCL_KRON_LOG1;
extern unsigned short BTOCL_KRON_ACCUM_SIGMA1;
extern unsigned short BTOCL_KRON_TILES1;
extern unsigned short BTOCL_CHOLUPDMAT_PURE;
extern unsigned short BTOCL_LTRI_LTBYL;
extern unsigned short BTOCL_WRITEUTOL_DIAG;
extern unsigned short BTOCL_LTRI_LTBYL_B;
extern unsigned short BTOCL_LTRI_LTBYL_B0;
extern unsigned short BTOCL_EXP_EIGEN;


cl_int btocl_load_all(int cont,int disc, int nos, int nsites);
// internal
int btocl_load_continuousKernels(BTOCL_RUNTIME* rt);
int btocl_load_discreteKernels(BTOCL_RUNTIME* rt);



#endif // if BTOCL defined

#endif