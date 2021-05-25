/*
*  BayesTriats 3.0
*
*  copyright 2017
*
*  Andrew Meade
*  School of Biological Sciences
*  University of Reading
*  Reading
*  Berkshire
*  RG6 6BX
*
* BayesTriats is free software: you can redistribute it and/or modify
* it under the terms of the GNU General Public License as published by
* the Free Software Foundation, either version 3 of the License, or
* (at your option) any later version.
* 
* This program is distributed in the hope that it will be useful,
* but WITHOUT ANY WARRANTY; without even the implied warranty of
* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
* GNU General Public License for more details.
* 
* You should have received a copy of the GNU General Public License
* along with this program.  If not, see <http://www.gnu.org/licenses/>
*
*/





#ifdef BTOCL

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <limits.h>

#include "typedef.h" 


#include "btocl_kernels_bayestraits.h"
#include "btocl_kernels_discrete.h"
#include "btocl_kernels_continuous.h"

// Declare kernel ids
// CL_USHRT_MAX
unsigned short BTOCL_CHOLUPDCOL_PURE = CL_USHRT_MAX;
unsigned short BTOCL_CHOLUPDMAT_P4 = CL_USHRT_MAX;
unsigned short BTOCL_CHOLUPDMAT_P8 = CL_USHRT_MAX;
unsigned short BTOCL_CHOLUPDMAT_P16 = CL_USHRT_MAX;
unsigned short BTOCL_CHOLUPDMAT_P32 = CL_USHRT_MAX;
unsigned short BTOCL_CHOLUPDMAT_P64 = CL_USHRT_MAX;
unsigned short BTOCL_CHOLUPDMAT_P128 = CL_USHRT_MAX;
unsigned short BTOCL_CHOLUPDCOL = CL_USHRT_MAX;
unsigned short BTOCL_CHOLUPDMAT = CL_USHRT_MAX;
unsigned short BTOCL_TRANSPOSE = CL_USHRT_MAX;
unsigned short BTOCL_LTRI_MBYL = CL_USHRT_MAX;
unsigned short BTOCL_LTRI_LBYM = CL_USHRT_MAX;
unsigned short BTOCL_EXPQT = CL_USHRT_MAX;
unsigned short BTOCL_EXPQT_ROWNOS2 = CL_USHRT_MAX;
unsigned short BTOCL_EXPQT_LOCAL = CL_USHRT_MAX;
unsigned short BTOCL_EXPQT_ROWNOS4 = CL_USHRT_MAX;
unsigned short BTOCL_PLH = CL_USHRT_MAX;
unsigned short BTOCL_PLHNODE = CL_USHRT_MAX;
unsigned short BTOCL_PLHROW = CL_USHRT_MAX;
unsigned short BTOCL_PLHROWG = CL_USHRT_MAX;
unsigned short BTOCL_PLHROWALLG = CL_USHRT_MAX;
unsigned short BTOCL_PLHROWGL = CL_USHRT_MAX;
unsigned short BTOCL_PLHNODE_NOS2 = CL_USHRT_MAX;
unsigned short BTOCL_PLHROWFULLG = CL_USHRT_MAX;
unsigned short BTOCL_PLHREDUCEGL = CL_USHRT_MAX;
unsigned short BTOCL_EXPQT_NOTEMP = CL_USHRT_MAX;
unsigned short BTOCL_EXPQT3 = CL_USHRT_MAX;
unsigned short BTOCL_EXPQT3_LOCAL = CL_USHRT_MAX;
unsigned short BTOCL_EXPQT3_LOCALERR = CL_USHRT_MAX;
unsigned short BTOCL_EXPQT3_NOS4 = CL_USHRT_MAX;
unsigned short BTOCL_EXPQT_LOCALERR = CL_USHRT_MAX;
unsigned short BTOCL_KRON_BASIC = CL_USHRT_MAX;
unsigned short BTOCL_KRON_ACCUM = CL_USHRT_MAX;
unsigned short BTOCL_KRON_LOG1 = CL_USHRT_MAX;
unsigned short BTOCL_KRON_ACCUM_SIGMA1 = CL_USHRT_MAX;
unsigned short BTOCL_KRON_TILES1 = CL_USHRT_MAX;
unsigned short BTOCL_CHOLUPDMAT_PURE = CL_USHRT_MAX;
unsigned short BTOCL_LTRI_LTBYL = CL_USHRT_MAX;
unsigned short BTOCL_WRITEUTOL_DIAG = CL_USHRT_MAX;
unsigned short BTOCL_LTRI_LTBYL_B = CL_USHRT_MAX;
unsigned short BTOCL_LTRI_LTBYL_B0 = CL_USHRT_MAX;
unsigned short BTOCL_EXP_EIGEN = CL_USHRT_MAX;


// New version that allows the use of compiler options and compilation of selected kernels
cl_int btocl_load_all(int continuous,int discrete, int nos, int nsites) {
  size_t log_size;
  char* program_log;
  int nos2;
  cl_int err;
  char options[64];
  char number[8];

  BTOCL_RUNTIME* rt = btocl_getruntime();



  options[0] = '\0';
  //nos = Trees->NoStates;   // BTOCL_NOS
  nos2 = nos*nos;            // BTOCL_NOS2
  //nsites = Trees->NoSites; // BTOCL_NSITES

  //if (Opt->ModelType == MT_CONTINUOUS) {
  if (continuous) {
	//printf("kernels continuous\n");
    err = btocl_load_continuousKernels(rt);
	//printf("...done\n");
    if (err != 0)
      return err;
  }

  //	if (Opt->ModelType == MT_DISCRETE) {
  if (discrete) {
	//printf("kernels discrete\n");
    err = btocl_load_discreteKernels(rt);
	//printf("....done\n");
    if (err != 0)
      return err;
  }

  if (discrete) {
    // load constants
    sprintf(number,"%d",nos);
    strcpy(options,"-D BT_NOS=");
    strcat(options,number);
    sprintf(number,"%d",nos2);
    strcat(options," -D BT_NOS2=");
    strcat(options,number);
    sprintf(number,"%d",nsites);
    strcat(options," -D BT_NSITES=");
    strcat(options,number);
    sprintf(number,"%d",nsites*nos);
    strcat(options," -D BT_NOSNSITES=");
    strcat(options,number);
#ifdef USE_BTEXP
	sprintf(number,"%d",1);
    strcat(options," -D USE_BTEXP=");
    strcat(options,number);
#endif

    //printf("%s \n",options);
    //exit(0);
  }

  //printf("loading kernrels\n");
  err = btocl_build_load_kernels(options);
  //printf("....done\n");

  return err;


}


int btocl_load_continuousKernels(BTOCL_RUNTIME* rt) {
  cl_int err;

  // kernel_index, kernel_name, program file name
//  err = load_kernel_file(&BTOCL_CHOLUPDCOL_PURE,"btocl_cholupdcol_pure",
//		   BTOCL_KERNELS_DIR "kernel_cholupdcol_pure.cl",rt);
  err = load_kernel_string(&BTOCL_CHOLUPDCOL_PURE,"btocl_cholupdcol_pure",
		   STRING_CHOLUPDCOL_PURE,rt);
  if (err != 0) return err;
  /*
  err = load_kernel_file(&BTOCL_CHOLUPDMAT_P4,"btocl_cholupdmat_p4",
		   BTOCL_KERNELS_DIR "kernel_cholupdmat_p4.cl",rt);
  if (err != 0) return err;
  err = load_kernel_file(&BTOCL_CHOLUPDMAT_P8,"btocl_cholupdmat_p8",
		   BTOCL_KERNELS_DIR "kernel_cholupdmat_p8.cl",rt);
  if (err != 0) return err;
  err = load_kernel_file(&BTOCL_CHOLUPDMAT_P16,"btocl_cholupdmat_p16",
		   BTOCL_KERNELS_DIR "kernel_cholupdmat_p16.cl",rt);
  if (err != 0) return err;
  err = load_kernel_file(&BTOCL_CHOLUPDMAT_P32,"btocl_cholupdmat_p32",
		   BTOCL_KERNELS_DIR "kernel_cholupdmat_p32.cl",rt);
  if (err != 0) return err;
  err = load_kernel_file(&BTOCL_CHOLUPDMAT_P64,"btocl_cholupdmat_p64",
		   BTOCL_KERNELS_DIR "kernel_cholupdmat_p64.cl",rt);
  if (err != 0) return err;
  err = load_kernel_file(&BTOCL_CHOLUPDMAT_P128,"btocl_cholupdmat_p128",
		   BTOCL_KERNELS_DIR "kernel_cholupdmat_p128.cl",rt);
  if (err != 0) return err;

  err = load_kernel_file(&BTOCL_CHOLUPDCOL,"btocl_cholupdcol",
		   BTOCL_KERNELS_DIR "kernel_cholupdcol.cl",rt);
  if (err != 0) return err;
  err = load_kernel_file(&BTOCL_CHOLUPDMAT,"btocl_cholupdmat",
		   BTOCL_KERNELS_DIR "kernel_cholupdmat.cl",rt);
  if (err != 0) return err;
  */
  // Chol - revisited
  //err = load_kernel_file(&BTOCL_CHOLUPDMAT_PURE,"btocl_cholupdmat_pure",
	//	   BTOCL_KERNELS_DIR "kernel_cholupdmat_pure.cl",rt);
  err = load_kernel_string(&BTOCL_CHOLUPDMAT_PURE,"btocl_cholupdmat_pure",
		   STRING_CHOLUPDMAT_PURE,rt);
  if (err != 0) return err;

  //err = load_kernel_file(&BTOCL_TRANSPOSE,"btocl_transpose",
  //		   BTOCL_KERNELS_DIR "kernel_transpose.cl",rt);
  //if (err != 0) return err;
  //err = load_kernel_file(&BTOCL_LTRI_LBYM,"btocl_ltri_LbyM",
	//	   BTOCL_KERNELS_DIR "kernel_ltri_LbyM.cl",rt);
  err = load_kernel_string(&BTOCL_LTRI_LBYM,"btocl_ltri_LbyM",
		   STRING_LTRI_LBYM,rt);
  if (err != 0) return err;

  //err = load_kernel_file(&BTOCL_LTRI_MBYL,"btocl_ltri_MbyL",
	//	   BTOCL_KERNELS_DIR "kernel_ltri_MbyL.cl",rt);
  err = load_kernel_string(&BTOCL_LTRI_MBYL,"btocl_ltri_MbyL",
		   STRING_LTRI_MBYL,rt);
  if (err != 0) return err;

  //err = load_kernel_file(&BTOCL_LTRI_LTBYL_B0,"btocl_ltri_LTbyL_b0",
	//	   BTOCL_KERNELS_DIR "kernel_ltri_LTbyL_b0.cl",rt);
  err = load_kernel_string(&BTOCL_LTRI_LTBYL_B0,"btocl_ltri_LTbyL_b0",
		   STRING_LTRI_LTBYL_B0,rt);
  if (err != 0) return err;
  //err = load_kernel_file(&BTOCL_LTRI_LTBYL_B,"btocl_ltri_LTbyL_b",
	//	   BTOCL_KERNELS_DIR "kernel_ltri_LTbyL_b.cl",rt);
  err = load_kernel_string(&BTOCL_LTRI_LTBYL_B,"btocl_ltri_LTbyL_b",
		   STRING_LTRI_LTBYL_B,rt);
  if (err != 0) return err;
  //err = load_kernel_file(&BTOCL_WRITEUTOL_DIAG,"btocl_writeUtoL_diag",
	//	   BTOCL_KERNELS_DIR "kernel_writeUtoL_diag.cl",rt);
   err = load_kernel_string(&BTOCL_WRITEUTOL_DIAG,"btocl_writeUtoL_diag",
		   STRING_WRITEUTOL_DIAG,rt);
  if (err != 0) return err;

  // KRONECKER PRODUCT
  //err = load_kernel_file(&BTOCL_KRON_TILES1,"btocl_kron_tiles1",
	//	   BTOCL_KERNELS_DIR "kernel_kron_tiles1.cl",rt);
  err = load_kernel_string(&BTOCL_KRON_TILES1,"btocl_kron_tiles1",
		   STRING_KRON_TILES1,rt);
  if (err != 0) return err;
   // err = load_kernel_file(&BTOCL_KRON_LOG1,"btocl_kron_log1",
	//	   BTOCL_KERNELS_DIR "kernel_kron_log1.cl",rt);
   err = load_kernel_string(&BTOCL_KRON_LOG1,"btocl_kron_log1",
		   STRING_KRON_LOG1,rt);
  if (err != 0) return err;
  //err = load_kernel_file(&BTOCL_KRON_ACCUM,"btocl_kron_accum",
	//	   BTOCL_KERNELS_DIR "kernel_kron_accum.cl",rt);
  err = load_kernel_string(&BTOCL_KRON_ACCUM,"btocl_kron_accum",
		   STRING_KRON_ACCUM,rt);
  if (err != 0) return err;

  /* not used
  err = load_kernel_file(&BTOCL_KRON_BASIC,"btocl_kron_basic",
		   BTOCL_KERNELS_DIR "kernel_kron_basic.cl",rt);
  if (err != 0) return err;
  err = load_kernel_file(&BTOCL_KRON_ACCUM_SIGMA1,"btocl_kron_accum_sigma1",
		   BTOCL_KERNELS_DIR "kernel_kron_accum_sigma1.cl",rt);
  if (err != 0) return err;
  */


  //err = load_kernel_file(&BTOCL_LTRI_LTBYL,"btocl_ltri_LTbyL",
//		   BTOCL_KERNELS_DIR "kernel_ltri_LTbyL.cl",rt);
  //if (err != 0) return err;


  return 0;

}

int btocl_load_discreteKernels(BTOCL_RUNTIME* rt) {
  cl_int err;
  // EXPQT

  //err = load_aux_file(BTOCL_KERNELS_DIR "aux_btexp.c",rt);
  err = load_aux_string(STRING_AUX_BEXP,rt);
  if (err != 0) return err;

  /*
  printf("loaded aux file\n");
  err = load_kernel_file(&BTOCL_EXPQT,"btocl_expqt",
		   BTOCL_KERNELS_DIR "kernel_expqt.cl",rt);
  if (err != 0) return err;
  printf("laoded btocl_expqt\n");

  err = load_kernel_file(&BTOCL_EXPQT_LOCAL,"btocl_expqt_local",
		   BTOCL_KERNELS_DIR "kernel_expqt_local.cl",rt);
  if (err != 0) return err;


  // EXPT revisited
  err = load_kernel_file(&BTOCL_EXPQT_NOTEMP,"btocl_expqt_notemp",
		   BTOCL_KERNELS_DIR "kernel_expqt_notemp.cl",rt);
  if (err != 0) return err;
  err = load_kernel_file(&BTOCL_EXPQT3,"btocl_expqt3",
		   BTOCL_KERNELS_DIR "kernel_expqt3.cl",rt);
  if (err != 0) return err;
  err = load_kernel_file(&BTOCL_EXPQT3_LOCAL,"btocl_expqt3_local",
		   BTOCL_KERNELS_DIR "kernel_expqt3_local.cl",rt);
  if (err != 0) return err;

  err = load_kernel_file(&BTOCL_EXPQT3_NOS4,"btocl_expqt3_nos4",
		   BTOCL_KERNELS_DIR "kernel_expqt3_nos4.cl",rt);
  if (err != 0) return err;
  err = load_kernel_file(&BTOCL_EXPQT_LOCALERR,"btocl_expqt_localerr",
		   BTOCL_KERNELS_DIR "kernel_expqt_localerr.cl",rt);
  if (err != 0) return err;

  // EXPQT NOS small
  err = load_kernel_file(&BTOCL_EXP_EIGEN,"btocl_exp_eigen",
		   BTOCL_KERNELS_DIR "kernel_exp_eigen.cl",rt);
  if (err != 0) return err;

  // PLH
  err = load_kernel_file(&BTOCL_PLH,"btocl_plh",
		   BTOCL_KERNELS_DIR "kernel_plh.cl",rt);
  if (err != 0) return err;
  err = load_kernel_file(&BTOCL_PLHNODE,"btocl_plhNode",
		   BTOCL_KERNELS_DIR "kernel_plhNode.cl",rt);
  if (err != 0) return err;
  err = load_kernel_file(&BTOCL_PLHROW,"btocl_plhRow",
		   BTOCL_KERNELS_DIR "kernel_plhRow.cl",rt);
  if (err != 0) return err;
  err = load_kernel_file(&BTOCL_PLHROWG,"btocl_plhRowG",
		   BTOCL_KERNELS_DIR "kernel_plhRowG.cl",rt);
  if (err != 0) return err;
  err = load_kernel_file(&BTOCL_PLHROWALLG,"btocl_plhRowAllG",
		   BTOCL_KERNELS_DIR "kernel_plhRowAllG.cl",rt);
  if (err != 0) return err;
  err = load_kernel_file(&BTOCL_PLHROWGL,"btocl_plhRowGL",
		   BTOCL_KERNELS_DIR "kernel_plhRowGL.cl",rt);
  if (err != 0) return err;
  err = load_kernel_file(&BTOCL_PLHNODE_NOS2,"btocl_plhNode_nos2",
		   BTOCL_KERNELS_DIR "kernel_plhNode_nos2.cl",rt);
  if (err != 0) return err;
  err = load_kernel_file(&BTOCL_PLHROWFULLG,"btocl_plhRowFullG",
		   BTOCL_KERNELS_DIR "kernel_plhRowFullG.cl",rt);
  if (err != 0) return err;
  err = load_kernel_file(&BTOCL_PLHREDUCEGL,"btocl_plhReduceGL",
		   BTOCL_KERNELS_DIR "kernel_plhReduceGL.cl",rt);
  if (err != 0) return err;
  */

  // Chosen kernels
  //err = load_kernel_file(&BTOCL_EXPQT_ROWNOS2,"btocl_expqt_rownos2",
	//	   BTOCL_KERNELS_DIR "kernel_expqt_rownos2.cl",rt);
  err = load_kernel_string(&BTOCL_EXPQT_ROWNOS2,"btocl_expqt_rownos2",
		   STRING_EXPQT_ROWNOS2,rt);
  if (err != 0) return err;
   // err = load_kernel_file(&BTOCL_EXPQT_ROWNOS4,"btocl_expqt_rownos4",
	//	   BTOCL_KERNELS_DIR "kernel_expqt_rownos4.cl",rt);
  err = load_kernel_string(&BTOCL_EXPQT_ROWNOS4,"btocl_expqt_rownos4",
		   STRING_EXPQT_ROWNOS4,rt);
  if (err != 0) return err;
  //err = load_kernel_file(&BTOCL_EXPQT3_LOCALERR,"btocl_expqt3_localerr",
	//	   BTOCL_KERNELS_DIR "kernel_expqt3_localerr.cl",rt);
	err = load_kernel_string(&BTOCL_EXPQT3_LOCALERR,"btocl_expqt3_localerr",
		   STRING_EXPQT3_LOCALERR,rt);
  if (err != 0) return err;


	return 0;
}



#endif // if BTOCL defined
