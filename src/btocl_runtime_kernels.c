#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#ifdef BTOCL

#include "btocl_runtime_kernels.h"


void initialise_kernel_info() {
	int i;
	BTOCL_RUNTIME* rt = btocl_getruntime();

	// Initialise kernel usage array
	for(i=0; i < NUM_KERNELS; i++) {
	  rt->kernel_use[i]=0;
	  rt->kernel_program_files[i] = NULL;
	  rt->kernel_names[i] = NULL;
	  // added - may have to go
	  rt->program_buffer[i] = NULL;
	  rt->program_size[i] = 0;
	  rt->kernel_num_used = 0;
	}

}

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
  //nos = Trees->NoOfStates;   // BTOCL_NOS
  nos2 = nos*nos;            // BTOCL_NOS2
  //nsites = Trees->NoOfSites; // BTOCL_NSITES
  
  //if (Opt->ModelType == MT_CONTINUOUS) {
  if (continuous) {
    err = btocl_load_continuousKernels(rt);
    if (err != 0)
      return err;
  }
  
  //	if (Opt->ModelType == MT_DISCRETE) {
  if (discrete) {
    err = btocl_load_discreteKernels(rt);
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
    
    //printf("%s \n",options);
    //exit(0);
  }
  
  printf("creating program\n\n");
  err = createProgramFromBuffers(rt->context,&rt->program,rt);
  if (err < 0) {
    printf("Couldn't create program\n");
    return err;
  }
  
  printf(".....building program\n");
  err = clBuildProgram(rt->program,1,&rt->device,options,NULL,NULL);
  if (err != CL_SUCCESS) {
    printf("Couldn't build program\n");
    clGetProgramBuildInfo(rt->program,rt->device,CL_PROGRAM_BUILD_LOG,0,NULL,&log_size);
    program_log=  (char*)calloc(log_size+1,sizeof(char));
    clGetProgramBuildInfo(rt->program,rt->device,CL_PROGRAM_BUILD_LOG,log_size,program_log,NULL);
    printf("Program Log:\n%s\n",program_log);
    free(program_log);
    exit(1);
  }
  printf("success...program built\n");
  
  
  err = loadKernels(rt);
  if (err != CL_SUCCESS) {
    printf("Problem while loading kernel\n");
    return err;
  }
  
  printf("OCL initialization successful\n");
  return CL_SUCCESS;
  

}


int btocl_load_continuousKernels(BTOCL_RUNTIME* rt) {
  cl_int err;
  
  // kernel_index, kernel_name, program file name
  err = load_kernel_file(BTOCL_CHOLUPDCOL_PURE,"btocl_cholupdcol_pure",
		   BTOCL_KERNELS "kernel_cholupdcol_pure.cl",rt);
  if (err != 0) return err;
  err = load_kernel_file(BTOCL_CHOLUPDMAT_P4,"btocl_cholupdmat_p4",
		   BTOCL_KERNELS "kernel_cholupdmat_p4.cl",rt);
  if (err != 0) return err;
  err = load_kernel_file(BTOCL_CHOLUPDMAT_P8,"btocl_cholupdmat_p8",
		   BTOCL_KERNELS "kernel_cholupdmat_p8.cl",rt);
  if (err != 0) return err;
  err = load_kernel_file(BTOCL_CHOLUPDMAT_P16,"btocl_cholupdmat_p16",
		   BTOCL_KERNELS "kernel_cholupdmat_p16.cl",rt);
  if (err != 0) return err;
  err = load_kernel_file(BTOCL_CHOLUPDMAT_P32,"btocl_cholupdmat_p32",
		   BTOCL_KERNELS "kernel_cholupdmat_p32.cl",rt);
  if (err != 0) return err;
  err = load_kernel_file(BTOCL_CHOLUPDMAT_P64,"btocl_cholupdmat_p64",
		   BTOCL_KERNELS "kernel_cholupdmat_p64.cl",rt);
  if (err != 0) return err;
  err = load_kernel_file(BTOCL_CHOLUPDMAT_P128,"btocl_cholupdmat_p128",
		   BTOCL_KERNELS "kernel_cholupdmat_p128.cl",rt);
  if (err != 0) return err;
  err = load_kernel_file(BTOCL_CHOLUPDCOL,"btocl_cholupdcol",
		   BTOCL_KERNELS "kernel_cholupdcol.cl",rt);
  if (err != 0) return err;
  err = load_kernel_file(BTOCL_CHOLUPDMAT,"btocl_cholupdmat",
		   BTOCL_KERNELS "kernel_cholupdmat.cl",rt);
  if (err != 0) return err;

  // Chol - revisited
  err = load_kernel_file(BTOCL_CHOLUPDMAT_PURE,"btocl_cholupdmat_pure",
		   BTOCL_KERNELS "kernel_cholupdmat_pure.cl",rt);
  if (err != 0) return err;
  err = load_kernel_file(BTOCL_TRANSPOSE,"btocl_transpose",
		   BTOCL_KERNELS "kernel_transpose.cl",rt);
  if (err != 0) return err;
  err = load_kernel_file(BTOCL_LTRI_MBYL,"btocl_ltri_MbyL",
		   BTOCL_KERNELS "kernel_ltri_MbyL.cl",rt);
  if (err != 0) return err;
  err = load_kernel_file(BTOCL_LTRI_LBYM,"btocl_ltri_LbyM",
		   BTOCL_KERNELS "kernel_ltri_LbyM.cl",rt);
  if (err != 0) return err;


  // KRONECKER PRODUCT
  err = load_kernel_file(BTOCL_KRON_BASIC,"btocl_kron_basic",
		   BTOCL_KERNELS "kernel_kron_basic.cl",rt);
  if (err != 0) return err;
  err = load_kernel_file(BTOCL_KRON_ACCUM,"btocl_kron_accum",
		   BTOCL_KERNELS "kernel_kron_accum.cl",rt);
  if (err != 0) return err;
  err = load_kernel_file(BTOCL_KRON_LOG1,"btocl_kron_log1",
		   BTOCL_KERNELS "kernel_kron_log1.cl",rt);
  if (err != 0) return err;
  err = load_kernel_file(BTOCL_LTRI_LTBYL,"btocl_ltri_LTbyL",
		   BTOCL_KERNELS "kernel_ltri_LTbyL.cl",rt);
  if (err != 0) return err;
  err = load_kernel_file(BTOCL_WRITEUTOL_DIAG,"btocl_writeUtoL_diag",
		   BTOCL_KERNELS "kernel_writeUtoL_diag.cl",rt);
  if (err != 0) return err;


}

int btocl_load_discreteKernels(BTOCL_RUNTIME* rt) {
  cl_int err;
  // EXPQT
  err = load_kernel_file(BTOCL_EXPQT,"btocl_expqt",
		   BTOCL_KERNELS "kernel_expqt.cl",rt);
  if (err != 0) return err;
  err = load_kernel_file(BTOCL_EXPQT_ROWNOS2,"btocl_expqt_rownos2",
		   BTOCL_KERNELS "kernel_expqt_rownos2.cl",rt);
  if (err != 0) return err;
  err = load_kernel_file(BTOCL_EXPQT_LOCAL,"btocl_expqt_local",
		   BTOCL_KERNELS "kernel_expqt_local.cl",rt);
  if (err != 0) return err;
  err = load_kernel_file(BTOCL_EXPQT_ROWNOS4,"btocl_expqt_rownos4",
		   BTOCL_KERNELS "kernel_expqt_rownos4.cl",rt);
  if (err != 0) return err;

  // EXPT revisited
  err = load_kernel_file(BTOCL_EXPQT_NOTEMP,"btocl_expqt_notemp",
		   BTOCL_KERNELS "kernel_expqt_notemp.cl",rt);
  if (err != 0) return err;
  err = load_kernel_file(BTOCL_EXPQT3,"btocl_expqt3",
		   BTOCL_KERNELS "kernel_expqt3.cl",rt);
  if (err != 0) return err;
  err = load_kernel_file(BTOCL_EXPQT3_LOCAL,"btocl_expqt3_local",
		   BTOCL_KERNELS "kernel_expqt3_local.cl",rt);
  if (err != 0) return err;
  err = load_kernel_file(BTOCL_EXPQT3_LOCALERR,"btocl_expqt3_localerr",
		   BTOCL_KERNELS "kernel_expqt3_localerr.cl",rt);
  if (err != 0) return err;
  err = load_kernel_file(BTOCL_EXPQT3_NOS4,"btocl_expqt3_nos4",
		   BTOCL_KERNELS "kernel_expqt3_nos4.cl",rt);  
  if (err != 0) return err;
  err = load_kernel_file(BTOCL_EXPQT_LOCALERR,"btocl_expqt_localerr",
		   BTOCL_KERNELS "kernel_expqt_localerr.cl",rt);
  if (err != 0) return err;

  // PLH
  err = load_kernel_file(BTOCL_PLH,"btocl_plh",
		   BTOCL_KERNELS "kernel_plh.cl",rt);
  if (err != 0) return err;
  err = load_kernel_file(BTOCL_PLHNODE,"btocl_plhNode",
		   BTOCL_KERNELS "kernel_plhNode.cl",rt);
  if (err != 0) return err;
  err = load_kernel_file(BTOCL_PLHROW,"btocl_plhRow",
		   BTOCL_KERNELS "kernel_plhRow.cl",rt);
  if (err != 0) return err;
  err = load_kernel_file(BTOCL_PLHROWG,"btocl_plhRowG",
		   BTOCL_KERNELS "kernel_plhRowG.cl",rt);
  if (err != 0) return err;
  err = load_kernel_file(BTOCL_PLHROWALLG,"btocl_plhRowAllG",
		   BTOCL_KERNELS "kernel_plhRowAllG.cl",rt);
  if (err != 0) return err;
  err = load_kernel_file(BTOCL_PLHROWGL,"btocl_plhRowGL",
		   BTOCL_KERNELS "kernel_plhRowGL.cl",rt);
  if (err != 0) return err;
  err = load_kernel_file(BTOCL_PLHNODE_NOS2,"btocl_plhNode_nos2",
		   BTOCL_KERNELS "kernel_plhNode_nos2.cl",rt);
  if (err != 0) return err;
  err = load_kernel_file(BTOCL_PLHROWFULLG,"btocl_plhRowFullG",
		   BTOCL_KERNELS "kernel_plhRowFullG.cl",rt);
  if (err != 0) return err;
  err = load_kernel_file(BTOCL_PLHREDUCEGL,"btocl_plhReduceGL",
		   BTOCL_KERNELS "kernel_plhReduceGL.cl",rt);
  if (err != 0) return err;

}



int load_kernel_file(unsigned short kernel_idx, char* kernel_name, char* kernel_file, BTOCL_RUNTIME* rt) {
  cl_ushort idx;
  size_t psize;
  char* pbuffer;
  char* str;
  FILE* program_handle;

  if (kernel_idx >= NUM_KERNELS) {
    printf("Kernel (%s,%s) index out of bounds: %d\n", kernel_name,kernel_file,kernel_idx);
    return kernel_idx;
  }


  str = (char*)malloc(strlen(kernel_file)+1);
  strcpy(str,kernel_file);
  rt->kernel_program_files[kernel_idx] = str;
  // Copy kernel name
  str = (char*)malloc(strlen(kernel_name)+1);
  strcpy(str,kernel_name);
  rt->kernel_names[kernel_idx] = str;
  // Set flag
  rt->kernel_use[kernel_idx] = 1;


  // Load File

  program_handle = fopen(kernel_file,"rb");
  if (program_handle == NULL) {
    printf("Couldn't open file: %s\n", kernel_file);
    return 1;
  }
  fseek(program_handle,0L,SEEK_END); // go to end of file
  psize = ftell(program_handle);
  rewind(program_handle);
  pbuffer=(char*)malloc(psize+1);
  fread(pbuffer,sizeof(char),psize,program_handle);
  fclose(program_handle);
  pbuffer[psize] = '\0';
  // Update rt - note that we don't use kernel_idx 
  idx = rt->kernel_num_used;  
  rt->program_buffer[idx] = pbuffer;
  rt->program_size[idx] = psize;
  rt->kernel_num_used++;  

  return 0; // success

}

cl_int createProgramFromBuffers(cl_context context, cl_program* program,BTOCL_RUNTIME* rt) {
  cl_int err;
  // Create program
  // Problem: all elements of program buffer must be initialised

  *program = clCreateProgramWithSource(context,rt->kernel_num_used,
				       (const char**)rt->program_buffer,
				       rt->program_size,
				       &err);
  if (err != CL_SUCCESS) {
    printf("Error while creating program from sources\n");
    return err;
  }

  return CL_SUCCESS;

}

cl_int freeProgramBuffers(BTOCL_RUNTIME* rt) {
  int i;
  for(i=0; i < NUM_KERNELS; i++) {
    if (rt->kernel_use[i])
      free(rt->program_buffer[i]);
  }
  return CL_SUCCESS;

}

// Re-used
cl_int loadKernels(BTOCL_RUNTIME* rt) {
  int i;
  cl_int err;

  // create kernel
  for(i=0; i < NUM_KERNELS; i++) {
    if (rt->kernel_use[i]) {
      printf("creating kernel %s ...\n",rt->kernel_names[i]);
      rt->kernels[i] = clCreateKernel(rt->program, rt->kernel_names[i],&err);
      if (err < 0) {
	printf("Couldn't create kernel:%s \n", rt->kernel_names[i]);
	return err;
      }
    }
  }
  //printf("kernel created\n");
  return CL_SUCCESS;
}

// Test - check later
cl_int createProgramVerbatim(cl_context context,cl_program *program) {

	cl_int err;
	const char *program_buffer =
	"__kernel                                                                    \n"
	" void matvec_test(__global float4* matrix,	 __global float4* result) {       \n"
	"    int i = get_global_id(0);                                                \n"
	"    float4 x = (float4)((i+5)*1.0f);                                         \n"
	"    result[i] = x;                                                           \n "
	" }                                                                           \n"
	;

	*program = clCreateProgramWithSource(context,1,(const char**)&program_buffer,NULL,&err);

	
	return err;
}









#endif // if BTOCL defined
