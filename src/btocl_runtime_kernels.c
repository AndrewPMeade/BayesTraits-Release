#ifdef BTOCL

#include <stdio.h>
#include <stdlib.h>
#include <string.h>


#include "btocl_runtime.h"

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


int load_kernel_file(cl_ushort *pkernel_idx,char* kernel_name, char* kernel_file, BTOCL_RUNTIME* rt) {
  cl_ushort kernel_idx;
  size_t psize;
  char* pbuffer;
  char* str;
  FILE* program_handle;

  // Next index available == number of allocated kernels and aux_files
  kernel_idx = rt->kernel_num_used;  
  if (kernel_idx >= NUM_KERNELS) {
    printf("Number of kernels exceeds allocated space (%d). Kernel %s in file %s\n",kernel_idx,kernel_name,kernel_file);
    return 0;  // error
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
 
  rt->program_buffer[kernel_idx] = pbuffer;
  rt->program_size[kernel_idx] = psize;
  rt->kernel_num_used++;  

  // return index in array - uniquely identifies kernel in array
  *pkernel_idx = kernel_idx;
  return 0;

}

int load_kernel_string(cl_ushort *pkernel_idx,char* kernel_name, const char* kernel_string, BTOCL_RUNTIME* rt) {
  cl_ushort kernel_idx;
  size_t psize;
  char* pbuffer;
  char* str;
  FILE* program_handle;
  char* kernel_file = "(string kernel)";

  // Next index available == number of allocated kernels and aux_files
  kernel_idx = rt->kernel_num_used;  
  if (kernel_idx >= NUM_KERNELS) {
    printf("Number of kernels exceeds allocated space (%d). Kernel %s in file %s\n",kernel_idx,kernel_name,kernel_file);
    return 0;  // error
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


 
 // use string instead of file
  pbuffer = (char*)malloc(strlen(kernel_string)+1);
  strcpy(pbuffer,kernel_string);
  rt->program_buffer[kernel_idx] = pbuffer;
  rt->program_size[kernel_idx] = strlen(kernel_string+1);
  rt->kernel_num_used++;  

  // return index in array - uniquely identifies kernel in array
  *pkernel_idx = kernel_idx;
  return 0;

}





int load_kernel_file_old(unsigned short kernel_idx, char* kernel_name, char* kernel_file, BTOCL_RUNTIME* rt) {
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

int load_aux_file(char* aux_file, BTOCL_RUNTIME* rt) {
	FILE* program_handle;
	cl_ushort idx;
	size_t psize;
	char* pbuffer;


	program_handle = fopen(aux_file,"rb");
	if (program_handle == NULL) {
		printf("Couldn't open file: %s\n", aux_file);
		return 1;
	}
	fseek(program_handle,0L,SEEK_END); // go to end of file
	psize = ftell(program_handle);
	rewind(program_handle);
	pbuffer=(char*)malloc(psize+1);
	fread(pbuffer,sizeof(char),psize,program_handle);
	fclose(program_handle);
	pbuffer[psize] = '\0';
	idx = rt->kernel_num_used;  
	rt->program_buffer[idx] = pbuffer;
	rt->program_size[idx] = psize;
	rt->kernel_use[idx] = 0;  // this is not a kernel
	rt->kernel_num_used++;  

	return 0; // success
}

int load_aux_string(const char* aux_string, BTOCL_RUNTIME* rt) {
	FILE* program_handle;
	cl_ushort idx;
	size_t psize;
	char* pbuffer;

	pbuffer = (char*)malloc(strlen(aux_string)+1);
	strcpy(pbuffer,aux_string);

	idx = rt->kernel_num_used;  
	rt->program_buffer[idx] = pbuffer;
	rt->program_size[idx] = psize;
	rt->kernel_use[idx] = 0;  // this is not a kernel
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
  for(i=0; i < rt->kernel_num_used; i++) {
		free(rt->program_buffer[i]);
  }
  return CL_SUCCESS;

}

// Re-used
cl_int loadKernels(BTOCL_RUNTIME* rt) {
  int i;
  cl_int err;

  // create kernel
  for(i=0; i < rt->kernel_num_used; i++) {
    if (rt->kernel_use[i]) { // Is this a kernel? or aux_file
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


cl_int btocl_build_load_kernels(const char* options) {

  size_t log_size;
  char* program_log;
  cl_int err;

  BTOCL_RUNTIME* rt = btocl_getruntime();

  //printf("creating program\n\n");
  err = createProgramFromBuffers(rt->context,&rt->program,rt);
  if (err < 0) {
    printf("Error: Couldn't create OpenCL program\n");
    //return err;
	exit(0);
  }
  
  // Program created, free proagram buffers
   freeProgramBuffers(rt);

  //printf(".....building program\n");
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
  //printf("success...program built\n");
  
  
  err = loadKernels(rt);
  if (err != CL_SUCCESS) {
    printf("Problem while loading kernel\n");
    return err;
  }
  
  printf("OCL initialization successful\n");
  return CL_SUCCESS;
  
}

// Release Program, kernel and related structures
cl_int btocl_clear_kernels(BTOCL_RUNTIME* rt) {
	int i;
	cl_int err;


	// release all kernels
	for(i=0; i < rt->kernel_num_used; i++) {
		if (rt->kernel_use[i]) { // Is this a kernel or aux file?
			err = clReleaseKernel(rt->kernels[i]);
			if (err != CL_SUCCESS) {
				printf("Couldn't release kernel:%s \n", rt->kernel_names[i]);
				return err;
			}
			free(rt->kernel_program_files[i]);
			rt->kernel_program_files[i] = NULL;
			free(rt->kernel_names[i]);
			rt->kernel_names[i] = NULL;
		}
		 rt->kernel_use[i]=0;
	}

	// Release Program
	if(rt->program != NULL) {
		err = clReleaseProgram(rt->program);
		rt->program = NULL;
		if (err != CL_SUCCESS) {
			printf("Error: couldn't release program\n");
			return err;
		}
	}

	 rt->kernel_num_used = 0;

	return CL_SUCCESS;
}



#endif // if BTOCL defined
