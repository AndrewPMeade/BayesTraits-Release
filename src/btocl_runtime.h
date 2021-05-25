#ifndef BTOCL_RUNTIME_H
#define BTOCL_RUNTIME_H

#ifdef BTOCL


#ifdef MAC
#include "OpenCL/cl.h"
#else
#include "CL/cl.h"
#endif

// OJO - update when needed
#define NUM_KERNELS 50


// The OCL runtime structure

typedef struct {
  char* kernel_dir;
  cl_platform_id* platforms;
  cl_uint platform;
  cl_uint num_platforms;

  cl_device_type device_type;
   cl_device_id device;
   cl_context context;
   cl_program program;

   cl_ushort kernel_use[NUM_KERNELS];
  cl_ushort kernel_num_used;
  char* kernel_program_files[NUM_KERNELS];
   cl_kernel kernels[NUM_KERNELS];
   char* kernel_names[NUM_KERNELS];
   cl_command_queue queue;
  // Program Buffer
  char* program_buffer[NUM_KERNELS];
  size_t program_size[NUM_KERNELS];
   

   cl_int err;
} BTOCL_RUNTIME;

BTOCL_RUNTIME* btocl_getruntime();
void btocl_free_runtime();

cl_int btocl_init_runtime(cl_device_type device_type);

cl_device_id * btocl_getDeviceID_ptr();
cl_context btocl_getContext();
cl_command_queue btocl_getCommandQueue();
cl_kernel btocl_getKernel(cl_ushort);
char* btocl_getKernelName(cl_ushort);

void btocl_printBufferInfo(cl_mem buffer);
void btocl_printRuntimeError(cl_int err);

#endif // if BTOCL defined

#endif
