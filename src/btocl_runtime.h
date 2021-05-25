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

  // Info selected GPU
   

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

// OpenCL functions plus error checking
//cl_mem clCreateBuffer (cl_context context, cl_mem_flags flags, size_t size, void *host_ptr, cl_int *errcode_ret)
cl_mem btocl_clCreateBuffer (const char* msg, cl_context context, cl_mem_flags flags, size_t size, void *host_ptr, cl_int *errcode_ret);
void btocl_printBufferInfo(cl_mem buffer);

// Error related
void btocl_printRuntimeError(cl_int err);
void btocl_checkErrorCode(cl_int error,const char* operation, const char* object);

// ------------------------------------------------------------------------------------

// Implemented by btocl_runtime_kernels.c

void initialise_kernel_info();

int load_kernel_file(cl_ushort* kernel_idx, char* kernel_name, char* kernel_file, BTOCL_RUNTIME* rt);
int load_kernel_string(cl_ushort* kernel_idx, char* kernel_name, const char* kernel_str, BTOCL_RUNTIME* rt);
int load_aux_file(char* kernel_file, BTOCL_RUNTIME* rt);
int load_aux_string(const char* aux_str, BTOCL_RUNTIME* rt);

cl_int createProgramFromBuffers(cl_context context, cl_program* program, BTOCL_RUNTIME* rt);
cl_int freeProgramBuffers(BTOCL_RUNTIME* rt);
cl_int loadKernels(BTOCL_RUNTIME* rt);
cl_int createProgramVerbatim(cl_context context,cl_program *program);

// new
cl_int btocl_build_load_kernels(const char* options);
cl_int btocl_clear_kernels(BTOCL_RUNTIME* rt);

#endif // if BTOCL defined

#endif
