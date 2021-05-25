#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#ifdef BTOCL

#include "btocl_runtime.h"

void printPlatformsInfo(cl_platform_id* platforms, cl_uint num_platforms);
void printDeviceInfo(cl_device_id device);
void printDeviceInfoAll(cl_device_id device);
cl_uint getFirstDevice(cl_platform_id* platforms, cl_uint num_platforms, 
						cl_device_type dt, cl_device_id *device);
cl_bool deviceSupportsFP64(cl_device_id *device);
cl_int deviceIsLittleEndian(cl_device_id *device, cl_bool *littleEndian);
void printDeviceVectorWidths(cl_device_id *device);


BTOCL_RUNTIME* rt = NULL;

/* Interface */

BTOCL_RUNTIME* btocl_getruntime() {
  return rt;
}

cl_int btocl_init_runtime(cl_device_type device_type) {
  cl_int err;
  size_t log_size;
  char* program_log;
  int i;
  
  if (device_type == CL_DEVICE_TYPE_CPU) {
	printf("Initializing OpenCL with CPU\n");
  } else if (device_type == CL_DEVICE_TYPE_GPU) {
	printf ("Initializing OpenCL with GPU\n");
  } else {
	printf("Incorrect device type\n");
	return 1;
  }
  
  if (rt != NULL) {
    return CL_SUCCESS;
  }

  rt = (BTOCL_RUNTIME*)malloc(sizeof(BTOCL_RUNTIME));
  
  rt->kernel_dir = NULL;
  rt->kernel_dir = getenv("BTOCL_KERNEL_DIR");

  /* Access the first installed platform */
  /* Changed: Max of 4 platforms considered */
   err = clGetPlatformIDs(4, NULL, &rt->num_platforms);
   if(err < 0) {
      perror("Couldn't find any platforms");
      return 1;
   }
   rt->platforms = (cl_platform_id*)malloc(sizeof(cl_platform_id) * rt->num_platforms);
   clGetPlatformIDs(rt->num_platforms, rt->platforms, NULL);
   
   printPlatformsInfo(rt->platforms,rt->num_platforms);
	
   //rt->gpuPlatform = getFirstDevice(rt->platforms,rt->num_platforms,device_type,&rt->device);
   rt->platform = getFirstDevice(rt->platforms,rt->num_platforms,device_type,&rt->device);
	
   if (rt->platform >= rt->num_platforms) {
     printf("Couldn't find device of indicated type\n");
     return 1;
   } else {
     printf("Found device type in Platform %d\n",rt->platform);
     printDeviceInfoAll(rt->device);
   }
	
   rt->context = clCreateContext(NULL,1,&rt->device,NULL,NULL,&err);
	
   // Create queue
   printf("Creating queue\n");
   rt->queue = clCreateCommandQueue(rt->context,rt->device,0,&err);
   if (err < 0) {
     printf("Couldn't create command queue\n");
     return err;
   }
   printf("queue created\n");
   
   // Prepares kernel structures
   initialise_kernel_info();
   
   return CL_SUCCESS;
}




void btocl_free_runtime () {
	int i;

	for (i=0; i < NUM_KERNELS; i++)
		if (rt->kernel_use[i])
			if (rt->kernels[i] != NULL) clReleaseKernel(rt->kernels[i]);
	if (rt->queue != NULL) clReleaseCommandQueue(rt->queue);
	if(rt->program != NULL) clReleaseProgram(rt->program);
	if(rt->context != NULL) clReleaseContext(rt->context);
	if (rt->platforms != NULL) (rt->platforms);	
}


/* ******** Auxiliary functions **********  */

void printPlatformsInfo(cl_platform_id* platforms, cl_uint num_platforms) {
	cl_uint i, j;
	size_t prop_size;
	cl_int err;
	char *plat_name, *plat_vendor;
	cl_uint num_devices;
	cl_device_id *devices;

	printf("Num platforms: %d\n",num_platforms);

	for(i=0; i < num_platforms; i++) {
		err = clGetPlatformInfo(platforms[i], CL_PLATFORM_NAME, 0, NULL, &prop_size);
		if(err < 0) {
			perror("Couldn't read platform name.");
			exit(1);
		}
		plat_name = (char*)malloc(prop_size);
		clGetPlatformInfo(platforms[i],CL_PLATFORM_NAME,prop_size, plat_name, NULL);

		err = clGetPlatformInfo(platforms[i], CL_PLATFORM_VENDOR, 0, NULL, &prop_size);
		if(err < 0) {
			perror("Couldn't read platform vendor.");
			exit(1);
		}
		plat_vendor = (char*)malloc(prop_size);
		clGetPlatformInfo(platforms[i],CL_PLATFORM_VENDOR,prop_size, plat_vendor, NULL);

		printf("Platform: %s, %s\n", plat_name,plat_vendor);

		free(plat_name);
		free(plat_vendor);
 
		num_devices=0;
		err = clGetDeviceIDs(platforms[i], CL_DEVICE_TYPE_ALL, 0, NULL, &num_devices);
		
		if(err < 0) {
			printf("error: ");
			switch (err) {
			case CL_INVALID_PLATFORM:
				printf("Platform invalid\n"); break;
			case CL_INVALID_DEVICE_TYPE:
				printf("invalid device type\n"); break;
			case CL_INVALID_VALUE:
				printf("Invalid value\n"); break;
			case CL_DEVICE_NOT_FOUND:
				printf("No devices found\n"); break;
			default:
				printf("Unknown error\n");
			}
		} else {
			printf("Number of devices: %d\n",num_devices);
		}
		
		devices = (cl_device_id*)malloc(sizeof(cl_device_id)*num_devices);
		err = clGetDeviceIDs(platforms[i], CL_DEVICE_TYPE_ALL, num_devices, devices, NULL);

		for(j=0; j < num_devices;j++) {
			//printDeviceInfo(devices[j]);
			printDeviceInfoAll(devices[j]);
		}

		
		free(devices);
		
	}

}

void printDeviceInfo(cl_device_id device) {
	size_t prop_size;
	char *dev_name;
	cl_device_type dev_type;
	clGetDeviceInfo(device,CL_DEVICE_NAME,0,NULL,&prop_size);
	dev_name = (char*)malloc(prop_size);
	clGetDeviceInfo(device,CL_DEVICE_NAME,prop_size,dev_name,NULL);
	clGetDeviceInfo(device,CL_DEVICE_TYPE,sizeof(dev_type),&dev_type,NULL);

	printf("-- Name: %s\n", dev_name);
	printf("   Device type: ");
	switch(dev_type) {
	case CL_DEVICE_TYPE_ACCELERATOR:  printf("Accelerator"); break;
	case CL_DEVICE_TYPE_GPU:  printf("GPU"); break;
	case CL_DEVICE_TYPE_CPU:  printf("CPU"); break;
//	case CL_DEVICE_TYPE_CUSTOM:  printf("Custom"); break;  -- not rsupport by NVIDIA headers 2008-2010
	case CL_DEVICE_TYPE_DEFAULT:  printf("default"); break;
	case CL_DEVICE_TYPE_ALL:  printf("all"); break;
	default: printf("Unknown");
	}
	printf("\n");
	free(dev_name);
	return;
}


void printDeviceInfoAll(cl_device_id device) {
	cl_uint dev_units, dev_max_dimensions;   // dev_ref_count=1;
	cl_ulong dev_global_mem, dev_local_mem, dev_max_alloc_size;
	size_t max_workgroup;
	cl_bool isLittleEndian, compiler_available;

	printDeviceInfo(device);

	dev_max_dimensions = 0;
	// clGetDeviceInfo(device,CL_DEVICE_REFERENCE_COUNT,sizeof(dev_ref_count),&dev_ref_count,NULL);
	// compiler couldn't find CL_DEVICE_REFERENCE_COUNT
	clGetDeviceInfo(device,CL_DEVICE_COMPILER_AVAILABLE,sizeof(compiler_available),&compiler_available,NULL);

	clGetDeviceInfo(device,CL_DEVICE_GLOBAL_MEM_SIZE,sizeof(dev_global_mem),&dev_global_mem,NULL);
	clGetDeviceInfo(device,CL_DEVICE_MAX_MEM_ALLOC_SIZE,sizeof(dev_max_alloc_size),&dev_max_alloc_size,NULL);
	clGetDeviceInfo(device,CL_DEVICE_LOCAL_MEM_SIZE,sizeof(dev_local_mem),&dev_local_mem,NULL);

	clGetDeviceInfo(device,CL_DEVICE_MAX_WORK_ITEM_DIMENSIONS,sizeof(dev_max_dimensions),&dev_max_dimensions,NULL);

	clGetDeviceInfo(device,CL_DEVICE_MAX_COMPUTE_UNITS,sizeof(dev_units),&dev_units,NULL);
	clGetDeviceInfo(device,CL_DEVICE_MAX_WORK_GROUP_SIZE,sizeof(max_workgroup),&max_workgroup,NULL);

	printf("  * Compute units: %d\n", dev_units);
	printf("  * Compiler available: %d\n", compiler_available);
//	printf("  * Reference Count: %d\n",dev_ref_count);
	printf("  * Global memory %lu bytes\n",(unsigned long)dev_global_mem);
	printf("  * Global memory Max Alloc Size %lu bytes\n",(unsigned long)dev_max_alloc_size);
	printf("  * Local memory %lu bytes\n",(unsigned long)dev_local_mem);
	printf("  * Max work item dimensions: %d\n",dev_max_dimensions);
	printf(" * Max work-group size: %lu\n",(unsigned long)max_workgroup);

	printf("  * FP64 ");
	if (deviceSupportsFP64(&device))
		printf("supported\n");
	else
		printf("NOT supported\n");

	if (deviceIsLittleEndian(&device,&isLittleEndian) != CL_SUCCESS)
		printf("Error while finding endianness\n");
	if (isLittleEndian)
		printf("* Device is Little Endian\n");
	else
		printf("* Device is Big Endian: %d\n",isLittleEndian);
	
	printDeviceVectorWidths(&device);
	
	printf("-------------------------------------------\n");

}

cl_uint getFirstDevice(cl_platform_id* platforms, cl_uint num_platforms, cl_device_type device_type, cl_device_id *device) {
	cl_uint i;
	cl_int err;
	cl_uint num_devices=0;
	for(i=0; i < num_platforms;i++) {
	  err = clGetDeviceIDs(platforms[i],device_type, 1, device, &num_devices);
	  printf("Devices found %d\n",num_devices);
	  if (err >= 0)
	    return i;
	}
	return num_platforms;
}



cl_int deviceIsLittleEndian(cl_device_id *device, cl_bool *littleEndian) {
	return clGetDeviceInfo(*device,CL_DEVICE_ENDIAN_LITTLE,sizeof(littleEndian),littleEndian,NULL);
}

cl_bool deviceSupportsFP64(cl_device_id *device) {
	char extensions[1024];
	size_t length;
	cl_int err;
	char* fp64;
	err = clGetDeviceInfo(*device,CL_DEVICE_EXTENSIONS,0,NULL,&length);
	if (err != CL_SUCCESS) {
		printf("Error while getting extensions\n");
		return CL_FALSE;
	} 
	
	clGetDeviceInfo(*device,CL_DEVICE_EXTENSIONS,sizeof(extensions),extensions,NULL);
	
	printf("Extensions: %s\n",extensions);
	fp64 = strstr(extensions,"fp64");
	
	return (fp64 != NULL);
}

void printDeviceVectorWidths(cl_device_id *device) {
	cl_uint width[7];
	int i;
	char* types[] = { "char", "short", "int", "long", "float", "half", "double" };
	clGetDeviceInfo(*device, CL_DEVICE_PREFERRED_VECTOR_WIDTH_CHAR,sizeof(cl_uint), &width[0], NULL);
	clGetDeviceInfo(*device, CL_DEVICE_PREFERRED_VECTOR_WIDTH_SHORT,sizeof(cl_uint), &width[1], NULL);
	clGetDeviceInfo(*device, CL_DEVICE_PREFERRED_VECTOR_WIDTH_INT,sizeof(cl_uint), &width[2], NULL);
	clGetDeviceInfo(*device, CL_DEVICE_PREFERRED_VECTOR_WIDTH_LONG,sizeof(cl_uint), &width[3], NULL);
	clGetDeviceInfo(*device, CL_DEVICE_PREFERRED_VECTOR_WIDTH_SHORT,sizeof(cl_uint), &width[4], NULL);
	clGetDeviceInfo(*device, CL_DEVICE_PREFERRED_VECTOR_WIDTH_HALF,sizeof(cl_uint), &width[5], NULL);
	clGetDeviceInfo(*device, CL_DEVICE_PREFERRED_VECTOR_WIDTH_LONG,sizeof(cl_uint), &width[6], NULL);
	printf("* Preferred vector widths \n");
	for(i = 0; i < 7; i++)
		printf("(%s: %u) ",types[i],width[i]);
	printf("\n");
}



cl_device_id * btocl_getDeviceID_ptr() {
	return &rt->device;
}

cl_context btocl_getContext() {
	return rt->context;
}

cl_command_queue btocl_getCommandQueue() {
	return rt->queue;
}

cl_kernel btocl_getKernel(cl_ushort index) {
  if (index >= NUM_KERNELS) {   
    rt->err = CL_INVALID_ARG_INDEX;
    return NULL;
  } else {
    rt->err = CL_SUCCESS;
    return rt->kernels[index];
  }
}

char* btocl_getKernelName(cl_ushort index) {
  if (index >= NUM_KERNELS) {
    rt->err = CL_INVALID_ARG_INDEX;
    return NULL;
  } else {
    rt->err = CL_SUCCESS;
    return rt->kernel_names[index];
  }
}

void btocl_printBufferInfo(cl_mem buffer) {
	size_t buffer_size;
	clGetMemObjectInfo(buffer, CL_MEM_SIZE,sizeof(buffer_size),&buffer_size,NULL);
	printf("Buffer Size: %lu\n",buffer_size);
	// add your favourite stats here
	
	return;
}

void btocl_printRuntimeError(cl_int err) {
	printf("Error (%d): ",err);
	switch (err) {
		case CL_INVALID_CONTEXT: printf("Invalid context\n"); break;
		case CL_INVALID_COMMAND_QUEUE: printf("invalid command queue\n"); break;
		case CL_INVALID_MEM_OBJECT: printf("invalid mem object\n");break;
		case CL_INVALID_VALUE: printf("invalid value\n"); break;
		case CL_INVALID_EVENT_WAIT_LIST: printf("event wait list\n"); break;
		case CL_MEM_OBJECT_ALLOCATION_FAILURE: printf("object allocation failure\n"); break;
		case CL_OUT_OF_HOST_MEMORY: printf("out of host memory\n"); break;
		case CL_OUT_OF_RESOURCES: printf("out of resources\n"); break;
		case CL_INVALID_OPERATION: printf("invalid operation\n"); break;
		default: printf("unknown\n");
	}
	return;
}

#endif // if BTOCL defined
