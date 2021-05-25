#ifndef BTOCL_KERNELS_DISCRETE_H
#define BTOCL_KERNELS_DISCRETE_H

#ifdef BTOCL


const char* STRING_AUX_BEXP = 
"     \n"
"#pragma OPENCL EXTENSION cl_khr_fp64: enable  \n"
"     \n"
"#ifdef USE_BTEXP  \n"
"#define FEXP btexp  \n"
"#else  \n"
"#define FEXP exp  \n"
"#endif   \n"
"  \n"
"#define C1 (*(double *)sc1)  \n"
"#define C2 (*(double *)sc2)  \n"
"#define P0  (*(double*)P)  \n"
"#define P1  (*((double*)P + 1))  \n"
"#define P2  (*((double*)P + 2))  \n"
"#define Q0  (*(double*)Q)  \n"
"#define Q1  (*((double*)Q + 1))  \n"
"#define Q2  (*((double*)Q + 2))  \n"
"#define Q3  (*((double*)Q + 3))  \n"
"  \n"
"double btexp(double x) {  \n"
"       double px,xx,poleval;  // the exponent of e**x  \n"
"      int n;  \n"
"       unsigned short sc1[] = {0x0000,0x0000,0x2e40,0x3fe6};  \n"
"      unsigned short sc2[] = {0xabca,0xcf79,0xf7d1,0x3eb7};  \n"
"       unsigned short P[] = {  \n"
"          0x4be8,0xd5e4,0x89cd,0x3f20,  \n"
"          0x2c7e,0x0cca,0x06d1,0x3f9f,  \n"
"          0x0000,0x0000,0x0000,0x3ff0,  \n"
"        };  \n"
"        unsigned short Q[] = {  \n"
"                 0x5fa0,0xbc36,0x2eb6,0x3ec9,  \n"
"                 0xb6c0,0xb508,0xae39,0x3f64,  \n"
"                 0xe074,0x9887,0x1709,0x3fcd,  \n"
"                 0x0000,0x0000,0x0000,0x4000,  \n"
"        };  \n"
"  \n"
"  \n"
"       // ********** E^x ************  \n"
"  \n"
"        px = floor(x/log(2.0) + 0.5);  \n"
"        n = px;  \n"
"  \n"
"        x -= px * C1;  \n"
"        x -= px * C2;  \n"
"        xx = x * x;  \n"
"  \n"
"        poleval = P0*xx;  \n"
"        poleval = (poleval+P1)*xx;  \n"
"        poleval += P2;  \n"
"  \n"
"        px = x * poleval;  \n"
"	  \n"
"        poleval = Q0*xx;  \n"
"        poleval = (poleval+Q1)*xx;  \n"
"        poleval = (poleval+Q2)*xx;  \n"
"        poleval += Q3;   \n"
"                         \n"
"        x = px / (poleval - px);  \n"
"        x = 1.0 + 2.0 * x;  \n"
"                         \n"
"        x = ldexp(x,n);  \n"
"                    \n"
"        return x;   \n"
"                    \n"
"}                   \n"
;


const char* STRING_EXPQT_ROWNOS2 = 
"  \n"
"#pragma OPENCL EXTENSION cl_khr_fp64: enable  \n"
"  \n"
"__kernel void btocl_expqt_rownos2(__global double2* V, __global double2* lambdas, __global double2* invV,  \n"
"                          __global double* vect_t, __global int* vect_id,  \n"
"						  __global double2* pmatrix, int nonodes, int num_nodes_wg,  \n"
"						  __global int* error  \n"
"						  ) {  \n"
"		  \n"	
"	int id, row;  \n"
"  \n"
"	double2 t01, l;  \n"
"	double s;  \n"
"	  \n"
"	id = get_global_id(0);  \n"
"	row = id % 2;  \n"
"	  \n"
"	id = id / 2;  // node id  \n"
"	if ((id >= nonodes)||   // dummy worker used to complete workgroup  \n"
"		(id == 0)) {        // skip root - it will set the error flag to true  \n"
"		return;			  \n"
"	}  \n"
"			  \n"
"	id = vect_id[id];  // real node id  \n"
"		  \n"
"	l = vect_t[id]*lambdas[0];  \n"
"	//l = FEXP(l);  \n"
"	l.s0 = FEXP(l.s0);  \n"
"	l.s1 = FEXP(l.s1);  \n"
"	t01 = l * V[row];  \n"
"	  \n"
"	// Two dot products (using invV transpose)  \n"
"	l.s0 = dot(t01,invV[0]); //  t0*invV[0] +   t1*invV[2];  \n"
"	l.s1 = dot(t01,invV[1]);   // t0*invV[1] +   t1*invV[3];  \n"
"	  \n"
"	pmatrix[2*id+row] = l;  // could be initial id if we remove vect_id  \n"
"	  \n"
"	if ((l.s0 < -0.000001) || (l.s1 < -0.000001)) {  \n"
"		error[0] = 1;  \n"
"		return;  \n"
"	}  \n"
"	s = 1 - (l.s0+l.s1);  \n"
"	if (s > 0.000001) {  \n"
"		error[0] = 1;  \n"
"		return;  \n"
"	} \n"
"	 \n"
"	return; \n"
"} \n"
;

const char* STRING_EXPQT_ROWNOS4 =
" \n"
"#pragma OPENCL EXTENSION cl_khr_fp64: enable \n"
" \n"
"__kernel void btocl_expqt_rownos4(__global double4* V, __global double4* lambdas, \n"
"						  __global double4* invV, \n"
"                         __global double* vect_t, __global int* vect_id, \n"
"						  __global double4* pmatrix, int nonodes, int num_nodes_wg, \n"
"						  __global int* error \n"
"						  ) { \n"
"			 \n"
"	int id, row; \n"
" \n"
"	double4 t03, l; \n"
"	double s; \n"
"		 \n"
"	id = get_global_id(0); \n"
"	row = id % 4; \n"
"	 \n"
"	id = id / 4;  // node id \n"
"	if (id >= nonodes) {  // dummy worker used to complete workgroup \n"
"		return;		 \n"	
"	} \n"
"	 \n"
"	id = vect_id[id];  // real node id \n"
"	 \n"
"	l = vect_t[id]*lambdas[0]; \n"
" \n"
"	l.s0 = FEXP(l.s0); \n"
"	l.s1 = FEXP(l.s1); \n"
"	l.s2 = FEXP(l.s2); \n"
"	l.s3 = FEXP(l.s3); \n"
"	t03 = l * V[row]; \n"
"	 \n"
"	// Four dot products (using transpose here) \n"
"	l.s0 = dot(t03,invV[0]);  \n"
"	l.s1 = dot(t03,invV[1]);   \n"
"	l.s2 = dot(t03,invV[2]); \n"
"	l.s3 = dot(t03,invV[3]); \n"
"	 \n"
"	pmatrix[4*id+row] = l;  // a single copy \n"
"		 \n"
"	if ((l.s0 < -0.000001) || (l.s1 < -0.000001) || \n"
"		(l.s2 < -0.000001) || (l.s3 < -0.000001)) { \n"
"		error[0] = 1; \n"
"		return; \n"
"	} \n"
"	s = 1 - (l.s0+l.s1+l.s2+l.s3); \n"
"	if (s > 0.000001) { \n"
"		error[0] = 1; \n"
"		return; \n"
"	} \n"
"	 \n"
"	return; \n"
"} \n"
;

const char* STRING_EXPQT3_LOCALERR = 
" \n"
"#pragma OPENCL EXTENSION cl_khr_fp64: enable \n"
" \n"
"__kernel void btocl_expqt3_localerr(__global double* V, __global double* lambdas, __global double* invV, \n"
 "                         __global double* vect_t, __global int* vect_id, \n"
"						  __global double* pmatrix, __global int* error, int num_cells, \n"
"						  __local double* local_temp) { \n"
"				 \n"		  
"	int id, local_idx; \n"
"	int i; \n"
"	double prodsum, t, temp_exp; \n"
"	int v_idx, inv_idx; \n"
"	int temp_start, temp_idx; \n"
" \n"
"	id =  get_global_id(0); \n"
"		 \n"
"	if (id >= num_cells)  // there are more workitems than cells due to extra rows \n"
"		return; \n"
"		 \n"
"	local_idx = id % BT_NOS2;  // position inside pmatrix (row,col)\n"
"	id = id / BT_NOS2;  // which pmatrix?\n"
"	id++;  // shift one, skip root node\n"
"	\n"
"	id = vect_id[id];  // real id\n"
"	\n"
"	// temp indexes - synchronised with local id\n"
"	temp_idx = get_local_id(0);\n"
"	temp_start = BT_NOS * (temp_idx / BT_NOS);\n"
"  \n"
"	i = local_idx % BT_NOS;  \n"
"	local_temp[temp_idx] = FEXP(lambdas[i]*vect_t[id])*V[local_idx]; \n"
"	 \n"
"	barrier(CLK_GLOBAL_MEM_FENCE); \n"
"			 \n"
"	// Matrix multiplication \n"
"	prodsum = 0; \n"
"	//v_idx = BT_NOS*(local_idx / BT_NOS);  // nos*row \n"
"	inv_idx =  local_idx % BT_NOS;    // col \n"
"	temp_idx = temp_start; \n"
"	for(i=0; i < BT_NOS; i++) {		 \n"
"		prodsum = fma(local_temp[temp_idx++],invV[inv_idx],prodsum);  \n"
"		//v_idx++; \n"
"		inv_idx +=BT_NOS; \n"
"	} \n"
"		 \n"
"	pmatrix[BT_NOS2*id+local_idx] = prodsum; \n"
"	 \n"
"	barrier(CLK_GLOBAL_MEM_FENCE); \n"
"	 \n"
"	local_temp[get_local_id(0)] = prodsum; \n"
"	 \n"
"	//barrier(CLK_GLOBAL_MEM_FENCE);  -- changed this due to problems in expqt_localerr \n"
"	barrier(CLK_LOCAL_MEM_FENCE); \n"
"	 \n"
"	// error checking \n"
"	// Need to pass extra argument for error message \n"
"	if (prodsum < -0.000001) { \n"
"		error[0] = 1; \n"
"		return; \n"
"	} \n"
"	// sum the errors \n"
"	temp_idx = temp_start; \n"
"	prodsum=0; \n"
"	for(i=0; i < BT_NOS; i++) \n"
"		prodsum += local_temp[temp_idx++];  \n"
"	prodsum = prodsum - 1.0; \n"
"	if (prodsum*prodsum > 0.00001) { \n"
"		error[0]= 1; \n"
"	} \n"
"	 \n"
"} \n"
;




#endif // if BTOCL defined

#endif