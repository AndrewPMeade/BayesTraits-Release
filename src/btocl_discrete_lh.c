/*
*  BayesTriats 4.0
*
*  copyright 2022
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
#include <math.h>
#include <string.h>

#include "TypeDef.h"
#include "Trees.h"
#include "Continuous.h"
#include "GenLib.h"
#include "Matrix.h"
// OpenCL headers



#include "btocl_discrete.h"
#include "btlin_basic.h"

typedef struct {
	// Maximum dimension is 3 - we mainly use 1
	size_t globalws[3];
	size_t localws[3];  // number of wokitems per dimension
	int dim;  //determines actual length of globalws and localws
	int usewg;
	int num_nodes_wg;  // number of nodes per workgroup
	int num_rows_wg;   // if workgroups must group by rows
	size_t num_wg;
	int check_error;  // indicates of kernel checks pmatrix error
	int transpose_inv;
	cl_ushort kernel_type;
	// Params info for PartialLh computation
	int action; // 0=Collect, 1=Collect/Compute 2=Compute (tips)
	int passedCut;
	int lastGroup;
	int num_nodes;
} KernelInfo_lh;

/* ******* Auxiliary functions ******* */
void calculate_height(TREE* Tree, NODE node,int current_height,int* height,int* max_nchildren, int parent, int posChild);
void btocl_computePartialLhKernel(RATES *Rates, TREES *Trees, OPTIONS *Opt, cl_ushort kernel_type);
void btocl_loadPartialLhBuffers(cl_command_queue queue, TREES* Trees, TREE* Tree, int NOS, int NoNodes);
void btocl_setPlhKernelInfo(TREES* Trees,TREE* Tree,int group,KernelInfo_lh* k);
int btocl_setPlhArgs(cl_kernel kernel, TREES* Trees, TREE* Tree, KernelInfo_lh k,  int group, int* err);
void printLhResult(double* pmatrix,double* plh,int nos);
void printTreeTraversal(TREE* Tree);
void printGroupPlh(double* plh, TREE* Tree,int group, int nos,int nsites);
void comparePlh(TREES* Trees, TREE* Tree);
void copyPlh(TREES* Trees, TREE* Tree,double* plh1, double* plh2);
void checkGroups(TREES* Trees, TREE* Tree);
void comparePlhDown(TREES* Trees, TREE* Tree);


// 	TREES: cl_mem      buffer_partialLh;
/*
Allocates OpenCl buffers and arrays with traversal information
*/
void btocl_AllocLhInfo(TREES* Trees) {
	int NOS,i;
	int Index, NIndex;
	cl_context context;
	int err;
	TREE* Tree;
	int group, node_idx, idx, children_idx, childrenIdx_idx;
	NODE N, C;
	int ntips, node_counter, groups_idx, groupsIdx_idx;
	int max_groupsize, max_ngroups, max_nnodes, max_nchildren;
	int *p_isTip;

	context = btocl_getContext();
	NOS = Trees->NoStates;
	max_nnodes = Trees->MaxNodes;


	// Allocate OpenCL buffer PartialLh and tree traversal information

	max_ngroups = 0;
	max_groupsize=0;
	max_nchildren = 0;
	for(Index=0;Index<Trees->NoTrees;Index++)
	{
		Tree = Trees->Tree[Index];
		// Some max info used to allocate OpenCl buffers
		if (max_ngroups < Tree->NoFGroups)
			max_ngroups = Tree->NoFGroups;
		if (Tree->NoFNodes[0] > max_groupsize)
			max_groupsize = Tree->NoFNodes[0];
		// Traversal group information
		Tree->groups = (int*)malloc(sizeof(int)*(Tree->NoNodes));
		Tree->groupsIdx = (int*)malloc(sizeof(int)*2*(Tree->NoFGroups + 1)); // one more group
		Tree->children = (int*)malloc(sizeof(int)*Tree->NoNodes);
		Tree->childrenIdx = (int*)malloc(sizeof(int)*2*Tree->NoNodes);
		Tree->parentInfo = (int*)malloc(sizeof(int)*2*max_nnodes);
		Tree->isTip = (int*)malloc(sizeof(int)*max_nnodes);

		// Initialize isTip
		p_isTip = Tree->isTip;
		for(i=0; i < max_nnodes; i++)
			p_isTip = 0;

		// Populate groupsIdx
		groups_idx = groupsIdx_idx=0;
		children_idx = childrenIdx_idx = 0;
		node_counter=0;

		// Calculate Height
		Tree->height = 0;
		Tree->max_nchildren = 0;
		// set tip info as well
		calculate_height(Tree, Tree->NodeList[0],0,&Tree->height,&Tree->max_nchildren,0,0);
		//printf("Height=%d\n",Tree->height);
		//printf("Num Groups=%d\n",Tree->NoFGroups);
		//printf("Max NChildren = %d\n",Tree->max_nchildren);

		if (Tree->max_nchildren > max_nchildren)
			max_nchildren = Tree->max_nchildren;

		// Find tips and include them in group [0]
		ntips = 0;
		Tree->groupsIdx[groupsIdx_idx++] = groups_idx;  // group start
		for(NIndex=0;NIndex < Tree->NoNodes;NIndex++)
		{
			if (Tree->NodeList[NIndex]->Tip) {
				ntips++;
				node_counter++; // for stats - add tip
				N = Tree->NodeList[NIndex];
				// note: N->NoNodes number of children
				Tree->groups[groups_idx++] = N->ID;
				childrenIdx_idx = 2*(N->ID);
				Tree->childrenIdx[childrenIdx_idx++] = -1;  // no children
				Tree->childrenIdx[childrenIdx_idx] = 0;    // no children, size 0
			}
		}
		Tree->groupsIdx[groupsIdx_idx++] = ntips; // group end

		// flat groups not including tips
		for(group=0;group < Tree->NoFGroups; group++)
		{
			Tree->groupsIdx[groupsIdx_idx++] = groups_idx;  // start
			Tree->groupsIdx[groupsIdx_idx++] = Tree->NoFNodes[group];  // size
			node_counter += Tree->NoFNodes[group];  // for stats
			//printf("[ Group %d NoFNodes %d : ",group,Tree->NoFNodes[group]);
			for(node_idx=0;node_idx < Tree->NoFNodes[group];node_idx++)
			{
				N = Tree->FNodes[group][node_idx]; // node belonging to group
				Tree->groups[groups_idx++] = N->ID;
				//printf("%d ",N->ID);
				// Traverse children
				//printf("<");
				// we will use the node id to access children data
				childrenIdx_idx = 2*(N->ID);
				Tree->childrenIdx[childrenIdx_idx++] = children_idx;  // first node
				Tree->childrenIdx[childrenIdx_idx] = N->NoNodes;  // number of nodes
				for(idx=0;idx < N->NoNodes; idx++) {
					C = N->NodeList[idx];
					Tree->children[children_idx++] = C->ID;
					//printf("%d ",C->ID);
				}

				//printf("> ");
			}

			//printf("] ");
		}
		//printf("\n");

		//printf("NoOfTips : %d\n",ntips);
		//printf("Nodes visited %d Tree NoNodes %d\n",node_counter,Tree->NoNodes);

		// Now lets test our traversal arrays
		//printTreeTraversal(Tree);


		//exit(0);

	}
	Trees->max_nchildren = max_nchildren;
	//printf("Max group size: %d\n",max_groupsize);
	//printf("No of Trees %d\n", Trees->NoTrees);

	//printf("Final MAx nchildren = %d\n",Trees->max_nchildren);


	// Allocate Open CL buffers
	// cl_mem      buffer_partialLh;
	Trees->buffer_partialLh  = clCreateBuffer(context, CL_MEM_READ_WRITE,
		sizeof(double)*(Trees->NoSites)*max_nnodes*NOS, NULL, &err);
	if (err != 0) {
		printf("Error allocating OpenCL buffer for partial likelihood\n");
		exit(0);
	}
	Trees->buffer_groups  = clCreateBuffer(context, CL_MEM_READ_ONLY,
		sizeof(int)*max_nnodes, NULL, &err);
	if (err != 0) {
			printf("Error allocating OpenCL buffer for group information\n");
			exit(0);
	}
	Trees->buffer_groupsIdx  = clCreateBuffer(context, CL_MEM_READ_ONLY,
		sizeof(int)*2*(max_ngroups+1), NULL, &err);
	if (err != 0) {
			printf("Error allocating OpenCL buffer for group information\n");
			exit(0);
	}
	// 		Tree->children = (int*)malloc(sizeof(int)*Tree->NoNodes);
	//	Tree->childrenIdx = (int*)malloc(sizeof(int)*2*Tree->NoNodes);
	Trees->buffer_children  = clCreateBuffer(context, CL_MEM_READ_ONLY,
		sizeof(int)*max_nnodes, NULL, &err);
	if (err != 0) {
			printf("Error allocating OpenCL buffer for children information\n");
			exit(0);
	}
	Trees->buffer_childrenIdx  = clCreateBuffer(context, CL_MEM_READ_ONLY,
		sizeof(int)*2*max_nnodes, NULL, &err);
	if (err != 0) {
			printf("Error allocating OpenCL buffer for children indexing\n");
			exit(0);
	}
	//	cl_mem		buffer_plhFactor; // MaxNodes*NoSites*NOS
	Trees->buffer_plhFactor  = clCreateBuffer(context, CL_MEM_READ_WRITE,
		sizeof(double)*(Trees->NoSites)*max_nnodes*NOS, NULL, &err);
	//Trees->buffer_plhFactor  = clCreateBuffer(context, CL_MEM_READ_WRITE,
	//	sizeof(double)*(Trees->NoSites)*max_nnodes*NOS*Trees->max_nchildren, NULL, &err);
	if (err != 0) {
			printf("Error allocating OpenCL buffer for children indexing\n");
			exit(0);
	}
	// debug
	Trees->buffer_debug_plhFactor  = clCreateBuffer(context, CL_MEM_READ_WRITE,
		sizeof(double)*(Trees->NoSites)*max_nnodes*NOS, NULL, &err);
	if (err != 0) {
			printf("Error allocating OpenCL buffer for children indexing\n");
			exit(0);
	}


	// New method
	Trees->buffer_parentInfo  = clCreateBuffer(context, CL_MEM_READ_ONLY,
		sizeof(int)*2*max_nnodes, NULL, &err);
	if (err != 0) {
			printf("Error allocating OpenCL buffer for Tree parent Info\n");
			exit(0);
	}
	Trees->buffer_isTip  = clCreateBuffer(context, CL_MEM_READ_ONLY,
		sizeof(int)*max_nnodes, NULL, &err);
	if (err != 0) {
			printf("Error allocating OpenCL buffer for Tree Tip info\n");
			exit(0);
	}

	// CPU buffer

	Trees->previous_plh = (double*)malloc(sizeof(double)*(Trees->NoSites)*max_nnodes*NOS);
	Trees->temp_plh = (double*)malloc(sizeof(double)*(Trees->NoSites)*max_nnodes*NOS);
	Trees->plhFactor = (double*)malloc(sizeof(double)*(Trees->NoSites)*max_nnodes*NOS);
	//Trees->plhFactor =
	//	(double*)malloc(sizeof(double)*(Trees->NoSites)*max_nnodes*NOS*Trees->max_nchildren);



	Trees->debug_plhFactor = (double*)malloc(sizeof(double)*(Trees->NoSites)*max_nnodes*NOS);
	//Trees->debug_plhFactor =
	//	(double*)malloc(sizeof(double)*(Trees->NoSites)*max_nnodes*NOS*Trees->max_nchildren);


	//exit(0);
}


void btocl_FreeLhInfo(TREES* Trees) {
	int i;
	TREE* Tree;


	// plh related
	free(Trees->plhFactor);


	// Must delete arrays for each tree
	// Make sure the tree has not been deleted before!!

	for(i=0; i < Trees->NoTrees; i++) {
		Tree = Trees->Tree[i];
		free(Tree->groups);
		free(Tree->groupsIdx);
		free(Tree->children);
		free(Tree->childrenIdx);
	}

	// Partial LH


	clReleaseMemObject(Trees->buffer_partialLh);
	clReleaseMemObject(Trees->buffer_groups);
	clReleaseMemObject(Trees->buffer_groupsIdx);
	clReleaseMemObject(Trees->buffer_children);
	clReleaseMemObject(Trees->buffer_childrenIdx);

}

// Recursive tree traversal function that sets and computes information for
// Tree: height, max_nchildren, isTip, parentInfo
void calculate_height(TREE* Tree, NODE node,int current_height,int* height,int* max_nchildren, int parent, int posChild) {
	int i, node_id;
	current_height++;
	node_id = node->ID;
	Tree->parentInfo[2*node_id] = parent;
	Tree->parentInfo[2*node_id+1] = posChild;
	if (node->Tip) {  // Tip
		if (node->NoNodes > 0) {
			printf("Error!! - tips with children\n");
			exit(0);
		}
		Tree->isTip[node_id] = 1;
		if (current_height > *height)
			*height = current_height;
		return;
	} else {
		if (node->NoNodes > *max_nchildren)
			*max_nchildren = node->NoNodes;
		for(i=0;i< node->NoNodes;i++)
			calculate_height(Tree,node->NodeList[i],current_height,height,max_nchildren,node_id,i);
	}
	return;
}



/* **************************************************************
Computation of Partial Likelihood vectors for all nodes in Tree
****************************************************************** */

/* Reference CPU version
*  it takes the group number as extra input
*  does NOT use new flat tree data structures
*/
void btocl_computePartialLhCPU(TREES *Trees, TREE* Tree, int group) {
	int g, gsize,node_idx, child_idx, nchildren, site;
	int row, col, nos;
	double *plh_node, *plh_node_start, *plh_child, *plh_child_start;
	double lh;
	NODE node, child;
	double *pmatrix, *pmatrix_child;

	nos = Trees->NoStates;
	pmatrix = Trees->PList[0]->me[0];

	// Original groups, tips not included
	for(g=group; g < Tree->NoFGroups; g++)
	{
		gsize = Tree->NoFNodes[g];
		//printf("gs ");
		// RunNodeGroup(GIndex, FALSE, Rates, Tree, Trees, Opt, SiteNo);
		//printf(" g %d",g);
		for(node_idx=0; node_idx < gsize; node_idx++)
		{
			node = Tree->FNodes[g][node_idx];
			nchildren = node->NoNodes;

			//printf("node id: \n",node->ID);
			for(site=0; site < Trees->NoSites; site++) {
				plh_node_start = node->Partial[site];
				plh_node = plh_node_start;
				for(col=0; col < nos; col++,plh_node++)
					*plh_node = 1.0;
				for(child_idx=0; child_idx < nchildren; child_idx++) {
					child = node->NodeList[child_idx];
					pmatrix_child = &pmatrix[nos*nos*(child->ID)];
					plh_child_start = child->Partial[site];
					plh_node = plh_node_start;
					for(row=0; row < nos; row++,plh_node++) {
						lh = 0.0;
						plh_child = plh_child_start;
						for(col=0; col < nos; col++) {
							lh += (*pmatrix_child)*(*plh_child);
							pmatrix_child++; plh_child++;
						}
						(*plh_node) *= lh;
					}
				}
			}
			//return; // after first node
		}
	}
}

/* Reference CPU version
*  it takes the group number as extra input
*  uses new flat tree data structures
*/
void btocl_computePartialLhCPUFlat(TREES *Trees, TREE* Tree, int group) {
	int g, gsize,node_idx, node_idx_start, child_idx, child_idx_start, nchildren, site, nodeid, childid;
	int row, col, nos;
	double *plh_node, *plh_node_start, *plh_child, *plh_child_start;
	double lh;
//	NODE node, child;
	double *pmatrix, *pmatrix_child, *plh;
	int nsites;
	int *children, *childrenIdx;

	if (group <= 0) {
		printf("ERROR while computing partial likelihood: group nust be > 0\n");
		exit(0);
	}
	// In order to compute whole tree, group = 1

	nsites = Trees->NoSites;
	nos = Trees->NoStates;
	pmatrix = Trees->PList[0]->me[0];
	plh = Tree->NodeList[0]->Partial[0];
	children = Tree->children;
	childrenIdx = Tree->childrenIdx;

	// Last group number is Tree->NoFGroups (for flat data structures)
	for(g=group; g <= Tree->NoFGroups; g++)
	{
		//printf("[Group %d: ",g);
		node_idx_start = Tree->groupsIdx[2*g];  // start
		gsize = Tree->groupsIdx[2*g+1];    // size
		//printf("start %d size %d ", node_idx,gsize);
		for(node_idx=node_idx_start; node_idx < node_idx_start+gsize; node_idx++)
		{
			// node = Tree->FNodes[g][node_idx];
			nodeid = Tree->groups[node_idx];  // current node in group
			// nchildren = node->NoNodes;
			child_idx_start = childrenIdx[2*nodeid];
			nchildren = childrenIdx[2*nodeid+1];

			//printf(" %d < ",nodeid);

			//printf("node id: \n",node->ID);
			for(site=0; site < nsites; site++) {

				// actually, the plh should be contiguos
				plh_node_start =  &plh[nsites*nos*nodeid+nos*site];
				plh_node = plh_node_start;
				for(col=0; col < nos; col++,plh_node++)
					*plh_node = 1.0;

				for(child_idx=child_idx_start; child_idx < child_idx_start+nchildren; child_idx++) {

					childid =  children[child_idx];
					//printf(" %d ",childid);
					pmatrix_child = &pmatrix[nos*nos*childid];

					plh_child_start = &plh[nsites*nos*childid+nos*site];
					plh_node = plh_node_start;
					for(row=0; row < nos; row++,plh_node++) {
						lh = 0.0;
						plh_child = plh_child_start;
						for(col=0; col < nos; col++) {
							lh += (*pmatrix_child)*(*plh_child);
							pmatrix_child++; plh_child++;
						}
						(*plh_node) *= lh;
					}
				}
			}
			//printf("> ",nodeid);
			//return; // after first node

		}
		//printf("]\n");
	}
}

/*  Correct - original version ***** */
//void btocl_computePartialLhOriginal(RATES* Rates, TREES *Trees, OPTIONS *Opt) {
//	int SiteNo;
//
//	for(SiteNo=0;SiteNo<Trees->NoSites;SiteNo++)
//	{
//			SumLhLiner(Rates, Trees, Opt, SiteNo);
//			//Err = SumLh(Rates, Trees, Opt, SiteNo);
//	}
//}



int btocl_computePartialLh(RATES* Rates, TREES *Trees, OPTIONS *Opt) {
	int NOS, NoSites;
	cl_ushort kernel_type;
	int err;


	double *plh, *pmatrix;
	TREE* Tree;


	NOS = Trees->NoStates;
	NoSites = Trees->NoSites;
	Tree = Trees->Tree[Rates->TreeNo];

	pmatrix = Trees->PList[0]->me[0];
	plh = Tree->NodeList[0]->Partial[0];
	err = 0;

	// original version
	//btocl_computePartialLhOriginal(Rates, Trees, Opt);
	//printLhResult(pmatrix,plh,NOS);
	//return;

	// new sequential version

	//btocl_computePartialLhCPU(Trees, Trees->Tree[Rates->TreeNo], 0);  // new cpu version
	//printLhResult(pmatrix,plh,NOS);
	//return;

	//printTreeTraversal(Tree);

	//btocl_computePartialLhCPUFlat(Trees, Trees->Tree[Rates->TreeNo], 1);  // new cpu version flat
	//return;


	// **** Testing!!!!
	//printf("***********   Correct results ***************\n");
	//plh = Tree->NodeList[0]->Partial[0];
	//printLhResult(pmatrix,plh,NOS);

	//exit(0);


	// Decides kernel type to use based on e.g. NOS
	switch(NOS) {
		case 2:
			//kernel_type = BTOCL_PLHNODE_NOS2; break;
			//kernel_type = BTOCL_PLHROWFULLG; break;
			kernel_type = BTOCL_PLHROWG; break;
			//kernel_type = BTOCL_PLHREDUCEGL; break;
		case 4:
			kernel_type = BTOCL_PLHROWG; break;
			//kernel_type = BTOCL_PLHREDUCEGL; break;

		default:
			kernel_type = BTOCL_PLHROWGL; break;
			//kernel_type = BTOCL_PLHREDUCEGL;
	}

	//printf("NOS %d Sites %d\n",NOS,NoSites);
	//printTreeTraversal(Tree);

	//printGroupPlh(plh, Tree,0,NOS,NoSites);

	//memcpy(Trees->previous_plh,plh,sizeof(plh));

	//printGroupPlh(Trees->previous_plh, Tree,0,NOS,NoSites);
	//exit(0);

	//checkGroups(Trees, Tree);

	//copyPlh(Trees,Tree,Trees->temp_plh,plh);
	//btocl_computePartialLhOriginal(Rates, Trees, Opt);
	//copyPlh(Trees,Tree,Trees->previous_plh,plh);
	//copyPlh(Trees,Tree,plh,Trees->temp_plh);

	//printf("Entering computePartialPlhKernel\n");
	btocl_computePartialLhKernel(Rates,Trees,Opt,kernel_type);

	//btocl_computePartialLhCPUFlat(Trees, Tree, 1);

	// GPU version has printLhResult
	//printLhResult(pmatrix,plh,NOS);


	//printGroupPlh(plh, Tree,0,NOS,NoSites);
	//comparePlhDown(Trees, Tree);
	//comparePlh(Trees, Tree);

	//exit(0);



	 //comparePlh(Trees, Tree);



	//exit(0);

	return err; // set correct value
}

// Prepares and executes call to GPU
void btocl_computePartialLhKernel(RATES *Rates, TREES *Trees, OPTIONS *Opt, cl_ushort kernel_type) {
	TREE *Tree;
	int NOS, MaxNodes, NoNodes, group;
	INVINFO* I;
	double *pmatrix, *plh;
//	int lastGroup;

	cl_command_queue queue;
	cl_kernel kernel;
	cl_context context;

	KernelInfo_lh k;

	int err;
	int argnum;

	//btdebug_enter("pmatrix");
	I = Trees->InvInfo[0];
	NOS = Trees->NoStates;
	MaxNodes = Trees->MaxNodes;

	Tree = Trees->Tree[Rates->TreeNo];
	NoNodes =  Tree->NoNodes;

	queue =  btocl_getCommandQueue();
	context = btocl_getContext();

	k.kernel_type = kernel_type;
	kernel = btocl_getKernel(kernel_type);
	if (kernel == NULL) {
		printf("Error: Couldn't load kernel %s\n",btocl_getKernelName(kernel_type));
		exit(1);
	}


	// load plh and group/tree info buffers - Done ONCE
	//printf("loading buffers\n");
	btocl_loadPartialLhBuffers(queue,Trees,Tree,NOS,Tree->NoNodes);
	//printf("finished loading buffers\n");

	//printf("Buffers loaded\n");

	// group 0 now contains the tips
	// groups:[0,1,2,....,Tree->NoFGroups]
	group=0;
	if (kernel_type == BTOCL_PLHROWFULLG) {
		btocl_setPlhKernelInfo(Trees,Tree,0,&k);
		argnum = btocl_setPlhArgs(kernel,Trees,Tree,k,0,&err);
		if (err < 0) {
			printf("Error. Couldn't set argument: %d\n",argnum-1);
		}

		//printf("nnodes %d nos %d nsites %d\n", NoNodes, NOS, Trees->NoSites);
		//printf(" no workitems %d\n", k.globalws[0]);
		clEnqueueNDRangeKernel(queue,kernel,k.dim,NULL,k.globalws,k.localws,0,NULL,NULL);
		// if (group <= Tree->NoFGroups)
		group = Tree->NoFGroups+1; // to skip CPU version
	} else if (kernel_type == BTOCL_PLHROWALLG) {  // Loop inside kernel
		btocl_setPlhKernelInfo(Trees,Tree,0,&k); // start from group 0
		argnum = btocl_setPlhArgs(kernel,Trees,Tree,k,0,&err);
		if (err < 0) {
			printf("Error. Couldn't set argument: %d\n",argnum-1);
		}
		if (k.usewg == FALSE) {
				//printf("group %d no workitems %d\n", group, k.globalws[0]);
				clEnqueueNDRangeKernel(queue,kernel,k.dim,NULL,k.globalws,NULL,0,NULL,NULL);
			} else {
				clEnqueueNDRangeKernel(queue,kernel,k.dim,NULL,k.globalws,k.localws,0,NULL,NULL);
		}
		group = k.lastGroup+1;  // next group to process
	} else
	{
		// print workgroup information:
		//for(group=0; group <= Tree->NoFGroups; group++)
		//{
		//	printf("Group: %d size %d  \n", group, Tree->groupsIdx[2*group+1]);
		//}

		for(group=0; group <= Tree->NoFGroups; group++)
		{

			//printf("Group: %d\n",group);

			//if ((group==0) &&(kernel_type==BTOCL_PLH)) {
			//	continue; // skip group zero
			//}

			// decide on workgroups and workitems
			//printf("setting kernel info\n");
			btocl_setPlhKernelInfo(Trees,Tree,group,&k);
			//printf("Finished setting kernel info\n");

			if (k.passedCut==TRUE) {
				//printf("TRUE!!");
				if (kernel_type == BTOCL_PLH)
					break;
				if (group==0)  // there's nothing to collect
					break;
			}

			//printf("before setPlhArgs w %d a %d", k.globalws[0],k.action);

			argnum = btocl_setPlhArgs(kernel,Trees,Tree,k,group,&err);
			if (err < 0) {
				printf("Error. Couldn't set argument: %d\n",argnum-1);
			}

			//printf("after setPlhARgs\n");

			//printf("gpuin--");
			if (k.usewg == FALSE) {
				//printf("group %d no workitems %d sites %d\n", group, k.globalws[0],Trees->NoSites);
				clEnqueueNDRangeKernel(queue,kernel,k.dim,NULL,k.globalws,NULL,0,NULL,NULL);
			} else {
				//printf("wgroup %d no workitems %d sites %d\n", group, k.globalws[0],Trees->NoSites);
				clEnqueueNDRangeKernel(queue,kernel,k.dim,NULL,k.globalws,k.localws,0,NULL,NULL);
			}

			if (k.passedCut) {
				group++; // this group is done try next one
				break; // finish loop for sure
			}
		}
	}
	// group: next group to process

	pmatrix = Trees->PList[0]->me[0];
	plh = Tree->NodeList[0]->Partial[0];

	clEnqueueReadBuffer(queue,Trees->buffer_partialLh,CL_TRUE,0,
		sizeof(double)*(Trees->NoSites)*NoNodes*NOS, plh,0,0,NULL);



	// debug
	//clEnqueueReadBuffer(queue,Trees->buffer_plhFactor,CL_TRUE,0,
	//	sizeof(double)*(Trees->NoSites)*NoNodes*NOS, Trees->plhFactor,0,0,NULL);
	// new clEnqueueReadBuffer(queue,Trees->buffer_plhFactor,CL_TRUE,0,
	//	sizeof(double)*(Trees->NoSites)*NoNodes*NOS*Trees->max_nchildren, Trees->plhFactor,0,0,NULL);

	//clEnqueueReadBuffer(queue,Trees->buffer_pmatrix,CL_TRUE,0,
	//			sizeof(double)*MaxNodes*NOS*NOS, Trees->check_pmatrix,0,0,NULL);
	// debug
	//clEnqueueReadBuffer(queue,Trees->buffer_debug_plhFactor,CL_TRUE,0,
	//	sizeof(double)*(Trees->NoSites)*NoNodes*NOS, Trees->debug_plhFactor,0,0,NULL);
	// new clEnqueueReadBuffer(queue,Trees->buffer_debug_plhFactor,CL_TRUE,0,
	//	sizeof(double)*(Trees->NoSites)*NoNodes*NOS*Trees->max_nchildren, Trees->debug_plhFactor,0,0,NULL);

	if (group <= Tree->NoFGroups) {
		//printf("CPU version from group %d \n",group);
		btocl_computePartialLhCPUFlat(Trees, Trees->Tree[Rates->TreeNo], group);
		//printf("out ");
	}

	//printSetPMatrixResults(Rates, Trees, Opt, RateMult, Kappa);

	//printLhResult(pmatrix,plh,NOS);



	//exit(0);

}


void btocl_loadPartialLhBuffers(cl_command_queue queue, TREES* Trees, TREE *Tree, int NOS, int NoNodes) {
	// PMatrix is already loaded in the GPU

	// PartialLh
	clEnqueueWriteBuffer(queue,Trees->buffer_partialLh,CL_TRUE,0,
		sizeof(double)*(Trees->NoSites)*NoNodes*NOS,
		Tree->NodeList[0]->Partial[0],0,0,NULL);
	// load groups, children, childrenIdx
	clEnqueueWriteBuffer(queue,Trees->buffer_groups,CL_TRUE,0,
		sizeof(int)*NoNodes, Tree->groups,0,0,NULL);
	//Tree->groupsIdx = (int*)malloc(sizeof(int)*2*(Tree->NoFGroups + 1)); // one more group
	clEnqueueWriteBuffer(queue,Trees->buffer_groupsIdx,CL_TRUE,0,
		sizeof(int)*2*(Tree->NoFGroups + 1), Tree->groupsIdx,0,0,NULL);
	clEnqueueWriteBuffer(queue,Trees->buffer_children,CL_TRUE,0,
		sizeof(int)*NoNodes, Tree->children,0,0,NULL);
	clEnqueueWriteBuffer(queue,Trees->buffer_childrenIdx,CL_TRUE,0,
		sizeof(int)*2*NoNodes, Tree->childrenIdx,0,0,NULL);
	// New version - read only buffers
	clEnqueueWriteBuffer(queue,Trees->buffer_parentInfo,CL_TRUE,0,
		sizeof(int)*2*NoNodes, Tree->parentInfo,0,0,NULL);
	clEnqueueWriteBuffer(queue,Trees->buffer_isTip,CL_TRUE,0,
		sizeof(int)*NoNodes, Tree->isTip,0,0,NULL);

}

// Could be merged with setExpqtKernelInfo
void btocl_setPlhKernelInfo(TREES* Trees,TREE* Tree,int group, KernelInfo_lh* k) {
	int NOS,NoNodes, nsites;
	int group_size;
	int wg_size;
	cl_ushort kernel_type;

	NOS = Trees->NoStates;
	NoNodes =  Tree->NoNodes;
	nsites = Trees->NoSites;
	// IMPORTANT - check before compiling

	k->dim = 1;
	k->passedCut = FALSE;

	kernel_type = k->kernel_type;
	if (kernel_type ==  BTOCL_PLH) {
		k->usewg = FALSE;  // No workgroups
		// use flat tree information
		//group_size = Tree->NoFNodes[group];
		group_size = Tree->groupsIdx[2*group+1];
		k->globalws[0] = (Trees->NoSites)*(group_size); // nodes per group
		if (k->globalws[0] < 20)
			k->passedCut = TRUE;

	} else if ((kernel_type == BTOCL_PLHNODE) || (kernel_type == BTOCL_PLHNODE_NOS2)) {
		k->usewg = FALSE;
		if (group==0) // tips, compute only
			k->action = 2;
		else
			k->action = 1;
		group_size = Tree->groupsIdx[2*group+1];
		k->globalws[0] = (Trees->NoSites)*(group_size); // nodes per group
		if (k->globalws[0] < 40) {  //
			k->passedCut = TRUE;
			k->action = 0; // collect
		}
	} else if (kernel_type == BTOCL_PLHROW) {
		k->usewg = FALSE;
		if (group==0) // tips, compute only
			k->action = 2;
		else
			k->action = 1;
		group_size = Tree->groupsIdx[2*group+1];
		k->globalws[0] = (Trees->NoSites)*(group_size)*NOS; // nodes per group
		if (k->globalws[0] < 800) {
			k->passedCut = TRUE;
			k->action = 0; // collect
		}

	} else if ((kernel_type ==  BTOCL_PLHROWG) // one call per group, using workgroups
			|| (kernel_type == BTOCL_PLHROWGL)) {
		//printf("Here!!\n");
		if (group==0) // tips, compute only
			k->action = 2;
		else
			k->action = 1;
		k->usewg = TRUE;
		// start from desired wg_size
		wg_size = 512;  // target. Final size will depend on nos*nsites
		k->num_nodes_wg = wg_size / (NOS*nsites);
		k->num_rows_wg = wg_size / NOS;
		//k->num_nodes_wg = 8;  // number of nodes per workgroup
		k->num_nodes = Tree->groupsIdx[2*group+1];  // number of nodes to compute

		//printf("Num nodes wg = %d nodes in group %d",k->num_nodes_wg, k->num_nodes);
		//printf("Num rows wg = %d nodes in group %d",k->num_rows_wg);

		k->num_wg = (k->num_nodes*nsites) / k->num_rows_wg;

		if ( (k->num_nodes*nsites) % k->num_rows_wg != 0)
			k->num_wg++;
		k->localws[0] = NOS*k->num_rows_wg; // workitems per workgroup
		k->globalws[0] = k->localws[0]*k->num_wg;
		// Figure out when to stop
		if (k->globalws[0] <1000) {
			k->passedCut = TRUE;
			k->action = 0; // collect
		}
	} else if (kernel_type == BTOCL_PLHREDUCEGL) {
		if (group==0) // tips, compute only
			k->action = 2;
		else
			k->action = 1;
		k->usewg = TRUE;
		// start from desired wg_size
		wg_size = 512;  // target. Final size will depend on nos*nsites
		k->num_nodes_wg = wg_size / (NOS*NOS*nsites); // pmatrix size mult by nsites
		//k->num_nodes_wg = 8;  // number of nodes per workgroup
		if (k->num_nodes_wg ==0) // This may blowup the workitem allocation per wg
			k->num_nodes_wg = 1;
		k->num_nodes = Tree->groupsIdx[2*group+1];  // number of nodes to compute
		k->num_wg = k->num_nodes / k->num_nodes_wg;
		if (k->num_nodes % k->num_nodes_wg != 0)  // not exact, then go one up
			k->num_wg++;
		k->localws[0] = NOS*NOS*nsites*k->num_nodes_wg; // workitems per workgroup
		k->globalws[0] = k->localws[0]*k->num_wg;
		// Figure out when to stop
		if (k->globalws[0] <1000) {
			k->passedCut = TRUE;
			k->action = 0; // collect
		}

	} else if (kernel_type == BTOCL_PLHROWALLG) {
		k->usewg = TRUE;
		// start from desired wg_size
		wg_size = 64;  // target. Final size will depend on nos*nsites
		k->num_nodes_wg = wg_size / (NOS*nsites);
		//k->num_nodes_wg = 8;  // number of nodes per workgroup
		k->num_nodes = Tree->groupsIdx[2*group+1];  // number of nodes to compute
		k->num_wg = k->num_nodes / k->num_nodes_wg;
		if (k->num_nodes % k->num_nodes_wg != 0)
			k->num_wg++;
		k->localws[0] = NOS*nsites*k->num_nodes_wg; // workitems per workgroup
		k->globalws[0] = k->localws[0]*k->num_wg;
		// Figure out when to stop
		k->lastGroup = (Tree->NoFGroups)/2+1;  // stop half-way and resume with CPU
		// k->lastGroup = (Tree->NoFGroups);      // do all groups

	} else if (kernel_type == BTOCL_PLHROWFULLG) {
		k->usewg = TRUE;
		// start from desired wg_size
		wg_size = 128;  // target. Final size will depend on nos*nsites
		k->num_nodes_wg = wg_size / (NOS*nsites);
		//k->num_nodes_wg = 8;  // number of nodes per workgroup
		if (k->num_nodes_wg==0)
			k->num_nodes_wg = 1;
		k->num_nodes = NoNodes;  // number of nodes to compute
		k->num_wg = k->num_nodes / k->num_nodes_wg;
		if (k->num_nodes % k->num_nodes_wg != 0)
			k->num_wg++;
		//printf("-> num_nodes_wg = %d\n",k->num_nodes_wg);
		k->localws[0] = NOS*nsites*k->num_nodes_wg; // workitems per workgroup
		//printf("Number workitems per workgroup= %d\n",k->localws[0]);
		k->globalws[0] = k->localws[0]*k->num_wg;  // total number of workitems

	} else {
		printf("PLH: Kernel type not supported\n");
		exit(0);
	}

}

/*
   Sets the kernel's actual parameters, depending on kernel type and group information.
*/
int btocl_setPlhArgs(cl_kernel kernel, TREES* Trees, TREE* Tree, KernelInfo_lh k, int group, int* err) {
	INVINFO* I;
	int NOS, NoNodes, nsites;
	int argnum, local_memSize;


	I = Trees->InvInfo[0];
	nsites = Trees->NoSites;
	NOS = Trees->NoStates;
	NoNodes =  Tree->NoNodes;

	*err = 0;
	argnum = 0;
	if ((*err = clSetKernelArg(kernel,argnum++,sizeof(cl_mem), &(Trees->buffer_pmatrix))) < 0)
		return argnum;
	if ((*err = clSetKernelArg(kernel,argnum++,sizeof(cl_mem), &(Trees->buffer_partialLh))) < 0)
		return argnum;
	if ((*err = clSetKernelArg(kernel,argnum++,sizeof(NOS), &NOS)) < 0)
			return argnum;
	if ((*err = clSetKernelArg(kernel,argnum++,sizeof(nsites), &nsites)) < 0)
			return argnum;

	if (k.kernel_type == BTOCL_PLHROWFULLG) {
	// parentInfo
	if ((*err = clSetKernelArg(kernel,argnum++,sizeof(cl_mem), &(Trees->buffer_parentInfo))) < 0)
		return argnum;
	// isTip
	if ((*err = clSetKernelArg(kernel,argnum++,sizeof(cl_mem), &(Trees->buffer_isTip))) < 0)
		return argnum;
	// plhFactor
	if ((*err = clSetKernelArg(kernel,argnum++,sizeof(cl_mem), &(Trees->buffer_plhFactor))) < 0)
		return argnum;
	// num_nodes
	if ((*err = clSetKernelArg(kernel,argnum++,sizeof(k.num_nodes), &k.num_nodes)) < 0)
		return argnum;
	// height
	if ((*err = clSetKernelArg(kernel,argnum++,sizeof(Tree->height), &Tree->height)) < 0)
		return argnum;
	// max_nchildren
	if ((*err = clSetKernelArg(kernel,argnum++,sizeof(Tree->max_nchildren), &Tree->max_nchildren)) < 0)
		return argnum;
	return argnum;  // that is, exit function
	}


	//printf("start %d size %d ",Tree->groupsIdx[2*group],Tree->groupsIdx[2*group+1]);
	if (k.kernel_type != BTOCL_PLHROWALLG) {  // gstart
		if ((*err = clSetKernelArg(kernel,argnum++,sizeof(int),
							&(Tree->groupsIdx[2*group]))) < 0)  // start index
			return argnum;
	}
	if ((*err = clSetKernelArg(kernel,argnum++,sizeof(cl_mem), &(Trees->buffer_groups))) < 0)
		return argnum;
	if (k.kernel_type == BTOCL_PLHROWALLG) {
		if ((*err = clSetKernelArg(kernel,argnum++,sizeof(cl_mem), &(Trees->buffer_groupsIdx))) < 0)
			return argnum;
		if ((*err = clSetKernelArg(kernel,argnum++,sizeof(k.lastGroup), &k.lastGroup)) < 0)
			return argnum;
	}
	if ((*err = clSetKernelArg(kernel,argnum++,sizeof(cl_mem), &(Trees->buffer_children))) < 0)
		return argnum;
	if ((*err = clSetKernelArg(kernel,argnum++,sizeof(cl_mem), &(Trees->buffer_childrenIdx))) < 0)
		return argnum;
	if (k.kernel_type == BTOCL_PLHROWALLG) {
		if ((*err = clSetKernelArg(kernel,argnum++,sizeof(cl_mem), &(Trees->buffer_plhFactor))) < 0)
			return argnum;
	}
	if ((k.kernel_type == BTOCL_PLHNODE)||(k.kernel_type == BTOCL_PLHNODE_NOS2) ||
		(k.kernel_type == BTOCL_PLHROW) ||
		(k.kernel_type == BTOCL_PLHROWG)||(k.kernel_type == BTOCL_PLHROWGL) ||
		(k.kernel_type == BTOCL_PLHREDUCEGL)) {
		if ((*err = clSetKernelArg(kernel,argnum++,sizeof(cl_mem), &(Trees->buffer_plhFactor))) < 0)
			return argnum;
		if ((*err = clSetKernelArg(kernel,argnum++,sizeof(k.action), &k.action)) < 0)
			return argnum;
	}
	if ((k.kernel_type == BTOCL_PLHROWG)||(k.kernel_type == BTOCL_PLHROWGL) ||
		(k.kernel_type == BTOCL_PLHREDUCEGL)) {  // pass fgroup size
		if ((*err = clSetKernelArg(kernel,argnum++,sizeof(k.num_nodes), &k.num_nodes)) < 0)
			return argnum;
	}
	// debug
	//if ((k.kernel_type == BTOCL_PLHNODE_NOS2)) {
	//	if ((*err = clSetKernelArg(kernel,argnum++,sizeof(cl_mem), &(Trees->buffer_debug_plhFactor))) < 0)
	//		return argnum;
	//}
	if (k.kernel_type == BTOCL_PLHROWGL) {
		local_memSize = sizeof(double)*NOS*k.num_rows_wg;  // changed
		if ((*err = clSetKernelArg(kernel,argnum++,local_memSize, NULL)) < 0)  // local_plh
				return argnum;
	}
	if (k.kernel_type == BTOCL_PLHREDUCEGL) {
		//if ((*err = clSetKernelArg(kernel,argnum++,sizeof(k.pow2), &k.pow2)) < 0)
		//	return argnum;
		local_memSize = sizeof(double)*nsites*NOS*NOS*k.num_nodes_wg;
		if ((*err = clSetKernelArg(kernel,argnum++,local_memSize, NULL)) < 0)  // dotprod
				return argnum;
	}
	return argnum;
}

void printLhResult(double* pmatrix,double* plh,int NOS) {
	int a,b,c;
	a=132;
	b = 133;
	c = 136;
	printf("Node a. LH and pamatrix\n");
	btlin_printR(plh+NOS*b,1,NOS);
	//printf("PMatrix\n");
	//btlin_printR(pmatrix+NOS*NOS*b,1,NOS);
	// Node c
	printf("Node b. LH and pamatrix\n");
	btlin_printR(plh+NOS*c,1,NOS);
	//printf("PMatrix\n");
	//btlin_printR(pmatrix+NOS*NOS*c,1,NOS);
	// Node a
	printf("Node A. Result LH \n");
	btlin_printR(plh+NOS*a,1,NOS);

}


void printTreeTraversal(TREE* Tree) {
		// print group information
	int start,end, start2, end2, group, node_id;
	for(group=0;group <= Tree->NoFGroups; group++)
	{
			start = Tree->groupsIdx[2*group];
			end = start + Tree->groupsIdx[2*group+1] - 1;
			printf("[ Group %d Size %d: ",group,end-start+1);
			while(start <= end) {
				node_id  = Tree->groups[start];
				printf(" %d",node_id);  // node id
				start2 = Tree->childrenIdx[2*node_id];
				end2 =   start2 + Tree->childrenIdx[2*node_id+1] - 1;
				printf("<");
				if (start2==-1)
					printf("Tip");
				else {
					while(start2 <= end2) {
						printf("%d ",Tree->children[start2]);
						start2++;
					}
				}
				printf(">");
				start++;
			}
			printf("]\n ");

		}
	printf("\n");
	// end print group info

}

void printGroupPlh(double* plh, TREE* Tree,int group, int nos,int nsites) {
	int start, end, node_id;
	start = Tree->groupsIdx[2*group];
	end = start + Tree->groupsIdx[2*group+1] - 1;
	while(start <= end) {
		node_id  = Tree->groups[start];
		printf("Node %d\n",node_id);
		btlin_printR(plh+nos*nsites*node_id,1,nos);
		start++;
	}
}

// Check correct order of groups in tree
void checkGroups(TREES* Trees, TREE* Tree) {
	int max_nnodes,i;
	int *nodes;
	int start,end, start2,end2, group, node_id, nos, child_id;
	max_nnodes = Trees->MaxNodes;
	nodes = (int*)malloc(sizeof(int)*max_nnodes);
	for(i=0; i < max_nnodes;i++) nodes[i] = 0;
	nos = Trees->NoStates;
	for(group=0;group <= Tree->NoFGroups; group++)
	{
		start = Tree->groupsIdx[2*group];
		end = start + Tree->groupsIdx[2*group+1] - 1;
		//printf("\n[ Group %d Size %d:\n ",group,end-start+1);

		while(start <= end) {
			node_id  = Tree->groups[start];
			nodes[node_id] = group;
			start++;
		}
	}
	// Now check that for all children of group g: group(child) < g
	for(group=0;group <= Tree->NoFGroups; group++)
	{
		start = Tree->groupsIdx[2*group];
		end = start + Tree->groupsIdx[2*group+1] - 1;
		//printf("\n[ Group %d Size %d:\n ",group,end-start+1);

		while(start <= end) {
			node_id  = Tree->groups[start];
			// visit children
			start2 = Tree->childrenIdx[2*node_id];
			end2 =   start2 + Tree->childrenIdx[2*node_id+1] - 1;

			if (start2!=-1) { // not a tip
				while(start2 <= end2) {
					child_id = Tree->children[start2];

					if (nodes[child_id] >= group) {
						printf("Problem group %d child %d\n",group,child_id);
						free(nodes);
						exit(0);
					}
					start2++;
				}
			}


			start++;
		}
	}

	printf("All nodes good\n");

}

void comparePlhDown(TREES* Trees, TREE* Tree) {
		// print group information
	int start,end, start2, end2, group, node_id, child_id, i, j, nos,nsites;
	double* plh, *plh1,*plh2, *plhf, *plhfdebug;

	double diff;
	int gsize,pos;
	double *pmatrix = Trees->PList[0]->me[0];

	plh = Tree->NodeList[0]->Partial[0];
	nsites = Trees->NoSites;
	nos = Trees->NoStates;

	for(group=0;group <= Tree->NoFGroups; group++)
	{

			start = Tree->groupsIdx[2*group];
			end = start + Tree->groupsIdx[2*group+1] - 1;
			//printf("\n[ Group %d Size %d:\n ",group,end-start+1);
			gsize = end-start+1;
			pos=0;
			while(start <= end) {
				node_id  = Tree->groups[start];
				//printf("<<Node %d:",node_id);  // node id
				// compare states
				plh1 = plh + nos*nsites*node_id;
				plh2 = Trees->previous_plh + nos*nsites*node_id;
				for(i=0; i < (nos)*(nsites); i++,plh1++,plh2++) {
					diff = fabs(*plh1 - *plh2);
					//printf("%d (%lf %lf) %lf ",i,*plh1,*plh2,diff);
					if (diff > 0.0001) {
						printf("nsites = %d\n",nsites);
						printf("\n[ Group %d size %d:\n ",group,gsize);
						printf("<< Node %d,%d:",node_id,pos);
						printf("%d (%lf %lf) %lf \n",i,*plh1,*plh2,diff);
						exit(0);
						// Track down error
						start2 = Tree->childrenIdx[2*node_id];
						end2 =   start2 + Tree->childrenIdx[2*node_id+1] - 1;
						printf("<");
						if (start2==-1)
						printf("Tip");
						else {
							printf("Children:\n");
							while(start2 <= end2) {
								child_id = Tree->children[start2];
								printf("----- Child %d plh:",child_id);
								// printf plh contents
								plh1 = plh + nos*nsites*child_id;
								plh2 = Trees->previous_plh + nos*nsites*child_id;
								for(j=0; j < (nos)*(nsites); j++,plh1++,plh2++) {
									printf("(%lf %lf)",*plh1,*plh2);
								}
								printf("\nplhFactor: ");
								plhf = Trees->plhFactor + nos*nsites*child_id;
								plhfdebug = Trees->debug_plhFactor + nos*nsites*child_id;
								for(j=0; j < (nos)*(nsites); j++,plhf++,plhfdebug++) {
									printf("%lf [%lf] ",*plhf,*plhfdebug);
								}
								printf("\nPMatrix\n");
								btlin_printR(pmatrix+nos*nos*child_id,nos,nos);
								//printf("GPU PMatrix\n");
								//btlin_printR(Trees->check_pmatrix+nos*nos*child_id,nos,nos);
								start2++;
								// debug of debug

							}
						}
						exit(0);
					}

				}
				//printf(">>");
				start++;
				pos++;
			}
		}
	printf("\nResult OK");
	// end print group info

}


void comparePlh(TREES* Trees, TREE* Tree) {
		// print group information
	int start,end, group, node_id, i, nos,nsites;
	double* plh, *plh1,*plh2;
	double diff;
	int counter,gsize,pos;

	counter = 0;
	plh = Tree->NodeList[0]->Partial[0];
	nsites = Trees->NoSites;
	nos = Trees->NoStates;

	for(group=0;group <= Tree->NoFGroups; group++)
	{

			start = Tree->groupsIdx[2*group];
			end = start + Tree->groupsIdx[2*group+1] - 1;
			//printf("\n[ Group %d Size %d:\n ",group,end-start+1);
			gsize = end-start+1;
			pos=0;
			while(start <= end) {
				node_id  = Tree->groups[start];
				//printf("<<Node %d:",node_id);  // node id
				// compare states
				plh1 = plh + nos*nsites*node_id;
				plh2 = Trees->previous_plh + nos*nsites*node_id;
				for(i=0; i < (nos)*(nsites); i++,plh1++,plh2++) {
					diff = fabs(*plh1 - *plh2);
					//printf("%d (%lf %lf) %lf ",i,*plh1,*plh2,diff);
					if (diff > 0.0001) {
						printf("\n[ Group %d size %d:\n ",group,gsize);
						printf("<< Node %d,%d:",node_id,pos);
						printf("%d (%lf %lf) %lf \n",i,*plh1,*plh2,diff);
						counter++;
					}
					if (counter >1)
						return;
				}
				//printf(">>");
				start++;
				pos++;
			}



		}
	printf("\n");
	// end print group info

}

// plh1 <-- plh2
void copyPlh(TREES* Trees,TREE* Tree,double* plh1, double* plh2) {
	int i;
	double *plh,*plh_previous;
	plh = plh2;
	plh_previous = plh1;
	for(i=0; i < Trees->NoSites * Trees->NoStates * Tree->NoNodes; i++) {
		*plh_previous = *plh;
		plh_previous++;
		plh++;
	}

}








#endif