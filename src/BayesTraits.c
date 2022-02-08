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

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "TypeDef.h"
#include "Trees.h"
#include "Data.h"
#include "Options.h"
#include "Rates.h"
#include "Likelihood.h"
#include "Priors.h"
#include "MCMC.h"
#include "Praxis.h"
#include "ML.h"
#include "GenLib.h"
#include "Continuous.h"
#include "Initialise.h"
#include "RandLib.h"
#include "BatchMode.h"
#include "AncStateV.h"

#ifdef	OPENMP_THR
	#include <omp.h>
#endif


#ifdef BTLAPACK
	#include "btlapack_interface.h"
#endif

#ifdef BTOCL
	#include "btocl_runtime.h"
	#include "btocl_runtime_kernels.h"
	#include "btdebug.h"
	#include "btocl_kernels_bayestraits.h"
#endif

// #include "btdebug.h"   -- igor

extern void BayesModeTest(OPTIONS *Opt, TREES *Trees);
extern void PMatrixTest(void);
/*

*/

#ifdef JNIRUN

int main(int argc, char** argv)
{
}

#else
/*
int		ValidTaxa(NODE N)
{
	int Index;
	NODE Ans, NN;

	if(N->Ans == NULL)
		return TRUE;

	Ans = N->Ans;

	for(Index=0;Index<Ans->NoNodes;Index++)
	{
		NN = Ans->NodeList[Index];
		if(NN != N)
		{
			if(NN->Tip == TRUE)
			{
				if(NN->Taxa->ConData[0] == N->Taxa->ConData[0])
					return FALSE;
			}
		}
	}

	return TRUE;
}

void	TestTreeLh(OPTIONS *Opt, TREES *Trees)
{
	TREE *Tree;
	int	Index;
	NODE N;

	Tree = Trees->Tree[0];

	for(Index=0;Index<Tree->NoNodes;Index++)
	{
		N = Tree->NodeList[Index];
		if(N->Tip == TRUE)
		{

			if(ValidTaxa(N) == FALSE)
			{
				printf("%s\t%f\tNo\n", N->Taxa->Name, N->Taxa->ConData[0]);
				N->Taxa->ConData[0] = -1;
			}
			else
				printf("%s\t%f\tYes\n", N->Taxa->Name, N->Taxa->ConData[0]);
		}
	}

	printf("====================\n");

	for(Index=0;Index<Trees->NoOfRemovedTaxa;Index++)
		printf("%s\t-\n", Trees->RemovedTaxa[Index]);
	exit(0);
}
*/

void GetTreeDataF(int argc, char** argv, char **TreeFN, char **DataFN)
{
	char Line[1024];

	if(argc == 3)
	{
		(*TreeFN) = StrMake(argv[1]);
		(*DataFN) = StrMake(argv[2]);
		return;
	}

	printf("BayesTraits take a tree file and a data file, it is run form the command line.\nPlease read the manual for more information.\n");
	printf("Press enter to leave.\n");
	fgets(&Line[0], 64, stdin);
	exit(0);
}
/*
PRIOR*		CreateWeibullPrior(char *Name, double Scale, double Exponent);
PRIOR*		CreateLogNormalPrior(char *Name, double Location, double Scale);
PRIOR*		CreateGammaPrior(char *Name, double Shape, double Scale);
PRIOR*		CreateUniformPrior(char *Name, double Min, double Max);
PRIOR*		CreateChiPrior(char *Name, double Mean);
PRIOR*		CreateExpPrior(char *Name, double Alpha);
PRIOR*		CreateSGammaPrior(char *Name, double Alpha, double Beta);

PRIOR*		CreateNormalPrior(char *Name, double Mean, double SD);
*/

// Full optermisation
//	cl /Ox /Oi /Ob2 /Ot /Oy /GL /w *.c ./MathLib/*.c

// gcc -O3 -fomit-frame-pointer -lm

// Big Lh + OpenMP
// gcc *.c -lm -O3 -DBIG_LH -lmpfr -lgmp -fomit-frame-pointer -static -DOPENMP_THR -fopenmp
// gcc *.c -lm -O3 -DBIG_LH -lmpfr -lgmp -fomit-frame-pointer -static -DOPENMP_THR -fopenmp -Dwarn _unused_result


// Threaded + quad math
// gcc *.c -O3 -fomit-frame-pointer -lgsl -DQUAD_DOUBLE -DOPENMP_THR -lquadmath -fopenmp


// ./Seq/francoisT.nex.trees ./Seq/francoisT.txt < ./Seq/in.txt > ./Seq/sout.txt
int main(int argc, char** argv)
{
	TREES*		Trees;
	OPTIONS*	Opt;
	char		*TreeFN, *DataFN;
	int			NoSites;

//	FatTailTest(argc, argv);

	DISPLAY_INFO;

	//btdebug_init();

	if(argc == 2)
	{
		#ifdef BTOCL
		btocl_init_runtime(CL_DEVICE_TYPE_GPU);
		#endif

		BatchRun(argv[1]);

		#ifdef BTOCL
		btocl_free_runtime();
		#endif

		return 0;
	}
		
	GetTreeDataF(argc, argv, &TreeFN, &DataFN);

	Trees  = LoadTrees(TreeFN);

	if(Trees->NoTrees == 0)
	{
		printf("Could not load any valid trees\n");
		exit(0);
	}

	LoadData(DataFN, Trees);

	Opt = SetUpOptions(Trees, TreeFN, DataFN);

	PrintOptions(stdout, Opt);
	
	GetOptions(Opt);
	CheckOptions(Opt);

//	MakeAncVMatrix(Opt, Trees);

	#ifdef BTOCL
	btocl_init_runtime(CL_DEVICE_TYPE_GPU);
	//btocl_load_all(Opt,Trees);
	if (btocl_load_all(Opt->ModelType == MT_CONTINUOUS,	Opt->ModelType == MT_DISCRETE, Trees->NoStates, Trees->NoSites) != 0)
	{
		printf("Error: Couldn't load OpenCL kernels\n");
		return 1;
	}
	#endif

	PreProcess(Opt, Trees);
		
	if(Opt->Analsis == ANALMCMC)
		MCMC(Opt, Trees);

	if(Opt->Analsis == ANALML)
		FindML(Opt, Trees);

	Finalise(Opt, Trees);

	NoSites = Trees->NoSites;
	FreeTrees(Trees, Opt);
	FreeOptions(Opt, NoSites);

	free(DataFN);
	free(TreeFN);

	#ifdef BTOCL
	btocl_free_runtime();
	#endif

	return 0;
}

#endif

