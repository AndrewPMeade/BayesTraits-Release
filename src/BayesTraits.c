#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "typedef.h"
#include "trees.h"
#include "data.h"
#include "options.h"
#include "rates.h"
#include "likelihood.h"
#include "priors.h"
#include "mcmc.h"
#include "praxis.h"
#include "ml.h"
#include "genlib.h"
#include "continuous.h"
#include "initialise.h"
#include "RandLib.h"
#include "BatchMode.h"


#include "mathlib.h"
//#include "./MathLib/mconf.h"

#ifdef	THREADED
	#include <omp.h>
#endif


#ifdef BTLAPACK
#include "btlapack_interface.h"
#endif

#ifdef BTOCL
#include "btocl_runtime.h"
#include "btocl_runtime_kernels.h"
#include "btdebug.h"
#endif

// #include "btdebug.h"   -- igor

extern void BayesModeTest(OPTIONS *Opt, TREES *Trees);
extern void PMatrixTest(void);
/*
	./Seq/AllDiono.trees ./Seq/AllDiono.txt < in.txt
	./Seq/allnodes.nex.tre ./Seq/InputData.txt < in.txt	 > sout.txt
	./Seq/a67.1000.trees ./Seq/descent.txt < a67nodes.txt > sout.txt
	./Seq/Primates.trees ./Seq/Primates-Data1.txt < in.txt > sout.txt
	./Seq/ExtantTree.tre ./Seq/ExtantData.txt < din.txt > sout.txt
	./Seq/AllDiono-Time.trees ./Seq/AllDiono.txt < din.txt > sout.txt
	./Seq/ATrees.trees ./Seq/Data.txt < in.txt > sout.txt
	./Seq/Tree1.trees ./Seq/body.txt < in.txt > sout.txt
	./Seq/MamTrees-50.trees ./Seq/MamData.txt < in.txt > sout.txt

	./Seq/Extra-ExtantTree.tre ./Seq/Extra-ExtantData.txt < in.txt > sout.txt
	./Seq/ExtantTree.tre ./Seq/ExtantData.txt < in.txt > sout.txt
	./Seq/Primates.trees ./Seq/Primates.txt < in.txt > sout.txt
	./Seq/TPrimates.trees ./Seq/TPrimates.txt < in.txt > sout.txt
	./Seq/treekjer_fbi.nex ./Seq/dat_kjer_fbi.tab < in.txt > sout.txt

	./Seq/Monarch.trees ./Seq/Monarch.txt < in.txt > sout.txt

	./Seq/MamTrees-50.trees ./Seq/MamData-Missing.txt < in.txt > sout.txt
	./Seq/MonarchRec.trees ./Seq/Monarch-T1-T2.txt < in2.txt > sout.txt 

	./Seq/MonarchRec.trees ./Seq/Monarch-T1.txt <in2.txt >sout.txt
	./Seq/MonarchRec.trees ./Seq/Monarch-T1.txt < in.txt > sout.txt
	./Seq/IE-M1P-RS.trees ./Seq/IE-MS.nex-0007.txt  < in.txt > sout.txt

	./Seq/IECon.trees ./Seq/
	./Seq/IECon.trees ./Seq/IE-MSAllCogs.txt < in.txt > sout.txt

	./Seq/CIE.trees ./Seq/IE-Lex.txt < in.txt > sout.txt

	./Seq/bayes3.trees ./Seq/SerSD.txt  < in.txt > sout.txt
	./Seq/TestTree.trees ./Seq/TestTree.txt < TestIn.txt > sout.txt

	./Seq/MamTrees-50.trees ./Seq/MamData.txt < in.txt > sout.txt

	./Seq/WConTree.trees ./Seq/Wres.txt < in.txt > sout.txt

	./Seq/SmallIE.trees ./Seq/SmallIELex.txt < in.txt > sout.txt

	./Seq/Primates.trees ./Seq/PrimatesEstData.txt < in.txt > sout.txt
	./Seq/9tree.trees ./Seq/9F.txt < in.txt > sout.txt
	./Seq/AustCon.trees ./Seq/0204.txt < in.txt > sout.txt

	./Seq/Mammal-ArtPrim.trees ./Seq/Mammal-ArtPrim.txt < in.txt > sout.txt

	./Seq/MamTrees-50.trees ./Seq/MamDataS1.txt < in.txt > sout.txt
	./Seq/LhInDep.trees ./Seq/LhInDep.txt < in.txt > sout.txt
	./Seq/ContrastTestLh.trees ./Seq/ContrastTestLh.txt < in.txt > sout.txt
	./Seq/MammalBig.trees ./Seq/MammalBig.txt < in.txt > sout.txt
	./Seq/MamBigTrim.trees ./Seq/MamBigTrim.txt < in.txt > sout.txt

	./Seq/MamTreesPoly-50.trees ./Seq/MamDataS1.txt < in.txt > sout.txt

	 ./Seq/Wing.trees ./Seq/Wing.txt < in.txt > sout.txt
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
		*TreeFN = StrMake(argv[1]);
		*DataFN = StrMake(argv[2]);
		return;
	}
	
	printf("BayesTraits take a tree file and a data file, it is run form the command line.\nPlease read the manual for more information.\n");
	printf("Press enter to leave.\n");
	fgets(&Line[0], 64, stdin);
	exit(0);
}

// Full optermisation
//	cl /Ox /Oi /Ob2 /Ot /Oy /GL /w *.c ./MathLib/*.c

// gcc -O3 -fomit-frame-pointer -lm 

// Big Lh + OpenMP
// gcc *.c -lm -O3 -DBIG_LH -lmpfr -lgmp -fomit-frame-pointer -static -DTHREADED -fopenmp
// gcc *.c -lm -O3 -DBIG_LH -lmpfr -lgmp -fomit-frame-pointer -static -DTHREADED -fopenmp -Dwarn _unused_result 
 
// ./Seq/MamTrees-1.trees ./Seq/MamData.txt < in.txt > sout.txt
// ./Seq/MamBrainBody.trees ./Seq/MamBrainBody.txt < in.txt. > sout.txt
// ./Seq/FritzMammalianSupertree.trees ./Seq/JustBM.txt < in.txt > sout.txt
// ./Seq/Primates.trees ./Seq/Primates.txt < in.txt > sout.txt

int main(int argc, char** argv)
{
	TREES*		Trees;
	OPTIONS*	Opt; 
	int			NoSites;
	char		*TreeFN, *DataFN; 

	DISPLAY_INFO;

	//btdebug_init();
	
	NoSites = 0;

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

	if(Trees->NoOfTrees == 0) 
	{
		printf("Could not load any valid trees\n");
		exit(0);
	}

	LoadData(DataFN, Trees);
	
	Opt = SetUpOptions(Trees, TreeFN, DataFN);

	PrintOptions(stdout, Opt);

	GetOptions(Opt);
	CheckOptions(Opt);
	
	#ifdef BTOCL
	btocl_init_runtime(CL_DEVICE_TYPE_GPU);
	//btocl_load_all(Opt,Trees);
	if (btocl_load_all(Opt->ModelType == MT_CONTINUOUS,	Opt->ModelType == MT_DISCRETE,
			Trees->NoOfStates, Trees->NoOfSites) != 0) {
		printf("Error: Couldn't load OpenCL kernels\n");
		return 1;
	}
	#endif
	
	PreProcess(Opt, Trees);

	if(Opt->Analsis == ANALMCMC)
		MCMC(Opt, Trees);

	if(Opt->Analsis == ANALML)
		FindML(Opt, Trees);

	NoSites = Trees->NoOfSites;
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

