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
#include "rand.h"
#include "priors.h"
#include "mcmc.h"
#include "praxis.h"
#include "ml.h"
#include "genlib.h"
#include "continuous.h"
#include "initialise.h"

#include "./MathLib/dcalc.h"
#include "./MathLib/mconf.h"

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
*/

#ifdef JNIRUN

int main(int argc, char** argv)
{
}

#else

int main(int argc, char** argv)
{
	TREES*		Trees=NULL;
	OPTIONS*	Opt=NULL; 

	SetSeed();

	if(argc != 3)
	{
		printf("The program takes 2 parmeters a tree file and a data file\n");
		exit(0);
	}

	Trees  = LoadTrees(argv[1]);

	if(Trees->NoOfTrees == 0) 
	{
		printf("Could not load any valid trees\n");
		exit(0);
	}

	LoadData(argv[2], Trees);
	
	Opt = SetUpOptions(Trees, argv[1], argv[2]);
	
	PrintOptions(stdout, Opt);

	GetOptions(Opt);
	
	PreProcess(Opt, Trees);

//	LhOverAllModels(Opt, Trees); 

	if(Opt->Analsis == ANALMCMC)
		MCMC(Opt, Trees);

	if(Opt->Analsis == ANALML)
		FindML(Opt, Trees);

	FreeTrees(Trees, Opt);
	FreeOptions(Opt);

	return 0;
}

#endif