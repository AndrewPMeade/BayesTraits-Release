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

#include "./MathLib/dcalc.h"
#include "./MathLib/mconf.h"

extern void BayesModeTest(OPTIONS *Opt, TREES *Trees);

extern void PMatrixTest(void);


OPTIONS*	SetUpOptions(TREES* Trees, char	*TreeFN, char	*DataFN)
{
	OPTIONS*	Opt=NULL;
	MODEL		Model;
	ANALSIS		Analsis;


	Model	= GetModel(Trees);
	Analsis = GetAnalsis(Trees);

	CheckDataWithModel(DataFN, Trees, Model);

	if((Model == DESCDEP) || (Model == DESCINDEP))
		SquashDep(Trees);

	if((Model == CONTINUOUSRR) || (Model == CONTINUOUSDIR) || (Model == CONTINUOUSREG))
		RemoveConMissingData(Trees);

	Opt = CreatOptions(Model, Analsis, Trees->NoOfStates, TreeFN, DataFN, Trees->SymbolList, Trees);

	return Opt;
}

void	PreProcess(OPTIONS *Opt, TREES* Trees)
{
	int		Index;

	Opt->LogFile		= OpenWrite(Opt->LogFN);
	Trees->UseCovarion	= Opt->UseCovarion;

	if(Opt->DataType == CONTINUOUS)
		InitContinus(Opt, Trees);
	else
	{
		if(Opt->UseCovarion == TRUE)
			Trees->NoOfStates = Trees->NoOfStates * 2;

		if((Opt->UseKappa == TRUE) && (Opt->FixKappa != -1))
		{
			for(Index=0;Index<Trees->NoOfTrees;Index++)
				TreeBLToPower(Trees, &Trees->Tree[Index], Opt->FixKappa);

			Opt->FixKappa = -1;
			Opt->UseKappa = FALSE;
		}

		AllocPartial(Trees, Opt->UseGamma);
		AllocLHInfo(Trees, Opt);

		SetFossiles(Trees, Opt);
	}
}

/*
	./Seq/AllDiono.trees ./Seq/AllDiono.txt < in.txt
	./Seq/allnodes.nex.tre ./Seq/InputData.txt < in.txt	 > sout.txt
	./Seq/a67.1000.trees ./Seq/descent.txt < a67nodes.txt > sout.txt
	./Seq/primeates.trees ./Seq/Primates-Data1.txt < in.txt > sout.txt
	./Seq/ExtantTree.tre ./Seq/ExtantData.txt < din.txt > sout.txt
	./Seq/AllDiono-Time.trees ./Seq/AllDiono.txt < din.txt > sout.txt
	./Seq/ATrees.trees ./Seq/Data.txt < in.txt > sout.txt
	./Seq/Tree1.trees ./Seq/body.txt < in.txt > sout.txt
	./Seq/MamTrees-1.trees ./Seq/MamData.txt < in.txt > sout.txt

	./Seq/Extra-ExtantTree.tre ./Seq/Extra-ExtantData.txt < in.txt > sout.txt
	./Seq/ExtantTree.tre ./Seq/ExtantData.txt < din.txt > sout.txt
	./Seq/Primates.trees ./Seq/Primates.txt < in.txt > sout.txt
	./Seq/treekjer_fbi.nex ./Seq/dat_kjer_fbi.tab < in.txt > sout.txt
*/

int main(int argc, char** argv)
{
	TREES*		Trees=NULL;
	OPTIONS*	Opt=NULL; 

	SetSeed();

	if(argc != 3)
	{
		printf("The program takes 2 paramiters a tree file and a data file\n");
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

	if(Opt->Analsis == ANALMCMC)
		MCMC(Opt, Trees);

	if(Opt->Analsis == ANALML)
		FindML(Opt, Trees);

	FreeTrees(Trees, Opt);

	FreeOptions(Opt);

	return 0;
}