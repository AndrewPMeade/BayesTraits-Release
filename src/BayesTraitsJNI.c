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

MODEL	GetModelFromNo(int ModelNo)
{
	switch(ModelNo)
	{
		case 0:	return 	MULTISTATE;
		case 1: return	DESCINDEP;
		case 2: return 	DESCDEP;
		case 3: return	CONTINUOUSRR;
		case 4: return	CONTINUOUSDIR;
		case 5: return 	CONTINUOUSREG;
	}

	exit(0);
}

ANALSIS	AnalsisFromNo(int AnalisisNo)
{
	if(AnalisisNo == 0)
		return ANALML;

	if(AnalisisNo == 1)
		return ANALMCMC;

	exit(0);
}

void JavaProgress(int No)
{
	
}

int mainJNI(int Size, char** RunP)
{
	TREES*		Trees=NULL;
	OPTIONS*	Opt=NULL; 
	char*		TreeFN;
	char*		DataFN;
	MODEL		Model;
	ANALSIS		Analysis;

	
	TreeFN	= RunP[0];
	DataFN	= RunP[1];
	Model	= GetModelFromNo(atoi(RunP[2]));
	Analysis= AnalsisFromNo(atoi(RunP[3]));
	SetSeed();

	Trees  = LoadTrees(TreeFN);
	LoadData(DataFN, Trees);
	
	Opt = CreatOptions(Model, Analysis, Trees->NoOfStates, TreeFN, DataFN, Trees->SymbolList, Trees);
	GetOptionsArry(Opt, Size-4, &RunP[4]);

	PreProcess(Opt, Trees);

	if(Opt->Analsis == ANALMCMC)
		MCMC(Opt, Trees);

	if(Opt->Analsis == ANALML)
		FindML(Opt, Trees);

	FreeTrees(Trees, Opt);
	FreeOptions(Opt);

	return 0;
}

int	main(int argc, char **argv)
{
	TEXTFILE*	TF;
	int			Ret;

	TF = LoadTextFile(argv[1], FALSE);

	Ret = mainJNI(TF->NoOfLines, TF->Data);

	FreeTextFile(TF);
	return Ret;
}