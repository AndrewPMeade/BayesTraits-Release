#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "typedef.h"
#include "initialise.h"
#include "typedef.h"
#include "trees.h"
#include "data.h"
#include "options.h"
#include "rates.h"
#include "likelihood.h"
#include "RandLib.h"
#include "priors.h"
#include "mcmc.h"
#include "praxis.h"
#include "ml.h"
#include "genlib.h"
#include "continuous.h"
#include "initialise.h"
#include "phyloplasty.h"
#include "BigLh.h"

OPTIONS*	SetUpOptions(TREES* Trees, char	*TreeFN, char *DataFN)
{
	OPTIONS*	Opt=NULL;
	MODEL		Model;
	ANALSIS		Analsis;

	Model	= GetModel(Trees);
	Analsis = GetAnalsis(Trees);
	CheckDataWithModel(DataFN, Trees, Model);

	Opt = CreatOptions(Model, Analsis, Trees->NoOfStates, TreeFN, DataFN, Trees->SymbolList, Trees);

	return Opt;
}

void	PreProcess(OPTIONS *Opt, TREES* Trees)
{
	int		Index;
	int		ID;

	FlattenRecNode(Opt);
	
	for(Index=0;Index<Trees->NoOfTrees;Index++)
	{
		ID = 0;
		SetNodeIDs(Trees->Tree[Index]->Root, &ID);
	}
	

	Opt->LogFile		= OpenWrite(Opt->LogFN);

	#ifdef JNIRUN
		Opt->LogFileRead = OpenRead(Opt->LogFN);
		Opt->LogFileBuffer = (char*)malloc(sizeof(char) * LOGFILEBUFFERSIZE);
		if(Opt->LogFileBuffer == NULL)
			MallocErr();
		Opt->PassedOut = (char**)malloc(sizeof(char*) * LOGFILEBUFFERSIZE);
		if(Opt->PassedOut == NULL)
			MallocErr();
	#endif

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
				TreeBLToPower(Trees, Trees->Tree[Index], Opt->FixKappa);

			Opt->FixKappa = -1;
			Opt->UseKappa = FALSE;
		}

		AllocPartial(Trees, Opt->UseGamma);
		AllocLHInfo(Trees, Opt);

		SetFossiles(Trees, Opt);

		SetNOSPerSite(Opt);

		InitTreeBigLh(Opt, Trees);
	}
}
