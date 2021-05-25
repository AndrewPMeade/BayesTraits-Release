#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "typedef.h"
#include "initialise.h"
#include "typedef.h"
#include "trees.h"
#include "data.h"
#include "options.h"
#include "Rates.h"
#include "likelihood.h"
#include "RandLib.h"
#include "priors.h"
#include "mcmc.h"
#include "praxis.h"
#include "ml.h"
#include "genlib.h"
#include "continuous.h"
#include "initialise.h"
#include "VarRates.h"
#include "BigLh.h"
#include "ptrees.h"
#include "Threaded.h"
#include "QuadDouble.h"
#include "contrasts.h"
#include "FatTail.h"

#ifdef BTOCL
#include "btocl_discrete.h"
#endif

OPTIONS*	SetUpOptions(TREES* Trees, char	*TreeFN, char *DataFN)
{
	OPTIONS		*Opt;
	MODEL		Model;
	ANALSIS		Analsis;


	Model	= GetModel(Trees);
	
	if(Model == M_FATTAIL)
		Analsis = ANALMCMC;
	else
		Analsis = GetAnalsis(Trees);

	CheckDataWithModel(DataFN, Trees, Model);

	Opt = CreatOptions(Model, Analsis, Trees->NoOfStates, TreeFN, DataFN, Trees->SymbolList, Trees);
	
	return Opt;
}

void	PreProcess(OPTIONS *Opt, TREES* Trees)
{
	int		Index;
	int		ID;

	if(Opt->ScaleTrees != -1.0)
		ScaleUserTrees(Trees, Opt->ScaleTrees);

	if(Opt->Model == M_CONTRAST_REG)
	{
		if(Opt->TestCorrel == FALSE)
			SetDataRegTC(Opt);
	}
	
	SetNoOfThreads(Opt->Cores);

	if(Opt->Stones != NULL)
		Opt->Stones->ItStart = Opt->Itters + 1;

	FlattenRecNode(Opt);
	
	for(Index=0;Index<Trees->NoOfTrees;Index++)
	{
		ID = 0;
		SetNodeIDs(Trees->Tree[Index]);
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

	SetPTrees(Opt, Trees);

	if(Opt->ModelType == MT_CONTINUOUS)
		InitContinus(Opt, Trees);

	if(Opt->ModelType == MT_CONTRAST)
		InitContrastAll(Opt, Trees);

	if(Opt->ModelType == MT_FATTAIL)
		InitFatTailTrees(Opt, Trees);

	if(Opt->ModelType == MT_DISCRETE)
	{
//		NormaliseTrees(Trees->NormConst, Trees);
		
		if(Opt->UseCovarion == TRUE)
			Trees->NoOfStates = Trees->NoOfStates * 2;

		if(Opt->Model == M_DESCCV)
			Trees->NoOfStates = Trees->NoOfStates * 2;

		if((Opt->UseKappa == TRUE) && (Opt->FixKappa != -1))
		{
			for(Index=0;Index<Trees->NoOfTrees;Index++)
				TreeBLToPower(Trees, Trees->Tree[Index], Opt->FixKappa);

			Opt->FixKappa = -1;
			Opt->UseKappa = FALSE;
		}

		AllocPartial(Opt, Trees, Opt->UseGamma);
		AllocLHInfo(Trees, Opt);

		SetFossiles(Trees, Opt);

		SetNOSPerSite(Opt);

		InitTreeBigLh(Opt, Trees);

#ifdef QUAD_DOUBLE
		InitQuadDoubleLh(Opt, Trees);
#endif

#ifdef BTOCL
		btocl_AllocPMatrixInfo(Trees);
		//btocl_AllocLhInfo(Trees);
#endif

	}

	if(FindNoEstDataPoint(Opt, Trees) > 0)
		Opt->EstData = TRUE;
	else
		Opt->EstData = FALSE;
		
	if(Opt->SaveTrees != NULL)
		SaveTrees(Opt->SaveTrees, Opt->Trees);
}
