#include <stdio.h>
#include <stdlib.h>

#include "GenLib.h"
#include "TypeDef.h"
#include "Output.h"
#include "Options.h"
#include "Rates.h"
#include "Priors.h"
#include "VarRates.h"
#include "RJDummy.h"
#include "Schedule.h"
#include "Stones.h"
#include "ModelFile.h"
#include "FatTail.h"
#include "IntraNode.h"
#include "Initialise.h"
#include "Trees.h"


void PrintMCMCStdOutHeader(OPTIONS *Opt)
{
	printf("Iteration\tLh\tLh Prior\tElapsed Seconds\tState\n");
	fflush(stdout);
}

void PrintMLStdOutHeader(void)
{
	printf("Tree No\tLh\tElapsed Seconds\n");
	fflush(stdout);
}

void	PrintPriorHeadder(FILE* Str, OPTIONS *Opt, RATES* Rates)
{
	int		PIndex;
	PRIOR	*Prior;

	for(PIndex=0;PIndex<Rates->NoPriors;PIndex++)
	{
		Prior = Rates->Priors[PIndex];

		if(Prior->UseHP == TRUE)
		{
			switch(Prior->Dist)
			{
				case PDIST_EXP:
					fprintf(Str, "%s - Mean\t", Prior->Name);
				break;

				case PDIST_GAMMA:
					fprintf(Str, "%s - Shape\t%s - Scale\t", Prior->Name, Prior->Name);
				break;

				case PDIST_LOGNORMAL:
					fprintf(Str, "%s - Mue\t%s - Sig\t", Prior->Name, Prior->Name);
				break;

				default:
					printf("%s - %s::%d Hyper Prior not supported.", Prior->Name, __FILE__, __LINE__);
					exit(0);
			}
		}
	}

	fprintf(Str, "\n");
}


FILE*	CreatStoneOuput(OPTIONS *Opt, STONES *Stones)
{
	FILE*	Ret;

	Ret = OpenWriteWithExt(Opt->BaseOutputFN, OUTPUT_EXT_STONES);

	OutputStoneHeadder(Ret, Stones);

	return Ret;
}

FILE*	SetScheduleFile(OPTIONS *Opt, SCHEDULE*	Shed)
{
	FILE *Ret;

	Ret = OpenWithExt(Opt->CheckPointAppendFiles, Opt->BaseOutputFN, OUTPUT_EXT_SCHEDULE);

	if(Opt->CheckPointAppendFiles == TRUE)
		return Ret;

	PrintShedHeadder(Opt, Shed, Ret);

	return Ret;
}


void	SetLogFile(OPTIONS *Opt, TREES *Trees, RATES *Rates)
{
	Opt->LogFile = OpenWithExt(Opt->CheckPointAppendFiles, Opt->BaseOutputFN, OUTPUT_EXT_LOG);

	if(Opt->CheckPointAppendFiles == TRUE)
		return;

	PrintOptions(Opt->LogFile, Opt, Trees);
	PrintRatesHeadder(Opt->LogFile, Opt, Trees);
	PrintPriorHeadder(Opt->LogFile, Opt, Rates);
	fflush(Opt->LogFile);
}

void	SetOutputFile(OPTIONS *Opt, TREES *Trees, RATES *Rates, SCHEDULE* Shed, STONES *Stones)
{
	SetLogFile(Opt, Trees, Rates);

	PrintOptions(stdout, Opt, Trees);

	if(Opt->SaveTrees == TRUE)
		InitialiseOutputTrees(Opt, Trees);

	if(Opt->Analsis == ANALYSIS_ML)
	{
		PrintMLStdOutHeader();
		return;
	}

	PrintMCMCStdOutHeader(Opt);

	if(UseNonParametricMethods(Opt) == TRUE)
		InitVarRatesFiles(Opt, Trees, Rates);

	if(Opt->RJDummy == TRUE)
		InitRJDummyFile(Opt, Trees);

	if(Opt->UseSchedule == TRUE)
		Opt->ShedFile = SetScheduleFile(Opt, Shed);

	if(Opt->SaveModels == TRUE)
		Opt->SaveModelFile = InitSaveModelFile(Opt->SaveModelsFN, Opt, Trees, Rates);

	if(Stones != NULL)
		Opt->StoneFile = CreatStoneOuput(Opt, Stones);

	if(Opt->ModelType == MT_FATTAIL)
	{
		InitFattailFile(Opt, Trees);

		if(Opt->UseIntraNode == TRUE)
			InitIntraNodeFile(Opt, Trees);
	}


}
