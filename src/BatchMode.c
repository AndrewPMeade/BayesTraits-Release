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

#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#include "GenLib.h"
#include "TypeDef.h"
#include "Trees.h"
#include "Data.h"
#include "Options.h"
#include "Initialise.h"
#include "MCMC.h"
#include "ML.h"
#include "BatchMode.h"

#ifdef BTOCL
#include "btocl_runtime.h"
#include "btocl_kernels_bayestraits.h"
#endif

char**	MakeCommands(char *Command, int *NoC)
{
	char **Ret;
	int Len;

	Len = (int)strlen(Command);
	*NoC = 0;

	Ret = (char**)malloc(sizeof(char*) * (Len + 1));
	if(Ret == NULL)
		MallocErr();

	*NoC = MakeArgvChar(Command, Ret, Len , ';');

	return Ret;
}

MODEL	GetBatchModel(TREES *Trees, char *MLine)
{
	MODEL Ret;
	int MNo;

	if(IsValidInt(MLine) == FALSE)
	{
		printf("%s is not a valid model number\n", MLine);
		exit(0);
	}

	MNo = atoi(MLine);
	Ret = IntToModel(MNo);

	if(ValidModelChoice(Trees, Ret) == FALSE)
	{
		printf("%s is not a valid model\n", MLine);
		exit(0);
	}

	return Ret;
}

ANALSIS GetBatchAnalsis(char *AType)
{
	int No;

	No = atoi(AType);

	if(No == 1)
		return ANALYSIS_ML;

	if(No == 2)
		return ANALYSIS_MCMC;

	printf("Unknown analysis type %s\n", AType);
	exit(0);
}

void	SetLogFName(OPTIONS *Opt, int BNo)
{
	char *Buffer;

	Buffer = (char*)SMalloc(sizeof(char) * 64);

	sprintf(Buffer, "BTBatchLog-%06d.txt", BNo);
	free(Opt->BaseOutputFN);

	Opt->BaseOutputFN = StrMake(Buffer);
	free(Buffer);
}

void	BatchRunLine(int BNo, char *TreeFN, char *DataFN, char **Coms, int NoComs)
{
	int			NoSites;
	TREES*		Trees;
	OPTIONS*	Opt;
	MODEL		Model;
	ANALSIS		Analsis;

	Trees  = LoadTrees(TreeFN);
	LoadData(DataFN, Trees);

	Model = GetBatchModel(Trees, Coms[0]);

	Analsis = GetBatchAnalsis(Coms[1]);

	CheckDataWithModel(DataFN, Trees, Model);

	PreProcessDataWithModel(Trees, Model);

	Opt = CreatOptions(Model, Analsis, Trees->NoStates, TreeFN, DataFN, Trees->SymbolList, Trees);
	SetLogFName(Opt, BNo);

	GetOptionsArry(Opt, NoComs-2, &Coms[2], Trees);

	PrintOptions(stdout, Opt, Trees);

	#ifdef BTOCL
	//printf("Loading kernels\n");
	if (btocl_load_all(Opt->ModelType == MT_CONTINUOUS,	Opt->ModelType == MT_DISCRETE, Trees->NoStates, Trees->NoSites) != 0)
	{
		printf("Error: Couldn't load OpenCL kernels\n");
		exit(0);
	}
	#endif

	PreProcess(Opt, Trees);

	if(Opt->Analsis == ANALYSIS_MCMC)
		MCMC(Opt, Trees);

	if(Opt->Analsis == ANALYSIS_ML)
		FindML(Opt, Trees);

	NoSites = Trees->NoSites;
	FreeTrees(Trees, Opt);
	FreeOptions(Opt, NoSites);

	#ifdef BTOCL
	//printf("Removing kernels\n");
	btocl_clear_kernels(btocl_getruntime());
	#endif

//	exit(0);
}

int		PassBatchLine(char *Line, char **TreeFN, char **DataFN, char **ComList)
{
	char	*FLine;
	char	*Buffer;

	FLine = Line;

	if(CountChar(Line, '\t') < 2)
	{
		printf("Error Passing Line %s\n", Line);
		return FALSE;
	}

	Buffer = (char*)SMalloc(sizeof(char) * BUFFERSIZE);


	*TreeFN = Line;

	while(*Line != '\t')
		Line++;
	*Line = '\0';
	Line++;

	*DataFN = Line;

	while(*Line != '\t')
		Line++;
	*Line = '\0';
	Line++;

	*ComList = Line;

	free(Buffer);

	return TRUE;
}


/*
void	BatchRun(char *BatchFN)
{
	TEXTFILE *TF;
	int		Index, Tokes, NoComs, BNo;
	char	**Buffer;
	char	**ComList;

	TF = LoadTextFile(BatchFN, FALSE);

	Buffer = (char**)malloc(sizeof(char*) * TF->MaxLine);
	if(Buffer == NULL)
		MallocErr();
	BNo = 1;
	for(Index=0;Index<TF->NoOfLines;Index++)
	{

		PassBatchLine(TF->Data[Index],
		Tokes = MakeArgvChar(TF->Data[Index], Buffer, TF->MaxLine, '\t');
		if(Tokes == 3)
		{
			ComList = MakeCommands(Buffer[2], &NoComs);
			BatchRunLine(Index, Buffer[0], Buffer[1], ComList, NoComs);
			free(ComList);
		}
	}

	free(Buffer);
	FreeTextFile(TF);
}

*/

void	RunLine(char *Line, int Index)
{
	char *TreeFN, *DataFN, *CList;
	char **ComList;
	int NoComs;

	if(PassBatchLine(Line, &TreeFN, &DataFN, &CList) == TRUE)
	{
		ComList = MakeCommands(CList, &NoComs);
		printf("\n\n\n\n\nBatch Command:\t%d\t%s\n", Index, Line);
		fflush(stdout);
		BatchRunLine(Index, TreeFN, DataFN,  ComList, NoComs);
		free(ComList);
	}
	
}

void	BatchRun(char *BatchFN)
{
	TEXTFILE *TF;
	int		Index, BNo;
	char	*TLine;
	char	**Buffer;

	TF = LoadTextFile(BatchFN, FALSE);

	Buffer = (char**)SMalloc(sizeof(char*) * TF->MaxLine);
	BNo = 1;

// this will not always work. Models which have dependencies, such as save / load models. 
#ifdef PBATCH
	#pragma omp parallel for private(TLine)
#endif
	for(Index=0;Index<TF->NoOfLines;Index++)
	{
		TLine = StrMake(TF->Data[Index]);
		printf("Running Line: %d\t%s\n", Index, TLine);
		RunLine(TLine, Index);
		free(TLine);
	}

	free(Buffer);
	FreeTextFile(TF);


	printf("Batch run complete.\n");
}