#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#include "genlib.h"
#include "typedef.h"
#include "trees.h"
#include "data.h"
#include "options.h"
#include "initialise.h"
#include "mcmc.h"
#include "ml.h"

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
		return ANALML;

	if(No == 2)
		return ANALMCMC;

	printf("Unknown analysis type %s\n", AType);
	exit(0);
}

void	SetLogFName(OPTIONS *Opt, int BNo)
{
	char *Buffer;

	Buffer = (char*)malloc(sizeof(char) * BUFFERSIZE);
	if(Buffer == NULL)
		MallocErr();

	sprintf(Buffer, "BTBatchLog-%06d.txt", BNo);
	free(Opt->LogFN);

	Opt->LogFN = StrMake(Buffer);
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

	Opt = CreatOptions(Model, Analsis, Trees->NoOfStates, TreeFN, DataFN, Trees->SymbolList, Trees);
	SetLogFName(Opt, BNo);
//	Opt->LogFN = StrMake(LogFN);

	GetOptionsArry(Opt, NoComs-2, &Coms[2]);

	PrintOptions(stdout, Opt);

	CheckOptions(Opt);
	
	#ifdef BTOCL
	//printf("Loading kernels\n");
	if (btocl_load_all(Opt->ModelType == MT_CONTINUOUS,	Opt->ModelType == MT_DISCRETE,
			Trees->NoOfStates, Trees->NoOfSites) != 0) {
		printf("Error: Couldn't load OpenCL kernels\n");
		exit(0);
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

	Buffer = (char*)malloc(sizeof(char) * BUFFERSIZE);
	if(Buffer == NULL)
		MallocErr();


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

void	BatchRun(char *BatchFN)
{
	TEXTFILE *TF;
	int		Index, NoComs, BNo;
	char	*TreeFN, *DataFN, *CList;
	char	*TLine;
	char	**Buffer;
	char	**ComList;

	TF = LoadTextFile(BatchFN, FALSE);

	Buffer = (char**)malloc(sizeof(char*) * TF->MaxLine);
	if(Buffer == NULL)
		MallocErr();
	BNo = 1;
	for(Index=0;Index<TF->NoOfLines;Index++)
	{
		TLine = StrMake(TF->Data[Index]);

		if(PassBatchLine(TF->Data[Index], &TreeFN, &DataFN, &CList) == TRUE)
		{
			ComList = MakeCommands(CList, &NoComs);
			printf("\n\n\n\n\nBatch Command:\t%d\t%s\n", Index, TLine);fflush(stdout);
			BatchRunLine(Index, TreeFN, DataFN,  ComList, NoComs);
			free(ComList);
		}
		free(TLine);
	}

	free(Buffer);
	FreeTextFile(TF);
}