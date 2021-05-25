#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

#include "genlib.h"
#include "typedef.h"
#include "modelfile.h"
#include "likelihood.h"
#include "rates.h"

int		GetParamPerModel(OPTIONS *Opt, TREES *Trees, RATES *Rates)
{
	int Ret;

	Ret = 0;

	if(Opt->DataType == CONTINUOUS)
	{
		Ret = Rates->NoOfRates;
		
		if(Opt->EstKappa == TRUE)
			Ret++;

		if(Opt->EstLambda == TRUE)
			Ret++;

		if(Opt->EstDelta == TRUE)
			Ret++;

		if(Opt->EstOU == TRUE)
			Ret++;

		return Ret;
	}
	else
	{
		Ret += Rates->NoOfFullRates;
		Ret += Trees->NoOfStates;
		if(Opt->UseCovarion == TRUE)
			Ret += 2;

		if(Opt->UseGamma == TRUE)
			Ret++;

		if(Opt->EstKappa == TRUE)
			Ret++;
	}

	return Ret;
}

void	MFWriteInt(int D, FILE *File)
{
	fwrite((void*)&D, sizeof(int), 1, File);
}

void	MFWriteDouble(double D, FILE *File)
{
	fwrite((void*)&D, sizeof(double), 1, File);
}

void	MFWriteDoubleArr(double *Arr, int No, FILE *File)
{ 
	fwrite((void*)Arr, sizeof(double), No, File);
}

FILE*	InitSaveModelFile(char *FName, OPTIONS *Opt, TREES *Trees, RATES *Rates)
{
	FILE *MFile;
	int		NoP;
		
	MFile = fopen(FName, "wb");

	if(MFile == NULL)
	{
		printf("Could not open model file %s for writting.\n", FName);
		exit(1);
	}
	
	NoP = GetParamPerModel(Opt, Trees, Rates);
	
	MFWriteInt(NoP, MFile);

	fflush(MFile);
//	exit(0);
	return MFile;
}

void	SaveModelFile(FILE *MFile, OPTIONS *Opt, TREES *Trees, RATES *Rates)
{
	if(Opt->DataType == CONTINUOUS)
	{
		MFWriteDoubleArr(Rates->Rates, Rates->NoOfRates, MFile);
	
		if(Opt->EstKappa == TRUE)
			MFWriteDouble(Rates->Kappa, MFile);

		if(Opt->EstLambda == TRUE)
			MFWriteDouble(Rates->Lambda, MFile);

		if(Opt->EstDelta == TRUE)
			MFWriteDouble(Rates->Delta, MFile);

		if(Opt->EstOU == TRUE)
			MFWriteDouble(Rates->OU, MFile);

	}
	else
	{
		MFWriteDoubleArr(Rates->FullRates, Rates->NoOfFullRates, MFile);
		MFWriteDoubleArr(Rates->Pis, Trees->NoOfStates, MFile);

		if(Opt->UseCovarion == TRUE)
		{
			MFWriteDouble(Rates->OffToOn, MFile);
			MFWriteDouble(Rates->OnToOff, MFile);
		}

		if(Opt->EstGamma == TRUE)
			MFWriteDouble(Rates->Gamma, MFile);
			
		if(Opt->EstKappa == TRUE)
			MFWriteDouble(Rates->Kappa, MFile);
	}

	fflush(MFile);
}

void		MapDescModelFile(OPTIONS *Opt, RATES *Rates)
{
	MODELFILE *MF;
	TREES *Trees;
	double *Data;
	int		Pos;

	Trees = Opt->Trees;
	MF = Rates->ModelFile;
	Data = MF->ModelP[Rates->ModelNo];

	memcpy((void*)Rates->FullRates, (void*)Data, sizeof(double) * Rates->NoOfFullRates);
	Pos = Rates->NoOfFullRates;
	memcpy((void*)Rates->Pis, (void*)&Data[Pos], sizeof(double) * Trees->NoOfStates);
	Pos += Trees->NoOfStates;

	if(Opt->UseCovarion == TRUE)
	{
		Rates->OffToOn = Data[Pos++];
		Rates->OnToOff = Data[Pos++];
	}

	if(Opt->EstGamma == TRUE)
		Rates->Gamma = Data[Pos++];

	if(Opt->EstKappa == TRUE)
		Rates->Kappa = Data[Pos++];
}

void	MapConModelFile(OPTIONS *Opt, RATES *Rates)
{
	MODELFILE *MF;
	TREES *Trees;
	double *Data;
	int		Pos;

	Trees = Opt->Trees;
	MF = Rates->ModelFile;
	Data = MF->ModelP[Rates->ModelNo];

	memcpy((void*)Rates->Rates, (void*)Data, sizeof(double) * Rates->NoOfRates);
	Pos = Rates->NoOfRates;

	if(Opt->EstKappa == TRUE)
		Rates->Kappa = Data[Pos++];

	if(Opt->EstLambda == TRUE)
		Rates->Lambda = Data[Pos++];

	if(Opt->EstDelta == TRUE)
		Rates->Delta = Data[Pos++];

	if(Opt->EstOU == TRUE)
		Rates->OU = Data[Pos++];

	MapMCMCConRates(Rates, Opt);
}

void		MapModelFile(OPTIONS *Opt, RATES *Rates)
{
	if(Opt->DataType == CONTINUOUS)
		MapConModelFile(Opt, Rates);
	else
		MapDescModelFile(Opt, Rates);
}

MODELFILE*		AllocModelFile(void)
{
	MODELFILE* Ret;

	Ret = (MODELFILE*)malloc(sizeof(MODELFILE));
	if(Ret == NULL)
		MallocErr();

	Ret->FName = NULL;
	Ret->ModelP = NULL;
	Ret->NoModels = 0;
	Ret->NoParam = 0;

	return Ret;
}

size_t			GetFileSize(FILE *F)
{
	size_t Ret, NR;
	char		*B;

	B = (char*)malloc(sizeof(char) * 1024);
	if(B == NULL)
		MallocErr();

	Ret = 0;
	do
	{
		NR = fread((void*)B, sizeof(char), 1024, F);
		Ret = Ret + NR;
	}while(NR == 1024);

	free(B);

	return Ret;
}

void	InitModelFile(MODELFILE* MF)
{
	FILE *In;
	size_t No;
	

	In = fopen(MF->FName, "rb");
	if(In == NULL)
	{
		printf("Could not open model file %s for reading.\n", MF->FName);
		exit(1);
	}

	No = fread(&MF->NoParam, sizeof(int), 1, In);
	if(No != 1)
	{ 
		printf("Error reading from model file %s.\n", MF->FName);
		exit(0);
	}

	No = GetFileSize(In);
	
	if(No % MF->NoParam != 0)
	{
		printf("Corrupt model file %s\n", MF->FName);
		exit(0);
	}

	MF->NoModels = (int)No / (MF->NoParam * sizeof(double));
	
	fclose(In);	
}

void	AlloModelMem(MODELFILE* MF)
{
	int Index;

	MF->ModelP = (double**)malloc(sizeof(double*) * MF->NoModels);
	if(MF->ModelP == NULL)
		MallocErr();

	for(Index=0;Index<MF->NoModels;Index++)
	{
		MF->ModelP[Index] = (double*)malloc(sizeof(double) * MF->NoParam);
		if(MF->ModelP[Index] == NULL)
			MallocErr();
	}
}

void	ReadModelFileData(MODELFILE* MF)
{
	FILE *In;
	int Index;

	In = fopen(MF->FName, "rb");
	if(In == NULL)
		MallocErr();
	fread((void*)&MF->NoParam, sizeof(int), 1, In);

	for(Index=0;Index<MF->NoModels;Index++)
		fread((void*)MF->ModelP[Index], sizeof(double), MF->NoParam, In);
	
	fclose(In);
}

void	PrintModelFile(MODELFILE *MF)
{
	int Index, x;

	printf("No Model\t%d\n", MF->NoModels);
	printf("No Param\t%d\n", MF->NoParam);
	for(Index=0;Index<MF->NoModels;Index++)
	{
		printf("%d\t", Index);
		for(x=0;x<MF->NoParam;x++)
			printf("%f\t", MF->ModelP[Index][x]);
		printf("\n");
	}
}

void			TestModelFile(OPTIONS *Opt, TREES *Trees, RATES *Rates)
{
	MODELFILE *MF;
	int Index;

	MF = Rates->ModelFile;

	printf("\n");
	for(Index=0;Index<MF->NoModels;Index++)
	{
		Rates->ModelNo = Index;
		Rates->Lh = Likelihood(Rates, Trees, Opt);

		printf("%d\t%f\n", Index, Rates->Lh);
	}

	exit(0);
}

MODELFILE*		LoadModelFile(char *FileName, OPTIONS *Opt, TREES *Trees, RATES *Rates)
{
	MODELFILE* Ret;

	Ret = AllocModelFile();
	Ret->FName = StrMake(FileName);
	
	InitModelFile(Ret);

	if(Ret->NoParam != GetParamPerModel(Opt, Trees, Rates))
	{
		printf("Number of parameter in the model file, does not match the number of parameter specified in the model.\n");
		exit(0);
	}

	AlloModelMem(Ret);

	ReadModelFileData(Ret);
//	PrintModelFile(Ret);
//	exit(0);
	return Ret;
}


void			FreeModelFile(MODELFILE *MF)
{
	int Index;

	for(Index=0;Index<MF->NoModels;Index++)
		free(MF->ModelP[Index]);
	free(MF->ModelP);
	free(MF);
}

void		ChangeModelFile(RATES *Rates, RANDSTATES *RS)
{
	MODELFILE *MF;

	MF = Rates->ModelFile;
	Rates->ModelNo = RandUSInt(RS) % MF->NoModels;
}