#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#include "typedef.h"
#include "data.h"
#include "genlib.h"
#include "trees.h"
#include "treenode.h"


void	PrintTaxaData(OPTIONS *Opt, TREES* Trees)
{
	int		TIndex;
	int		SIndex;
	TAXA	*Taxa;

	for(TIndex=0;TIndex<Trees->NoOfTaxa;TIndex++)
	{
		Taxa = &Trees->Taxa[TIndex];

		printf("%s\t", Taxa->Name);
		for(SIndex=0;SIndex<Trees->NoOfSites;SIndex++)
		{
			if(Opt->DataType == DISCRETE)
				printf("%d\t", Taxa->DesDataChar[SIndex]);
			else
				printf("%f\t", Taxa->ConData[SIndex]);
		}
		printf("\n");
	}
}

void	PrintData(TREES* Trees)
{
	int		NIndex;
	int		SiteIndex,StateIndex;
	TAXA	*T;
	NODE	N;


	printf("Symbol List: %s\n", Trees->SymbolList);

	for(NIndex=0;NIndex<Trees->NoOfNodes;NIndex++)
	{
		N =&Trees->Tree[0].NodeList[NIndex];
			
		if(N->Tip == TRUE)
		{
			T = N->Taxa;

			printf("%s\t", T->Name);
			for(SiteIndex=0;SiteIndex<Trees->NoOfSites;SiteIndex++)
			{
				for(StateIndex=0;StateIndex<Trees->NoOfStates;StateIndex++)
				{
					printf("%1.0f", N->Partial[SiteIndex][StateIndex]);
				}
			}

			printf("\n");
		}

		
	}
}

TAXA*	FindTaxaFromName(char *Name, TREES* Trees)
{
	int	Index;

	for(Index=0;Index<Trees->NoOfTaxa;Index++)
	{
		if(strcmp(Trees->Taxa[Index].Name, Name) == 0)
			return &Trees->Taxa[Index];
	}

	return NULL;
}

int		ValidDouble(char *Str)
{
	int	Point=FALSE;
	int	Len,Index;

	Len = strlen(Str);

	for(Index=0;Index<Len;Index++)
	{
		if(Str[Index] == '.')
			if(Point==FALSE)
				Point = TRUE;
			else
				return FALSE;
		else
		{
			if(!((Str[Index] >= '0') && (Str[Index] <= '9')))
				return FALSE;
		}
	}

	return TRUE;
}

void	AllocEstDataInfo(TAXA *Taxa, int NoSites)
{
	int	Index;

	if(Taxa->EstDataP != NULL)
		return;

	Taxa->EstData = FALSE;
	Taxa->EstDataP = (char*)malloc(sizeof(char) * NoSites);
	if(Taxa->EstDataP == NULL)
		MallocErr();

	for(Index=0;Index<NoSites;Index++)
		Taxa->EstDataP[Index] = FALSE;
}


int		HadEstData(char *Site)
{
	
	while(*Site != '\0')
	{
		if(*Site == '?')
			return TRUE;
		Site++;
	}

	return FALSE;
}

void	AddDesTaxaData(int Tokes, char** Passed, TREES* Trees)
{
	TAXA	*Taxa=NULL;
	int		Index;

	Taxa = FindTaxaFromName(Passed[0], Trees);
	if(Taxa == NULL)
	{
		printf("Could not find a matching taxa name for data point %s\n", Passed[0]);
		return;
	}

	if(Taxa->DesDataChar != NULL)
	{
		printf("Warrning: Taxa %s all ready had data, and will be over written\n", Taxa->Name);
		free(Taxa->DesDataChar);
	}

	Taxa->DesDataChar = (char**)malloc(sizeof(char*)*Trees->NoOfSites);
	
	AllocEstDataInfo(Taxa, Trees->NoOfSites);

	for(Index=1;Index<Trees->NoOfSites+1;Index++)
	{
		if(HadEstData(Passed[Index]) == FALSE)
		{
			Taxa->DesDataChar[Index-1] = (char*)malloc(sizeof(char)*strlen(Passed[Index])+1);
			if(Taxa->DesDataChar[Index-1] == NULL)
				MallocErr();
			strcpy(Taxa->DesDataChar[Index-1], Passed[Index]);
		}
		else
		{
			Taxa->DesDataChar[Index-1] = (char*)malloc(sizeof(char)*2);
			Taxa->DesDataChar[Index-1][0] = '-';
			Taxa->DesDataChar[Index-1][1] = '\0';

			Taxa->EstDataP[Index-1] = TRUE;
			Taxa->EstData = TRUE;
		}
	}
}

int		EstData(TREES *Trees)
{
	TAXA	*Taxa;
	int	TIndex;
	int	SIndex;

	for(TIndex=0;TIndex<Trees->NoOfTaxa;TIndex++)
	{
		Taxa  = &Trees->Taxa[TIndex];

		for(SIndex=0;SIndex<Trees->NoOfSites;SIndex++)
			if(Taxa->EstDataP[SIndex] == TRUE)
				return TRUE;
	}

	return FALSE;
}

void	AddContinuousTaxaData(int Tokes, char** Passed, TREES* Trees)
{
	TAXA	*Taxa=NULL;
	int		Index;

	Taxa = FindTaxaFromName(Passed[0], Trees);
	if(Taxa == NULL)
	{
		printf("Could not find a matching taxa name for data point %s\n", Passed[0]);
		return;
	}

	if(Taxa->ConData != NULL)
	{
		printf("Warrning: Taxa %s all ready had data, and will be over written\n", Taxa->Name);
		free(Taxa->ConData);
	}

	Taxa->ConData		= (double*)malloc(sizeof(double)*Trees->NoOfSites);
	if(Taxa->ConData == NULL)
		MallocErr();

	AllocEstDataInfo(Taxa, Trees->NoOfSites);

	for(Index=1;Index<Trees->NoOfSites+1;Index++)
	{
		Taxa->ConData[Index-1] = atof(Passed[Index]);
	
		if((Taxa->ConData[Index-1] == 0.0) && (ValidDouble(Passed[Index])== FALSE))
		{
			if((strcmp(Passed[Index], "*") == 0) || (strcmp(Passed[Index], "-") == 0))
				Taxa->Exclude = TRUE;
			else
			{
				if(strcmp(Passed[Index], "?") == 0)
				{
					Taxa->EstDataP[Index-1] = TRUE;
					Taxa->EstData = TRUE;
				}
				else
					Trees->ValidCData = FALSE;
			}
		}
	}
}

void	LoadTaxaData(char* FileName, TREES* Trees)
{
	char*	Buffer=NULL;
	char**	Passed=NULL;
	int		Tokes;
	int		Line=0;
	TEXTFILE	*DataFile;

	Buffer = (char*)malloc(sizeof(char)*BUFFERSIZE);
	Passed = (char**)malloc(sizeof(char**)*BUFFERSIZE);
	if((Buffer == NULL) || (Passed == NULL))
		MallocErr();

	Trees->NoOfSites = -1;
	DataFile = LoadTextFile(FileName, TRUE);

	for(Line=0;Line<DataFile->NoOfLines;Line++)
	{
		strcpy(&Buffer[0], DataFile->Data[Line]);
		Tokes = MakeArgv(&Buffer[0], Passed, BUFFERSIZE);

		if(Tokes > 1)
		{
			if(Trees->NoOfSites == -1)
				Trees->NoOfSites = Tokes - 1;

			if(Tokes - 1  != Trees->NoOfSites)
			{
				printf("Line %d has %d sites but was expecting %d\n", Line, Tokes - 1, Trees->NoOfSites);
				exit(0);
			}

			AddDesTaxaData(Tokes, Passed, Trees);
			AddContinuousTaxaData(Tokes, Passed, Trees);
		}
	}

	free(Buffer);
	free(Passed);

	FreeTextFile(DataFile);
}

int		IsSymbolInList(char Symbol, char* List)
{
	int Index;

	if(Symbol == UNKNOWNSTATE)
		return TRUE	;

	for(Index=strlen(List);Index>=0;Index--)
		if(List[Index] == Symbol)
			return TRUE;

	return FALSE;
}

void	BildSymbolList(TREES *Trees)
{
	char*	Temp=NULL;
	int		TIndex,SIndex,TokeIndex;

//	Temp = (char*)malloc(sizeof(char) * (Trees->NoOfTaxa * Trees->NoOfTaxa) + 1);
	Temp = (char*)malloc(sizeof(char) * BUFFERSIZE);
	if(Temp == NULL)
		MallocErr();
	Temp[0]='\0';

	for(TIndex=0;TIndex<Trees->NoOfTaxa;TIndex++)
	{
		for(SIndex=0;SIndex<Trees->NoOfSites;SIndex++)
		{
			for(TokeIndex=0;TokeIndex<(int)strlen(Trees->Taxa[TIndex].DesDataChar[SIndex]);TokeIndex++)
			{
				if(IsSymbolInList(Trees->Taxa[TIndex].DesDataChar[SIndex][TokeIndex], Temp) == FALSE)
				{
					Temp[strlen(Temp)+1] = '\0';
					Temp[strlen(Temp)] = Trees->Taxa[TIndex].DesDataChar[SIndex][TokeIndex];
				}
			}
		}
	}

	/* To add a hiden state */
/*	Temp[2] = '\0';
	if(Temp[0] == '0')
		Temp[1] = '1';
	else
		Temp[1] = '0';
*/
	Trees->SymbolList = (char*)malloc(sizeof(char)*strlen(Temp)+1);
	
	if(Trees->SymbolList == NULL)
		MallocErr();
	strcpy(Trees->SymbolList, Temp);

	Trees->NoOfStates = strlen(Trees->SymbolList);

	free(Temp);
}

int CompChars(char *char1, char *char2)
{
  if (*char1 <  *char2) 
	  return -1;
  if (*char1 == *char2) 
	  return  0;
  return 1;
}

void	FindSiteSymbols(TREES *Trees, int SiteNo)
{
	int		TIndex;
	TAXA	*Taxa;
	char	*Buffer;
	char	*DP;
	int		DPLen;
	int		DPIndex;

	Buffer = (char*)malloc(sizeof(char) * BUFFERSIZE);
	if(Buffer == NULL)
		MallocErr();
	Buffer[0] = '\0';

	for(TIndex=0;TIndex<Trees->NoOfTaxa;TIndex++)
	{
		Taxa = &Trees->Taxa[TIndex];
		DP = Taxa->DesDataChar[SiteNo];
		DPLen = strlen(DP);
		for(DPIndex=0;DPIndex<DPLen;DPIndex++)
		{
			if(IsSymbolInList(DP[DPIndex], Buffer) == FALSE)
			{
				Buffer[strlen(Buffer)+1] = '\0';
				Buffer[strlen(Buffer)] = DP[DPIndex];
			}
		}
	}

	qsort(Buffer, strlen(Buffer), sizeof(char), (void *)CompChars);
	Trees->SiteSymbols[SiteNo] = StrMake(Buffer);
	Trees->NOSList[SiteNo] = strlen(Trees->SiteSymbols[SiteNo]);

	free(Buffer);	
}

int		ValidDescDataStr(char* Str)
{
	while(*Str != '\0')
	{
		if(!((*Str == '0') || (*Str == '1') || (*Str == '-')))
			return FALSE;
		Str++;
	}

	return TRUE;
}

void	CheckDescData(TREES* Trees)
{
	int		TIndex;
	TAXA*	Taxa;
	int		SIndex;

	for(TIndex=0;TIndex<Trees->NoOfTaxa;TIndex++)
	{
		Taxa = &Trees->Taxa[TIndex];
		for(SIndex=0;SIndex<Trees->NoOfSites;SIndex++)
		{
			if(ValidDescDataStr(Taxa->DesDataChar[SIndex]) == FALSE)
			{
				printf("Taxa %s has invalid discrete data for site %d.\nOnly 0,1 and - are valid discrete data character.\n", Taxa->Name, SIndex+1);
				exit(0);
			}
		}
	}
}

void	CheckDataWithModel(char* FileName, TREES *Trees, MODEL Model)
{
	FILE*	ErrFile;
	char	ErrFileName[1024];

	if(Model == MULTISTATE)
	{
		qsort(Trees->SymbolList, Trees->NoOfStates, sizeof(char), (void *)CompChars);

		if(strlen(Trees->SymbolList) == 1)
		{
			sprintf(&ErrFileName[0], "%s.%s", FileName, LOGFILEEXT);
			ErrFile = OpenWrite(&ErrFileName[0]);

			fprintf(ErrFile, "There has to be more then one state in file %s\n", FileName);
			fclose(ErrFile);
			printf("There has to be more then one state in file %s\n", FileName);
			exit(0);
		}
	}
	else
	{
		if((Model == DESCDEP) || (Model == DESCINDEP))
		{
			CheckDescData(Trees);
		}
	}

	if((Model == DESCDEP) || (Model == DESCINDEP) || (Model == MULTISTATE))
		SetMinBL(Trees);
}

void	LoadData(char* FileName, TREES *Trees)
{
	int		Index;

	LoadTaxaData(FileName, Trees);

	for(Index=0;Index<Trees->NoOfTaxa;Index++)
	{
		if(Trees->Taxa[Index].DesDataChar == NULL)
		{
			printf("Could not load data for taxa %s\n", Trees->Taxa[Index].Name);
			exit(0);
		}
	}

	BildSymbolList(Trees);

	return;
}

void	FreeTaxa(TAXA *Taxa, int NoOfSites)
{
	int		SIndex;

	if(Taxa->DesDataChar != NULL)
	{
		for(SIndex=0;SIndex<NoOfSites;SIndex++)
			free(Taxa->DesDataChar[SIndex]);
		free(Taxa->DesDataChar);
	}

	if(Taxa->ConData != NULL)
		free(Taxa->ConData);
		
	if(Taxa->EstDataP != NULL)
		free(Taxa->EstDataP);

	free(Taxa->Name);
}

void	FreeData(OPTIONS *Opt)
{
	int		Index;
	TAXA	*Taxa;
	int		NOS;
	TREES	*Trees;

	Trees = Opt->Trees;

	NOS = Trees->NoOfSites;
	if(Opt->Model == CONTINUOUSREG)
		NOS++;

	for(Index=0;Index<Opt->Trees->NoOfTaxa;Index++)
	{
		Taxa = &Opt->Trees->Taxa[Index];
		FreeTaxa(Taxa, NOS);
	}
}

//char*	SetDescUnknownStates(char** Sites)
char*	SetDescUnknownStates(char S1, char S2)
{
	char	*Ret=NULL;
	
	if((S1 == UNKNOWNSTATE) && (S2 == UNKNOWNSTATE))
	{
		Ret = (char*)malloc(sizeof(char)*5);
		if(Ret == NULL)
			MallocErr();

		Ret[0] = '0';
		Ret[1] = '1';
		Ret[2] = '2';
		Ret[3] = '3';
		Ret[4] = '\0';
		return Ret;
	}

	Ret = (char*)malloc(sizeof(char)*3);
	if(Ret == NULL)
		MallocErr();

	Ret[2] = '\0';

	if((S1 == '0') && (S2 == UNKNOWNSTATE))
	{
		Ret[0] = '0';
		Ret[1] = '1';
		return Ret;
	}

	if((S1 == '1') && (S2 == UNKNOWNSTATE))
	{
		Ret[0] = '2';
		Ret[1] = '3';
		return Ret;
	}

	if((S1 == UNKNOWNSTATE) && (S2 == '0'))
	{
		Ret[0] = '0';
		Ret[1] = '2';
		return Ret;
	}

	if((S1 == UNKNOWNSTATE) && (S2 == '1'))
	{
		Ret[0] = '1';
		Ret[1] = '3';
		return Ret;
	}

	return Ret;
}

void	SquashDep(TREES	*Trees)
{
	int		TIndex;
	TAXA	*Taxa;
	char	*TempS;
	int		Seen;

	for(TIndex=0;TIndex<Trees->NoOfTaxa;TIndex++)
	{
		Taxa = &Trees->Taxa[TIndex];
		Seen = FALSE;
		
		if(Taxa->EstData == TRUE)
		{
			Taxa->RealData = (char*)malloc(sizeof(char) * 3);
			if(Taxa->RealData == NULL)
				MallocErr();
			Taxa->RealData[0] = Taxa->DesDataChar[0][0];
			Taxa->RealData[1] = Taxa->DesDataChar[1][0];
			Taxa->RealData[2] = '\0';
		}
		
		if((Taxa->DesDataChar[0][0] == '0') && (Taxa->DesDataChar[1][0] == '0') && (Seen == FALSE))
		{
			Taxa->DesDataChar[0][0] = '0';
			Taxa->DesDataChar[0][1] = '\0';
			Seen = TRUE;
		}
		
		if((Taxa->DesDataChar[0][0] == '0') && (Taxa->DesDataChar[1][0]== '1') && (Seen == FALSE))
		{
			Taxa->DesDataChar[0][0] = '1';
			Taxa->DesDataChar[0][1] = '\0';
			Seen = TRUE;
		}

		if((Taxa->DesDataChar[0][0] == '1') && (Taxa->DesDataChar[1][0]== '0') && (Seen == FALSE))
		{
			Taxa->DesDataChar[0][0] = '2';
			Taxa->DesDataChar[0][1] = '\0';
			Seen = TRUE;
		}

		if((Taxa->DesDataChar[0][0] == '1') && (Taxa->DesDataChar[1][0]== '1') && (Seen == FALSE))
		{
			Taxa->DesDataChar[0][0] = '3';
			Taxa->DesDataChar[0][1] = '\0';
			Seen = TRUE;
		}

		if((Taxa->DesDataChar[0][0] == UNKNOWNSTATE) || (Taxa->DesDataChar[1][0] == UNKNOWNSTATE))
		{
			TempS = SetDescUnknownStates(Taxa->DesDataChar[0][0], Taxa->DesDataChar[1][0]);

			free(Taxa->DesDataChar[0]);
			Taxa->DesDataChar[0] = TempS;
			Seen = TRUE;
		}

		free(Taxa->DesDataChar[1]);
		Taxa->DesDataChar[1] = NULL;

/*		printf("%s\t%c\n", Taxa->Name, Taxa->DesDataChar[0][0]); */
	}

	Trees->NoOfStates = 4;
	Trees->NoOfSites = 1;
	free(Trees->SymbolList);
	
	Trees->SymbolList = (char*)malloc(sizeof(char)*strlen("0123")+1);
	if(Trees->SymbolList == NULL)
		MallocErr();
	strcpy(Trees->SymbolList, "0123");

}

void	RemoveConMissingData(TREES* Trees)
{
	int		Index;

	FreePartitions(Trees);

	for(Index=0;Index<Trees->NoOfTaxa;Index++)
	{
		if(Trees->Taxa[Index].Exclude == TRUE)
		{
			RemoveTaxa(NULL, Trees, Trees->Taxa[Index].Name);
			Index=-1;
		}
	}

	SetPartitions(Trees);
}

int		NoOfNodesBelow(NODE N)
{
	int Ret=0;

	while(N->Ans != NULL)
	{
		if(N->Length != 0)
			Ret++;
		N = N->Ans;
	}

	return Ret - 1;
}

double	RootToTipLen(NODE N)
{
	double	Ret;

	Ret = 0;

	do
	{
		Ret += N->Length;
		N = N->Ans;
	} while(N->Ans != NULL);

	return Ret;
}

void	SetTreeAsData(OPTIONS *Opt, TREES *Trees, int TreeNo)
{
	int		NIndex;
	TAXA	*Taxa;
	NODE	N;

	for(NIndex=0;NIndex<Trees->NoOfNodes;NIndex++)
	{
		N = &Trees->Tree[TreeNo].NodeList[NIndex];

		if(N->Tip == TRUE)
		{
			Taxa = N->Taxa;
			Taxa->ConData[0] = (double)NoOfNodesBelow(N);
			if(Opt->NodeBLData == TRUE)
				Taxa->ConData[1] = RootToTipLen(N);
		}		
	}

/*
	printf("\n\n\n");
	printf("\n\nMy Data\n");
	for(NIndex=0;NIndex<Trees->NoOfNodes;NIndex++)
	{
		N = &Trees->Tree[TreeNo].NodeList[NIndex];

		if(N->Tip == TRUE)
		{
			Taxa = N->Taxa;
			printf("%s\t%f\t%f\n", Taxa->Name, Taxa->ConData[0], Taxa->ConData[1]);
		}		
	}
	exit(0);
*/
}

int		TaxaInList(char* Taxa, char** List, int ListSize)
{
	int	Index;

	for(Index=0;Index<ListSize;Index++)
		if(strcmp(Taxa, List[Index]) == 0)
			return TRUE;

	return FALSE;
}

void	LoadVarDataTaxa(OPTIONS *Opt, TEXTFILE* File)
{
	char	**Passed;
	int		Tokes;
	TREES*	Trees;
	int		TIndex;

	Trees = Opt->Trees;

	Passed = (char**)malloc(sizeof(char*) * File->MaxLine);
	if(Passed == NULL)
		MallocErr();

	Tokes = MakeArgv(File->Data[0], Passed, File->MaxLine);

	if(Tokes != Trees->NoOfTaxa)
	{
		printf("Error Number of taxa in %s (%d) does not match number of taxa in the trees (%d)\n", File->FileName, Tokes, Trees->NoOfTaxa);
		exit(0);
	}


	for(TIndex=0;TIndex<Trees->NoOfTaxa;TIndex++)
	{

		if(TaxaInList(Trees->Taxa[TIndex].Name, Passed, Trees->NoOfTaxa) == FALSE)
		{
			printf("Error could not find matching coloum for taxa %s\n", Trees->Taxa[TIndex].Name);
			exit(0);
		}
	}

	Opt->VarData->TaxaNames = (char**)malloc(sizeof(char*) * Trees->NoOfTaxa);
	if(Opt->VarData->TaxaNames == NULL)
		MallocErr();

	for(TIndex=0;TIndex<Trees->NoOfTaxa;TIndex++)
		Opt->VarData->TaxaNames[TIndex] = StrMake(Passed[TIndex]);

	free(Passed);
}

int		FindAndLoadVarData(OPTIONS *Opt, TEXTFILE *File, VARDATA *VarData)
{
	char	**Passed;
	int		Tokes;
	TREES*	Trees;
	int		Ret=0;
	int		Index;
	int		Hits;
	int		TIndex;
	char	*Buffer;
	
	Trees = Opt->Trees;

	Passed = (char**)malloc(sizeof(char*) * File->MaxLine);
	Buffer = (char*)malloc(sizeof(char) * File->MaxLine);
	if((Passed == NULL) || (Buffer == NULL))
		MallocErr();


	Hits = 0;
	for(Index=1;Index<File->NoOfLines;Index++)
	{
		strcpy(Buffer, File->Data[Index]);
		Tokes = MakeArgv(Buffer, Passed, File->MaxLine);
		if(Tokes == Trees->NoOfTaxa)
		{
			if(VarData != NULL)
			{
				
				VarData->Data[Hits] = (double*)malloc(sizeof(double) * Trees->NoOfTaxa);
				if(VarData->Data[Hits] == NULL)
					MallocErr();

				for(TIndex=0;TIndex<Trees->NoOfTaxa;TIndex++)
					VarData->Data[Hits][TIndex] = atof(Passed[TIndex]);
								
			}
			Hits++;
		}
	}

	free(Passed);
	free(Buffer);
	return Hits;
}

int		GetVarDataTaxaPos(char* Taxa, VARDATA* VarData, int No)
{
	int	Index;

	for(Index=0;Index<No;Index++)
	{
		if(strcmp(Taxa, VarData->TaxaNames[Index])  == 0)
			return Index;
	}

	return -1;
}


void	OrderVarData(OPTIONS *Opt)
{
	int			Pos;
	int			TIndex;
	double**	Data;
	TREES*		Trees;
	int			Index;
	char**		TempC;


	Trees = Opt->Trees;
	Data = (double**)malloc(sizeof(double) * Opt->VarData->NoPoints);
	if(Data == NULL)
		MallocErr();

	for(TIndex=0;TIndex<Opt->VarData->NoPoints;TIndex++)
	{
		Data[TIndex] = (double*)malloc(sizeof(double) * Trees->NoOfTaxa);
		if(Data[TIndex] == NULL)
			MallocErr();
	}
	TempC = (char**)malloc(sizeof(char*) * Trees->NoOfTaxa);
	if(TempC == NULL)
		MallocErr();
	
	for(TIndex=0;TIndex<Trees->NoOfTaxa;TIndex++)
	{
		Pos = GetVarDataTaxaPos(Trees->Taxa[TIndex].Name, Opt->VarData, Trees->NoOfTaxa);
		for(Index=0;Index<Opt->VarData->NoPoints;Index++)
			Data[Index][TIndex] = Opt->VarData->Data[Index][Pos];
		TempC[TIndex] = Opt->VarData->TaxaNames[Pos];
	}

	for(Index=0;Index<Opt->VarData->NoPoints;Index++)
		free(Opt->VarData->Data[Index]);
	free(Opt->VarData->Data);

	free(Opt->VarData->TaxaNames);

	Opt->VarData->TaxaNames = TempC;
	Opt->VarData->Data = Data;
}


void	LoadVarData(OPTIONS* Opt)
{
	VARDATA		*Ret;
	TEXTFILE	*File;
/*	int			TIndex,Index; */

	File = LoadTextFile(Opt->VarDataFile, FALSE);

	Ret = (VARDATA*)malloc(sizeof(VARDATA));
	if(Ret == NULL)
		MallocErr();

	Opt->VarData = Ret;

	LoadVarDataTaxa(Opt, File);

	Ret->NoPoints = FindAndLoadVarData(Opt, File, NULL);

	Ret->Data = (double**)malloc(sizeof(double*) * Ret->NoPoints);
	if(Ret->Data == NULL)
		MallocErr();

	FindAndLoadVarData(Opt, File, Opt->VarData);

	FreeTextFile(File);

	OrderVarData(Opt);
	
/*	for(TIndex=0;TIndex<Opt->Trees->NoOfTaxa;TIndex++)
		printf("%s\t", Ret->TaxaNames[TIndex]);

	printf("\n");
	for(Index=0;Index<Ret->NoPoints;Index++)
	{

		for(TIndex=0;TIndex<Opt->Trees->NoOfTaxa;TIndex++)
			printf("%f\t", Ret->Data[Index][TIndex]);
		printf("\n");
	}

	exit(0);
*/

}


void		FreeVarData(OPTIONS *Opt, VARDATA *VarData)
{
	int	Index;

	for(Index=0;Index<Opt->Trees->NoOfTaxa;Index++)
		free(VarData->TaxaNames[Index]);

	free(VarData->TaxaNames);

	for(Index=0;Index<VarData->NoPoints;Index++)
		free(VarData->Data[Index]);

	free(VarData->Data);

	free(VarData);	
}

void	SetVarData(TREES* Trees, VARDATA *VarData, int Site)
{
	int	TIndex;

	for(TIndex=0;TIndex<Trees->NoOfTaxa;TIndex++)
		Trees->Taxa[TIndex].ConData[0] = VarData->Data[Site][TIndex];
}



int		FreeTaxaNo(int No, TREES* Trees)
{
	int	Index;

	for(Index=0;Index<Trees->NoOfTaxa;Index++)
	{
		if(Trees->Taxa[Index].No == No)
			return FALSE;
	}

	return TRUE;
}

int		GetFreeTaxaNo(TREES* Trees)
{
	int	No;

	No = Trees->NoOfTaxa;
	while(FreeTaxaNo(No, Trees)==FALSE)
		No--;

	return No;
}
/*
	Taxa
	int		No;
	char*	Name;
	char**	DesDataChar;
	double*	ConData;
	int		Exclude;
	double	Dependant;

	int		EstData;
	char	*EstDataP;
	int		EstDepData;
*/

void	SetNewConTaxa(TAXA* Taxa, char* Name, TREES* Trees)
{
	int	Index;

	Taxa->Name			= StrMake(Name);
	Taxa->No			= GetFreeTaxaNo(Trees);
	Taxa->DesDataChar	= NULL;

	Taxa->ConData		= (double*)malloc(sizeof(double) * Trees->NoOfSites);
	if(Taxa->ConData == NULL)
		MallocErr();

	Taxa->Exclude		= FALSE;
	Taxa->EstDepData	= FALSE;

	Taxa->EstData		= TRUE;
	Taxa->EstDataP		= (char*)malloc(sizeof(char) * Trees->NoOfSites);
	if(Taxa->EstDataP == NULL)
		MallocErr();
	
	for(Index=0;Index<Trees->NoOfSites;Index++)
	{
		Taxa->EstDataP[Index] = TRUE;	
		Taxa->ConData[Index] = 0;
	}
}

void	AddNewConTaxa(TREES* Trees, char* Name)
{
	TAXA	*NewTaxa;
	TAXA	*NT;

	NewTaxa = (TAXA*)malloc(sizeof(TAXA) * (Trees->NoOfTaxa + 1));
	if(NewTaxa == NULL)
		MallocErr();
	memset(NewTaxa, '\0', sizeof(TAXA) * (Trees->NoOfTaxa + 1));

	NewTaxa = memcpy(NewTaxa, Trees->Taxa, sizeof(TAXA) * Trees->NoOfTaxa);
 
	NT = &NewTaxa[Trees->NoOfTaxa];

	SetNewConTaxa(NT, Name, Trees);
	
	free(Trees->Taxa);
	Trees->Taxa = NewTaxa;
	Trees->NoOfTaxa++;
}

void	AddRecNodes(OPTIONS *Opt, TREES *Trees)
{
	int		Index;
	RECNODE	RNode;

	/* Add Rec Taxa */
	for(Index=0;Index<Opt->NoOfRecNodes;Index++)
	{
		RNode = Opt->RecNodeList[Index];
		AddNewConTaxa(Trees, RNode->Name);
		AddNewRecNode(Trees, RNode);
	}
}