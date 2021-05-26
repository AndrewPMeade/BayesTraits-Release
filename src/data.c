/*
*  BayesTriats 3.0
*
*  copyright 2017
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

#include "TypeDef.h"
#include "Data.h"
#include "GenLib.h"
#include "Trees.h"
#include "Part.h"
#include "Options.h"
#include "Geo.h"
#include "FatTail.h"
#include "Rates.h"

void	PrintTaxaData(OPTIONS *Opt, TREES* Trees)
{
	int		TIndex;
	int		SIndex;
	TAXA	*Taxa;

	for(TIndex=0;TIndex<Trees->NoTaxa;TIndex++)
	{
		Taxa = Trees->Taxa[TIndex];

		printf("%s\t", Taxa->Name);
		for(SIndex=0;SIndex<Trees->NoSites;SIndex++)
		{
			if(Opt->DataType == DISCRETE)
				printf("%s\t", Taxa->DesDataChar[SIndex]);
			else
				printf("%f\t", Taxa->ConData[SIndex]);
		}
		printf("\n");
	}
}

void	PrintDataDesc(TREES* Trees)
{
	int		NIndex;
	int		SiteIndex,StateIndex;
	TAXA	*T;
	TREE	*Tree;
	NODE	N;

	Tree = Trees->Tree[0];

	printf("Symbol List: %s\n", Trees->SymbolList);

	for(NIndex=0;NIndex<Tree->NoNodes;NIndex++)
	{
		N = Tree->NodeList[NIndex];

		if(N->Tip == TRUE)
		{
			T = N->Taxa;

			printf("%s\t", T->Name);
			for(SiteIndex=0;SiteIndex<Trees->NoSites;SiteIndex++)
			{

				for(StateIndex=0;StateIndex<Trees->NoStates;StateIndex++)
				{
					printf("%1.0f", N->Partial[SiteIndex][StateIndex]);
				}
			}

			printf("\n");
		}
	}
}

void	PrintDataCon(TREES* Trees, OPTIONS *Opt)
{
	int		TIndex;
	int		SIndex;
	TAXA	*Taxa;

	for(TIndex=0;TIndex<Trees->NoTaxa;TIndex++)
	{
		Taxa = Trees->Taxa[TIndex];

		printf("%s\t", Taxa->Name);
		if(strcmp(Taxa->Name, "donii") == 0)
			printf("Yetp\n");
		for(SIndex=0;SIndex<Trees->NoSites;SIndex++)
			printf("%f\t", Taxa->ConData[SIndex]);
		if(Opt->Model == M_CONTINUOUS_REG)
			printf("|\t%f", Taxa->Dependant);
		printf("\n");
	}
}

void	PrintData(TREES* Trees, OPTIONS *Opt)
{
	if(Opt->DataType == CONTINUOUS)
	{
		PrintDataCon(Trees, Opt);
		return;
	}

	PrintDataDesc(Trees);
}


TAXA*	FindTaxaFromName(char *Name, TREES* Trees)
{
	int	Index;

	for(Index=0;Index<Trees->NoTaxa;Index++)
	{
		if(strcmp(Trees->Taxa[Index]->Name, Name) == 0)
			return Trees->Taxa[Index];
	}

	return NULL;
}

int		ValidDouble(char *Str)
{
	int	Point=FALSE;
	int	Len,Index;

	Len = (int)strlen(Str);

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
	Taxa->EstDataP = (char*)SMalloc(sizeof(char) * NoSites);

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


	if(Taxa->DesDataChar != NULL)
	{
		printf("Warrning: Taxa %s all ready had data, and will be over written\n", Taxa->Name);
		free(Taxa->DesDataChar);
	}

	Taxa->DesDataChar = (char**)SMalloc(sizeof(char*)*Trees->NoSites);



	for(Index=1;Index<Trees->NoSites+1;Index++)
	{
		if(HadEstData(Passed[Index]) == FALSE)
			Taxa->DesDataChar[Index-1] = StrMake(Passed[Index]);
		else
		{
			Taxa->DesDataChar[Index-1] = (char*)SMalloc(sizeof(char)*2);
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

	for(TIndex=0;TIndex<Trees->NoTaxa;TIndex++)
	{
		Taxa  = Trees->Taxa[TIndex];

		if(Taxa->EstData == TRUE)
			return TRUE;
	}
	return FALSE;
}


void	AddContinuousTaxaData(int Tokes, char** Passed, TREES* Trees)
{
	TAXA	*Taxa;
	int		Index;

	Taxa = FindTaxaFromName(Passed[0], Trees);

	if(Taxa->ConData != NULL)
	{
		printf("Warrning: Taxa %s all ready had data, and will be over written\n", Taxa->Name);
		free(Taxa->ConData);
		fflush(stdout);
	}

	Taxa->ConData		= (double*)SMalloc(sizeof(double)*Trees->NoSites);

	Taxa->Exclude = FALSE;

	for(Index=1;Index<Trees->NoSites+1;Index++)
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
	char*		Buffer;
	char**		Passed;
	int			Tokes;
	int			Line;
	TEXTFILE	*DataFile;
	size_t		MSize;
	TAXA		*Taxa;


	Trees->NoSites = -1;
	DataFile = LoadTextFile(FileName, TRUE);

	MSize = DataFile->MaxLine + 1;

	Buffer = (char*)SMalloc(sizeof(char) * MSize);
	Passed = (char**)SMalloc(sizeof(char**) * MSize);

	for(Line=0;Line<DataFile->NoOfLines;Line++)
	{
		strcpy(&Buffer[0], DataFile->Data[Line]);
		Tokes = MakeArgv(&Buffer[0], Passed, (int)MSize);

		if(Tokes > 1)
		{
			if(Trees->NoSites == -1)
				Trees->NoSites = Tokes - 1;

			if(Tokes - 1  != Trees->NoSites)
			{
				
				printf("Line %d has %d traits but was expecting %d\n", Line+1, Tokes - 1, Trees->NoSites);
				printf("Error occurred on line.\n%s", DataFile->Data[Line]);
				exit(0);
			}


			Taxa = FindTaxaFromName(Passed[0], Trees);
			if(Taxa != NULL)
			{
				AllocEstDataInfo(Taxa, Trees->NoSites);
				AddDesTaxaData(Tokes, Passed, Trees);
				AddContinuousTaxaData(Tokes, Passed, Trees);
			}
			else
				printf("Could not find a matching taxa name for data point %s\n", Passed[0]);
		}
	}

	free(Buffer);
	free(Passed);

	Trees->NoUserSites = Trees->NoSites;

	FreeTextFile(DataFile);
}

int		IsSymbolInList(char Symbol, char* List)
{
	int Index;

	if(Symbol == UNKNOWNSTATE)
		return TRUE	;

	for(Index=(int)strlen(List);Index>=0;Index--)
		if(List[Index] == Symbol)
			return TRUE;

	return FALSE;
}

void	BildSymbolList(TREES *Trees)
{
	char*	Temp;
	int		TIndex,SIndex,TokeIndex;

	Temp = (char*)SMalloc(sizeof(char) * BUFFERSIZE);
	for(TIndex=0;TIndex<BUFFERSIZE;TIndex++)
		Temp[TIndex] = '\0';

	for(TIndex=0;TIndex<Trees->NoTaxa;TIndex++)
	{
		for(SIndex=0;SIndex<Trees->NoSites;SIndex++)
		{
			for(TokeIndex=0;TokeIndex<(int)strlen(Trees->Taxa[TIndex]->DesDataChar[SIndex]);TokeIndex++)
			{
				if(IsSymbolInList(Trees->Taxa[TIndex]->DesDataChar[SIndex][TokeIndex], Temp) == FALSE)
					Temp[strlen(Temp)] = Trees->Taxa[TIndex]->DesDataChar[SIndex][TokeIndex];
			}
		}
	}

	/* To add a hiden state */
/*	if(strlen(Temp) == 1)
	{
		if(IsSymbolInList('0', Temp) == FALSE)
			Temp[strlen(Temp)] = '0';
		else
			Temp[strlen(Temp)] = '1';
	}
	*/
	Trees->SymbolList = StrMake(Temp);

	Trees->NoStates = (int)strlen(Trees->SymbolList);

	free(Temp);
}

//int CompChars(char *char1, char *char2)
int CompChars(const void* c1, const void *c2)
{
	char *char1, *char2;

	char1 = (char*)c1;
	char2 = (char*)c2;

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

	for(TIndex=0;TIndex<Trees->NoTaxa;TIndex++)
	{
		Taxa = Trees->Taxa[TIndex];
		DP = Taxa->DesDataChar[SiteNo];
		DPLen = (int)strlen(DP);
		for(DPIndex=0;DPIndex<DPLen;DPIndex++)
		{
			if(IsSymbolInList(DP[DPIndex], Buffer) == FALSE)
			{
				Buffer[strlen(Buffer)+1] = '\0';
				Buffer[strlen(Buffer)] = DP[DPIndex];
			}
		}
	}

//	qsort(Buffer, strlen(Buffer), sizeof(char), (void *)CompChars);
	qsort(Buffer, strlen(Buffer), sizeof(char), CompChars);
	Trees->SiteSymbols[SiteNo] = StrMake(Buffer);
	Trees->NOSList[SiteNo] = (int)strlen(Trees->SiteSymbols[SiteNo]);

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

	for(TIndex=0;TIndex<Trees->NoTaxa;TIndex++)
	{
		Taxa = Trees->Taxa[TIndex];
		for(SIndex=0;SIndex<Trees->NoSites;SIndex++)
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
	if(ValidModelChoice(Trees, Model) == FALSE)
		exit(1);

	if(Model == M_MULTISTATE)
	{
		qsort(Trees->SymbolList, Trees->NoStates, sizeof(char), CompChars);

		if(strlen(Trees->SymbolList) == 1)
		{
			printf("There has to be more then one state in file %s\n", FileName);
			exit(0);
		}
	}
	else
	{
		if(Model == M_DISC_DEP || Model == M_DISC_INDEP)
			CheckDescData(Trees);
	}

	if(Model == M_GEO)
		ValidGeoData(Trees);

	if(GetModelType(Model) == MT_DISCRETE)
		SetMinBL(Trees);

	if(GetModelType(Model) == MT_FATTAIL)
		CheckFatTailBL(Trees);
}

void	PreProcessDataWithModel(TREES *Trees, MODEL Model)
{
	if(ModelDep(Model) == TRUE)
		SquashDep(Trees);

	if(GetModelType(Model) == MT_CONTINUOUS || GetModelType(Model) == MT_CONTRAST || GetModelType(Model) == MT_FATTAIL)
		RemoveConMissingData(Trees);

	if(Model == M_GEO)
		PreProcessGeoData(Trees);
}

void	LoadData(char* FileName, TREES *Trees)
{
	int		Index;

	LoadTaxaData(FileName, Trees);

	for(Index=0;Index<Trees->NoTaxa;Index++)
	{
		if(Trees->Taxa[Index]->DesDataChar == NULL)
		{
			printf("Could not load data for taxa %s\n", Trees->Taxa[Index]->Name);
			exit(0);
		}
	}
/*
	for(Index=0;Index<Trees->NoTaxa;Index++)
	{
		printf("SData:\t%s\t%s\t%s\n", Trees->Taxa[Index]->Name, Trees->Taxa[Index]->DesDataChar[0], Trees->Taxa[Index]->DesDataChar[1]);
	}
	exit(0);
*/
	BildSymbolList(Trees);

	return;
}

void	FreeTaxa(TAXA *Taxa, int NoSites)
{
	int		SIndex;

	if(Taxa->DesDataChar != NULL)
	{
		for(SIndex=0;SIndex<NoSites;SIndex++)
			free(Taxa->DesDataChar[SIndex]);
		free(Taxa->DesDataChar);
	}

	if(Taxa->ConData != NULL)
		free(Taxa->ConData);

	if(Taxa->EstDataP != NULL)
		free(Taxa->EstDataP);

	if(Taxa->RealData != NULL)
		free(Taxa->RealData);

	free(Taxa->Name);

	free(Taxa);
}

void	FreeData(OPTIONS *Opt)
{
	int		Index;
	TAXA	*Taxa;
	int		NOS;
	TREES	*Trees;

	Trees = Opt->Trees;

	NOS = Trees->NoSites;
	if(Opt->Model == M_CONTINUOUS_REG)
		NOS++;

	for(Index=0;Index<Opt->Trees->NoTaxa;Index++)
	{
		Taxa = Opt->Trees->Taxa[Index];
		FreeTaxa(Taxa, NOS);
	}
}

/* char*	SetDescUnknownStates(char** Sites) */
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

int		Dep01Site(char *Site)
{
	int S0, S1, SM, Index, Len;

	S0 = S1 = SM = FALSE;
	Len = (int)strlen(Site);
	for(Index=0;Index<Len;Index++)
	{
		if(Site[Index] == '0')
			S0 = TRUE;

		if(Site[Index] == '1')
			S1 = TRUE;

		if(Site[Index] == UNKNOWNSTATE)
			SM = TRUE;
	}

	if((S0 == TRUE) && (S1 == TRUE))
		return TRUE;

	if(SM == TRUE)
		return TRUE;

	return FALSE;
}

void	SetDep01Unknown(TAXA *Taxa)
{
	if(Taxa->EstData == TRUE)
		return;

	if(Dep01Site(Taxa->DesDataChar[0]) == TRUE)
	{
		Taxa->DesDataChar[0][0] = UNKNOWNSTATE;
		Taxa->DesDataChar[0][1] = '\0';
	}

	if(Dep01Site(Taxa->DesDataChar[1]) == TRUE)
	{
		Taxa->DesDataChar[1][0] = UNKNOWNSTATE;
		Taxa->DesDataChar[1][1] = '\0';
	}
}

void	SquashDep(TREES	*Trees)
{
	int		TIndex;
	TAXA	*Taxa;
	char	*TempS;
	int		Seen;

	for(TIndex=0;TIndex<Trees->NoTaxa;TIndex++)
	{
		Taxa = Trees->Taxa[TIndex];
		Seen = FALSE;

		if(Taxa->EstData == TRUE)
		{
			Taxa->RealData = (char*)SMalloc(sizeof(char) * 3);

			Taxa->RealData[0] = Taxa->DesDataChar[0][0];
			Taxa->RealData[1] = Taxa->DesDataChar[1][0];
			Taxa->RealData[2] = '\0';
		}

		SetDep01Unknown(Taxa);

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

		if(		(Taxa->DesDataChar[0][0] == UNKNOWNSTATE) 
			||	(Taxa->DesDataChar[1][0] == UNKNOWNSTATE))
		{
			TempS = SetDescUnknownStates(Taxa->DesDataChar[0][0], Taxa->DesDataChar[1][0]);

			free(Taxa->DesDataChar[0]);
			Taxa->DesDataChar[0] = TempS;
			Seen = TRUE;
		}

		free(Taxa->DesDataChar[1]);
		Taxa->DesDataChar[1] = NULL;
	}
	Trees->NoStates = 4;
	Trees->NoSites = 1;
	free(Trees->SymbolList);


	Trees->SymbolList = StrMake("0123");
}

void	RemoveConMissingData(TREES* Trees)
{
	int		Index;

	FreeTreeParts(Trees);

	for(Index=0;Index<Trees->NoTaxa;Index++)
	{
		if(Trees->Taxa[Index]->Exclude == TRUE)
		{
			RemoveTaxa(Trees, Trees->Taxa[Index]->Name);
			Index=-1;
		}

		if(Trees->NoTaxa == 1)
		{
			printf("Deleting taxa, 1 remaining taxa.\n");
			exit(1);
		}
	}

	SetParts(Trees);
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
	TREE	*Tree;
	TAXA	*Taxa;
	NODE	N;

	Tree = Trees->Tree[TreeNo];

	for(NIndex=0;NIndex<Tree->NoNodes;NIndex++)
	{
		N = Tree->NodeList[NIndex];

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



int		FreeTaxaNo(int No, TREES* Trees)
{
	int	Index;

	for(Index=0;Index<Trees->NoTaxa;Index++)
	{
		if(Trees->Taxa[Index]->No == No)
			return FALSE;
	}

	return TRUE;
}

int		GetFreeTaxaNo(TREES* Trees)
{
	int	No;

	No = Trees->NoTaxa;
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

void	SetNewConTaxaData(TAXA *Taxa, RECNODE *RNode, TREES* Trees)
{
	int	Index;

	Taxa->ConData		= (double*)SMalloc(sizeof(double) * Trees->NoSites);

/* Will have to chagne */

	Taxa->EstDataP		= (char*)SMalloc(sizeof(char) * Trees->NoSites);

	Taxa->EstData		= FALSE;

	for(Index=0;Index<Trees->NoSites;Index++)
	{
		if(strcmp(RNode->ConData[Index], "?") == 0)
		{
			Taxa->EstData = TRUE;
			Taxa->EstDataP[Index] = TRUE;
			Taxa->ConData[Index] = 0;
		}
		else
		{
			Taxa->EstDataP[Index] = FALSE;
			Taxa->ConData[Index] = atof(RNode->ConData[Index]);
		}
	}

}

TAXA*	SetNewConTaxa(RECNODE *RNode, TREES* Trees)
{
	TAXA *Ret;
	Ret = (TAXA*)SMalloc(sizeof(TAXA));

	Ret->Name			= StrMake(RNode->Name);
	Ret->No				= GetFreeTaxaNo(Trees);
	Ret->DesDataChar	= NULL;

	Ret->Exclude		= FALSE;
	Ret->EstDepData		= FALSE;
	Ret->RealData		= NULL;


	SetNewConTaxaData(Ret, RNode, Trees);

	return Ret;
}

void	AddNewConTaxa(TREES* Trees, RECNODE	*RNode)
{
	TAXA	**NewTaxa;

	NewTaxa = (TAXA**)SMalloc(sizeof(TAXA*) * (Trees->NoTaxa + 1));

	memcpy(NewTaxa, Trees->Taxa, sizeof(TAXA*) * Trees->NoTaxa);

	NewTaxa[Trees->NoTaxa] = SetNewConTaxa(RNode, Trees);

	free(Trees->Taxa);
	Trees->Taxa = NewTaxa;
	Trees->NoTaxa++;
}

void	AddRecNodes(OPTIONS *Opt, TREES *Trees)
{
	int		Index;
	RECNODE	*RNode;

	for(Index=0;Index<Opt->NoOfRecNodes;Index++)
	{
		RNode = Opt->RecNodeList[Index];
		AddNewConTaxa(Trees, RNode);
		AddNewRecNode(Trees, RNode);
	}

	FreeTreeParts(Trees);
	SetParts(Trees);

	SetTreesDistToRoot(Trees);
}


void		SetDataRegTC(OPTIONS *Opt)
{
	int TIndex, SIndex;
	TREES *Trees;
	TAXA	*Taxa;

	Trees = Opt->Trees;


	for(TIndex=0;TIndex<Trees->NoTaxa;TIndex++)
	{
		Taxa = Trees->Taxa[TIndex];
		for(SIndex=1;SIndex<Trees->NoSites;SIndex++)
		{
			Taxa->ConData[SIndex] = 1;
		}
	}
}
