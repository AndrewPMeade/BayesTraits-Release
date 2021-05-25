#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>

#include "typedef.h"
#include "genlib.h"
#include "part.h"
#include "tag.h"

TAG*	GetTagFromName(OPTIONS *Opt, char *Name)
{
	int Index;

	for(Index=0;Index<Opt->NoTags;Index++)
		if(strcmp(Name, Opt->TagList[Index]->Name) == 0)
			return Opt->TagList[Index];

	return NULL;
}

TAG*	CreateTag(TREES *Trees, char *Name, int NoTaxa, char **TaxaNames)
{
	TAG *Ret;
	int Index;

	Ret = (TAG*)malloc(sizeof(TAG));
	if(Ret == NULL)
		MallocErr();

	Ret->Name = StrMake(Name);
	Ret->NoTaxa = NoTaxa;

	Ret->Taxa = (char**)malloc(sizeof(char*) * NoTaxa);
	if(Ret->Taxa == NULL)
		MallocErr();

	for(Index=0;Index<NoTaxa;Index++)
		Ret->Taxa[Index] = StrMake(TaxaNames[Index]);

	Ret->Part = CreatePart(Trees, Ret->NoTaxa, Ret->Taxa);

	Ret->NodeList = (NODE*)malloc(sizeof(NODE) * Trees->NoOfTrees);
	if(Ret->NodeList == NULL)
		MallocErr();

	for(Index=0;Index<Trees->NoOfTrees;Index++)
		Ret->NodeList[Index] = PartGetMRCA(Trees->Tree[Index], Ret->Part);		
	
	return Ret;
}

void	AddTag(OPTIONS *Opt, int Tokes, char **Passed)
{
	
}