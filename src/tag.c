#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>

#include "typedef.h"
#include "genlib.h"
#include "part.h"
#include "Tag.h"

TAG*	GetTagFromName(OPTIONS *Opt, char *Name)
{
	int Index;

	for(Index=0;Index<Opt->NoTags;Index++)
		if(strcmp(Name, Opt->TagList[Index]->Name) == 0)
			return Opt->TagList[Index];

	return NULL;
}


void	FreeTag(TAG *Tag)
{
	int Index;

	free(Tag->Name);

	for(Index=0;Index<Tag->NoTaxa;Index++)
		free(Tag->Taxa[Index]);
	free(Tag->Taxa);

	free(Tag->NodeList);
	FreePart(Tag->Part);
	free(Tag);
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

	Ret->NodeList = (NODE*)malloc(sizeof(NODE) * Trees->NoTrees);
	if(Ret->NodeList == NULL)
		MallocErr();

	for(Index=0;Index<Trees->NoTrees;Index++)
		Ret->NodeList[Index] = PartGetMRCA(Trees->Tree[Index], Ret->Part);		
	
	return Ret;
}

void	AddTag(OPTIONS *Opt, int Tokes, char **Passed)
{
	TAG*	Tag;
	TREES*	Trees;
	char*	Name;
	
	Trees = Opt->Trees;

	if(Tokes < 3)
	{
		printf("A tag must have a name and list of one or more taxa.\n");
		exit(0);
	}

	Name = Passed[1];
	
	Tag = GetTagFromName(Opt, Name);

	if(Tag != NULL)
	{
		printf("Tag %s already exists.\n", Name);
		exit(0);
	}
		
	Tag = CreateTag(Trees, Name, Tokes-2, &Passed[2]);
	
	Opt->TagList = (TAG**)AddToList(&Opt->NoTags, (void**)Opt->TagList, Tag);
}

void	PrintTag(FILE *Out, TAG *Tag)
{
	int Index;

	fprintf(Out, "\t%s\t%d\t", Tag->Name, Tag->NoTaxa);

	for(Index=0;Index<Tag->NoTaxa;Index++)
		fprintf(Out, "%s ", Tag->Taxa[Index]);
	fprintf(Out, "\n");
}

void	PrintTags(FILE *Out, OPTIONS *Opt)
{
	int Index;

	if(Opt->NoTags < 1)
		return;

	fprintf(Out, "Tags:\t%d\n", Opt->NoTags);
	
	for(Index=0;Index<Opt->NoTags;Index++)
		PrintTag(Out, Opt->TagList[Index]);
}