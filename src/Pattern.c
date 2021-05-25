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



#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "Pattern.h"
#include "genlib.h"
#include "typedef.h"
#include "Tag.h"

PATTERN*	AllocPattern(void)
{
	PATTERN* Ret;

	Ret = (PATTERN*)SMalloc(sizeof(PATTERN));

	Ret->NoTags = 0;
	Ret->Name = NULL;
	Ret->TagList = NULL;
	
	return Ret;
}

void		FreePattern(PATTERN *Pattern)
{
	free(Pattern->Name);
	free(Pattern->TagList);
	free(Pattern);
}

PATTERN*	GetPatternFromName(char *Name, int NoP, PATTERN **PList)
{
	int Index;

	for(Index=0;Index<NoP;Index++)
	{
		if(StrICmp(Name, PList[Index]->Name) == 0)
			return PList[Index];
	}

	return NULL;
}

void	PrintPatterns(FILE *Str, int NoPatterns, PATTERN **PList)
{
	int Index, TIndex;
	PATTERN *P;

	if(NoPatterns == 0)
		return;

	fprintf(Str, "Patterns %d:\n", NoPatterns);

	for(Index=0;Index<NoPatterns;Index++)
	{
		P = PList[Index];

		fprintf(Str, "\t%s\t", P->Name);

		for(TIndex=0;TIndex<P->NoTags;TIndex++)
			fprintf(Str, "%s\t", P->TagList[TIndex]->Name);

		fprintf(Str, "\n");
	}
}

void	AddPattern(OPTIONS *Opt, char *Name, int NoTags, char **TagNameList)
{
	PATTERN *NPat;
	TAG		*Tag;
	int Index;

	NPat = GetPatternFromName(Name, Opt->NoPatterns, Opt->PatternList);

	if(NPat != NULL)
	{
		printf("Pattern name %s allready in use.\n", Name);
		exit(1);
	}

	NPat = AllocPattern();

	NPat->Name = StrMake(Name);
	NPat->NoTags = NoTags;
	NPat->TagList = (TAG**)SMalloc(sizeof(TAG*) * NoTags);

	for(Index=0;Index<NoTags;Index++)
	{
		Tag = GetTagFromName(Opt, TagNameList[Index]);
		NPat->TagList[Index] = Tag;
	}

	Opt->PatternList = (PATTERN**)AddToList(&Opt->NoPatterns, (void**)Opt->PatternList, (void*)NPat);
}

void	DelOldRateNames(OPTIONS *Opt)
{
	int Index;

	for(Index=0;Index<Opt->NoOfRates;Index++)
		free(Opt->RateName[Index]);

	free(Opt->RateName);
}

char*	MakeLocalRateName(char *DName, char *PName)
{
	char *Ret;
	size_t Size;

	Size = strlen(DName) + strlen(PName) + 2;

	Ret = (char*)SMalloc(sizeof(char) * Size);

	sprintf(Ret, "%s-%s", DName, PName);

	return Ret;
}

void	SetPatternRateNames(OPTIONS *Opt)
{
	int Pos, Index, PIndex;
	char *PName;

	DelOldRateNames(Opt);
	
	Opt->NoOfRates = Opt->DefNoRates * (1 + Opt->NoPatterns);
	Opt->RateName = (char**)SMalloc(sizeof(char*) * Opt->NoOfRates);

	Pos = 0;
	for(Index=0;Index<Opt->DefNoRates;Index++,Pos++)
		Opt->RateName[Pos] = StrMake(Opt->DefRateNames[Index]);

	for(PIndex=0;PIndex<Opt->NoPatterns;PIndex++)
	{
		PName = Opt->PatternList[PIndex]->Name;
		for(Index=0;Index<Opt->DefNoRates;Index++,Pos++)
			Opt->RateName[Pos] = MakeLocalRateName(Opt->DefRateNames[Index], PName);
	}
}

void	RecSetNodePNo(NODE N, int PNo)
{
	int Index;

	N->PatternNo = PNo;
	if(N->Tip == TRUE)
		return;

	for(Index=0;Index<N->NoNodes;Index++)
		RecSetNodePNo(N->NodeList[Index], PNo);
}

void	SetPatternNoTree(OPTIONS *Opt, TREES *Trees, int TreeNo)
{
	int PIndex, TIndex;
	PATTERN *Pattern;
	TREE *Tree;
	NODE N;

	Tree = Trees->Tree[TreeNo];

	RecSetNodePNo(Tree->Root, 0);

	for(PIndex=0;PIndex<Opt->NoPatterns;PIndex++)
	{
		Pattern = Opt->PatternList[PIndex];
		for(TIndex=0;TIndex<Pattern->NoTags;TIndex++)
		{
			N = Pattern->TagList[TIndex]->NodeList[TreeNo];
			RecSetNodePNo(N, PIndex+1);
		}
	}
}

void	SetPatternNo(OPTIONS *Opt, TREES *Trees)
{
	int TIndex;

	for(TIndex=0;TIndex<Trees->NoTrees;TIndex++)
		SetPatternNoTree(Opt, Trees, TIndex);
}