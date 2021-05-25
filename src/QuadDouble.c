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
#include <stdio.h>

#include "TypeDef.h"
#include "QuadDouble.h"
#include "GenLib.h"


#ifndef QUAD_DOUBLE
void	InitQuadDoubleLh(OPTIONS *Opt, TREES *Trees) {}
void	FreeQuadLh(OPTIONS *Opt, TREES *Trees) {}

void	NodeLhQuadDouble(NODE N, OPTIONS *Opt, TREES *Trees, int SiteNo) {}


double	CombineQuadDoubleLh(RATES* Rates, TREES *Trees, OPTIONS *Opt, int SiteNo, int NOS) {return -1;}

#else

//#ifdef QUAD_DOUBLE

void	AllocNodeQuadMem(NODE N, OPTIONS *Opt, TREES *Trees)
{
	int SIndex, Index;

	N->BigPartial = (QDOUBLE**)malloc(sizeof(QDOUBLE*) * Trees->NoOfSites);
	if(N->BigPartial == NULL)
		MallocErr();

	for(SIndex=0;SIndex<Trees->NoOfSites;SIndex++)
	{
		N->BigPartial[SIndex] = (QDOUBLE*)malloc(sizeof(QDOUBLE) * Trees->NoStates);
		if(N->BigPartial[SIndex] == NULL)
			MallocErr();

		for(Index=0;Index<Trees->NoStates;Index++)
			N->BigPartial[SIndex][Index] = N->Partial[SIndex][Index];
	}
}

void	InitQuadDoubleLh(OPTIONS *Opt, TREES *Trees)
{
	int TIndex, NIndex;
	TREE *Tree;
	NODE N;

	for(TIndex=0;TIndex<Trees->NoTrees;TIndex++)
	{
		Tree = Trees->Tree[TIndex];
		for(NIndex=0;NIndex<Tree->NoNodes;NIndex++)
		{
			N = Tree->NodeList[NIndex];
			AllocNodeQuadMem(N, Opt, Trees);
		}
	}
}

void	FreeQuadLh(OPTIONS *Opt, TREES *Trees)
{
	int TIndex, NIndex, Index;
	NODE N;
	TREE *T;

	for(TIndex=0;TIndex<Trees->NoTrees;TIndex++)
	{
		T = Trees->Tree[TIndex];
		for(NIndex=0;NIndex<T->NoNodes;NIndex++)
		{
			N = T->NodeList[NIndex];
			for(Index=0;Index<Trees->NoOfSites;Index++)
				free(N->BigPartial[Index]);
			free(N->BigPartial);
		}
	}
}

void	NodeLhQuadDouble(NODE N, OPTIONS *Opt, TREES *Trees, int SiteNo)
{
	int		Inner, Outter, NIndex;
	QDOUBLE	Lh;
	QDOUBLE **TBigLh;
	double **Mat;

	for(Outter=0;Outter<Trees->NoStates;Outter++)
	{
		N->BigPartial[SiteNo][Outter] = 1.0;

		for(NIndex=0;NIndex<N->NoNodes;NIndex++)
		{
			Mat = Trees->PList[N->NodeList[NIndex]->ID]->me;
			TBigLh = N->NodeList[NIndex]->BigPartial;

			Lh = 0;
			for(Inner=0;Inner<Trees->NoStates;Inner++)
				Lh += TBigLh[SiteNo][Inner] * (QDOUBLE)Mat[Outter][Inner];

			N->BigPartial[SiteNo][Outter] *= Lh;
		}
	}

	if(N->FossilMask != NULL)
		FossilLh(N, Opt, Trees, SiteNo);
}

double	CombineQuadDoubleLh(RATES* Rates, TREES *Trees, OPTIONS *Opt, int SiteNo, int NOS)
{
	int Index;
	double Ret;
	QDOUBLE Sum;
	TREE *Tree;

	Tree = Trees->Tree[Rates->TreeNo];

	Sum = 0;

	for(Index=0;Index<NOS;Index++)
		Sum += Tree->Root->BigPartial[SiteNo][Index] * (QDOUBLE)Rates->Pis[Index];

	for(Index=0;Index<NOS;Index++)
	{
		Tree->Root->Partial[SiteNo][Index] = Tree->Root->BigPartial[SiteNo][Index] / Sum;
	}

	Ret = (double)logq(Sum);


	return Ret;
}

void	FossilDepLhQuadDobule(NODE N, int SiteNo, int s00, int s01, int s10, int s11)
{
	if(s00 == 0)
		N->BigPartial[SiteNo][0] = 0;

	if(s01 == 0)
		N->BigPartial[SiteNo][1] = 0;

	if(s10 == 0)
		N->BigPartial[SiteNo][2] = 0;

	if(s11 == 0)
		N->BigPartial[SiteNo][3] = 0;
}

void	SetQuadDoubleNodeRec(NODE N, int NOS, int NoOfSites, RATES *Rates, OPTIONS *Opt)
{
	int SIndex, Index;
	QDOUBLE Sum;

	for(SIndex=0;SIndex<NoOfSites;SIndex++)
	{
		Sum = 0;
		for(Index=0;Index<NOS;Index++)
			Sum += N->BigPartial[SIndex][Index];

		for(Index=0;Index<NOS;Index++)
			N->Partial[SIndex][Index] =  N->BigPartial[SIndex][Index] / Sum;
	}
}

#endif

