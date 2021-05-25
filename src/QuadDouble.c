#include <stdio.h>
#include <stdlib.h>
#include <stdio.h>

#include "typedef.h"
#include "QuadDouble.h"
#include "genlib.h"


#ifndef QUAD_DOUBLE
void	InitQuadDoubleLh(OPTIONS *Opt, TREES *Trees) {}
void	FreeQuadLh(OPTIONS *Opt, TREES *Trees) {}

void	NodeLhQuadDouble(NODE N, TREES *Trees, int SiteNo) {}
		

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
		N->BigPartial[SIndex] = (QDOUBLE*)malloc(sizeof(QDOUBLE) * Trees->NoOfStates);
		if(N->BigPartial[SIndex] == NULL)
			MallocErr();

		for(Index=0;Index<Trees->NoOfStates;Index++)
			N->BigPartial[SIndex][Index] = N->Partial[SIndex][Index];
	}
}

void	InitQuadDoubleLh(OPTIONS *Opt, TREES *Trees)
{
	int TIndex, NIndex;
	TREE *Tree;
	NODE N;

	for(TIndex=0;TIndex<Trees->NoOfTrees;TIndex++)
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
	
	for(TIndex=0;TIndex<Trees->NoOfTrees;TIndex++)
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

void	NodeLhQuadDouble(NODE N, TREES *Trees, int SiteNo)
{
	int		Inner, Outter, NIndex;
	QDOUBLE	Lh;
	QDOUBLE **TBigLh;
	double **Mat;

	for(Outter=0;Outter<Trees->NoOfStates;Outter++)
	{
		N->BigPartial[SiteNo][Outter] = 1.0;

		for(NIndex=0;NIndex<N->NoNodes;NIndex++)
		{
			Mat = Trees->PList[N->NodeList[NIndex]->ID]->me;
			TBigLh = N->NodeList[NIndex]->BigPartial;
			
			Lh = 0;
			for(Inner=0;Inner<Trees->NoOfStates;Inner++)
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

