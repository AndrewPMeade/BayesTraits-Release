#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

#include "typedef.h"
#include "genlib.h"
#include "BigLh.h"
#include "likelihood.h"

#ifndef BIG_LH
void	InitTreeBigLh(OPTIONS *Opt, TREES *Trees) { }
void	FreeTreeBigLh(OPTIONS *Opt, TREES *Trees) { }

void	LhBigLh(NODE N, TREES *Trees, int Pre, int SiteNo) { }
double	CombineBigLh(RATES* Rates, TREES *Trees, OPTIONS *Opt, int SiteNo, int NOS) {return 0.0;}

void	SetBigLhNodeRec(NODE N, int NOS, int NoOfSites, RATES *Rates, OPTIONS *Opt) { }

void	FossilLhBig(NODE N, TREES *Trees, int SiteNo) { }
void	FossilDepLhBig(NODE N, int SiteNo, int s00, int s01, int s10, int s11) { }
#else

void	AllocNodeMemBigLh(NODE N, OPTIONS *Opt, TREES *Trees)
{
	int i,j;

	N->BigPartial = (mpfr_t**)malloc(sizeof(mpfr_t*) * Trees->NoOfSites);
	if(N->BigPartial == NULL)
		MallocErr();

	for(i=0;i<Trees->NoOfSites;i++)
	{
		N->BigPartial[i] = (mpfr_t*)malloc(sizeof(mpfr_t) * Trees->NoOfStates);
		if(N->BigPartial[i] == NULL)
			MallocErr();
		for(j=0;j<Trees->NoOfStates;j++)
		{
			mpfr_init2(N->BigPartial[i][j], Opt->Precision);
			mpfr_set_d(N->BigPartial[i][j], 0.0, DEF_ROUND);
		}

	}


	mpfr_init2(N->t1, Opt->Precision);
	mpfr_init2(N->t2, Opt->Precision);
	mpfr_init2(N->t3, Opt->Precision);

	mpfr_set_d(N->t1, 0.0, DEF_ROUND);
	mpfr_set_d(N->t2, 0.0, DEF_ROUND);
	mpfr_set_d(N->t3, 0.0, DEF_ROUND);

//	mpfr_t		**BigPartial;

}

void	AllocMemory(OPTIONS *Opt, TREES *Trees)
{
	int TIndex, NIndex;
	NODE N;
	TREE *T;

	for(TIndex=0;TIndex<Trees->NoOfTrees;TIndex++)
	{
		T = Trees->Tree[TIndex];
		for(NIndex=0;NIndex<T->NoNodes;NIndex++)
		{
			N = T->NodeList[NIndex];
			AllocNodeMemBigLh(N, Opt, Trees);
		}
	}
}

void	SetTipDataNodeBigLh(NODE N, OPTIONS *Opt, TREES *Trees)
{
	int i,j;

	for(i=0;i<Trees->NoOfSites;i++)
	{
		for(j=0;j<Trees->NoOfStates;j++)
			mpfr_set_d(N->BigPartial[i][j], N->Partial[i][j], DEF_ROUND);
	}
}



void	SetTipDataBigLh(OPTIONS *Opt, TREES *Trees)
{
	int TIndex, NIndex;
	NODE N;
	TREE *T;

	for(TIndex=0;TIndex<Trees->NoOfTrees;TIndex++)
	{
		T = Trees->Tree[TIndex];
		for(NIndex=0;NIndex<T->NoNodes;NIndex++)
		{
			N = T->NodeList[NIndex];
			SetTipDataNodeBigLh(N, Opt, Trees);
		}
	}

}

void	InitTreeBigLh(OPTIONS *Opt, TREES *Trees)
{
	AllocMemory(Opt, Trees);

	SetTipDataBigLh(Opt, Trees);
}

void	FreeNodeMemBigLh(NODE N, OPTIONS *Opt, TREES *Trees)
{
	int i,j;

	for(i=0;i<Trees->NoOfSites;i++)
	{
		for(j=0;j<Trees->NoOfStates;j++)
			mpfr_clear(N->BigPartial[i][j]);
		free(N->BigPartial[i]);
	}

	free(N->BigPartial);

	mpfr_clears(N->t1, N->t2, N->t3, NULL);
}

void	FreeTreeBigLh(OPTIONS *Opt, TREES *Trees)
{
	int TIndex, NIndex;
	NODE N;
	TREE *T;

	for(TIndex=0;TIndex<Trees->NoOfTrees;TIndex++)
	{
		T = Trees->Tree[TIndex];
		for(NIndex=0;NIndex<T->NoNodes;NIndex++)
		{
			N = T->NodeList[NIndex];
			FreeNodeMemBigLh(N, Opt, Trees);
		}
	}
}

void SetBigLHZero(mpfr_t State)
{
	mpfr_set_d(State, 0.0, DEF_ROUND);
}

void	FossilLhBig(NODE N, TREES *Trees, int SiteNo)
{
	int Index;

	for(Index=0;Index<Trees->NoOfStates;Index++)
	{
		if(Index != N->FossilState)
			SetBigLHZero(N->BigPartial[SiteNo][Index]);
	}
}

void	FossilDepLhBig(NODE N, int SiteNo, int s00, int s01, int s10, int s11)
{
	if(s00 == 0)
		SetBigLHZero(N->BigPartial[SiteNo][0]);

	if(s01 == 0)
		SetBigLHZero(N->BigPartial[SiteNo][1]);

	if(s10 == 0)
		SetBigLHZero(N->BigPartial[SiteNo][2]);

	if(s11 == 0)
		SetBigLHZero(N->BigPartial[SiteNo][3]);
}



void	LhBigLh(NODE N, TREES *Trees, int Pre, int SiteNo)
{
	int Inner, Outter, NIndex;
	double Lh;
	double **Mat;
	mpfr_t **Partail;

	// Set N->BigPartial to know value
	for(Outter=0;Outter<Trees->NoOfStates;Outter++)
	{
		mpfr_set_d(N->BigPartial[SiteNo][Outter], 1.0, DEF_ROUND);
		for(NIndex=0;NIndex<N->NoNodes;NIndex++)
		{
			Mat = Trees->PList[N->NodeList[NIndex]->ID]->me;
			Partail = N->NodeList[NIndex]->BigPartial;

			mpfr_set_d(N->t1, 0.0, DEF_ROUND);
			for(Inner=0;Inner<Trees->NoOfStates;Inner++)
			{
				mpfr_mul_d(N->t2, Partail[SiteNo][Inner], Mat[Outter][Inner], DEF_ROUND);
				mpfr_add(N->t3, N->t1, N->t2, DEF_ROUND);
				mpfr_set(N->t1, N->t3, DEF_ROUND);
			}

			mpfr_mul(N->t2, N->t1, N->BigPartial[SiteNo][Outter], DEF_ROUND);
			mpfr_set(N->BigPartial[SiteNo][Outter], N->t2, DEF_ROUND);
		}
	}

	if(N->FossilState != -1)
		FossilLh(N, Trees, SiteNo);
}

double CombineBigLh(RATES* Rates, TREES *Trees, OPTIONS *Opt, int SiteNo, int NOS)
{
	int Index;
	double Ret;
	NODE Root;
	TREE *Tree;

	Ret = 0;

	Tree = Trees->Tree[Rates->TreeNo];
	Root = Tree->Root;

	mpfr_set_d(Root->t1, 0.0, DEF_ROUND);
	mpfr_set_d(Root->t2, 0.0, DEF_ROUND);
	mpfr_set_d(Root->t3, 0.0, DEF_ROUND);

	for(Index=0;Index<NOS;Index++)
	{
		mpfr_mul_d(Root->t1, Root->BigPartial[SiteNo][Index], Rates->Pis[Index], DEF_ROUND);
		mpfr_add(Root->t3, Root->t1, Root->t2, DEF_ROUND);
		mpfr_set(Root->t2, Root->t3, DEF_ROUND);
	}

	mpfr_log(Root->t1, Root->t2, DEF_ROUND);
	Ret = mpfr_get_d(Root->t1, DEF_ROUND);

	for(Index=0;Index<NOS;Index++)
	{
		mpfr_div(Root->t3, Root->BigPartial[SiteNo][Index], Root->t2, DEF_ROUND);
		Root->Partial[SiteNo][Index] = mpfr_get_d(Root->t3, DEF_ROUND);
	}

	return Ret;
}

void	SetBigLhNodeRec(NODE N, int NOS, int NoOfSites, RATES *Rates, OPTIONS *Opt)
{
	int SIndex, Index;

	for(SIndex=0;SIndex<NoOfSites;SIndex++)
	{
		mpfr_set_d(N->t1, 0.0, DEF_ROUND);
		for(Index=0;Index<NOS;Index++)
		{
			mpfr_add(N->t2, N->t1, N->BigPartial[SIndex][Index], DEF_ROUND);
			mpfr_set(N->t1, N->t2, DEF_ROUND);
		}

		for(Index=0;Index<NOS;Index++)
		{
			mpfr_div(N->t2, N->BigPartial[SIndex][Index], N->t1, DEF_ROUND);
			N->Partial[SIndex][Index] =  mpfr_get_d(N->t2, DEF_ROUND);
		}
	}
}
#endif

