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
#include <math.h>

#include "typedef.h"
#include "genlib.h"
#include "RJDummy.h"
#include "likelihood.h"
#include "options.h"
#include "contrasts.h"

int			GetMaxDummy(OPTIONS *Opt, TREES *Trees)
{
	return Trees->Tree[0]->NoNodes * 2;
	return 10;
}


void		FreeDummyCode(DUMMYCODE* DC)
{
	free(DC->Beta);
	free(DC);
}

DUMMYCODE*	AllocDummyCode(NODE N)
{
	DUMMYCODE *Ret;

	Ret = (DUMMYCODE*)malloc(sizeof(DUMMYCODE));
	if(Ret == NULL)
		MallocErr();

	Ret->Beta = (double*)malloc(sizeof(double) * 2);
	if(Ret->Beta == NULL)
		MallocErr();

	Ret->Beta[0] = -1;
	Ret->Beta[1] = -1;

	Ret->Node = N;

	return Ret;
}

DUMMYCODE*	CloneDummyCode(DUMMYCODE* DC)
{
	DUMMYCODE	*Ret;
	
	Ret = AllocDummyCode(DC->Node);

	Ret->Type = DC->Type;
		
	Ret->Beta[0] = DC->Beta[0];
	Ret->Beta[1] = DC->Beta[1];

	Ret->Iteration = DC->Iteration;

	return Ret;
}

RJDUMMY*	CreatRJDummyCode(OPTIONS *Opt, TREES *Trees)
{
	RJDUMMY *Ret;

	Ret = (RJDUMMY*)malloc(sizeof(RJDUMMY));
	if(Ret == NULL)
		MallocErr();

	Ret->NoDummyCode	= 0;
	
	Ret->NoMaxDummy		= GetMaxDummy(Opt, Trees);
	
	Ret->DummyList		= (DUMMYCODE**)malloc(sizeof(DUMMYCODE*) * Ret->NoMaxDummy);
	if(Ret->DummyList == NULL)
		MallocErr();

	Ret->DummyBeta		= (double*)malloc(sizeof(double) * Ret->NoMaxDummy * 2);
	if(Ret->DummyBeta == NULL)
		MallocErr();

	return Ret;
}

void		FreeRJDummyCode(RJDUMMY *RJDummy)
{
	int Index;

	for(Index=0;Index<RJDummy->NoDummyCode;Index++)
		FreeDummyCode(RJDummy->DummyList[Index]);
	
	free(RJDummy->DummyList);

	free(RJDummy->DummyBeta);

	free(RJDummy);
}

void		BuildDummyCodeBeta(RJDUMMY *RJDummy)
{
	int Index;
	DUMMYCODE *DC;

	for(Index=0;Index<RJDummy->NoDummyCode;Index++)
	{
		DC = RJDummy->DummyList[Index];

		RJDummy->DummyBeta[Index] = DC->Beta[0];
	}
}

int			ValidDummyNode(TREE *Tree, NODE N)
{
	if(N == Tree->Root)
		return FALSE;

	if(N->Part->NoTaxa < 5)
		return FALSE;

//	if(N->Tip == TRUE)
//		return FALSE;

	return TRUE;
}

NODE		GetDummyCodeNode(RANDSTATES *RS, TREE *Tree)
{
	int Pos;
	NODE N;

	do
	{
		Pos = RandUSInt(RS) % Tree->NoNodes;
		N = Tree->NodeList[Pos];
	} while(ValidDummyNode(Tree, N) == FALSE);

	return N;
}

int CompRJDummy(const void *av, const void *bv)
{
	DUMMYCODE **a, **b;

	a = (DUMMYCODE**)av;
	b = (DUMMYCODE**)bv;

	return (*b)->Node->Part->NoTaxa - (*a)->Node->Part->NoTaxa;
}


void		PrintDummyCodeData(OPTIONS *Opt, TREES *Trees, RATES *Rates)
{
	int NoSites, Index, TIndex;
	TAXA *Taxa;
	NODE N;
	TREE *Tree;

	MapDummyCodes(Trees, Rates);
//	ClearDummyCode(Trees, Rates);

	printf("Output\t");
	for(Index=0;Index<Rates->RJDummy->NoDummyCode;Index++)
		printf("%d\t", Rates->RJDummy->DummyList[Index]->Node->Part->NoTaxa);
	printf("\n\n");

	NoSites = Rates->Contrast->NoSites;

	Tree = Trees->Tree[Rates->TreeNo];

	for(TIndex=0;TIndex<Tree->NoNodes;TIndex++)
	{
		N = Tree->NodeList[TIndex];
		if(N->Tip == TRUE)
		{
			Taxa = N->Taxa;
			printf("%s\t", Taxa->Name);
			for(Index=0;Index<NoSites;Index++)
				printf("%f\t", N->ConData->Contrast[0]->Data[Index]);
			printf("\n");
		}
	}

	exit(0);
}

void		CheckDummyList(RJDUMMY *RJDummy)
{
	int Index;
	DUMMYCODE *DC1, *DC2;

	for(Index=1;Index<RJDummy->NoDummyCode;Index++)
	{
		DC1 = RJDummy->DummyList[Index-1];
		DC2 = RJDummy->DummyList[Index];

		if(DC1->Node->Part->NoTaxa < DC2->Node->Part->NoTaxa)
		{
			printf("err\n");
			exit(0);
		}
	}
}

void		SortDummyCodes(RJDUMMY *RJDummy)
{
	qsort(RJDummy->DummyList, RJDummy->NoDummyCode, sizeof(DUMMYCODE*), CompRJDummy);
}

void		AddRJDummyCode(long long Itter, OPTIONS *Opt, TREES *Trees, RATES *Rates, NODE N)
{
	CONTRASTR	*CRates;
	TREE		*Tree;
	RJDUMMY		*RJDummy;
	DUMMYCODE	*DC;

	RJDummy = Rates->RJDummy;
	Tree	= Trees->Tree[Rates->TreeNo];

	CRates = Rates->Contrast;

	DC = AllocDummyCode(N);

	RJDummy->DummyList[RJDummy->NoDummyCode] = DC;

	DC->Iteration = Itter;
	DC->Beta[0] = (RandDouble(Rates->RS) * 10) - 5;
	DC->Type = RJDUMMY_INTER;
	
	if(RandDouble(Rates->RS) < 0.5)
	{
		DC->Type = RJDUMMY_INTER_SLOPE;
		DC->Beta[1] = (RandDouble(Rates->RS) * 10) - 5;
		CRates->NoSites++;
	}
	
	CRates->NoSites++;
	RJDummy->NoDummyCode++;

	Rates->LnHastings	= 0;
	Rates->LnJacobion	= 0;
	Rates->LhPrior		= 0;

	SortDummyCodes(RJDummy);
//	CheckDummyList(RJDummy);
//	PrintDummyCodeData(Opt, Trees, Rates);
}


void		DelDummyCode(RATES *Rates, RJDUMMY *RJDummy, int Pos)
{
	CONTRASTR		*CRates;
	int				Index;
	RJDUMMY_TYPE	Type;

	CRates = Rates->Contrast;

	Type = RJDummy->DummyList[Pos]->Type;

	FreeDummyCode(RJDummy->DummyList[Pos]);
	
	for(Index=Pos;Index<RJDummy->NoDummyCode-1;Index++)
		RJDummy->DummyList[Index] = RJDummy->DummyList[Index+1];
	
	RJDummy->NoDummyCode--;
	CRates->NoSites--;

	if(Type == RJDUMMY_INTER_SLOPE)
		CRates->NoSites--;

	SortDummyCodes(RJDummy);
}

void		DeleteRJDummyCode(OPTIONS *Opt, TREES *Trees, RATES *Rates, int Pos)
{
	CONTRASTR	*CRates;
	RJDUMMY		*RJDummy;


	RJDummy = Rates->RJDummy;

	CRates = Rates->Contrast;

	DelDummyCode(Rates, RJDummy, Pos);

	Rates->LnHastings	= 0;
	Rates->LnJacobion	= 0;
	Rates->LhPrior		= 0;
			
//	PrintDummyCodeData(Opt, Trees, Rates);
}

int			GetPosFromNode(NODE N, RJDUMMY *RJDummy)
{
	int Index;

	for(Index=0;Index<RJDummy->NoDummyCode;Index++)
	{
		if(RJDummy->DummyList[Index]->Node == N)
			return Index;
	}

	return -1;
}


void		RJDummyMove(long long Itter, OPTIONS *Opt, TREES *Trees, RATES *Rates)
{
	RJDUMMY		*RJDummy;
	NODE		N;
	TREE		*Tree;
	int			Pos;

	Tree = Trees->Tree[Rates->TreeNo];

	RJDummy = Rates->RJDummy;

	N = GetDummyCodeNode(Rates->RS, Tree);

	Pos = GetPosFromNode(N, RJDummy);
	if(Pos == -1)
		AddRJDummyCode(Itter, Opt, Trees, Rates, N);
	else
		DeleteRJDummyCode(Opt, Trees, Rates, Pos);

	
}


/*
void		RJDummyMove(OPTIONS *Opt, TREES *Trees, RATES *Rates)
{
	RJDUMMY		*RJDummy;
	RJDummy = Rates->RJDummy;

	if(RJDummy->NoDummyCode == 0)
	{
		AddRJDummyCode(Opt, Trees, Rates);
		return;
	}

	if(RJDummy->NoDummyCode == RJDummy->NoMaxDummy)
	{
		DeleteRJDummyCode(Opt, Trees, Rates);
		return;
	}

	if(RandDouble(Rates->RS) < 0.5)
		AddRJDummyCode(Opt, Trees, Rates);
	else
		DeleteRJDummyCode(Opt, Trees, Rates);
}
*/

void		RJDummyMoveNode( OPTIONS *Opt, TREES *Trees, RATES *Rates)
{
	TREE		*Tree;
	NODE		N;
	RJDUMMY		*RJDummy;
	DUMMYCODE	*DC;
	int			Pos;

//	return;
	Tree = Trees->Tree[Rates->TreeNo];

	RJDummy = Rates->RJDummy;
	
	Pos = RandUSInt(Rates->RS) % RJDummy->NoDummyCode;
	DC = RJDummy->DummyList[Pos];


	N = DC->Node;

	if((RandDouble(Rates->RS) < 0.33) || (N->Tip == TRUE))
		N = N->Ans;
	else
	{
		Pos = RandUSInt(Rates->RS) % N->NoNodes;
		N = N->NodeList[Pos];
	}
		
	if(ValidDummyNode(Tree, N) == TRUE)
		DC->Node = N;

	SortDummyCodes(RJDummy);
}

void	RJDummyBlank(RJDUMMY *RJDummy)
{
	int Index;

	for(Index=0;Index<RJDummy->NoDummyCode;Index++)
		FreeDummyCode(RJDummy->DummyList[Index]);

	RJDummy->NoDummyCode = 0;
}

void	RJDummyCopy(RATES *A, RATES *B)
{
	RJDUMMY	*DA, *DB;
	int Index;

	DA = A->RJDummy;
	DB = B->RJDummy;
	
	RJDummyBlank(DA);

	for(Index=0;Index<DB->NoDummyCode;Index++)
		DA->DummyList[Index] = CloneDummyCode(DB->DummyList[Index]);

	DA->NoDummyCode = DB->NoDummyCode;

	A->Contrast->NoSites = B->Contrast->NoSites;
}

void		SetDummyCodeNode(int TNoSites, DUMMYCODE *DC, NODE N, int Pos)
{
	int Index;
	CONDATA *CD;

	if(N->Tip == TRUE)
	{
		CD = N->ConData;

		for(Index=TNoSites;Index<Pos;Index++)
			CD->Contrast[0]->Data[Index] = 0;
		
		CD->Contrast[0]->Data[Pos] = 1;
		if(DC->Type == RJDUMMY_INTER_SLOPE)
			CD->Contrast[0]->Data[Pos] = CD->Contrast[0]->Data[1];
		
		return;
	}

	for(Index=0;Index<N->NoNodes;Index++)
		SetDummyCodeNode(TNoSites, DC, N->NodeList[Index], Pos);
}

void		ClearDummyCode(TREES *Trees, RATES *Rates)
{
	int			Index, X;
	CONDATA		*CD;
	CONTRAST	*C;
	TREE		*Tree;
	NODE		N;

	Tree = Trees->Tree[Rates->TreeNo];

	for(Index=0;Index<Tree->NoNodes;Index++)
	{
		N = Tree->NodeList[Index];
		if(N->Tip == TRUE)
		{
			CD = N->ConData;
			C = CD->Contrast[0];

			for(X=Trees->NoSites;X<Rates->Contrast->NoSites;X++)
				C->Data[X] = 0;
		}
	}
}

int			TaxaHasCrossProduct(CONTRAST *C, int NoRSites, int TotalSites)
{
	int Index;

	for(Index=NoRSites;Index<TotalSites;Index++)
	{
		if((C->Data[Index] != 0) && (C->Data[Index] != 1))
			return TRUE;
	}

	return FALSE;
}

void		ZeroImpliedDummyCode(TREES *Trees, TREE *Tree, int NoRSites, int TotalSites)
{
	int Index;
	NODE N;
	CONTRAST *C;

//	return;

	for(Index=0;Index<Tree->NoNodes;Index++)
	{
		N = Tree->NodeList[Index];

		if(N->Tip == TRUE)
		{
			C = N->ConData->Contrast[0];
			if(TaxaHasCrossProduct(C, NoRSites, TotalSites) == TRUE)
				C->Data[1] = 0;
			else
				C->Data[1] = DataToZScore(N->Taxa->ConData[1], Trees->PMean[1], Trees->PSD[1]);
		}
	}
}

void		MapDummyCodes(TREES *Trees, RATES *Rates)
{
	RJDUMMY		*RJDummy;
	DUMMYCODE	*DC;
	int			Index, Pos;
	TREE		*Tree;
	
	RJDummy = Rates->RJDummy;

	ClearDummyCode(Trees, Rates);

	Pos = Trees->NoSites;
	for(Index=0;Index<RJDummy->NoDummyCode;Index++)
	{		
		DC = RJDummy->DummyList[Index];

		SetDummyCodeNode(Trees->NoSites, DC, DC->Node, Pos);
	
		if(DC->Type == RJDUMMY_INTER_SLOPE)
			Pos += 2;
		else
			Pos++;
	}

	Tree = Trees->Tree[Rates->TreeNo];

	ZeroImpliedDummyCode(Trees, Tree, Trees->NoSites, Rates->Contrast->NoSites);
}


void	InitRJDummyFile(OPTIONS *Opt)
{
	Opt->RJDummyLog = OpenWriteWithExt(Opt->BaseOutputFN, OUTPUT_EXT_DUMMY_CODE);

	PrintOptions(Opt->RJDummyLog, Opt);

	fprintf(Opt->RJDummyLog, "Iteration\tLh\tNo Dummy Codes\tBeta\tNo Taxa\tTaxa List");


	fprintf(Opt->RJDummyLog, "\n");

	fflush(Opt->RJDummyLog);
}

void	PrintRJDummyTaxaList(FILE *Out, TREES *Trees, DUMMYCODE *DC)
{
	int Index;
	PART *P;
	TAXA *Taxa;

	P = DC->Node->Part;

	for(Index=0;Index<P->NoTaxa;Index++)
	{
		Taxa = Trees->Taxa[P->Taxa[Index]];
		fprintf(Out, "%s", Taxa->Name);
		if(Index<P->NoTaxa-1)
			fprintf(Out, ",");
	}

	fprintf(Out, "\t");
}

void	PrintRJDummy(long long Itter, OPTIONS *Opt, TREES *Trees, RATES *Rates)
{
	RJDUMMY		*RJDummy;
	DUMMYCODE	*DC;
	int			Index, X;
	NODE		N;
	
	RJDummy = Rates->RJDummy;

	fprintf(Opt->RJDummyLog, "%lld\t%f\t%d\t", Itter, Rates->Lh, RJDummy->NoDummyCode);

	for(Index=0;Index<RJDummy->NoDummyCode;Index++)
	{
		DC = RJDummy->DummyList[Index];

		if(DC->Type == RJDUMMY_INTER)
		{
			fprintf(Opt->RJDummyLog, "Intercept\t");
			fprintf(Opt->RJDummyLog, "%f\t", DC->Beta[0]);
		}
		else
		{
			fprintf(Opt->RJDummyLog, "Intercept+Slope\t");
			fprintf(Opt->RJDummyLog, "%f,", DC->Beta[0]);
			fprintf(Opt->RJDummyLog, "%f\t", DC->Beta[1]);
		}
		
		fprintf(Opt->RJDummyLog, "%llu\t", DC->Iteration);
		
		fprintf(Opt->RJDummyLog, "%d\t", DC->Node->Part->NoTaxa);
		PrintRJDummyTaxaList(Opt->RJDummyLog, Trees, DC);
	}

	fprintf(Opt->RJDummyLog, "\n");

	for(Index=0;Index<RJDummy->NoDummyCode;Index++)
	{
		DC = RJDummy->DummyList[Index];
		fprintf(Opt->RJDummyLog, "#FF0000\tNoName\t1\t");
		N = DC->Node;
		for(X=0;X<N->Part->NoTaxa;X++)
		{
			fprintf(Opt->RJDummyLog, "%s\t", Trees->Taxa[N->Part->Taxa[X]]->Name);
		}

		fprintf(Opt->RJDummyLog, "\n");
	}


	fflush(Opt->RJDummyLog);
}

void	RJDummyChange(OPTIONS *Opt, TREES *Trees, RATES *Rates)
{
	RJDUMMY		*RJDummy;
	DUMMYCODE	*DC;
	int			Pos;
	double		Dev;

	RJDummy = Rates->RJDummy;

	Pos = RandUSInt(Rates->RS) % RJDummy->NoDummyCode;

	DC = RJDummy->DummyList[Pos];

	Dev = Opt->RJDummyBetaDev;

	if((DC->Type == RJDUMMY_INTER_SLOPE) && (RandDouble(Rates->RS) < 0.5))
		DC->Beta[1] += (RandDouble(Rates->RS) * Dev) - (Dev / 2.0);
	else
		DC->Beta[0] += (RandDouble(Rates->RS) * Dev) - (Dev / 2.0);
}

void	TestDummyCodeSig(OPTIONS *Opt, TREES *Trees, RATES *Rates)
{
	RJDUMMY	*RJDummy;
	double	Lh;

	RJDummy = Rates->RJDummy;

	while(RJDummy->NoDummyCode)
	{
		Lh = Likelihood(Rates, Trees, Opt);

		printf("%d\t%f\t", RJDummy->NoDummyCode, Lh);
		DelDummyCode(Rates, RJDummy, 0);
		
		Lh = Likelihood(Rates, Trees, Opt);
		printf("%f\n", Lh);
	}

	exit(0);
}

