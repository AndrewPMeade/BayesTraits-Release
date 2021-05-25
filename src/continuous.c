#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

//#define MATHMAT

#include "typedef.h"
#include "trees.h"
#include "genlib.h"
#include "data.h"
#include "likelihood.h"
#include "matrix.h"
#include "linalg.h"
#include "rand.h"
#include "rates.h"
#include "ckappa.h"
#include "contrasts.h"


void	InitEstData(OPTIONS *Opt, TREES *Trees)
{
	int		*TempEst;
	int		TIndex;
	int		SIndex;
	TAXA	*Taxa;

	TempEst = (int*)malloc(sizeof(int) * (Trees->NoOfSites+1));
	if(TempEst == NULL)
		MallocErr();

	for(SIndex=0;SIndex<Trees->NoOfSites+1;SIndex++)
		TempEst[SIndex] = FALSE;
	Opt->NoEstDataSite = 0;

	for(TIndex=0;TIndex<Trees->NoOfTaxa;TIndex++)
	{
		Taxa = &Trees->Taxa[TIndex];


		for(SIndex=0;SIndex<Trees->NoOfSites;SIndex++)
		{
			if(Taxa->EstDataP[SIndex] == TRUE)
			{
				if(TempEst[SIndex] == FALSE)
				{
					TempEst[SIndex] = TRUE;
					Opt->NoEstDataSite++;
				}
			}
		}

		if(Taxa->EstDepData == TRUE)
		{
			if(TempEst[Trees->NoOfSites] == FALSE)
			{
				TempEst[Trees->NoOfSites] = TRUE;
				Opt->NoEstDataSite++;
			}
		}
	}

	if(Opt->NoEstDataSite == 0)
	{
		free(TempEst);
		Opt->EstDataSites = NULL;
		return;
	}

	Opt->EstDataSites = (int*)malloc(sizeof(int) * Opt->NoEstDataSite);
	if(Opt->EstDataSites == NULL)
		MallocErr();

	TIndex=0;
	for(SIndex=0;SIndex<Trees->NoOfSites;SIndex++)
	{
		if(TempEst[SIndex] == TRUE)
			Opt->EstDataSites[TIndex++] = SIndex;
	}

	if(TempEst[Trees->NoOfSites] == TRUE)
		Opt->EstDataSites[TIndex] = -1;

	free(TempEst);
}

void	RemoveDependantData(OPTIONS *Opt, TREES *Trees)
{
	int		Index;
	int		TIndex;
	TAXA*	Taxa;
	int		DepNo;

	DepNo = Opt->DependantSite;

	for(TIndex=0;TIndex<Trees->NoOfTaxa;TIndex++)
	{
		Taxa = &Trees->Taxa[TIndex];

		if(Taxa->EstDataP[DepNo] == TRUE)
		{
			Taxa->EstDepData = TRUE;
			for(Index=DepNo+1;Index<Trees->NoOfSites;Index++)
				Taxa->EstDataP[Index-1] = Taxa->EstDataP[Index];
		}

		Taxa->Dependant = Taxa->ConData[DepNo];

		for(Index=DepNo+1;Index<Trees->NoOfSites;Index++)
			Taxa->ConData[Index-1] = Taxa->ConData[Index];
	}

	Trees->NoOfSites--;
}

NODE	TaxaToNode(TREES* Trees, TREE *Tree, TAXA *Taxa)
{
	int	NIndex;

	for(NIndex=0;NIndex<Tree->NoNodes;NIndex++)
	{
		if(Tree->NodeList[NIndex]->Tip == TRUE)
		{
			if(strcmp(Tree->NodeList[NIndex]->Taxa->Name, Taxa->Name)==0)
				return Tree->NodeList[NIndex];
		}
	}

	return NULL;
}

double	GetTinyBL(TREES *Trees)
{
	int		TIndex;
	int		NIndex;
	NODE	Node;
	double	Ret;

	Ret = 100000;

	for(TIndex=0;TIndex<Trees->NoOfTrees;TIndex++)
	{
		for(NIndex=0;NIndex<Trees->Tree[TIndex].NoNodes;NIndex++)
		{
			Node = Trees->Tree[TIndex].NodeList[NIndex];
			if((Node != Trees->Tree[TIndex].Root) &&
				(Node->Length != 0))
			{
				if(Node->Length < Ret)
					Ret = Node->Length;
			}
		}
	}

	return Ret;
}

void	CheckZeroTaxaBL(TREES *Trees)
{
	int		TIndex;
	int		NIndex;
	NODE	Node;
	double	TinyBL;

	/* Get the smallest branch in the tree */
	TinyBL = GetTinyBL(Trees);

	/* Make it even smaller */
	TinyBL = TinyBL * 0.001;

	for(TIndex=0;TIndex<Trees->NoOfTrees;TIndex++)
	{
		for(NIndex=0;NIndex<Trees->Tree[TIndex].NoNodes;NIndex++)
		{
			Node = Trees->Tree[TIndex].NodeList[NIndex];

			if((Node->Tip == TRUE) && (Node != Trees->Tree[TIndex].Root))
			{
				if(Node->Length < TinyBL)
					Node->Length = TinyBL;
			}
		}
	}
}

double	DistToRoot(NODE N, NODE Root)
{
	if(N == Root)
		return 0;

	return N->Length  + (DistToRoot(N->Ans, Root));
}

double	FindCoVar(TREES* Trees, TREE *Tree, int T1, int T2)
{
	int		Index;
	double	Ret;
	NODE	N1;
	NODE	N2;

	Ret = 0;

    for(Index=0;Index<Tree->NoNodes;Index++)
		Tree->NodeList[Index]->Visited = FALSE;

	N1 = TaxaToNode(Trees, Tree, &Trees->Taxa[T1]);
	N2 = TaxaToNode(Trees, Tree, &Trees->Taxa[T2]);

	while(N1!=Tree->Root)
	{
		N1 = N1->Ans;
		N1->Visited = TRUE;
	}

	do
	{
		N2 = N2->Ans;
	} while(N2->Visited != TRUE);

	Ret = DistToRoot(N2, Tree->Root);

	return Ret;
}

void	CalcZ(TREES* Trees, TREE *Tree, OPTIONS *Opt)
{
	int	SIndex, TIndex;
	int	ZPos;

	ZPos = 0;

	if(Opt->Model == CONTINUOUSREG)
	{
		for(SIndex=0;SIndex<Trees->NoOfTaxa;SIndex++)
			Tree->ConVars->Z[SIndex] = Trees->Taxa[SIndex].Dependant;
	}
	else
	{
		for(SIndex=0;SIndex<Trees->NoOfSites;SIndex++)
		{
			for(TIndex=0;TIndex<Trees->NoOfTaxa;TIndex++,ZPos++)
				Tree->ConVars->Z[ZPos] = Trees->Taxa[TIndex].ConData[SIndex];
		}
	}
}

void	CalcPVarCoVar(TREES* Trees, TREE *Tree)
{
	int		x,y;
	double	CoVar;
	NODE	N;
/*	double	*WV1, *WV2;
	MATRIX	*TMat;
	int	Ret;
*/
	for(x=0;x<Trees->NoOfTaxa;x++)
	{
		for(y=x;y<Trees->NoOfTaxa;y++)
		{
			if(x!=y)
			{
				CoVar = FindCoVar(Trees, Tree, x,y);
				Tree->ConVars->V->me[x][y] = CoVar;
				Tree->ConVars->V->me[y][x] = CoVar;
			}
			else
			{
				N = TaxaToNode(Trees, Tree, &Trees->Taxa[x]);
				CoVar = DistToRoot(N, Tree->Root);
				Tree->ConVars->V->me[x][x] = CoVar;
			}
		}
	}

	#ifdef IDMATRIX
		SetIdentityMatrix(Tree->ConVars->V);
	#endif

/*
	PrintMathematicaMatrix(Tree->VarCoVar, "Var Co var for a tree", stdout);
	PrintMatrix(Tree->VarCoVar, "Var Co var for a tree", stdout);

	WV1 = (double*)malloc(sizeof(double)*Tree->VarCoVar->NoOfCols);
	WV2 = (double*)malloc(sizeof(double)*Tree->VarCoVar->NoOfCols);

	TMat = AllocMatrix(Tree->VarCoVar->NoOfRows, Tree->VarCoVar->NoOfCols);

	Ret = InvertMatrix(Tree->VarCoVar->me, Tree->VarCoVar->NoOfCols, WV1, WV2, TMat->me);
	Ret = InvertMatrix(Tree->VarCoVar->me, Tree->VarCoVar->NoOfCols, WV1, (int*)WV2, TMat->me);

	printf("Ret\t%d\n", Ret);
	fflush(stdout);

	free(WV1);
	free(WV2);
	FreeMatrix(TMat);
	PrintMathematicaMatrix(TMat, "Inv vec", stdout);
*/
}

double	FindSum(TAXA* Taxa, int NoOfTaxa, int TraitNo)
{
	double	Ret=0;
	int		Index;

	for(Index=0;Index<NoOfTaxa;Index++)
		Ret += Taxa[Index].ConData[TraitNo];

	return Ret;
}

double	FindSumSqu(TAXA* Taxa, int NoOfTaxa, int TraitNo)
{
	double	Ret=0;
	int		Index;

	for(Index=0;Index<NoOfTaxa;Index++)
		Ret += (Taxa[Index].ConData[TraitNo] * Taxa[Index].ConData[TraitNo]);

	return Ret;
}

double	CalcVar(TREES* Trees, int TraitNo)
{
	double	Ret=0;
	double	Sum=0;
	double	SumSqu=0;

	Sum		= FindSum(Trees->Taxa, Trees->NoOfTaxa, TraitNo);
	SumSqu	= FindSumSqu(Trees->Taxa, Trees->NoOfTaxa, TraitNo);

	Ret = (Sum * Sum) / (double)Trees->NoOfTaxa;
	Ret = SumSqu - Ret;
	Ret = Ret / (double)(Trees->NoOfTaxa - 1);

	return Ret;
}

double	CalcCoVar(TREES* Trees, int T1, int T2)
{
	double	Ret=0;
	double	Sum1, Sum2;
	int		Index;

	Sum1 = FindSum(Trees->Taxa, Trees->NoOfTaxa, T1);
	Sum2 = FindSum(Trees->Taxa, Trees->NoOfTaxa, T2);

	for(Index=0;Index<Trees->NoOfTaxa;Index++)
		Ret = Ret + (Trees->Taxa[Index].ConData[T1] * Trees->Taxa[Index].ConData[T2]);

	Ret = Ret - ((Sum1 * Sum2) / (double)Trees->NoOfTaxa);
/*	Ret = Ret / (double)(Trees->NoOfTaxa - 1); */
	Ret = Ret / (double)(Trees->NoOfTaxa);

	return Ret;
}

void	CalcDVarCoVar(TREES* Trees, TREE *Tree)
{
	int		x,y;
	double	Val;

	for(x=0;x<Trees->NoOfSites;x++)
	{
		for(y=x;y<Trees->NoOfSites;y++)
		{
			if(x!=y)
			{
				Val = CalcCoVar(Trees, x, y);
				Tree->ConVars->Sigma->me[x][y] = Val;
				Tree->ConVars->Sigma->me[y][x] = Val;
			}
			else
			{
				Val = CalcVar(Trees, x);
				Tree->ConVars->Sigma->me[x][x] = Val;
			}
		}
	}
}

void	FreeConVar(CONVAR* ConVar, int NoTaxa)
{
	FreeMatrix(ConVar->V);
	FreeMatrix(ConVar->InvV);
	FreeMatrix(ConVar->Sigma);
	FreeMatrix(ConVar->InvSigma);
	FreeMatrix(ConVar->KProd);
	FreeMatrix(ConVar->InvKProd);

	if(ConVar->TrueV != NULL)
		FreeMatrix(ConVar->TrueV);

	if(ConVar->DepVect != NULL)
		free(ConVar->DepVect);

	free(ConVar->Alpha);
	free(ConVar->Z);
	free(ConVar->ZA);
	free(ConVar->ZATemp);

	free(ConVar->TVect1);
	free(ConVar->TVect2);
	free(ConVar->TVect3);
	free(ConVar->SVect);

	if(ConVar->TVT != NULL)
		FreeMatrix(ConVar->TVT);

	if(ConVar->TVTTemp != NULL)
		FreeMatrix(ConVar->TVTTemp);

	if(ConVar->InvXVX != NULL)
		FreeMatrix(ConVar->InvXVX);

	if(ConVar->Beta	!= NULL)
		free(ConVar->Beta);

	if(ConVar->TaxaDist != NULL)
		FreeCKappaTree(ConVar, NoTaxa);

	free(ConVar);
}

CONVAR*	AllocConVar(OPTIONS *Opt, TREES* Trees)
{
	CONVAR* Ret=NULL;
	int		Lager;

	Ret = (CONVAR*)malloc(sizeof(CONVAR));
	if(Ret==NULL)
		MallocErr();

	Ret->TVT		=	NULL;

	Ret->TVTTemp	=	NULL;
	Ret->InvXVX		=	NULL;
	Ret->Beta		=	NULL;

	if(Opt->Model == CONTINUOUSDIR)
	{
		Ret->TVT	=	AllocMatrix(2, 2);
		Ret->TVTTemp=	AllocMatrix(2, Trees->NoOfTaxa);
	}

	if(Opt->Model == CONTINUOUSREG)
	{
		if(Opt->AlphaZero == FALSE)
		{
			Ret->TVT	=	AllocMatrix(Trees->NoOfTaxa, Trees->NoOfSites + 1);
			Ret->TVTTemp=	AllocMatrix(Trees->NoOfTaxa, Trees->NoOfSites + 1);
		}
		else
		{
			Ret->TVT	=	AllocMatrix(Trees->NoOfTaxa, Trees->NoOfSites);
			Ret->TVTTemp=	AllocMatrix(Trees->NoOfTaxa, Trees->NoOfSites);
		}

		if(Opt->Analsis == ANALML)
			Ret->InvXVX	= AllocMatrix(Trees->NoOfSites+1, Trees->NoOfSites+1);
		else
			Ret->InvXVX = NULL;
	}

	Ret->V		= AllocMatrix(Trees->NoOfTaxa, Trees->NoOfTaxa);
	Ret->InvV	= AllocMatrix(Trees->NoOfTaxa, Trees->NoOfTaxa);

	if(Opt->InvertV == TRUE)
		Ret->TrueV = AllocMatrix(Trees->NoOfTaxa, Trees->NoOfTaxa);
	else
		Ret->TrueV = NULL;

	Ret->TaxaDist = NULL;


	if(Opt->Model == CONTINUOUSREG)
	{
		Ret->KProd		= AllocMatrix(Trees->NoOfTaxa, Trees->NoOfTaxa);
		Ret->InvKProd	= AllocMatrix(Trees->NoOfTaxa, Trees->NoOfTaxa);
		Ret->Sigma		= AllocMatrix(1, 1);
		Ret->InvSigma	= AllocMatrix(1, 1);
	}
	else
	{
		Ret->KProd		= AllocMatrix(Trees->NoOfSites * Trees->NoOfTaxa, Trees->NoOfSites * Trees->NoOfTaxa);
		Ret->InvKProd	= AllocMatrix(Trees->NoOfSites * Trees->NoOfTaxa, Trees->NoOfSites * Trees->NoOfTaxa);
		Ret->Sigma		= AllocMatrix(Trees->NoOfSites, Trees->NoOfSites);
		Ret->InvSigma	= AllocMatrix(Trees->NoOfSites, Trees->NoOfSites);
	}

	Ret->Beta	=	NULL;

	if((Opt->Model == CONTINUOUSDIR) || (Opt->Model == CONTINUOUSRR))
		Ret->Alpha	=	(double*)malloc(sizeof(double) * Trees->NoOfSites);
	else
		Ret->Alpha	=	(double*)malloc(sizeof(double) * 1);

	if(Ret->Alpha == NULL)
		MallocErr();

	if((Opt->Model == CONTINUOUSDIR) || (Opt->Model == CONTINUOUSREG))
	{
		Ret->Beta = (double*)malloc(sizeof(double) * Trees->NoOfSites);
		if(Ret->Beta == NULL)
			MallocErr();
	}

	if(Opt->Model == CONTINUOUSREG)
	{
		Ret->Z		=	(double*)malloc(sizeof(double) * Trees->NoOfTaxa);
		Ret->ZA		=	(double*)malloc(sizeof(double) * Trees->NoOfTaxa);
		Ret->ZATemp	=	(double*)malloc(sizeof(double) * Trees->NoOfTaxa);
		Ret->DepVect=	(double*)malloc(sizeof(double) * Trees->NoOfTaxa);
		if(Ret->DepVect == NULL)
			MallocErr();
	}
	else
	{
		Ret->Z		=	(double*)malloc(sizeof(double) * Trees->NoOfSites * Trees->NoOfTaxa);
		Ret->ZA		=	(double*)malloc(sizeof(double) * Trees->NoOfSites * Trees->NoOfTaxa);
		Ret->ZATemp	=	(double*)malloc(sizeof(double) * Trees->NoOfSites * Trees->NoOfTaxa);
		Ret->DepVect=	NULL;
	}

	if(	(Ret->Z		== NULL) ||
		(Ret->ZA	== NULL) ||
		(Ret->ZATemp== NULL))
		MallocErr();

	if(Trees->NoOfTaxa > Trees->NoOfSites)
		Lager = Trees->NoOfTaxa;
	else
		Lager = Trees->NoOfSites;

	Ret->TVect1	=	(double*)malloc(sizeof(double) * Lager);
	Ret->TVect2	=	(double*)malloc(sizeof(double) * Lager);
	Ret->TVect3	=	(double*)malloc(sizeof(double) * Lager);
	Ret->SVect	=	(double*)malloc(sizeof(double) * Lager);

	if( (Ret->TVect1 == NULL) ||
		(Ret->TVect2 == NULL) ||
		(Ret->TVect3 == NULL) ||
		(Ret->SVect  == NULL))
		MallocErr();

	Ret->LogDetOfSigma	= 1;
	Ret->LogDetOfV		= 1;

	return Ret;
}

void	FindInvV(TREES *Trees, TREE* Tree)
{
	int		Err;
	TEMPCONVAR*	TempCon;

	TempCon = Trees->TempConVars;

	CopyMatrix(TempCon->TMat, Tree->ConVars->V);

	Err = InvertMatrixAndDet(TempCon->TMat->me, Trees->NoOfTaxa, TempCon->T1, TempCon->T2, Tree->ConVars->InvV->me, &Tree->ConVars->LogDetOfV);

	if(Err == ERROR)
	{
		printf("V Matrix inverstion error in %s %d\n", __FILE__, __LINE__);
		PrintMathematicaMatrix(Tree->ConVars->V, "V=", stdout);
		exit(0);
	}
}

void	FindTVT(TREES* Trees, TREE *Tree, int AlphaZero)
{
	int			Row, Col, Index;
	double		Temp;
	MATRIX*		TMat;
	CONVAR*		ConVars;

	ConVars = Tree->ConVars;

	TMat = AllocMatrix(2, Trees->NoOfTaxa);

	for(Col=0;Col<Trees->NoOfTaxa;Col++)
	{
		if(AlphaZero == FALSE)
			ConVars->TVTTemp->me[0][Col] = 1;
		else
			ConVars->TVTTemp->me[0][Col] = 0;

		ConVars->TVTTemp->me[1][Col] = Tree->ConVars->V->me[Col][Col];
	}

	for(Row=0;Row<2;Row++)
	{
		for(Col=0;Col<Trees->NoOfTaxa;Col++)
		{
			Temp = 0;
			for(Index=0;Index<Trees->NoOfTaxa;Index++)
			{
				Temp += ConVars->TVTTemp->me[Row][Index] * ConVars->InvV->me[Index][Col];
			}
			TMat->me[Row][Col] = Temp;
		}
	}

	Temp = 0;
	for(Index=0;Index<Trees->NoOfTaxa;Index++)
		Temp += ConVars->TVTTemp->me[0][Index] * TMat->me[0][Index];

	ConVars->TVT->me[0][0] = Temp;

	Temp = 0;
	for(Index=0;Index<Trees->NoOfTaxa;Index++)
		Temp += ConVars->TVTTemp->me[1][Index] * TMat->me[0][Index];

	ConVars->TVT->me[0][1] = Temp;

	Temp = 0;
	for(Index=0;Index<Trees->NoOfTaxa;Index++)
		Temp += ConVars->TVTTemp->me[0][Index] * TMat->me[1][Index];

	ConVars->TVT->me[1][0] = Temp;

	Temp = 0;
	for(Index=0;Index<Trees->NoOfTaxa;Index++)
		Temp += ConVars->TVTTemp->me[1][Index] * TMat->me[1][Index];

	ConVars->TVT->me[1][1] = Temp;

	FreeMatrix(TMat);

	if(AlphaZero == TRUE)
	{
		ConVars->TVT->me[1][1] = 1.0 / ConVars->TVT->me[1][1];
		return;
	}

	Temp = ConVars->TVT->me[0][0];
	ConVars->TVT->me[0][0] = ConVars->TVT->me[1][1];
	ConVars->TVT->me[1][1] = Temp;

	Temp = 1.0 / (ConVars->TVT->me[0][0] * ConVars->TVT->me[1][1] - ConVars->TVT->me[1][0] * ConVars->TVT->me[0][1]);

	ConVars->TVT->me[0][0] = Temp * ConVars->TVT->me[0][0];
	ConVars->TVT->me[0][1] = Temp * -ConVars->TVT->me[0][1];
	ConVars->TVT->me[1][0] = Temp * -ConVars->TVT->me[1][0];
	ConVars->TVT->me[1][1] = Temp * ConVars->TVT->me[1][1];
}

void	FindMLRagVals(TREES* Trees, TREE *Tree, OPTIONS *Opt)
{
	TEMPCONVAR*	TempCon;
	int		x,y;
	CONVAR	*CV;


	CV = Tree->ConVars;
	TempCon = Trees->TempConVars;

	for(x=0;x<Trees->NoOfTaxa;x++)
		TempCon->Y[x] = Trees->Taxa[x].Dependant;


	if(Opt->AlphaZero == FALSE)
	{
		for(x=0;x<Trees->NoOfTaxa;x++)
			TempCon->X->me[x][0] = 1;

		for(x=1;x<Trees->NoOfSites+1;x++)
		{
			for(y=0;y<Trees->NoOfTaxa;y++)
				TempCon->X->me[y][x] = Trees->Taxa[y].ConData[x-1];
		}
	}
	else
	{
		for(x=0;x<Trees->NoOfSites;x++)
		{
			for(y=0;y<Trees->NoOfTaxa;y++)
				TempCon->X->me[y][x] = Trees->Taxa[y].ConData[x];
		}
	}

	/* Calc X'.InvV.Y */
	MatrixByVectMult(CV->InvV, TempCon->Y, CV->TVect1);
	VectByMatrixMult(CV->TVect1, TempCon->X, CV->TVect2);
	/* Now X'.InvV.Y is in TVect2 */

	/* Calc d */
	MatrixMult(CV->InvV, TempCon->X, CV->TVTTemp);

	Transpose(TempCon->X, TempCon->TranX);
	MatrixMult(TempCon->TranX, CV->TVTTemp, TempCon->NX);
	InvertMatrix(TempCon->NX->me, TempCon->NX->NoOfCols, CV->TVect1, (int*)CV->TVect3, CV->InvXVX->me);
	MatrixByVectMult(CV->InvXVX, CV->TVect2, CV->TVect1);

	if(Opt->AlphaZero == FALSE)
	{
		CV->Alpha[0] = CV->TVect1[0];
		for(x=1;x<Trees->NoOfSites+1;x++)
			CV->Beta[x-1] = CV->TVect1[x];
	}
	else
	{
		CV->Alpha[0] = 0;
		for(x=0;x<Trees->NoOfSites;x++)
			CV->Beta[x] = CV->TVect1[x];
	}

/*
	printf("ML Vals\n");
	for(x=0;x<Trees->NoOfSites+1;x++)
		printf("%f\t", CV->TVect1[x]);
	printf("\nML Done\n");
*/

/*	To Check With Mathematica.
	PrintMathematicaMatrix(X, "X = ", stdout);
	PrintMathematicaMatrix(CV->InvV, "InvV = ", stdout);
	PrintMathematicaVect(Y, Trees->NoOfTaxa, "Y = ", stdout);
	printf("B1 = Inverse[Transpose[X].InvV.X]\n");
	printf("B2 = Transpose[X].InvV.Y\n");
	printf("B = B1.B2\n");

	exit(0);
*/
}

void	MLFindAlphaBeta(TREES* Trees, TREE *Tree, int Site, int AlphaZero)
{
	int			Row, Col, Index;
	double		Temp;
	double		TempTVX[2];
	MATRIX*		TMat;
	CONVAR*		ConVars;

	ConVars = Tree->ConVars;


	TMat = AllocMatrix(2, Trees->NoOfTaxa);

	for(Row=0;Row<2;Row++)
	{
		for(Col=0;Col<Trees->NoOfTaxa;Col++)
		{
			Temp = 0;
			for(Index=0;Index<Trees->NoOfTaxa;Index++)
			{
				Temp += ConVars->TVTTemp->me[Row][Index] * ConVars->InvV->me[Index][Col];
			}
			TMat->me[Row][Col] = Temp;
		}
	}

	Temp = 0;
	for(Index=0;Index<Trees->NoOfTaxa;Index++)
		Temp += Trees->Taxa[Index].ConData[Site] * TMat->me[0][Index];

	TempTVX[0] = Temp;

	Temp = 0;
	for(Index=0;Index<Trees->NoOfTaxa;Index++)
		Temp += Trees->Taxa[Index].ConData[Site] * TMat->me[1][Index];

	TempTVX[1] = Temp;

	if(AlphaZero == TRUE)
	{
		Tree->ConVars->Alpha[Site] = 0;
		Tree->ConVars->Beta[Site] = ConVars->TVT->me[1][1] * TempTVX[1];

	}
	else
	{
		Tree->ConVars->Alpha[Site]	= (ConVars->TVT->me[0][0] * TempTVX[0]) + (ConVars->TVT->me[0][1] * TempTVX[1]);
		Tree->ConVars->Beta[Site]	= (ConVars->TVT->me[1][0] * TempTVX[0]) + (ConVars->TVT->me[1][1] * TempTVX[1]);
	}

	FreeMatrix(TMat);
}

double	MLFindAlphaReg(TREES* Trees, TREE *Tree, double *Data)
{
	double	P1=0;
	double	P2;
	double	ColTemp;
	int		x,y;

	for(x=0;x<Trees->NoOfTaxa;x++)
	{
		for(y=x+1;y<Trees->NoOfTaxa;y++)
			P1 = P1 + (2 * Tree->ConVars->InvV->me[x][y]);
	}

	for(x=0;x<Trees->NoOfTaxa;x++)
		P1 = P1 + Tree->ConVars->InvV->me[x][x];

	P1 = 1 / P1;

	P2 = 0;
	for(y=0;y<Trees->NoOfTaxa;y++)
	{
		ColTemp = 0;
		for(x=0;x<Trees->NoOfTaxa;x++)
			ColTemp = ColTemp + Tree->ConVars->InvV->me[x][y];
		P2 += ColTemp * Data[y];
	}

	return P1 * P2;
}

/*

double	MLFindAlphaReg(TREES* Trees, TREE *Tree)
{
	double	P1=0;
	double	P2;
	double	ColTemp;
	int		x,y;

	for(x=0;x<Trees->NoOfTaxa;x++)
	{
		for(y=x+1;y<Trees->NoOfTaxa;y++)
			P1 = P1 + (2 * Tree->ConVars->InvV->me[x][y]);
	}

	for(x=0;x<Trees->NoOfTaxa;x++)
		P1 = P1 + Tree->ConVars->InvV->me[x][x];

	P1 = 1 / P1;

	P2 = 0;
	for(y=0;y<Trees->NoOfTaxa;y++)
	{
		ColTemp = 0;
		for(x=0;x<Trees->NoOfTaxa;x++)
			ColTemp = ColTemp + Tree->ConVars->InvV->me[x][y];
		P2 += ColTemp * Trees->Taxa[y].Dependant;
	}

	return P1 * P2;
}
*/
double	MLFindAlphaMean(TREES* Trees, TREE *Tree, int Site)
{
	double	P1=0;
	double	P2;
	double	ColTemp;
	int		x,y;

	for(x=0;x<Trees->NoOfTaxa;x++)
	{
		for(y=x+1;y<Trees->NoOfTaxa;y++)
			P1 = P1 + (2 * Tree->ConVars->InvV->me[x][y]);
	}

	for(x=0;x<Trees->NoOfTaxa;x++)
		P1 = P1 + Tree->ConVars->InvV->me[x][x];

	P1 = 1 / P1;

	P2 = 0;
	for(y=0;y<Trees->NoOfTaxa;y++)
	{
		ColTemp = 0;
		for(x=0;x<Trees->NoOfTaxa;x++)
			ColTemp = ColTemp + Tree->ConVars->InvV->me[x][y];
		P2 += ColTemp * Trees->Taxa[y].ConData[Site];
	}

	return P1 * P2;
}

double	FindMLVarMatic(TREES* Trees, TREE *Tree)
{
	double	Ret;
	int		x,y;
	CONVAR	*CV;

	CV = Tree->ConVars;

	for(y=0;y<Trees->NoOfTaxa;y++)
	{
		Ret = 0;
		for(x=0;x<Trees->NoOfTaxa;x++)
			Ret = Ret + (CV->TVect1[x] * CV->InvV->me[x][y]);

		CV->TVect2[y] = Ret;
	}

	Ret = 0;
	for(x=0;x<Trees->NoOfTaxa;x++)
		Ret = Ret + (CV->TVect2[x] * CV->TVect3[x]);

/* Old values */
/*	return Ret * ((double)1/(Trees->NoOfTaxa));  */

	/* Least Squas results */
	return Ret * ((double)1/(Trees->NoOfTaxa - (Trees->NoOfSites+1)));
}

double	FindMLVar(TREES* Trees, TREE *Tree, int Site1, double Alpha1, double Beta1, int Site2, double Alpha2, double Beta2)
{
	int		x;
	CONVAR	*CV;

	CV = Tree->ConVars;

	for(x=0;x<Trees->NoOfTaxa;x++)
	{
		CV->TVect1[x] = Trees->Taxa[x].ConData[Site1]  - (Alpha1 + (Beta1 * CV->V->me[x][x]));
		CV->TVect3[x] = Trees->Taxa[x].ConData[Site2]  - (Alpha2 + (Beta2 * CV->V->me[x][x]));

		CV->TVect2[x] = 0;
	}

	return FindMLVarMatic(Trees, Tree);
}


double	FindMLRegVar(TREES* Trees, TREE *Tree)
{
	int		TIndex;
	int		x;
	CONVAR	*CV;
	double	Reg;

	CV = Tree->ConVars;

	for(TIndex=0;TIndex<Trees->NoOfTaxa;TIndex++)
	{
		Reg = CV->Alpha[0];
		for(x=0;x<Trees->NoOfSites;x++)
			Reg += CV->Beta[x] * Trees->Taxa[TIndex].ConData[x];

		CV->TVect1[TIndex] = Trees->Taxa[TIndex].Dependant  - Reg;
		CV->TVect3[TIndex] = CV->TVect1[TIndex];
		CV->TVect2[TIndex] = 0;
	}

	Reg = FindMLVarMatic(Trees, Tree);

	return Reg;
}

void	CalcSigma(OPTIONS *Opt, TREES* Trees, TREE *Tree, double* Means, double* Beta)
{
	int		x,y;
	double	Val;

	if((Opt->Analsis == ANALML) && (Opt->Model == CONTINUOUSDIR))
		FindTVT(Trees, Tree, Opt->AlphaZero);

	if((Opt->Analsis == ANALML) && (Opt->Model == CONTINUOUSREG))
		FindMLRagVals(Trees, Tree, Opt);

	for(x=0;x<Trees->NoOfSites;x++)
	{
		if(Opt->Analsis == ANALML)
		{
			if(Opt->Model == CONTINUOUSDIR)
				MLFindAlphaBeta(Trees, Tree, x, Opt->AlphaZero);

			if(Opt->Model == CONTINUOUSRR)
			{
				if(Opt->AlphaZero == FALSE)
					Tree->ConVars->Alpha[x] = MLFindAlphaMean(Trees, Tree, x);
				else
					Tree->ConVars->Alpha[x] = 0.0;
			}
		}
		else
		{
			if((Opt->Model == CONTINUOUSRR) || (Opt->Model == CONTINUOUSDIR))
				Tree->ConVars->Alpha[x] = Means[x];
			else
				Tree->ConVars->Alpha[0] = Means[0];

			if((Opt->Model == CONTINUOUSDIR) || (Opt->Model == CONTINUOUSREG))
				Tree->ConVars->Beta[x] = Beta[x];
		}
	}

	if(Opt->Model == CONTINUOUSREG)
	{
		Tree->ConVars->Sigma->me[0][0] = FindMLRegVar(Trees, Tree);
		return;
	}

	for(x=0;x<Trees->NoOfSites;x++)
	{
		for(y=x;y<Trees->NoOfSites;y++)
		{
			if(Opt->TestCorrel == TRUE)
			{
				if(Opt->Model == CONTINUOUSRR)
					Val = FindMLVar(Trees, Tree, x, Tree->ConVars->Alpha[x], 0, y, Tree->ConVars->Alpha[y], 0);
				else
					Val = FindMLVar(Trees, Tree, x, Tree->ConVars->Alpha[x], Tree->ConVars->Beta[x], y, Tree->ConVars->Alpha[y], Tree->ConVars->Beta[y]);


				Tree->ConVars->Sigma->me[x][y] = Val;
				Tree->ConVars->Sigma->me[y][x] = Val;
			}
			else
			{
				if(x==y)
				{
					if(Opt->Model == CONTINUOUSRR)
						Tree->ConVars->Sigma->me[x][x] = FindMLVar(Trees, Tree, x, Tree->ConVars->Alpha[x], 0, y, Tree->ConVars->Alpha[y], 0);
					else
						Tree->ConVars->Sigma->me[x][x] = FindMLVar(Trees, Tree, x, Tree->ConVars->Alpha[x], Tree->ConVars->Beta[x], y, Tree->ConVars->Alpha[y], Tree->ConVars->Beta[y]);

				}
				else
				{
					Tree->ConVars->Sigma->me[x][y] = 0;
					Tree->ConVars->Sigma->me[y][x] = 0;
				}
			}
		}
	}
}

double	CalcZAReg(TREES* Trees, TREE* Tree, int TaxaNo)
{
	double	Ret;
	int		Index;
	CONVAR	*CV;
	TAXA	*Taxa;

	CV	= Tree->ConVars;
	Taxa= &Trees->Taxa[TaxaNo];

	Ret = CV->Alpha[0];
	for(Index=0;Index<Trees->NoOfSites;Index++)
		Ret += CV->Beta[Index] * Taxa->ConData[Index];

	return Ret;
}

void	RegCalcZAlpha(TREES* Trees, TREE *Tree)
{
	int	TIndex;
	CONVAR	*CV;

	CV = Tree->ConVars;

	for(TIndex=0;TIndex<Trees->NoOfTaxa;TIndex++)
		CV->ZA[TIndex] = CV->Z[TIndex] - CalcZAReg(Trees, Tree, TIndex);

}

void	CalcZAlpha(TREES* Trees, TREE *Tree, MODEL Model)
{
	int	SIndex, TIndex;
	int	ZPos;

	ZPos = 0;

	if(Model == CONTINUOUSREG)
	{
		RegCalcZAlpha(Trees, Tree);
		return;
	}

	for(SIndex=0;SIndex<Trees->NoOfSites;SIndex++)
	{
		for(TIndex=0;TIndex<Trees->NoOfTaxa;TIndex++,ZPos++)
		{
			switch(Model)
			{
				case CONTINUOUSRR:
					Tree->ConVars->ZA[ZPos] = Tree->ConVars->Z[ZPos] - Tree->ConVars->Alpha[SIndex];
				break;

				case CONTINUOUSDIR:
					Tree->ConVars->ZA[ZPos] = Tree->ConVars->Z[ZPos] - (Tree->ConVars->Alpha[SIndex] + (Tree->ConVars->Beta[SIndex] * Tree->ConVars->V->me[TIndex][TIndex]));
				break;
			}
		}
	}
}


double	FindDet(TREES* Trees, TREE *Tree, OPTIONS *Opt)
{
	if(Opt->Model == CONTINUOUSREG)
		return (Tree->ConVars->LogDetOfV) + ((double)Trees->NoOfTaxa * Tree->ConVars->LogDetOfSigma);

	return ((double)Trees->NoOfSites * Tree->ConVars->LogDetOfV) + ((double)Trees->NoOfTaxa * Tree->ConVars->LogDetOfSigma);
}

/*
<< LinearAlgebra`MatrixManipulation`;
BlockMatrix[Outer[Times, S, V]]
*/

void	CalcDelta(MATRIX *V, double Delta)
{
	int x,y;

	for(x=0;x<V->NoOfCols;x++)
	{
		for(y=0;y<V->NoOfRows;y++)
		{
			V->me[x][y]= pow(V->me[x][y], Delta);
		}
	}

}

void	CalcLabda(MATRIX *V, double Labda)
{
	int x,y;

	for(x=0;x<V->NoOfCols;x++)
	{
		for(y=x;y<V->NoOfRows;y++)
		{
			if(x!=y)
			{
				V->me[x][y] = V->me[x][y] * Labda;
				V->me[y][x] = V->me[x][y];
			}
		}
	}
}
#ifdef akdlkljkajlk
MATRIX*	FindRegVar(TREES *Trees, RATES* Rates)
{
	MATRIX	*Ret;
	static MATRIX	*XT=NULL;
	static MATRIX	*X=NULL;
	static MATRIX	*TempV1;
	static MATRIX	*TempV2;
	TAXA	*Taxa;
	TREE	*Tree;
	CONVAR	*CV;
	int		x,y;

	Tree = &Trees->Tree[Rates->TreeNo];
	CV = Tree->ConVars;

	if(X == NULL)
	{
		X		= AllocMatrix(Trees->NoOfTaxa, Trees->NoOfSites+1);
		XT		= AllocMatrix(Trees->NoOfSites+1, Trees->NoOfTaxa);
		TempV1	= AllocMatrix(Trees->NoOfSites+1, Trees->NoOfTaxa);
		TempV2	= AllocMatrix(Trees->NoOfSites+1, Trees->NoOfSites+1);

		for(x=0;x<Trees->NoOfTaxa;x++)
		{
			Taxa = &Trees->Taxa[x];
			X->me[x][0] = 1;
			for(y=0;y<Trees->NoOfSites;y++)
				X->me[x][y+1] = Taxa->ConData[y];
		}

		Transpose(X, XT);
	}

	Ret = AllocMatrix(Trees->NoOfSites+1, Trees->NoOfSites+1);

	/* Do
		Sig*Inverse[Transpose[X].InvV.X]
	*/
	MatrixMult(XT, CV->InvV, TempV1);
	MatrixMult(TempV1, X, TempV2);

	InvertMatrix(TempV2->me, Trees->NoOfSites+1, CV->TVect1,(int*)CV->TVect2, Ret->me);

	ScaleMatrix(Ret, CV->Sigma->me[0][0]);

/*	printf("\n");
	PrintMathematicaMatrix(X, "X=", stdout);
	PrintMathematicaMatrix(CV->InvV, "InvV=", stdout);
	printf("Sig=%f;\n", CV->Sigma->me[0][0]);

	printf("Sqrt[Sig*Inverse[Transpose[X].InvV.X]]\n");
	exit(0);
*/

	return Ret;
}
#endif

MATRIX*	FindRegVar(TREES *Trees, RATES* Rates, int AlphaZero)
{
	MATRIX	*Ret;
	TEMPCONVAR	*TempCon;

	/*	static MATRIX	*XT=NULL;
	static MATRIX	*X=NULL;
	static MATRIX	*TempV1;
	static MATRIX	*TempV2; */
	TREE	*Tree;
	CONVAR	*CV;
	int		Size;

	TempCon = Trees->TempConVars;

	Tree = &Trees->Tree[Rates->TreeNo];
	CV = Tree->ConVars;

	if(AlphaZero == FALSE)
		Size = Trees->NoOfSites+1;
	else
		Size = Trees->NoOfSites;

	Ret = AllocMatrix(Size, Size);

	/* Do
		Sig*Inverse[Transpose[X].InvV.X]
	*/
	MatrixMult(TempCon->XT, CV->InvV, TempCon->TempV1);
	MatrixMult(TempCon->TempV1, TempCon->RVX, TempCon->TempV2);

	InvertMatrix(TempCon->TempV2->me, Size, CV->TVect1,(int*)CV->TVect2, Ret->me);

	ScaleMatrix(Ret, CV->Sigma->me[0][0]);

/*	printf("\n");
	PrintMathematicaMatrix(X, "X=", stdout);
	PrintMathematicaMatrix(CV->InvV, "InvV=", stdout);
	printf("Sig=%f;\n", CV->Sigma->me[0][0]);

	printf("Sqrt[Sig*Inverse[Transpose[X].InvV.X]]\n");
	exit(0);
*/

	return Ret;
}


void	SetEstData(TREES *Trees, RATES* Rates)
{
	int		TIndex;
	TAXA*	Taxa;
	int		SIndex;
	int		RIndex;

	RIndex=0;
	for(TIndex=0;TIndex<Trees->NoOfTaxa;TIndex++)
	{
		Taxa = &Trees->Taxa[TIndex];

		if(Taxa->EstData == TRUE)
		{
			for(SIndex=0;SIndex<Trees->NoOfSites;SIndex++)
			{
				if(Taxa->EstDataP[SIndex] == TRUE)
					Taxa->ConData[SIndex] = Rates->EstData[RIndex++];
			}

			if(Taxa->EstDepData == TRUE)
				Taxa->Dependant = Rates->EstData[RIndex++];
		}
	}
}

void	PrintMathmatCode(void)
{
/*
	printf("Part3 = -0.5*ZA.Inverse[BlockMatrix[Outer[Times, Sigma, (V^Delta)]]].ZA\n");
	printf("Part2 = Log[Det[BlockMatrix[Outer[Times, (Sigma), (V^Delta)]]]^-0.5]\n");
	printf("Part1 = 1.8378770664093453 ((-(NoOfTaxa*NoOfSites))/2)\n");
	printf("Lh = Part1 + Part2 + Part3\n");
	printf("Lh = (Log[2 Pi] ((-(NoOfTaxa*NoOfSites))/2)) + (Log[Det[BlockMatrix[Outer[Times, Sigma, (V^Delta)]]]^-0.5]) + (-0.5*ZA.Inverse[BlockMatrix[Outer[Times, Sigma, (V^Delta)]]].ZA)\n");
*/

	printf("Part3 = -0.5*ZA.Inverse[ArrayFlatten[Outer[Times, Sigma, V]]].ZA;\n");
	printf("Part2 = Log[Det[ArrayFlatten[Outer[Times, (Sigma), V]]]^-0.5];\n");
	printf("Part1 = 1.8378770664093453 ((-(NoOfTaxa*NoOfSites))/2);\n");
	printf("Lh = Part1 + Part2 + Part3;\n");
	printf("Lh = (Log[2 Pi] ((-(NoOfTaxa*NoOfSites))/2)) + (Log[Det[ArrayFlatten[Outer[Times, Sigma, V]]]^-0.5]) + (-0.5*ZA.Inverse[ArrayFlatten[Outer[Times, Sigma, V]]].ZA)\n");

}

double	LHRandWalk(OPTIONS *Opt, TREES* Trees, RATES* Rates)
{
	double	Val;
	MATRIX	*TMat;
	int		Index;
	int		Len;
	double	Det;
	double	Ret;
	TREE	*Tree;
	CONVAR	*CV;

	Tree = &Trees->Tree[Rates->TreeNo];
	CV = Tree->ConVars;

	if(Rates->UseEstData == TRUE)
		SetEstData(Trees, Rates);

	if(Opt->UseVarData == TRUE)
		SetVarData(Trees, Opt->VarData, Rates->VarDataSite);

	if((Rates->UseEstData == TRUE) || (Opt->UseVarData == TRUE))
		CalcZ(Trees, Tree, Opt);

#ifdef MATHMAT

	printf("(* **************** Start Tree %d *********************** *)\n", Rates->TreeNo);
	printf("NoOfSites = %d;\n", Trees->NoOfSites);
	printf("NoOfTaxa = %d;\n", Trees->NoOfTaxa);

	if(Tree->ConVars->TrueV != NULL)
		PrintMathematicaTFMatrix(Tree->ConVars->TrueV, "V = ", stdout);
	else
		PrintMathematicaTFMatrix(Tree->ConVars->V, "V = ", stdout);

//	PrintMatrix(Tree->ConVars->V, "V = ", stdout);
#endif

	if(Opt->InvertV	== TRUE)
	{
		if(Opt->EstKappa == TRUE)
			MakeKappaV(Trees, Tree, Rates->Kappa);
		else
			CopyMatrix(Tree->ConVars->V, Tree->ConVars->TrueV);

		if(Opt->EstDelta == TRUE)
			CalcDelta(Tree->ConVars->V, Rates->Delta);

		if(Opt->FixDelta != -1)
			CalcDelta(Tree->ConVars->V, Opt->FixDelta);

		if(Opt->EstLambda == TRUE)
			CalcLabda(Tree->ConVars->V, Rates->Lambda);

		if(Opt->FixLambda != -1)
			CalcLabda(Tree->ConVars->V, Opt->FixLambda);

		FindInvV(Trees, Tree);
	}

	CalcSigma(Opt, Trees, Tree, Rates->Means, Rates->Beta);

#ifdef MATHMAT
	PrintMathematicaMatrix(Tree->ConVars->Sigma, "Sigma = ", stdout);
#endif

//	if(TMat == NULL)
		TMat = AllocMatrix(Trees->NoOfSites, Trees->NoOfSites);

	if(Opt->Model == CONTINUOUSREG)
	{
		CV->InvSigma->me[0][0] = 1 / CV->Sigma->me[0][0];
		CV->LogDetOfSigma = log(CV->Sigma->me[0][0]);
	}
	else
	{
		if(TMat == NULL)
			TMat = AllocMatrix(Trees->NoOfSites, Trees->NoOfSites);
		CopyMatrix(TMat, CV->Sigma);
		InvertMatrixAndDet(TMat->me, Trees->NoOfSites, CV->TVect1, (int*)CV->TVect2, CV->InvSigma->me, &CV->LogDetOfSigma);
	}

	CalcZAlpha(Trees, Tree, Opt->Model);

#ifdef MATHMAT
	PrintMathematicaVect(Tree->ConVars->ZA, Trees->NoOfSites * Trees->NoOfTaxa, "ZA = ", stdout);
#endif

	KroneckerProduct(Tree->ConVars->InvSigma, Tree->ConVars->InvV, Tree->ConVars->InvKProd);

	VectByMatrixMult(Tree->ConVars->ZA, Tree->ConVars->InvKProd, Tree->ConVars->ZATemp);

	if(Opt->Model == CONTINUOUSREG)
		Len = Trees->NoOfTaxa;
	else
		Len = Trees->NoOfSites * Trees->NoOfTaxa;

	Val = 0;
	for(Index=0;Index<Len;Index++)
		Val += Tree->ConVars->ZATemp[Index] * Tree->ConVars->ZA[Index];

	Val = (-0.5) * Val;

	Det = FindDet(Trees, Tree, Opt);

	Det = -0.5 * Det;

	if(Opt->Model == CONTINUOUSREG)
		Ret = -(double)(Trees->NoOfTaxa);
	else
		Ret = -(double)(Trees->NoOfSites * Trees->NoOfTaxa);

	Ret = Ret / 2.0;

	Ret = Ret * 1.837877066409345483560659472811;

#ifdef MATHMAT

	PrintMathmatCode();

	printf("(* Part3 = %f *)\n", Val);
	printf("(* Part2 = %f *)\n", Det);
	printf("(* Part1 = %f *)\n", Ret);
	printf("(* LH = %f *)\n", Ret + Det + Val);

	printf("(* **************** End Tree %d *********************** *)\n", Rates->TreeNo);
#endif

	Ret = Ret + Det + Val;

	if(Ret == Ret + 1)
		return ERRLH;

	FreeMatrix(TMat);

	Rates->Lh = Ret;

	return Ret;
}

void	TreeBLToPower(TREES *Trees, TREE *Tree, double Power)
{
	int	Index;
	NODE N;

	for(Index=0;Index<Tree->NoNodes;Index++)
	{
		N = Tree->NodeList[Index];
		if(N != Tree->Root)
			N->Length = pow(N->Length, Power);
	}
}


void	InitContinusTree(OPTIONS *Opt, TREES* Trees, int TreeNo)
{
	int		Index;
	CONVAR	*CV;

	if(Opt->Model == CONTRASTM)
	{
		InitContrastTree(Opt, Trees, TreeNo);
		return;
	}

	Trees->Tree[TreeNo].ConVars = AllocConVar(Opt, Trees);
	CV = Trees->Tree[TreeNo].ConVars;

	if((Opt->NodeData == TRUE) || (Opt->NodeBLData == TRUE))
		SetTreeAsData(Opt, Trees, TreeNo);

	if(Opt->FixKappa != -1)
		TreeBLToPower(Trees, &Trees->Tree[TreeNo], Opt->FixKappa);

	CalcPVarCoVar(Trees, &Trees->Tree[TreeNo]);

	if(Opt->InvertV == TRUE)
		CopyMatrix(Trees->Tree[TreeNo].ConVars->TrueV, Trees->Tree[TreeNo].ConVars->V);

	if(Opt->FixDelta != -1)
		CalcDelta(CV->V, Opt->FixDelta);

	if(Opt->FixLambda != -1)
		CalcLabda(CV->V, Opt->FixLambda);

	CalcZ(Trees, &Trees->Tree[TreeNo], Opt);
	FindInvV(Trees, &Trees->Tree[TreeNo]);

	if(Opt->EstKappa == TRUE)
	{
		InitCKappaTree(Trees, &Trees->Tree[TreeNo]);
		KappaVarCoVar(Trees, &Trees->Tree[TreeNo]);
	}

	if(Opt->Model == CONTINUOUSREG)
	{
		for(Index=0;Index<Trees->NoOfTaxa;Index++)
			CV->DepVect[Index] = Trees->Taxa[Index].Dependant;
	}
}

void		FreeTempConVars(TEMPCONVAR* TempCon)
{
	free(TempCon->T1);
	free(TempCon->T2);
	FreeMatrix(TempCon->TMat);

	FreeMatrix(TempCon->X);
	FreeMatrix(TempCon->TranX);
	FreeMatrix(TempCon->NX);
	free(TempCon->Y);

	FreeMatrix(TempCon->RVX);
	FreeMatrix(TempCon->XT);
	FreeMatrix(TempCon->TempV1);
	FreeMatrix(TempCon->TempV2);

	free(TempCon);
}

TEMPCONVAR* AllocTempConVars(OPTIONS *Opt, TREES* Trees)
{
	TEMPCONVAR*	Ret;
	int			Size;
	TAXA		*Taxa;
	int			x,y;

	Ret = (TEMPCONVAR*) malloc(sizeof(TEMPCONVAR));
	if(Ret == NULL)
		MallocErr();

	/* Statics form FindInvV */
	Ret->T1 = (double*)malloc(sizeof(double)*Trees->NoOfTaxa);
	Ret->T2 = (int*)malloc(sizeof(int)*Trees->NoOfTaxa);
	Ret->TMat = AllocMatrix(Trees->NoOfTaxa, Trees->NoOfTaxa);
	if((Ret->T1==NULL) || (Ret->T2==NULL) || (Ret->TMat == NULL))
			MallocErr();

	/* Statics from FindMLRagVals */
	if(Opt->AlphaZero == FALSE)
	{
		Ret->X		= AllocMatrix(Trees->NoOfTaxa, Trees->NoOfSites + 1);
		Ret->TranX	= AllocMatrix(Trees->NoOfSites + 1, Trees->NoOfTaxa);
		Ret->NX		= AllocMatrix(Trees->NoOfSites + 1, Trees->NoOfSites + 1);
	}
	else
	{
		Ret->X		= AllocMatrix(Trees->NoOfTaxa, Trees->NoOfSites);
		Ret->TranX	= AllocMatrix(Trees->NoOfSites, Trees->NoOfTaxa);
		Ret->NX		= AllocMatrix(Trees->NoOfSites, Trees->NoOfSites);
	}

	Ret->Y	= (double*)malloc(sizeof(double) * Trees->NoOfTaxa);

	/* Statics from FindRegVar */
	if(Opt->AlphaZero == FALSE)
		Size = Trees->NoOfSites+1;
	else
		Size = Trees->NoOfSites;

	Ret->RVX	= AllocMatrix(Trees->NoOfTaxa, Size);
	Ret->XT		= AllocMatrix(Size, Trees->NoOfTaxa);
	Ret->TempV1	= AllocMatrix(Size, Trees->NoOfTaxa);
	Ret->TempV2	= AllocMatrix(Size, Size);

	if(Opt->AlphaZero == FALSE)
	{
		for(x=0;x<Trees->NoOfTaxa;x++)
		{
			Taxa = &Trees->Taxa[x];
			Ret->RVX->me[x][0] = 1;
			for(y=0;y<Trees->NoOfSites;y++)
				Ret->RVX->me[x][y+1] = Taxa->ConData[y];
		}
	}
	else
	{
		for(x=0;x<Trees->NoOfTaxa;x++)
		{
			Taxa = &Trees->Taxa[x];
			for(y=0;y<Trees->NoOfSites;y++)
				Ret->RVX->me[x][y] = Taxa->ConData[y];
		}
	}

	Transpose(Ret->RVX, Ret->XT);

	return Ret;
}

void	InitContinus(OPTIONS *Opt, TREES* Trees)
{
	int		TIndex;

	CheckZeroTaxaBL(Trees);

	if(Opt->Model == CONTRASTM)
	{
		if(Opt->Analsis == ANALMCMC)
			InitContrastAll(Opt, Trees);
		return;
	}

	if(Opt->Model == CONTINUOUSREG)
		RemoveDependantData(Opt, Trees);


	AddRecNodes(Opt, Trees);

	Trees->TempConVars = AllocTempConVars(Opt, Trees);

	InitEstData(Opt, Trees);

	if(Opt->UseVarData == TRUE)
		LoadVarData(Opt);

	Opt->InvertV = FALSE;
	if(	(Opt->EstDelta == TRUE) ||
		(Opt->EstKappa == TRUE) ||
		(Opt->EstLambda== TRUE) ||
		(Opt->UsePhyloPlasty == TRUE))
		Opt->InvertV = TRUE;

	if(Opt->Analsis == ANALMCMC)
	{
		for(TIndex=0;TIndex<Trees->NoOfTrees;TIndex++)
			InitContinusTree(Opt, Trees, TIndex);
	}
}
