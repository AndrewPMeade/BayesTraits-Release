#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "typedef.h"
#include "rates.h"
#include "genlib.h"
#include "rand.h"
#include "trees.h"
#include "continuous.h"
#include "revjump.h"
#include "priors.h"
#include "likelihood.h"
#include "data.h"
#include "matrix.h"

double**	LoadModelFile(char* FileName, int *NoModels, RATES* Rates);
void		SetFixedModel(RATES *Rates);

int		FindNoOfRates(OPTIONS *Opt)
{
	int	Index;
	int	Ret=0;

	if(Opt->UseRModel == TRUE)
	{
		if(Opt->RModelP == -1)
			Ret++;
	}
	else
	{
		for(Index=0;Index<Opt->NoOfRates;Index++)
			if(Opt->ResTypes[Index] == RESNONE)
				Ret++;	
	}

	if((Opt->UseCovarion == TRUE) && (Opt->Analsis == ANALML))
		Ret+=2;

	if((Opt->EstKappa == TRUE) && (Opt->Analsis == ANALML))
		Ret++;

	if((Opt->EstGamma == TRUE) && (Opt->Analsis == ANALML))
		Ret++;

	return Ret;
}

int		FindRatePos(int Rate, OPTIONS *Opt)
{
	int	Pos;

	if((Opt->ResTypes[Rate] == RESNONE) ||(Opt->ResTypes[Rate] == RESCONST))
		return Rate;

	Pos = Rate;
	do
	{
		Pos = Opt->ResNo[Pos];
		if((Opt->ResTypes[Pos] == RESNONE) || (Opt->ResTypes[Pos] == RESCONST))
			return Pos;

	}while(1==1);
	return Pos;
}

double	FindRateVal(int Pos, RATES *Rates, OPTIONS *Opt)
{
	int	RateIndex=0;
	int	OptIndex=0;

	if(Opt->ResTypes[Pos] == RESCONST)
		return Opt->ResConst[Pos];

	for(;;)
	{
		if(OptIndex==Pos)
			return Rates->Rates[RateIndex];

		
		if(Opt->ResTypes[OptIndex] == RESNONE)
			RateIndex++;
		OptIndex++;
	}
	return -1;
}

void	MapMCMCConRates(RATES* Rates, OPTIONS *Opt)
{
	int	Index;

	if(Opt->Model == CONTINUOUSREG)
		Rates->Means[0] = Rates->Rates[0];
	else
	{
		for(Index=0;Index<Opt->Trees->NoOfSites;Index++)
				Rates->Means[Index] = Rates->Rates[Index];
	}

	if(Opt->Model == CONTINUOUSRR)
		return;

	if(Opt->Model == CONTINUOUSDIR)
	{
		for(;Index<Rates->NoOfRates;Index++)
			Rates->Beta[Index - Opt->Trees->NoOfSites] = Rates->Rates[Index];
		return;
	}

	for(Index=1;Index<Rates->NoOfRates;Index++)
		Rates->Beta[Index - 1] = Rates->Rates[Index];
}

void	MapRates(RATES* Rates, OPTIONS *Opt)
{
	int	Index;
	int	Pos;

	if(Opt->DataType == CONTINUOUS)
	{
		if(Opt->Analsis == ANALMCMC)
			MapMCMCConRates(Rates, Opt);
	
		return;
	}

	if(Opt->UseRJMCMC == TRUE)
	{
		MapRJRates(Rates->Rates, Rates->MappingVect, Rates->NoOfFullRates, Rates->FullRates);
		return;
	}

	for(Index=0;Index<Rates->NoOfRates;Index++)
	{
		if((Rates->Rates[Index] < MINRATE) || (IsNum(Rates->Rates[Index]) == FALSE))
			Rates->Rates[Index] = MINRATE;

		if(Rates->Rates[Index] > MAXRATE)
			Rates->Rates[Index] = MAXRATE;
	}

	if(Opt->UseRModel == TRUE)
	{
		if(Opt->RModelP != -1)
			Rates->FullRates[0] = Opt->RModelP;
		else
			Rates->FullRates[0] = Rates->Rates[0];
	}
	else
	{
		for(Index=0;Index<Rates->NoOfFullRates;Index++)
		{
			Pos = FindRatePos(Index, Opt);
			Rates->FullRates[Index] = FindRateVal(Pos, Rates, Opt);

			if(Rates->FullRates[Index] < MINRATE)
				Rates->FullRates[Index] = MINRATE;

			if(Rates->FullRates[Index] > MAXRATE)
				Rates->FullRates[Index] = MAXRATE;
		}
	}
	

	Pos = Rates->NoOfRates;
	if((Opt->UseCovarion == TRUE) && (Opt->Analsis == ANALML))
		Pos = Pos - 2;

	if((Opt->EstKappa == TRUE) && (Opt->Analsis == ANALML))
		Pos = Pos - 1;

	if((Opt->EstGamma == TRUE) && (Opt->Analsis == ANALML))
		Pos = Pos - 1;

	if((Opt->UseCovarion == TRUE) && (Opt->Analsis == ANALML))
	{
		Rates->OnToOff = Rates->Rates[Pos++];
		Rates->OffToOn = Rates->Rates[Pos++];
	/*	Rates->OffToOn = Rates->OnToOff; */
	}

	if((Opt->EstKappa == TRUE) && (Opt->Analsis == ANALML))
	{
		Rates->Kappa = Rates->Rates[Pos++];
	
		if(Rates->Kappa < 0)
			Rates->Kappa = 0.0000001;

		if(Rates->Kappa > 5)
			Rates->Kappa = 5;
	}

	if((Opt->EstGamma == TRUE) && (Opt->Analsis == ANALML))
	{
		Rates->Gamma = Rates->Rates[Pos++];
	
		if(Rates->Gamma < GAMMAMIN)
			Rates->Gamma = GAMMAMIN;

		if(Rates->Gamma > GAMMAMAX)
			Rates->Gamma = GAMMAMAX;
	}

}

void	FindEmpPis(RATES *Rates, OPTIONS *Opt)
{
	TREES	*Trees;
	double	*TempPis;
	int		State;
	double	Weight;
	double	Total;
	int		SIndex,TIndex,SymbolIndex;
	TAXA	*Taxa;

	Trees = Opt->Trees;

	TempPis = (double*)malloc(sizeof(double)*Trees->NoOfStates);
	if(TempPis == NULL)
		MallocErr();
	for(SIndex=0;SIndex<Trees->NoOfStates;SIndex++)
		TempPis[SIndex] = 0;

	for(SIndex=0;SIndex<Trees->NoOfSites;SIndex++)
	{
		for(TIndex=0;TIndex<Trees->NoOfTaxa;TIndex++)
		{
			Taxa = &Trees->Taxa[TIndex];

			if(SiteHadUnKnownState(Taxa->DesDataChar[SIndex]) == FALSE)
			{
				Weight = (double)1/(double)strlen(Taxa->DesDataChar[SIndex]);

				for(SymbolIndex=0;SymbolIndex<(int)strlen(Taxa->DesDataChar[SIndex]);SymbolIndex++)
				{
					State = SymbolToPos(Taxa->DesDataChar[SIndex][SymbolIndex], Trees->SymbolList);
					TempPis[State] += Weight;
				}
			}
			else
			{
				Weight = (double)1/(double)Trees->NoOfStates;
				for(SymbolIndex=0;SymbolIndex<Trees->NoOfStates;SymbolIndex++)
					TempPis[SymbolIndex] += Weight;
			}
		}
	}

	Total = (double)Trees->NoOfSites * (double)Trees->NoOfTaxa;

	for(SIndex=0;SIndex<Trees->NoOfStates;SIndex++)
		Rates->Pis[SIndex] = TempPis[SIndex] / Total;

	free(TempPis);
}

void	SetPiValues(RATES *Rates, OPTIONS *Opt)
{
	int		Index;
	TREES	*Trees;

	Trees = Opt->Trees;

	if((Opt->PiTypes == PIUNI) || (Opt->PiTypes == PIEST))
	{
		for(Index=0;Index<Trees->NoOfStates;Index++)
			Rates->Pis[Index] = (double)1/Trees->NoOfStates;

		return;
	}

	if(Opt->PiTypes == PINONE)
	{
		for(Index=0;Index<Trees->NoOfStates;Index++)
			Rates->Pis[Index] = 1;

		return;
	}

	FindEmpPis(Rates, Opt);
}

void	CreatCRates(OPTIONS *Opt, RATES *Rates)
{
	int		Index;
	TREES	*Trees;
	TAXA	*Taxa;

	Rates->Delta = 1;
	Rates->Kappa = 1;
	Rates->Lambda= 1;

	Rates->Prios = NULL;

	if(Opt->Analsis == ANALMCMC)
	{
		Rates->NoOfRates = Opt->Trees->NoOfSites;
		
		if(Opt->Model == CONTINUOUSDIR)
			Rates->NoOfRates = Rates->NoOfRates  * 2;

		if(Opt->Model == CONTINUOUSREG)
			Rates->NoOfRates = Opt->Trees->NoOfSites + 1;

		Rates->NoOfFullRates = Rates->NoOfRates;

		Rates->Rates = (double*)malloc(sizeof(double) * Rates->NoOfRates);
		if(Rates->Rates == NULL)
			MallocErr();

		if(Opt->Model == CONTINUOUSREG)
			Rates->Means = (double*)malloc(sizeof(double) * 1);
		else
			Rates->Means = (double*)malloc(sizeof(double) * Opt->Trees->NoOfSites);
		
		if(Rates->Means == NULL)
			MallocErr();

		if((Opt->Model == CONTINUOUSDIR) || (Opt->Model == CONTINUOUSREG))
		{
			Rates->Beta = (double*)malloc(sizeof(double) * Opt->Trees->NoOfSites);
			if(Rates->Beta== NULL)
				MallocErr();
		}
		else
			Rates->Beta = NULL;


		if(Opt->UseVarData == TRUE)
			Rates->VarDataSite = 0;

	}
	else
	{
		Rates->NoOfRates = 0;
		Rates->Means = NULL;

		if(Opt->EstDelta == TRUE)
			Rates->NoOfRates++;

		if(Opt->EstKappa == TRUE)
			Rates->NoOfRates++;

		if(Opt->EstLambda == TRUE)
			Rates->NoOfRates++;

		if(Opt->NoOfRates > 0)
		{
			Rates->Rates = (double*)malloc(sizeof(double) * Rates->NoOfRates);
			if(Rates==NULL)
				MallocErr();
			for(Index=0;Index<Rates->NoOfRates;Index++)
				Rates->Rates[Index] = 1;
		}
	}

	Trees = Opt->Trees;

	Rates->NoEstData = 0;

	if(Opt->Model == CONTINUOUSREG)
	{
		for(Index=0;Index<Trees->NoOfTaxa;Index++)
		{
			Taxa = &Trees->Taxa[Index];
			if(Taxa->EstData == TRUE)
				Rates->NoEstData++;
		}

		if(Rates->NoEstData > 0)
		{
			Rates->EstData = (double*)malloc(sizeof(double) * Rates->NoEstData);
			if(Rates->EstData == NULL)
				MallocErr();

			for(Index=0;Index<Rates->NoEstData;Index++)
				Rates->EstData[Index] = 0;
		}
		else
			Rates->EstData = NULL;	
	}

	if(Opt->UseModelFile == TRUE)
	{
		Rates->FixedModels = LoadModelFile(Opt->ModelFile, &Rates->NoOfModels, Rates);
		SetFixedModel(Rates);		
	}
}

RATES*	CreatRates(OPTIONS *Opt)
{
	RATES*	Ret=NULL;
	int		Index;

	Ret = (RATES*)malloc(sizeof(RATES));
	if(Ret==NULL)
		MallocErr();

	Ret->NoOfFullRates	= Opt->NoOfRates;

	if(Opt->UseRModel == TRUE)
		Ret->NoOfFullRates = 1;

	Ret->NoOfRates		= FindNoOfRates(Opt);

	Ret->TreeNo			= 0;
	Ret->Prios			= NULL;
	Ret->Root			= NULL;
	Ret->Rates			= NULL;
	Ret->Pis			= NULL;
	Ret->FullRates		= NULL;
	Ret->Means			= NULL;
	Ret->Beta			= NULL;
	Ret->MappingVect	= NULL;
	Ret->LhPrior		= 0;
	Ret->LogHRatio		= 0;
	Ret->LogJacobion	= 0;

	Ret->Gamma			= -1;
	Ret->GammaCats		= 1;
	Ret->GammaMults		= NULL;
	Ret->LastGamma		= -1;
	Ret->GammaPrior		= NULL;
	Ret->HMeanCount		= 0;
	Ret->HMeanSum		= 0.0;
	Ret->NoEstData		= 0;
	Ret->EstData		= NULL;
	Ret->NoOfModels		= 0;
	Ret->FixedModels	= NULL;
	Ret->VarDataSite	= -1;

	if(Opt->UseGamma == TRUE)
	{
		Ret->GammaMults= (double*)malloc(sizeof(double) * Opt->GammaCats);
		if(Ret->GammaMults == NULL)
			MallocErr();
		Ret->GammaCats = Opt->GammaCats;
		
		if(Opt->EstGamma == FALSE)
			Ret->Gamma = Opt->FixGamma;
		else
			Ret->Gamma = 1;
	}

	if(Opt->UseKappa == FALSE)
		Ret->Kappa = -1;
	else
	{
		if(Opt->EstKappa == FALSE)
			Ret->Kappa = Opt->FixKappa;
		else
			Ret->Kappa = 1;
	}
	

	if(Opt->DataType == CONTINUOUS)
	{
		CreatCRates(Opt, Ret);
		return Ret;
	}

	if(Opt->Model == DESCINDEP)
		if(Opt->UseCovarion == FALSE)
			Ret->Root = (double*)malloc(sizeof(double)*4);
		else
			Ret->Root = (double*)malloc(sizeof(double)*8);
	
	if(Ret->NoOfRates > 0)
	{
		Ret->Rates = (double*)malloc(sizeof(double)*Ret->NoOfRates);
		if(Ret->Rates == NULL)
			MallocErr();
	}

	Ret->FullRates = (double*)malloc(sizeof(double)*Ret->NoOfFullRates);
	if(Ret->FullRates == NULL)
		MallocErr();

	for(Index=0;Index<Ret->NoOfRates;Index++)
		Ret->Rates[Index] = 1;

	Ret->Pis = (double*)malloc(sizeof(double)*Opt->Trees->NoOfStates);

	if(Ret->Pis == NULL)
		MallocErr();

	SetPiValues(Ret, Opt);

	if(Opt->UseCovarion == TRUE)
	{
		Ret->OffToOn = 1;
		Ret->OnToOff = 1;
		Ret->CoVarPis[0] = 1;
		Ret->CoVarPis[1] = 1;
	}

	if(Opt->UseRJMCMC == TRUE)
	{
		Ret->MappingVect = (int*)malloc(sizeof(int) * Ret->NoOfFullRates);
		if(Ret->MappingVect == NULL)
			MallocErr();
		for(Index=0;Index<Ret->NoOfFullRates;Index++)
			Ret->MappingVect[Index] = Index;
	}
	else
		Ret->MappingVect = NULL;

	MapRates(Ret, Opt);

	return Ret;
}

void	PrintConRegVarCoVarHeadder(FILE* Str, int NoOfSites, int DepSiteNo)
{
	int	x,y;

	NoOfSites++;

	for(x=0;x<NoOfSites;x++)
	{
		if(x!=DepSiteNo)
			fprintf(Str, "Trait %d Var\t", x+1);
	}

	for(x=0;x<NoOfSites;x++)
		for(y=x+1;y<NoOfSites;y++)
			if((y != DepSiteNo) && (x != DepSiteNo))
				fprintf(Str, "Trait %d,%d CoVar\t", x+1,y+1);
	
}

void	PrintRatesHeadderCon(FILE* Str, OPTIONS *Opt)
{
	int		Index;
	int		x,y;
	TREES	*Trees;

	if(Opt->Analsis == ANALMCMC)
		fprintf(Str, "Iteration\t");

	fprintf(Str, "Tree No\t");
	fprintf(Str, "Lh\t");

	if(Opt->Analsis == ANALMCMC)
		fprintf(Str, "HMean\t");

	if((Opt->Model == CONTINUOUSDIR) || (Opt->Model == CONTINUOUSRR))
	{
		for(Index=0;Index<Opt->Trees->NoOfSites;Index++)
		{
			fprintf(Str, "Alpha Trait %d\t", Index+1);
		}
	}
	else
		fprintf(Str, "Alpha\t");

	if(Opt->Model == CONTINUOUSDIR)
	{
		for(Index=0;Index<Opt->Trees->NoOfSites;Index++)
			fprintf(Str, "Beta Trait %d\t", Index+1);
	}

	if(Opt->Model == CONTINUOUSREG)
	{
		for(Index=0;Index<Opt->Trees->NoOfSites+1;Index++)
		{
			if(Index != Opt->DependantSite)
				fprintf(Str, "Beta Trait %d\t", Index+1);
		}

		fprintf(Str, "Var\t");
		fprintf(Str, "R^2\t");
		fprintf(Str, "Error Ratio\t");

		PrintConRegVarCoVarHeadder(Str, Opt->Trees->NoOfSites, Opt->DependantSite);
	}
	else
	{
		for(Index=0;Index<Opt->Trees->NoOfSites;Index++)
			fprintf(Str, "Trait %d Var\t", Index+1);

		for(x=0;x<Opt->Trees->NoOfSites;x++)
			for(y=x+1;y<Opt->Trees->NoOfSites;y++)
				fprintf(Str, "Trait %d %d Co-Var\t", x+1, y+1);
	}		

	if(Opt->UseVarData == TRUE)
		fprintf(Str, "Var Data Site\t");

	if((Opt->EstKappa == TRUE) || (Opt->FixKappa != -1))
		fprintf(Str, "Kappa\t");

	if((Opt->EstDelta == TRUE) || (Opt->FixDelta != -1))
		fprintf(Str, "Delta\t");

	if((Opt->EstLambda == TRUE) || (Opt->FixLambda != -1))
		fprintf(Str, "Lambda\t");

	if(Opt->NodeBLData == TRUE)
	{
		fprintf(Str, "Slope Nodes\tSlope  Root to Tip\t");
		fprintf(Str, "Min Nodes\tMax Nodes\t");
	}

	Trees = Opt->Trees;
	for(Index=0;Index<Trees->NoOfTaxa;Index++)
		if(Trees->Taxa[Index].EstData == TRUE)
			fprintf(Str, "Est %s\t", Trees->Taxa[Index].Name);

	if(Opt->Analsis != ANALMCMC)
		fprintf(Str, "\n");
}

void	PrintRecNodeHeadder(FILE* Str, OPTIONS *Opt, char* Name, int SiteNo)
{
	int	Index;
	int	NOS;

	if(Opt->Model == DESCINDEP)
	{
/*		if(Opt->UseCovarion == TRUE)
		{
			fprintf(Str, "%s - Status T1 P(0)\t", Name);
			fprintf(Str, "%s - Status T1 P(1)\t", Name);
			fprintf(Str, "%s - Status T2 P(0)\t", Name);
			fprintf(Str, "%s - Status T2 P(1)\t", Name);
		}
*/
		fprintf(Str, "%s - T1 - P(0)\t", Name);
		fprintf(Str, "%s - T1 - P(1)\t", Name);
		fprintf(Str, "%s - T2 - P(0)\t", Name);
		fprintf(Str, "%s - T2 - P(1)\t", Name);

		return;
	}

	if(Opt->Model == DESCDEP)
	{
/*		if(Opt->UseCovarion == TRUE)
		{
			fprintf(Str, "%s - Status P(0,0)\t", Name);
			fprintf(Str, "%s - Status P(0,1)\t", Name);
			fprintf(Str, "%s - Status P(1,0)\t", Name);
			fprintf(Str, "%s - Status P(1,1)\t", Name);
		}
*/
		fprintf(Str, "%s - P(0,0)\t", Name);
		fprintf(Str, "%s - P(0,1)\t", Name);
		fprintf(Str, "%s - P(1,0)\t", Name);
		fprintf(Str, "%s - P(1,1)\t", Name);
		return;
	}

	if(Opt->UseCovarion == TRUE)
		NOS = (Opt->Trees->NoOfStates / 2);
	else
		NOS = Opt->Trees->NoOfStates;

	for(Index=0;Index<NOS;Index++)
	{
		if(SiteNo != -1)
			fprintf(Str, "%s - S(%d) - P(%c)\t", Name, SiteNo,Opt->Trees->SymbolList[Index]);
		else
			fprintf(Str, "%s P(%c)\t", Name, Opt->Trees->SymbolList[Index]);
	}
}

void	PrintRatesHeadder(FILE* Str, OPTIONS *Opt)
{
	int			Index;
	int			SiteIndex;
	RECNODE		RNode=NULL;
	
	if(Opt->DataType == CONTINUOUS)
	{
		PrintRatesHeadderCon(Str, Opt);
		return;
	}	

	if(Opt->Analsis == ANALMCMC)
		fprintf(Str, "Iteration\tLh\tHarmonic Mean\tTree No\t");
	else
		fprintf(Str, "Tree No\tLh\t");

	if(Opt->UseRJMCMC == TRUE)
	{
		fprintf(Str, "No Off Paramiters\t");
		fprintf(Str, "Model string\t");
		if(Opt->Model == DESCDEP)
			fprintf(Str, "Dep / InDep\t");
	}

	if(Opt->UseRModel == FALSE)
	{
		for(Index=0;Index<Opt->NoOfRates;Index++)
			fprintf(Str, "%s\t", Opt->RateName[Index]);
	}
	else
		fprintf(Str, "R Model\t");

	if(Opt->UseCovarion == TRUE)
		fprintf(Str, "Covar On to Off\t Covar Off to On\t");

	if(Opt->UseKappa == TRUE)
		fprintf(Str, "Kappa\t");

	if(Opt->UseGamma == TRUE)
		fprintf(Str, "Gamma\t");

	for(SiteIndex=0;SiteIndex<Opt->Trees->NoOfSites;SiteIndex++)
	{
		if(Opt->Trees->NoOfSites == 1)
			PrintRecNodeHeadder(Str, Opt, "Root", -1);
		else
			PrintRecNodeHeadder(Str, Opt, "Root", SiteIndex);
	}

	for(Index=0;Index<Opt->NoOfRecNodes;Index++)
	{
		RNode = Opt->RecNodeList[Index];

		for(SiteIndex=0;SiteIndex<Opt->Trees->NoOfSites;SiteIndex++)
		{
			if(Opt->Trees->NoOfSites == 1)
				PrintRecNodeHeadder(Str, Opt, RNode->Name, -1);
			else
				PrintRecNodeHeadder(Str, Opt, RNode->Name, SiteIndex);
		}
	}

	if(Opt->Analsis == ANALML)
		fprintf(Str, "\n");
}

double	TransVarCoVar(int N, double x)
{
	return (x * (double)N) / (double)(N-1);
}

double	FindMean(double *List, int WSize)
{
	double	Ret=0;
	int		Index;

	for(Index=0;Index<WSize;Index++)
		Ret += List[Index];

	return Ret / WSize;
}

double	FindCorrelation(double *ListX, double *ListY, int WSize)
{
	double	Top;
	double	BotX;
	double	BotY;
	double	XMean;
	double	YMean;
	int		Index;

	XMean = FindMean(ListX, WSize);
	YMean = FindMean(ListY, WSize);

	Top = 0;
	for(Index=0;Index<WSize;Index++)
		Top += (ListX[Index] - XMean) * (ListY[Index] - YMean);

	BotX = 0;
	BotY = 0;
	for(Index=0;Index<WSize;Index++)
	{
		BotX += (ListX[Index] - XMean) * (ListX[Index] - XMean);
		BotY += (ListY[Index] - YMean) * (ListY[Index] - YMean);
	}

	BotX = sqrt(BotX);
	BotY = sqrt(BotY);

	return Top / (BotX * BotY);
}

double	FindYEst(double Alpha, double *Beta, double *Sites, int NoSites)
{
	double	Ret;
	int		Index;

	Ret = Alpha;
	for(Index=0;Index<NoSites;Index++)
		Ret += (Beta[Index] * Sites[Index]);

	return Ret;
}

double	FindERatio(RATES* Rates, OPTIONS *Opt)
{
	TREES	*Trees;
	TREE	*Tree;
	TAXA	*Taxa;
	CONVAR	*CV;
	double	*Y;
	double	*YP;
	double	*TempV;
	int		Index;
	double	Ret;
	double	SSy;
	double	SSe;
	double	Mean;

	Trees	= Opt->Trees;
	Tree	= &Trees->Tree[Rates->TreeNo];
	CV		= Tree->ConVars;
	
	Y		= (double*)malloc(sizeof(double) * Trees->NoOfTaxa);
	YP		= (double*)malloc(sizeof(double) * Trees->NoOfTaxa);
	TempV	= (double*)malloc(sizeof(double) * Trees->NoOfTaxa);

	if((Y == NULL) || (YP == NULL) || (TempV == NULL))
		MallocErr();

	for(Index=0;Index<Trees->NoOfTaxa;Index++)
	{
		Taxa = &Trees->Taxa[Index];

		Y[Index] = Taxa->Dependant;
		YP[Index] = Taxa->Dependant - FindYEst(CV->Alpha[0], CV->Beta, Taxa->ConData, Trees->NoOfSites);
	}

	Mean = MLFindAlphaReg(Trees, Tree, Y);
	for(Index=0;Index<Trees->NoOfTaxa;Index++)
		Y[Index] -= Mean;

	VectByMatrixMult(Y, Tree->ConVars->InvV, TempV);
	SSy = VectByVectMult(Y, TempV, Trees->NoOfTaxa);

	VectByMatrixMult(YP, Tree->ConVars->InvV, TempV);
	SSe = VectByVectMult(YP, TempV, Trees->NoOfTaxa);

	Ret = (SSy - SSe) / SSy;

 	free(Y);
	free(YP);
	free(TempV);

	return Ret;
}


double	FindRSquared(RATES* Rates, OPTIONS *Opt)
{
	TREES	*Trees;
	TREE	*Tree;
	TAXA	*Taxa;
	CONVAR	*CV;
	double	*Y;
	double	*YP;
	double	*TempV;
	int		Index;
	double	Ret;
	double	MeanY;
	double	MeanYP;
	double	T,B1,B2;

	Trees	= Opt->Trees;
	Tree	= &Trees->Tree[Rates->TreeNo];
	CV		= Tree->ConVars;
	
	Y		= (double*)malloc(sizeof(double) * Trees->NoOfTaxa);
	YP		= (double*)malloc(sizeof(double) * Trees->NoOfTaxa);
	TempV	= (double*)malloc(sizeof(double) * Trees->NoOfTaxa);

	if((Y == NULL) || (YP == NULL) || (TempV == NULL))
		MallocErr();

	for(Index=0;Index<Trees->NoOfTaxa;Index++)
	{
		Taxa = &Trees->Taxa[Index];

		Y[Index] = Taxa->Dependant;
		YP[Index] = FindYEst(CV->Alpha[0], CV->Beta, Taxa->ConData, Trees->NoOfSites);
	}

	MeanY = MLFindAlphaReg(Trees, Tree, Y);
	MeanYP= MLFindAlphaReg(Trees, Tree, YP);
	for(Index=0;Index<Trees->NoOfTaxa;Index++)
	{
		Y[Index] -= MeanY;
		YP[Index] -= MeanYP;
	}


	VectByMatrixMult(Y, Tree->ConVars->InvV, TempV);
	T = VectByVectMult(YP, TempV, Trees->NoOfTaxa);
	T = T * T;

	VectByMatrixMult(Y, Tree->ConVars->InvV, TempV);
	B1 = VectByVectMult(Y, TempV, Trees->NoOfTaxa);

	VectByMatrixMult(YP, Tree->ConVars->InvV, TempV);
	B2 = VectByVectMult(YP, TempV, Trees->NoOfTaxa);

	Ret = T / (B1 * B2);
	
	
 	free(Y);
	free(YP);
	free(TempV);

	return Ret;
}


void	PrintRegVarCoVar(FILE* Str, RATES *Rates, OPTIONS *Opt)
{
	int		x,y;
	int		x1,y1;
	CONVAR	*CV;
	TREES	*Trees;
	double	SigScalar;
	int		NOS;

	NOS = Opt->Trees->NoOfSites+1;

	Trees = Opt->Trees;

	CV = Trees->Tree[Rates->TreeNo].ConVars;
	SigScalar = CV->Sigma->me[0][0];

/*	for(x=0;x<NoOfSites;x++)
		for(y=x+1;y<NoOfSites;y++)
			if((y != DepSiteNo) && (x != DepSiteNo))
				fprintf(Str, "Trait %d,%d CoVar\t", x+1,y+1);
*/

	for(x=0;x<NOS;x++)
		for(y=x+1;y<NOS;y++)
			if((y != Opt->DependantSite) && (x != Opt->DependantSite))
			{
				y1 = y;
				if(y1 >= Opt->DependantSite)
					y1--;

				x1 = x;
				if(x1 >= Opt->DependantSite)
					x1--;

				fprintf(Str, "%f\t", CV->InvXVX->me[x1][y1]);
			}
	
}

void	PrintRatesCon(FILE* Str, RATES* Rates, OPTIONS *Opt)
{
	int		Index;
	int		x,y;
	CONVAR	*ConVar;
	int		MinNodes, MaxNodes;
	TAXA	*Taxa;
	double	HMean;
	
	ConVar = Opt->Trees->Tree[Rates->TreeNo].ConVars;

	fprintf(Str, "%d\t", Rates->TreeNo+1);
	fprintf(Str, "%f\t", Rates->Lh);

	if(Opt->Analsis == ANALMCMC)
	{
		HMean = log(Rates->HMeanCount / Rates->HMeanSum);
		fprintf(Str, "%f\t", HMean);
	}

	if((Opt->Model == CONTINUOUSRR) || (Opt->Model == CONTINUOUSDIR))
	{
		for(Index=0;Index<Opt->Trees->NoOfSites;Index++)
		{
			fprintf(Str, "%f\t", ConVar->Alpha[Index]);
			if(Opt->Model == CONTINUOUSDIR)
				fprintf(Str, "%f\t", ConVar->Beta[Index]);
		}	
			for(Index=0;Index<Opt->Trees->NoOfSites;Index++)
				fprintf(Str, "%f\t", TransVarCoVar(Opt->Trees->NoOfTaxa, ConVar->Sigma->me[Index][Index]));
		
			for(x=0;x<Opt->Trees->NoOfSites;x++)
				for(y=x+1;y<Opt->Trees->NoOfSites;y++)
					fprintf(Str, "%f\t", TransVarCoVar(Opt->Trees->NoOfTaxa, ConVar->Sigma->me[x][y]));
	}
	else
	{
		fprintf(Str, "%f\t", ConVar->Alpha[0]);
		for(Index=0;Index<Opt->Trees->NoOfSites;Index++)
			fprintf(Str, "%f\t", ConVar->Beta[Index]);

		fprintf(Str, "%f\t", ConVar->Sigma->me[0][0]);
		fprintf(Str, "%f\t", FindRSquared(Rates, Opt));
		fprintf(Str, "%f\t", FindERatio(Rates, Opt));

		PrintRegVarCoVar(Str, Rates, Opt);
	}

	if(Opt->UseVarData == TRUE)
		fprintf(Str, "%d\t", Rates->VarDataSite);

	if(Opt->EstKappa == TRUE)
		fprintf(Str, "%f\t", Rates->Kappa);

	if(Opt->FixKappa != -1)
		fprintf(Str, "%f\t", Opt->FixKappa);

	if(Opt->EstDelta == TRUE)
		fprintf(Str, "%f\t", Rates->Delta);

	if(Opt->FixDelta != -1)
		fprintf(Str, "%f\t", Opt->FixDelta);

	if(Opt->EstLambda == TRUE)
		fprintf(Str, "%f\t", Rates->Lambda);

	if(Opt->FixLambda != -1)
		fprintf(Str, "%f\t", Opt->FixLambda);

	if(Opt->NodeBLData == TRUE)
	{
		fprintf(Str, "%f\t", ConVar->Sigma->me[0][1] / ConVar->Sigma->me[0][0]);
		fprintf(Str, "%f\t", ConVar->Sigma->me[0][1] / ConVar->Sigma->me[1][1]);
		
		MinNodes = MaxNodes = (int)Opt->Trees->Taxa[0].ConData[0];
		for(Index=1;Index<Opt->Trees->NoOfTaxa;Index++)
		{
			Taxa = &Opt->Trees->Taxa[Index];
			if((int)Taxa->ConData[0] > MaxNodes)
				MaxNodes = (int)Taxa->ConData[0];

			if((int)Taxa->ConData[0] < MinNodes)
				MinNodes = (int)Taxa->ConData[0];
		}

		fprintf(Str, "%d\t%d\t", MinNodes, MaxNodes);
	}

	for(Index=0;Index<Rates->NoEstData;Index++)
		fprintf(Str, "%f\t", Rates->EstData[Index]);
}

void	PrintNodeRec(FILE *Str, NODE Node, int NOS, int NoOfSites, RATES* Rates, OPTIONS *Opt)
{
	double	Tot=0;
	int		Index;
	int		SiteIndex;
	int		TrueNOS;

	if(Node == NULL)
	{
		for(SiteIndex=0;SiteIndex<NoOfSites;SiteIndex++)
			for(Index=0;Index<NOS;Index++)
				fprintf(Str, "--\t");
		return;
	}

	for(SiteIndex=0;SiteIndex<NoOfSites;SiteIndex++)
	{
		Tot=0;

		for(Index=0;Index<NOS;Index++)
			Tot += (Node->Partial[SiteIndex][Index] * Rates->Pis[Index]);

		if(Opt->Model == DESCINDEP)
		{
			if(Opt->UseCovarion == FALSE)
			{
				fprintf(Str, "%f\t", (Node->Partial[0][0] + Node->Partial[0][1])/Tot);
				fprintf(Str, "%f\t", (Node->Partial[0][2] + Node->Partial[0][3])/Tot);
				fprintf(Str, "%f\t", (Node->Partial[0][0] + Node->Partial[0][2])/Tot);
				fprintf(Str, "%f\t", (Node->Partial[0][1] + Node->Partial[0][3])/Tot);
			}
			else
			{
				fprintf(Str, "%f\t", (Node->Partial[0][0] + Node->Partial[0][1] + Node->Partial[0][4] + Node->Partial[0][5])/Tot);
				fprintf(Str, "%f\t", (Node->Partial[0][2] + Node->Partial[0][3] + Node->Partial[0][6] + Node->Partial[0][7])/Tot);
				fprintf(Str, "%f\t", (Node->Partial[0][0] + Node->Partial[0][2] + Node->Partial[0][4] + Node->Partial[0][6])/Tot);
				fprintf(Str, "%f\t", (Node->Partial[0][1] + Node->Partial[0][3] + Node->Partial[0][5] + Node->Partial[0][7])/Tot);
			}
		}
		if(Opt->Model == DESCDEP)
		{
			if(Opt->UseCovarion == FALSE)
			{
				fprintf(Str, "%f\t", (Node->Partial[0][0])/Tot);
				fprintf(Str, "%f\t", (Node->Partial[0][1])/Tot);
				fprintf(Str, "%f\t", (Node->Partial[0][2])/Tot);
				fprintf(Str, "%f\t", (Node->Partial[0][3])/Tot);
			}
			else
			{
				fprintf(Str, "%f\t", (Node->Partial[0][0] + Node->Partial[0][4])/Tot);
				fprintf(Str, "%f\t", (Node->Partial[0][1] + Node->Partial[0][5])/Tot);
				fprintf(Str, "%f\t", (Node->Partial[0][2] + Node->Partial[0][6])/Tot);
				fprintf(Str, "%f\t", (Node->Partial[0][3] + Node->Partial[0][7])/Tot);
			}
		}

		if(Opt->Model == MULTISTATE)
		{
			if(Opt->UseCovarion == FALSE)
			{
				for(Index=0;Index<NOS;Index++)
					fprintf(Str, "%f\t", (Node->Partial[SiteIndex][Index] *  Rates->Pis[Index]) / Tot);
			}
			else
			{
				TrueNOS = NOS / 2;
				for(Index=0;Index<TrueNOS;Index++)
					fprintf(Str, "%f\t", ((Node->Partial[SiteIndex][Index] *  Rates->Pis[Index]) +
										  (Node->Partial[SiteIndex][Index+TrueNOS] *  Rates->Pis[Index+TrueNOS]))/ Tot);
			}

		}
	}
}

char	RJModelType(int *ModelStr)
{

	if(	(ModelStr[0] == ModelStr[5]) &&
		(ModelStr[1] == ModelStr[3]) &&
		(ModelStr[2] == ModelStr[7]) &&
		(ModelStr[4] == ModelStr[6])
		)
		return 'I';
	else
		return 'D';
}

void	PrintRates(FILE* Str, RATES* Rates, OPTIONS *Opt)
{
	int		Index;
	double	HMean;

	if(Opt->DataType == CONTINUOUS)
	{
		PrintRatesCon(Str, Rates, Opt);
		return;
	}

	if(Opt->Analsis == ANALMCMC)
	{
		HMean = log(Rates->HMeanCount / Rates->HMeanSum);
	/*	fprintf(Str, "%f\t%f\t%f\t%d\t", Rates->Lh, Rates->Lh + Rates->LhPrior, HMean, Rates->TreeNo+1); */
		fprintf(Str, "%f\t%f\t%d\t", Rates->Lh, HMean, Rates->TreeNo+1);
	}
	else
		fprintf(Str, "%d\t%f\t", Rates->TreeNo+1,Rates->Lh);

	if(Opt->UseRJMCMC == TRUE)
	{
		fprintf(Str, "%d\t'", NoOfPramGroups(Rates, NULL, NULL));
		for(Index=0;Index<Rates->NoOfFullRates;Index++)
		{
			if(Rates->MappingVect[Index] == ZERORATENO)
				fprintf(Str, "Z");
			else
			{
				if(Rates->MappingVect[Index] <= 9)
					fprintf(Str, "%d", Rates->MappingVect[Index]);
				else
					fprintf(Str, "%c", Rates->MappingVect[Index] + 'A');
			}
		}
		fprintf(Str, "\t");

		if(Opt->Model == DESCDEP)
			fprintf(Str, "%c\t", RJModelType(Rates->MappingVect));
	}

	if(Opt->UseRModel == TRUE)
		fprintf(Str, "%f\t", Rates->FullRates[0]);
	else
	{
		for(Index=0;Index<Opt->NoOfRates;Index++)
			fprintf(Str, "%f\t", Rates->FullRates[Index]);
	}

	if(Opt->UseCovarion == TRUE)
		fprintf(Str, "%f\t%f\t", Rates->OffToOn, Rates->OnToOff);

	if(Opt->UseKappa == TRUE)
		fprintf(Str, "%f\t", Rates->Kappa);

	if(Opt->UseGamma == TRUE)
		fprintf(Str, "%f\t", Rates->Gamma);

	PrintNodeRec(Str, Opt->Trees->Tree[Rates->TreeNo].Root, Opt->Trees->NoOfStates, Opt->Trees->NoOfSites, Rates, Opt);

	for(Index=0;Index<Opt->NoOfRecNodes;Index++)
		PrintNodeRec(Str, Opt->RecNodeList[Index]->TreeNodes[Rates->TreeNo], Opt->Trees->NoOfStates, Opt->Trees->NoOfSites, Rates, Opt);
}

void	CopyRJRtaes(RATES *A, RATES *B, OPTIONS *Opt)
{
	if(A->NoOfRates != B->NoOfRates)
	{
		free(A->Rates);
		A->Rates = (double*)malloc(sizeof(double) * B->NoOfRates);
		if(A->Rates == NULL)
			MallocErr();
		A->NoOfRates = B->NoOfRates;
	}

	memcpy(A->Rates, B->Rates, sizeof(double)*B->NoOfRates);
	memcpy(A->MappingVect, B->MappingVect, sizeof(int)*B->NoOfFullRates);
}

void	CopyRates(RATES *A, RATES *B, OPTIONS *Opt)
{
	int	Index;

	A->Delta = B->Delta;
	A->Kappa = B->Kappa;
	A->Lambda= B->Lambda;

	A->HMeanCount	= B->HMeanCount;
	A->HMeanSum		= B->HMeanSum;

	if(Opt->UseRJMCMC == FALSE)
	{
		if(A->Rates!=NULL)
		{
			for(Index=0;Index<A->NoOfRates;Index++)
				A->Rates[Index] = B->Rates[Index];
		}
	}
	else
		CopyRJRtaes(A, B, Opt);

	A->TreeNo = B->TreeNo;
	A->Lh	  = B->Lh;
	A->LhPrior= B->LhPrior;

	if(Opt->UseCovarion==TRUE)
	{
		A->OffToOn	= B->OffToOn;
		A->OnToOff	= B->OnToOff;
	}

	if(Opt->Analsis == ANALMCMC)
	{
		A->NoOfPriors = B->NoOfPriors;
		CopyRatePriors(A->Prios, B->Prios, B->NoOfPriors);
	}

	if(Opt->UseGamma == TRUE)
	{
		A->Gamma	= B->Gamma;
		A->GammaCats= B->GammaCats;
		A->LastGamma= B->LastGamma;
		memcpy(A->GammaMults, B->GammaMults, sizeof(double) * A->GammaCats);

		CopyPrior(A->GammaPrior, B->GammaPrior);
	}

	if(B->NoEstData > 0)
	{
		A->NoEstData = B->NoEstData;
		for(Index=0;Index<B->NoEstData;Index++)
			A->EstData[Index] = B->EstData[Index];
	}

	A->VarDataSite = B->VarDataSite;
} 


/* MrBays Rare Changes */
double ChangeRates(double Rate, double dev)
{
	int		Exit;
	double	Ret;

	if(Rate >= MAXRATE)
		return  MAXRATE;
	do
	{
		Exit = TRUE;
		Ret = (GenRand() * dev) - (dev / 2.0); 
		Ret += Rate;

		if(Ret > MAXRATE)
			Exit = FALSE;

		if(Ret < MINRATE)
			Exit = FALSE;

	} while(Exit == FALSE);

	return Ret;
}

double	MultePram(double Val, double Min, double Max, double Dev)
{
	double	Ret;
	int		Exit;

	do
	{
		Exit = TRUE;

		Ret = (GenRand() * Dev) - (Dev / 2.0); 
		Ret += Val;

		if(Ret > Max)
			Exit = FALSE;

		if(Ret < Min)
			Exit = FALSE;

	} while(Exit == FALSE);

	return Ret;		
}


void	MutateRatesOld(OPTIONS* Opt, RATES* Rates)
{
	int		Index;

	if(Rates->Rates != NULL)
		for(Index=0;Index<Rates->NoOfRates;Index++)
			Rates->Rates[Index] = ChangeRates(Rates->Rates[Index], Opt->RateDev);

	if(Opt->DataType == CONTINUOUS)
	{
		for(Index=0;Index<Opt->Trees->NoOfSites;Index++)
			Rates->Means[Index] += (GenRand() * Opt->RateDev) - (Opt->RateDev / 2.0);

		if(Opt->EstDelta == TRUE)
			Rates->Delta = MultePram(Rates->Delta, 0.000001, 3.0, Opt->RateDev);

		if(Opt->EstKappa == TRUE)
			Rates->Kappa = MultePram(Rates->Kappa, 0.000001, 3.0, Opt->RateDev);

		if(Opt->EstLambda == TRUE)
			Rates->Lambda = MultePram(Rates->Lambda, 0.000001, 1, Opt->RateDev);
	}
	else
		if(Opt->EstKappa == TRUE)
			Rates->Kappa = MultePram(Rates->Kappa, 0.000001, 3.0, Opt->RateDev);

	if(Opt->UseCovarion == TRUE)
	{
		Rates->OffToOn = ChangeRates(Rates->OffToOn, Opt->RateDev);
		Rates->OnToOff = ChangeRates(Rates->OnToOff, Opt->RateDev);
		
	}
	

	Rates->TreeNo = IntRand() % Opt->Trees->NoOfTrees;
}


int	PickACat(double *Vect, int Size)
{
	double	Val;
	int		Index;

	Val = GenRand();
	for(Index=0;Index<Size;Index++)
		if(Val<Vect[Index])
			return Index;

	return -1;
}


void	MutateRates(OPTIONS* Opt, RATES* Rates, SCHEDULE*	Shed)
{
	int		Index;
	int		NoOfGroups;

	Shed->Op = PickACat(Shed->OptFreq, Shed->NoOfOpts);

	switch(Shed->Op)
	{
		case(SRATES):

			if(Opt->UseModelFile == TRUE)
			{
				SetFixedModel(Rates);
				break;
			}
			if(Opt->DataType == DISCRETE)
			{
				for(Index=0;Index<Rates->NoOfRates;Index++)
					Rates->Rates[Index] = ChangeRates(Rates->Rates[Index], Opt->RateDev);
			}
			else
			{
				for(Index=0;Index<Rates->NoOfRates;Index++)
					Rates->Rates[Index] += (GenRand() * Opt->RateDev) - (Opt->RateDev / 2.0);

				if(Opt->AlphaZero == TRUE)
				{
					for(Index=0;Index<Opt->Trees->NoOfSites;Index++)
						Rates->Rates[Index] = 0;
				}
			}
		break;

		case(SCV):
			Rates->OffToOn = ChangeRates(Rates->OffToOn, Opt->RateDev);
			Rates->OnToOff = ChangeRates(Rates->OnToOff, Opt->RateDev);

			Rates->OffToOn = Rates->OnToOff;
		break;

		case(SKAPPA):
			Rates->Kappa = MultePram(Rates->Kappa, 0.000001, 5.0, Opt->RateDev);
		break;

		case(SDELTA):
			Rates->Delta = MultePram(Rates->Delta, 0.000001, 3.0, Opt->RateDev);
		break;

		case(SLABDA):
			Rates->Lambda = MultePram(Rates->Lambda, 0.000001, 1.0, Opt->RateDev);
		break;

		case(SJUMP):
			NoOfGroups = NoOfPramGroups(Rates, NULL, NULL);

			if(GenRand() < 0.75)
			{
				switch(NoOfGroups)
				{
					case 1:
						RJSplit(Rates, Opt);
					break;

					default:
						if(NoOfGroups == Rates->NoOfFullRates)
						{
							RJMerge(Rates, Opt);
							break;
						}

						if(GenRand() < 0.5)
							RJSplit(Rates, Opt);
						else
							RJMerge(Rates, Opt);
					break;
				}
			}
			else
			{
				if(GenRand() < 0.5)
					RJAugment(Rates, Opt);
				else
					RJReduce(Rates, Opt);
			}
		break;

		case(SPPROR):
			MutatePriorsNormal(Rates->Prios, Rates->NoOfPriors, Opt->HPDev);
		break;

		case(SESTDATA):
			for(Index=0;Index<Rates->NoEstData;Index++)
				Rates->EstData[Index] += (GenRand() * Opt->EstDataDev) - (Opt->EstDataDev / 2.0);
		break;

		case(SVARDATA):
			Rates->VarDataSite = rand() % Opt->VarData->NoPoints;
		break;

	}

	Shed->Tryed[Shed->Op]++;

	Rates->TreeNo = IntRand() % Opt->Trees->NoOfTrees;
}

void	FreeRates(RATES *Rates)
{
	if(Rates->Rates != NULL)
		free(Rates->Rates);

	if(Rates->Root != NULL)
		free(Rates->Root);

	if(Rates->Pis != NULL)
		free(Rates->Pis);

	if(Rates->FullRates != NULL)
		free(Rates->FullRates);

	if(Rates->Means != NULL)
		free(Rates->Means);

	if(Rates->Beta	!= NULL)
		free(Rates->Beta);

	if(Rates->MappingVect != NULL)
		free(Rates->MappingVect);

	if(Rates->GammaMults != NULL)
		free(Rates->GammaMults);

	free(Rates);
}

void	UpDataSummaryNo(SUMMARYNO *SumNo, double No)
{
	SumNo->Sum		+= No;
	SumNo->SumSqrs	+= No * No;
	SumNo->N		+= 1;
}

double	GetSummaryAve(SUMMARYNO *SumNo)
{
	return (SumNo->Sum / SumNo->N);
}

double	GetSummaryVar(SUMMARYNO *SumNo)
{
	double	Ret;

	Ret = (SumNo->Sum * SumNo->Sum) / (double)SumNo->N;
	Ret = SumNo->SumSqrs - Ret;

	return Ret / (SumNo->N - 1);	
}

void	UpDataSummary(SUMMARY *Summary, RATES* Rates, OPTIONS *Opt)
{
	int		Index;
	double	Pct;
	double	*RootP=NULL;

	RootP = Opt->Trees->Tree[Rates->TreeNo].Root->Partial[0];

	UpDataSummaryNo(&Summary->Lh, Rates->Lh);

	for(Index=0;Index<Opt->NoOfRates;Index++)
		UpDataSummaryNo(&Summary->Rates[Index], Rates->Rates[Index]);

	for(Index=0;Index<Opt->Trees->NoOfStates;Index++)
	{
		Pct = GetStateProbPct(Index, Opt->Trees->NoOfStates, RootP);
		UpDataSummaryNo(&Summary->Root[Index], Pct);
	}
}

void	InitSummaryNo(SUMMARYNO *SumNo)
{
	SumNo->N		=	0;
	SumNo->Sum		=	0;
	SumNo->SumSqrs	=	0;
}

void	FreeSummary(SUMMARY*	Summary)
{
	free(Summary->Rates);
	free(Summary->Root);
}

SUMMARY*	CreatSummary(OPTIONS *Opt)
{
	SUMMARY* Ret=NULL;
	int		 Index;

	Ret = (SUMMARY*)malloc(sizeof(SUMMARY));
	if(Ret==NULL)
		MallocErr();

	InitSummaryNo(&Ret->Lh);
	
	Ret->Rates = (SUMMARYNO*)malloc(sizeof(SUMMARYNO)*Opt->NoOfRates);
	if(Ret->Rates == NULL)
		MallocErr();
	for(Index=0;Index<Opt->NoOfRates;Index++)
		InitSummaryNo(&Ret->Rates[Index]);


	Ret->Root = (SUMMARYNO*)malloc(sizeof(SUMMARYNO)*Opt->Trees->NoOfStates);
	if(Ret->Root== NULL)
		MallocErr();
	for(Index=0;Index<Opt->Trees->NoOfStates;Index++)
		InitSummaryNo(&Ret->Root[Index]);

	return Ret;
}

void	PrintSummaryHeadder(FILE* Str, SUMMARY	*Summary, OPTIONS *Opt)
{
	int	Index;

	fprintf(Str, "Tree File\tData File\tNOS\t");

	fprintf(Str, "Lh - Ave\tLh - Var\t");


	fprintf(Str, "%s - Ave\t%s - Var\t", Opt->RateName[0], Opt->RateName[0]);
	fprintf(Str, "%s - Ave\t%s - Var\t", Opt->RateName[1], Opt->RateName[1]);

	
	for(Index=0;Index<Opt->Trees->NoOfStates;Index++)
		fprintf(Str, "%d - Ave\t%d - Var\t", Index, Index);

	fprintf(Str, "\n");
}

void	PrintSummary(FILE* Str, SUMMARY	*Summary, OPTIONS *Opt)
{
	int	Index;

	fprintf(Str, "%s\t%s\t%d\t", Opt->TreeFN, Opt->DataFN, Opt->Trees->NoOfStates);

	fprintf(Str, "%f\t%f\t", GetSummaryAve(&Summary->Lh), GetSummaryVar(&Summary->Lh));

	fprintf(Str, "%f\t%f\t", GetSummaryAve(&Summary->Rates[0]), GetSummaryVar(&Summary->Rates[0]));
	fprintf(Str, "%f\t%f\t", GetSummaryAve(&Summary->Rates[1]), GetSummaryVar(&Summary->Rates[1]));

	for(Index=0;Index<Opt->Trees->NoOfStates;Index++)
		fprintf(Str, "%f\t%f\t", GetSummaryAve(&Summary->Root[Index]), GetSummaryVar(&Summary->Root[Index]));

	fprintf(Str, "\n");
}

void	BlankSchedule(SCHEDULE*	Shed)
{
	int	Index;

	for(Index=0;Index<Shed->NoOfOpts;Index++)
	{
		Shed->Accepted[Index] = 0;
		Shed->Tryed[Index] = 0;
	}
}

void		ScaleVect(double *Vect, int VectSize)
{
	double	SF=0;
	int		Index=0;

	for(Index=0;Index<VectSize;Index++)
		SF += Vect[Index];

	SF = 1 / SF;

	for(Index=0;Index<VectSize;Index++)
		Vect[Index] *= SF;

	SF = 0;
	for(Index=0;Index<VectSize;Index++)
	{
		Vect[Index] = Vect[Index] + SF;
		SF = Vect[Index];
	}
}

int		UsingHP(OPTIONS *Opt)
{
	int Index;

	for(Index=0;Index<Opt->NoOfRates;Index++)
	{
		if(Opt->Priors[Index]->UseHP == TRUE)
			return TRUE;
	}

	if(Opt->RJPrior->UseHP == TRUE)
		return TRUE;

	return FALSE;
}

void	SetSchedule(SCHEDULE*	Shed, OPTIONS *Opt)
{
	double	Left = 1.0;
	int		Index;
	int			Rates;


	for(Index=0;Index<Shed->NoOfOpts;Index++)
		Shed->OptFreq[Index] = 0.0;

	if(Opt->UseCovarion == TRUE)
	{
		Shed->OptFreq[1] = 0.2;
		Left = Left - Shed->OptFreq[1];
	}

	if(Opt->EstKappa == TRUE)
	{
		Shed->OptFreq[2] = 0.1;
		Left = Left - Shed->OptFreq[2];
	}

	if(Opt->EstDelta == TRUE)
	{
		Shed->OptFreq[3] = 0.1;
		Left = Left - Shed->OptFreq[3];
	}

	if(Opt->EstLambda == TRUE)
	{
		Shed->OptFreq[4] = 0.1;
		Left = Left - Shed->OptFreq[4];
	}

	if(Opt->UseRJMCMC == TRUE)
	{
		Shed->OptFreq[5] = 0.1;
		Left = Left - Shed->OptFreq[5];
	}

	if(UsingHP(Opt) == TRUE)
	{
		Shed->OptFreq[6] = 0.1;
		Left = Left - Shed->OptFreq[6];
	}

	if(EstData(Opt->Trees) == TRUE)
	{
		Shed->OptFreq[7] = 0.3;
		Left = Left - Shed->OptFreq[7];
	}

	if(Opt->UseVarData == TRUE)
	{
		Shed->OptFreq[8] = 0.1;
		Left = Left - Shed->OptFreq[8];
	}
	

	Rates = 0;
	if(Opt->DataType == CONTINUOUS)
		Rates = Opt->Trees->NoOfSites;
	else
		for(Index=0;Index<Opt->NoOfRates;Index++)
			if(Opt->ResTypes[Index] == RESNONE)
				Rates++;

	if(Rates == 0)
		Shed->OptFreq[0] = 0;
	else
		Shed->OptFreq[0] = Left;

	if(Opt->UseModelFile == TRUE)
		Shed->OptFreq[0] = 0.3;
	
	ScaleVect(Shed->OptFreq, Shed->NoOfOpts);	
}

SCHEDULE*	CreatSchedule(OPTIONS *Opt)
{
	SCHEDULE*	Ret=NULL;	
	
	Ret = (SCHEDULE*)malloc(sizeof(SCHEDULE));
	if(Ret==NULL)
		MallocErr();

	Ret->NoOfOpts = NOOFOPERATORS;
	BlankSchedule(Ret);

	SetSchedule(Ret, Opt);

	return Ret;
}

void	PrintShedHeadder(OPTIONS* Opt, SCHEDULE* Shed, FILE* Str)
{
	int	Index;
	double	Last=0;
	
	for(Index=0;Index<Shed->NoOfOpts;Index++)
	{
		fprintf(Str, "%s\t%2.2f%%\n", SHEDOP[Index], (Shed->OptFreq[Index]-Last)*100);
		Last = Shed->OptFreq[Index];
	}
	
	
	for(Index=0;Index<Shed->NoOfOpts;Index++)
		fprintf(Str, "Pct %s taken\t(Tried)\t", SHEDOP[Index]);
	fprintf(Str, "\n");	

	fflush(Str);

}

void	PrintShed(OPTIONS* Opt, SCHEDULE* Shed, FILE* Str)
{
	int	Index;
	
	for(Index=0;Index<Shed->NoOfOpts;Index++)
		if(Shed->Tryed[Index] != 0)
			fprintf(Str, "%f\t%d\t", (double)Shed->Accepted[Index] / (double)Shed->Tryed[Index], Shed->Tryed[Index]);
		else
			fprintf(Str, "0.0\t0\t");
	fprintf(Str, "\n");
}

double**	LoadModelFile(char* FileName, int *NoModels, RATES* Rates)
{
	double	**Ret;
	char*	Buffer;
	char**	Passed;
	int		Tokes;
	int		Index;
	int		RIndex;
	TEXTFILE	*Data;
	
	Data = LoadTextFile(FileName, FALSE);

	Buffer = (char*)malloc(sizeof(char) * (Data->MaxLine + 1));
	Passed = (char**)malloc(sizeof(char*) * (Data->MaxLine + 1));
	if((Buffer == NULL) || (Passed == NULL))
		MallocErr();

	(*NoModels) = 0;
	for(Index=0;Index<Data->NoOfLines;Index++)
	{
		strcpy(Buffer, Data->Data[Index]);
		Tokes = MakeArgv(Buffer, Passed, Data->MaxLine);
		if(Tokes != Rates->NoOfRates)
			printf("Expecting %d but found %d values, reading file %s line %d", Rates->NoOfRates, Tokes, FileName, Index+1);
		else
			(*NoModels)++;
	}

	Ret = AllocMatMem(*NoModels, Rates->NoOfRates);

	(*NoModels) = 0;
	for(Index=0;Index<Data->NoOfLines;Index++)
	{
		strcpy(Buffer, Data->Data[Index]);
		Tokes = MakeArgv(Buffer, Passed, Data->MaxLine);
		if(Tokes == Rates->NoOfRates)
		{
			for(RIndex=0;RIndex<Rates->NoOfRates;RIndex++)
				Ret[(*NoModels)][RIndex] = atof(Passed[RIndex]);
			(*NoModels)++;
		}
	}

/*	printf("\nMy Models\n");

	for(Index=0;Index<*NoModels;Index++)
	{
		for(RIndex=0;RIndex<Rates->NoOfRates;RIndex++)
			printf("%f\t", Ret[Index][RIndex]);
		printf("\n");
	}
*/
	free(Buffer);
	free(Passed);
	FreeTextFile(Data);

	return Ret;
}

void	SetFixedModel(RATES *Rates)
{
	int	No;

	No = rand() % Rates->NoOfModels;

	memcpy(Rates->Rates, Rates->FixedModels[No], sizeof(double) * Rates->NoOfRates);
}