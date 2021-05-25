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



#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "ContrastRegOutput.h"
#include "ContrastRegCalcReg.h"
#include "genlib.h"
#include "typedef.h"
#include "contrasts.h"

void			PrintRegConPost(REG_CON_POST* RegConPost)
{
	int SIndex, CIndex;

	for(CIndex=0;CIndex<RegConPost->NoContrasts;CIndex++)
	{
		printf("%d\t%f\t%f\t|\t", CIndex, RegConPost->Predicated[CIndex], RegConPost->Residuals[CIndex]);

		for(SIndex=0;SIndex<RegConPost->NoSites;SIndex++)
			printf("%f\t", RegConPost->Contrasts[CIndex][SIndex]);
		//	printf("%f\t", RegConPost->ContrastsVar[CIndex][SIndex]);
		printf("\n");
	}
}

void			SetRegConPost(TREE *Tree, REG_CON_POST* RegConPost)
{
	int Index, Pos, CIndex, SIndex;
	NODE N;

	Pos = 0;
	for(Index=0;Index<Tree->NoNodes;Index++)
	{
		N = Tree->NodeList[Index];
		if(N->Tip == FALSE)
		{
			for(CIndex=0;CIndex<N->ConData->NoContrast;CIndex++)
			{
				for(SIndex=0;SIndex<RegConPost->NoSites;SIndex++)
				{
					RegConPost->Contrasts[Pos][SIndex] = N->ConData->Contrast[CIndex]->Cont[SIndex] / sqrt(N->ConData->Contrast[CIndex]->Var[SIndex]);
				//	RegConPost->Contrasts[Pos][SIndex] = N->ConData->Contrast[CIndex]->Cont[SIndex];
				//	RegConPost->ContrastsVar[Pos][SIndex] = N->ConData->Contrast[CIndex]->Cont[SIndex] / sqrt(N->ConData->Contrast[CIndex]->Var[SIndex]);
				}

				Pos++;
			}
		}
	}
}

void			SetRegConPostPreRes(double *Beta, REG_CON_POST* RegConPost)
{
	int Index, SIndex;
	double Pre;

	for(Index=0;Index<RegConPost->NoContrasts;Index++)
	{
		Pre = 0;

		for(SIndex=0;SIndex<RegConPost->NoSites-1;SIndex++)
			Pre += Beta[SIndex] * RegConPost->Contrasts[Index][SIndex+1];

		RegConPost->Predicated[Index] = Pre;
		RegConPost->Residuals[Index] = RegConPost->Contrasts[Index][0] - Pre;
	}
}

REG_CON_POST*	CreatRegConPost(TREES *Trees, RATES *Rates)
{
	REG_CON_POST	*Ret;
	TREE			*Tree;
	CONTRASTR		*CR; 
	int				Index;

	CR					=	Rates->Contrast;
	Tree				=	Trees->Tree[Rates->TreeNo];

	Ret					=	(REG_CON_POST*)SMalloc(sizeof(REG_CON_POST));
	
	Ret->NoContrasts	=	Tree->NoContrast;
	Ret->NoSites		=	Trees->NoSites;

	Ret->Predicated		= (double*)SMalloc(sizeof(double) * Ret->NoContrasts);
	Ret->Residuals		= (double*)SMalloc(sizeof(double) * Ret->NoContrasts);
	Ret->Contrasts		= (double**)SMalloc(sizeof(double*) * Ret->NoContrasts);


	Ret->SSE	= (double*)SMalloc(sizeof(double) * (Ret->NoSites-1));
	Ret->R2		= (double*)SMalloc(sizeof(double) * (Ret->NoSites-1));
	Ret->SSRes	= (double*)SMalloc(sizeof(double) * (Ret->NoSites-1));
	Ret->SSTotal= (double*)SMalloc(sizeof(double) * (Ret->NoSites-1));
	Ret->SE		= (double*)SMalloc(sizeof(double) * (Ret->NoSites-1));

	for(Index=0;Index<Ret->NoContrasts;Index++)
		Ret->Contrasts[Index] = (double*)SMalloc(sizeof(double) * Ret->NoSites);

	SetRegConPost(Tree, Ret);
	SetRegConPostPreRes(Rates->Contrast->RegBeta, Ret);
	
	return Ret;
}

void			FreeRegConPost(REG_CON_POST* RegConPost)
{
	int Index;

	free(RegConPost->Predicated);
	free(RegConPost->Residuals);
	
	for(Index=0;Index<RegConPost->NoContrasts;Index++)
		free(RegConPost->Contrasts[Index]);

	free(RegConPost->Contrasts);
	
	free(RegConPost->SSE);
	free(RegConPost->R2);
	free(RegConPost->SSRes);
	free(RegConPost->SSTotal);
	free(RegConPost->SE);

	free(RegConPost);
}

double	GetRegConPostSiteMean(int SiteNo, REG_CON_POST* RegConPost)
{
	double Ret;
	int Index;

	Ret = 0;
	for(Index=0;Index<RegConPost->NoContrasts;Index++)
		Ret += RegConPost->Contrasts[Index][SiteNo];

	Ret = Ret / RegConPost->NoContrasts;

	return Ret;
}

void	CalcSSERegConPost(REG_CON_POST* RegConPost)
{
	double Mean, SS;
	int Index, CIndex;

		
	for(Index=1;Index<RegConPost->NoSites;Index++)
	{
		Mean = GetRegConPostSiteMean(Index, RegConPost);
		SS = 0;
		for(CIndex=0;CIndex<RegConPost->NoContrasts;CIndex++)
			SS += (RegConPost->Contrasts[CIndex][Index] - Mean) * (RegConPost->Contrasts[CIndex][Index] - Mean);
		RegConPost->SSE[Index-1] = SS;
	}
}

void	SetSSData(MATRIX *M, int SiteNo, int NoSites, int NoCon, double **Data)
{
	
	int		Index, x, Pos;

	for(Index=0;Index<NoCon;Index++)
	{
		M->me[Index][0] = Data[Index][SiteNo];
		Pos = 1;
		for(x=1;x<NoSites;x++)
		{
			if(x != SiteNo)
			{
				M->me[Index][Pos] = Data[Index][x];
				Pos++;
			}
		}
	}
}

double	CaclSSSiteResid(MATRIX *M, double *Beta)
{
	double Ret, Pre;
	int Index, x;


	Ret = 0;
	for(Index=0;Index<M->NoOfRows;Index++)
	{
		Pre = 0;

		for(x=1;x<M->NoOfCols;x++)
			Pre += Beta[x-1] * M->me[Index][x];

		Ret += (Pre - M->me[Index][0]) * (Pre - M->me[Index][0]);
	}

	return Ret;
}

double	CaclSSSiteTotal(MATRIX *M)
{
	double Ret, Mean;
	int Index;

	Mean = 0;
	for(Index=0;Index<M->NoOfRows;Index++)
		Mean += M->me[Index][0];
	Mean = Mean / M->NoOfRows;

	Ret = 0;
	for(Index=0;Index<M->NoOfRows;Index++)
		Ret += (M->me[Index][0] - Mean) * (M->me[Index][0] - Mean);
	
	return Ret;
}

double	GetSiteR2(int SiteNo, REG_CON_POST* RegConPost)
{
	double Ret;

	Ret = RegConPost->SSRes[SiteNo] / RegConPost->SSTotal[SiteNo];

	Ret = 1.0 - Ret;
	return Ret;
}

double	CalcGlobalR2(REG_CON_POST* RegConPost)
{
	double Ret, Mean, SSR, SST;
	int	Index;

	SSR = 0;
	for(Index=0;Index<RegConPost->NoContrasts;Index++)
		SSR += RegConPost->Residuals[Index] * RegConPost->Residuals[Index];

	Mean = 0;
	for(Index=0;Index<RegConPost->NoContrasts;Index++)
		Mean += RegConPost->Contrasts[Index][0];
	Mean = Mean / RegConPost->NoContrasts;

	SST = 0;
	for(Index=0;Index<RegConPost->NoContrasts;Index++)
		SST += (Mean - RegConPost->Contrasts[Index][0]) * (Mean - RegConPost->Contrasts[Index][0]);


	Ret = 1.0 - (SSR / SST);
	return Ret;
}


double	GetSSTotal(REG_CON_POST* RegConPost, int SiteNo)
{
	int Index;
	double Ret, Mean;

	Mean = 0;
	for(Index=0;Index<RegConPost->NoContrasts;Index++)
		Mean += RegConPost->Contrasts[Index][SiteNo];
	Mean = Mean / RegConPost->NoContrasts;

	Ret = 0;
	for(Index=0;Index<RegConPost->NoContrasts;Index++)
		Ret += (RegConPost->Contrasts[Index][SiteNo] - Mean) * (RegConPost->Contrasts[Index][SiteNo] - Mean);

	return Ret;
}

double	CaclSSResidual(REG_CON_POST* RegConPost)
{
	int Index;
	double Ret;

	Ret = 0;

	for(Index=0;Index<RegConPost->NoContrasts;Index++)
		Ret += RegConPost->Residuals[Index] * RegConPost->Residuals[Index];
	

	return Ret;
}

void	CalcSSResiduals(REG_CON_POST* RegConPost, int TestCorrel)
{
	int SIndex;
	MATRIX *M;
	double	*Beta;
	double	T;
	
	if(RegConPost->NoSites == 2)
	{
	//	RegConPost->SSRes[0] = CaclSSResidual(RegConPost);
		RegConPost->SSRes[0] = 0;
		RegConPost->SSTotal[0] = GetSSTotal(RegConPost, 1);
		


		T = 1.0 - (RegConPost->SSRes[0] / RegConPost->SSTotal[0]);
		return;
	}

	for(SIndex=1;SIndex<RegConPost->NoSites;SIndex++)
	{
		M = AllocMatrix(RegConPost->NoContrasts, RegConPost->NoSites - 1);

//		SetSSData(M, SIndex, RegConPost->NoSites, RegConPost->NoContrasts, RegConPost->ContrastsVar);
		SetSSData(M, SIndex, RegConPost->NoSites, RegConPost->NoContrasts, RegConPost->Contrasts);
		Beta = ContrastMultiReg(M, TestCorrel);
	
		SetSSData(M, SIndex, RegConPost->NoSites, RegConPost->NoContrasts, RegConPost->Contrasts);

		RegConPost->SSRes[SIndex-1] = CaclSSSiteResid(M, Beta);
		RegConPost->SSTotal[SIndex-1] = CaclSSSiteTotal(M);
		
		free(Beta);
		FreeMatrix(M);
	}


}

void	CalcMSE(REG_CON_POST*	RegConPost)
{
	int Index;
	double Ret;

	Ret = 0;
	for(Index=0;Index<RegConPost->NoContrasts;Index++)
	{
		Ret += RegConPost->Residuals[Index] * RegConPost->Residuals[Index];
	}

	Ret = Ret / (RegConPost->NoContrasts - RegConPost->NoSites);

	RegConPost->MSE = Ret;
}

void	CalcSEMultiReg(REG_CON_POST*	RegConPost)
{
	int Index;
	double SE, R2, SSE, MSE;

	MSE = RegConPost->MSE;

	if(RegConPost->NoSites == 2)
	{
		RegConPost->SE[0] = sqrt(MSE) / sqrt(RegConPost->SSTotal[0]);
//		R2 = CalcGlobalR2(RegConPost);
		return;
	}

	
	for(Index=0;Index<RegConPost->NoSites-1;Index++)
	{
		R2 = GetSiteR2(Index, RegConPost);
		SSE = RegConPost->SSTotal[Index];

		SE = SSE * (1.0 - R2);
		SE = sqrt(MSE / SE);
		RegConPost->SE[Index] = SE;
	}
}


void	ProcRegConPost(REG_CON_POST *RegConPost, int TestCorrel)
{
	CalcMSE(RegConPost);

	CalcSSERegConPost(RegConPost);

	CalcSSResiduals(RegConPost, TestCorrel);

	CalcSEMultiReg(RegConPost);
	
	RegConPost->GR2 = CalcGlobalR2(RegConPost);
}


void	OutputConReg(FILE *Str, OPTIONS *Opt, TREES *Trees, RATES *Rates)
{
	int Index;
	REG_CON_POST*	RegConPost;

	if(Opt->RJDummy == TRUE)
		fprintf(Str, "%d\t", Rates->RJDummy->NoDummyCode);

	fprintf(Str, "%0.12f\t", Rates->Contrast->RegAlpha);

	for(Index=0;Index<Trees->NoSites-1;Index++)
		fprintf(Str, "%0.12f\t", Rates->Contrast->RegBeta[Index]);

	RegConPost = CreatRegConPost(Trees, Rates);
	ProcRegConPost(RegConPost, Opt->TestCorrel);
	
	fprintf(Str, "%0.12f\t%0.12f\t", Rates->Contrast->GlobalVar, RegConPost->GR2);

	fprintf(Str, "N\\A\t");

	for(Index=0;Index<Trees->NoSites-1;Index++)
		fprintf(Str, "%0.12f\t", RegConPost->SE[Index]);

//	PrintRegConPost(RegConPost);
//	exit(0);


	FreeRegConPost(RegConPost);
}