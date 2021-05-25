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



#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#include "TypeDef.h"
#include "GenLib.h"
#include "Options.h"
#include "Data.h"
#include "Trees.h"
#include "Priors.h"
#include "RandLib.h"
#include "Threaded.h"
#include "Part.h"
#include "Rates.h"
#include "Stones.h"
#include "RJLocalScalar.h"
#include "Tag.h"
#include "VarRates.h"
#include "LocalTransform.h"
#include "DistData.h"
#include "Schedule.h"
#include "TimeSlices.h"
#include "NLOptBT.h"
#include "Pattern.h"

#define	RATEOUTPUTLEN	33
#define	RIGHT_INDENT	4

void	FreeRecNodes(OPTIONS *Opt, int NoSites);
void	SetLocalTransformPrior(OPTIONS *Opt, TRANSFORM_TYPE	Type);
void	AddConAnsStatePrior(OPTIONS *Opt, int SiteNo);

char*	FormatStr(char* RateName)
{
	char* Ret, *Temp;
	int	Index;
	
	Temp = (char*)SMalloc(sizeof(char) * BUFFERSIZE);

	for(Index=0;Index<RATEOUTPUTLEN;Index++)
		Temp[Index] = ' ';
	Temp[RATEOUTPUTLEN] = '\0';

	sprintf(&Temp[RIGHT_INDENT], "%s", RateName);
	Temp[RIGHT_INDENT+strlen(RateName)] = ' ';
	
	Ret = StrMake(Temp);

	free(Temp);

	return Ret;
}

void	PrintOptRes(FILE* Str, OPTIONS *Opt)
{
	int		Index;
	char*	FRateName;

	if((Opt->AnalyticalP == TRUE) || (Opt->UseRModel == TRUE) || (Opt->NOSPerSite == TRUE))
		return;

	fprintf(Str, "Restrictions:\n");
	for(Index=0;Index<Opt->NoOfRates;Index++)
	{
		FRateName = FormatStr(Opt->RateName[Index]);
		fprintf(Str, "%s", FRateName);

		free(FRateName);

		switch(Opt->ResTypes[Index])
		{
			case RESNONE:
				if(Opt->UseRJMCMC == TRUE)
					fprintf(Str, "RJ MCMC\n");
				else
					fprintf(Str, "None\n");
			break;
			
			case RESCONST:
				fprintf(Str, "%f\n", Opt->ResConst[Index]);
			break;

			case RESRATE:
				fprintf(Str, "%s\n", Opt->RateName[Opt->ResNo[Index]]);
			break;
		}

	}
}

void	PrintPriorVals(FILE *Str, PRIOR *P)
{
	int		VIndex;

	fprintf(Str, "        %s - ", P->Name);
	
	if(P->UseHP == TRUE)
		fprintf(Str, "Hyper Prior ");

	fprintf(Str, "%s ", DISTNAMES[(int)P->Dist]);
		

	if(P->UseHP == FALSE)
	{
		for(VIndex=0;VIndex<DISTPRAMS[P->Dist];VIndex++)
			fprintf(Str, "%2.2f ", P->DistVals[VIndex]);
	}
	else
	{
		for(VIndex=0;VIndex<DISTPRAMS[P->Dist];VIndex++)
			fprintf(Str, "(%2.2f,%2.2f) ", P->HP[VIndex*2], P->HP[(VIndex*2)+1]);
	}
	fprintf(Str, "\n");
}

int		RateNameToPos(char *Name, OPTIONS *Opt)
{
	int Index;

	for(Index=0;Index<Opt->NoOfRates;Index++)
	{
		if(strcmp(Name, Opt->RateName[Index]) == 0)
			return Index;
	}

	return -1;
}

void	PrintPriorOpt(FILE* Str, OPTIONS *Opt)
{
	int		Index;
	PRIOR	*Prior;
	int		RatePos;

	fprintf(Str, "Prior Information:\n");
	if(Opt->ModelType == MT_DISCRETE)
		fprintf(Str, "    Prior Categories:            %d\n", Opt->PriorCats);

	fprintf(Str, "    Priors\n");

	for(Index=0;Index<Opt->NoAllPriors;Index++)
	{
		Prior = Opt->AllPriors[Index];
		RatePos = RateNameToPos(Prior->Name, Opt);

		if(RatePos == -1)
			PrintPriorVals(Str, Prior);
		else
		{
			if(Opt->ResTypes[RatePos] == RESNONE && Opt->UseRJMCMC == FALSE)
				PrintPriorVals(Str, Prior);
		}
	}
}

void	PrintConPar(FILE* Str, int InUse, int Est, double Const)
{
	if(Est == TRUE)
	{
		fprintf(Str, "Estimate\n");
		return;
	}

	if(Const != -1)
	{
		fprintf(Str, "%f\n", Const);
		return;
	}

	fprintf(Str, "Not in use\n");
}

void	PrintEstData(FILE *Str, OPTIONS *Opt)
{
	int		TIndex;
	int		SIndex;
	TREES	*Trees;
	TAXA	*Taxa;

	Trees = Opt->Trees;
		
	if(EstData(Trees) == FALSE)
		return;

	fprintf(Str, "Estimating values for taxa and Sites\n");

	for(TIndex=0;TIndex<Trees->NoTaxa;TIndex++)
	{
		Taxa = Trees->Taxa[TIndex];
		if(Taxa->EstData == TRUE)
		{
			fprintf(Str, "\t");
			PrintFixSize(Taxa->Name, 20, Str);
			for(SIndex=0;SIndex<Trees->NoSites;SIndex++)
			{
				if(Taxa->EstDataP[SIndex] == TRUE)
					fprintf(Str, " %d ", SIndex+1);
			}
			fprintf(Str, "\n");
		}
	}
}

void	PrintHetMapPart(FILE *Str, TREES *Trees, PART *Part)
{
	int ID, Index;
	TAXA *Taxa;

	for(Index=0;Index<Part->NoTaxa;Index++)
	{
		ID = Part->Taxa[Index];
		Taxa = Trees->Taxa[ID];
		fprintf(Str, "%s\t", Taxa->Name);
	}
}


void	PrintHetMap(FILE *Str, OPTIONS *Opt, TREES *Trees)
{
	int Index;
	TREE *Tree;
	NODE N;

	Tree = Trees->Tree[0];

	fprintf(Str, "Hetro Model Key:\n");

	for(Index=1;Index<Tree->NoNodes;Index++)
	{
		N = Tree->NodeList[Index];

		fprintf(Str, "\tHNode\t%d\t", Index);
		PrintHetMapPart(Str, Trees, N->Part);

		fprintf(Str, "\n");
	}
}

void	PrintRecNodes(FILE* Str, OPTIONS *Opt)
{
	int Index, RIndex;
	RECNODE *RNode;

	if(Opt->NoOfRecNodes == 0)
		return;

	fprintf(Str, "Node reconstruction / fossilisation:\n");

	for(RIndex=0;RIndex<Opt->NoOfRecNodes;RIndex++)
	{

		RNode = Opt->RecNodeList[RIndex];
		if(RNode->NodeType == MRCA)
			fprintf(Str, "\tMRCA %s %s\n", RNode->Name, RNode->Tag->Name);

		if(RNode->NodeType == NODEREC)
			fprintf(Str, "\tNode %s %s\n", RNode->Name, RNode->Tag->Name);

		if(RNode->NodeType == FOSSIL)
		{
			if(Opt->Model == M_MULTISTATE)
			{
				fprintf(Str, "\tFossil %s %s ", RNode->Name, RNode->Tag->Name);
				for(Index=0;Index<RNode->NoFossilStates;Index++)
					fprintf(Str, "%c ", Opt->Trees->SymbolList[RNode->FossilStates[Index]]);
				fprintf(Str, "\n");
			}
			else
				fprintf(Str, "\tFossil %s %s (%d)\n", RNode->Name, RNode->Tag->Name, RNode->FossilStates[0]);
		}
	}
}

void	PrintRJLocalTrans(FILE* Str, OPTIONS *Opt)
{
	fprintf(Str, "Min Trans Taxa No:               %d\n", Opt->MinTransTaxaNo);
	
	fprintf(Str, "RJ Local Branch:                 ");
	if(Opt->UseRJLocalScalar[VR_BL] == TRUE)
		fprintf(Str, "True\n");
	else
		fprintf(Str, "False\n");	

	fprintf(Str, "RJ Local Node:                   ");
	if(Opt->UseRJLocalScalar[VR_NODE] == TRUE)
		fprintf(Str, "True\n");
	else
		fprintf(Str, "False\n");	

	fprintf(Str, "RJ Local Kappa:                  ");
	if(Opt->UseRJLocalScalar[VR_KAPPA] == TRUE)
		fprintf(Str, "True\n");
	else
		fprintf(Str, "False\n");
			
	fprintf(Str, "RJ Local Lambda:                 ");
	if(Opt->UseRJLocalScalar[VR_LAMBDA] == TRUE)
		fprintf(Str, "True\n");
	else
		fprintf(Str, "False\n");
		
	fprintf(Str, "RJ Local Delta:                  ");
	if(Opt->UseRJLocalScalar[VR_DELTA] == TRUE)
		fprintf(Str, "True\n");
	else
		fprintf(Str, "False\n");
		
	fprintf(Str, "RJ Local OU:                     ");
	if(Opt->UseRJLocalScalar[VR_OU] == TRUE)
		fprintf(Str, "True\n");
	else
		fprintf(Str, "False\n");	

	fprintf(Str, "RJ Local Landscape BL:           ");
	if(Opt->UseRJLocalScalar[VR_LS_BL] == TRUE)
		fprintf(Str, "True\n");
	else
		fprintf(Str, "False\n");	

}

void	PrintEmpPis(FILE* Str, OPTIONS *Opt)
{
	int NOS, Index;
	double *Pis;
	
	NOS = Opt->Trees->NoStates;

	Pis = GetEmpPis(Opt);

	fprintf(Str, "(");
	for(Index=0;Index<NOS-1;Index++)
		fprintf(Str, "%f,", Pis[Index]);
	fprintf(Str, "%f)", Pis[Index]);

	free(Pis);
}

void	PrintOptions(FILE* Str, OPTIONS *Opt)
{
	int		Index, NOS;
	

	fprintf(Str, "Options:\n");

	fprintf(Str, "Model:                           %s\n", MODELNAMES[Opt->Model]);
	
	fprintf(Str, "Tree File Name:                  %s\n", Opt->TreeFN);
	fprintf(Str, "Data File Name:                  %s\n", Opt->DataFN);
	fprintf(Str, "Log File Name:                   %s%s\n", Opt->BaseOutputFN, OUTPUT_EXT_LOG);

	fprintf(Str, "Save Initial Trees:              ");
	if(Opt->SaveInitialTrees == NULL)
		fprintf(Str, "False\n");
	else
		fprintf(Str, "%s\n", Opt->SaveInitialTrees);

	fprintf(Str, "Save Trees:                      ");
	if(Opt->SaveTrees == FALSE)
		fprintf(Str, "False\n");
	else
		fprintf(Str, "True\n");		

	fprintf(Str, "Summary:                         ");
	if(Opt->Summary == FALSE)
		fprintf(Str, "False\n");
	else
		fprintf(Str, "True\n");
	
	fprintf(Str, "Seed:                            %lu\n", Opt->Seed);

	if(Opt->FatTailNormal == TRUE)
		fprintf(Str, "Fat Tail Normal:             True\n");

	if(Opt->MakeUM == TRUE)
		fprintf(Str, "Make UM:                     True\n");

	if(Opt->Analsis == ANALML)
	{
		fprintf(Str, "Analsis Type:                    Maximum Likelihood\n" );
		fprintf(Str, "ML Attempt Per Tree:             %d\n", Opt->MLTries);
		fprintf(Str, "ML Max Evaluations:              %d\n", Opt->MLMaxEVals);
		fprintf(Str, "ML Tolerance:                    %f\n", Opt->MLTol);
		fprintf(Str, "ML Algorithm:                    %s\n", Opt->MLAlg);
		fprintf(Str, "Rate Range:                      %f - %f\n", Opt->RateMin, Opt->RateMax);
	}
	
	fprintf(Str, "Precision:                       %d bits\n", Opt->Precision);
	fprintf(Str, "Cores:                           %d\n", Opt->Cores);

	if(Opt->Analsis == ANALMCMC)
	{
		fprintf(Str, "Analysis Type:                   MCMC\n" );
		fprintf(Str, "Sample Period:                   %d\n", Opt->Sample);
		fprintf(Str, "Iterations:                      %lld\n", Opt->Itters);
		fprintf(Str, "Burn in:                         %lld\n", Opt->BurnIn);

		fprintf(Str, "MCMC ML Start:                   ");
		
		if(Opt->MCMCMLStart == FALSE)
			fprintf(Str, "False\n");
		else
			fprintf(Str, "True\n");

		if(Opt->UseRJMCMC == TRUE)
		{
			fprintf(Str, "Use RJ MCMC:                     True\n");
			if(Opt->CapRJRatesNo != -1)
				fprintf(Str, "Cap RJ Rate Number:              %d\n", Opt->CapRJRatesNo);

			fprintf(Str, "Zero Rates:                      ");
			if(Opt->RJZero == TRUE)
				fprintf(Str, "Yes\n");
			else
				fprintf(Str, "No\n");
		}
			
		if(Opt->UseSchedule	== TRUE)
			fprintf(Str, "Schedule File:                   %s%s\n", Opt->BaseOutputFN, OUTPUT_EXT_SCHEDULE);

		fprintf(Str, "Rate Dev:                        AutoTune\n");

		if(Opt->Model == M_FATTAIL)
			fprintf(Str, "No Slice Sample Steps:           %d\n", Opt->NoSliceSampleSteps);
	}

	if(Opt->RJDummy == TRUE)
		fprintf(Str, "RJDummy Codeing:                 True\n");
		

	if(Opt->NOSPerSite == TRUE)
		fprintf(Str, "Fit no of states per Site:       Yes\n");
	fprintf(Str, "No of Rates:                     %d\n", Opt->NoOfRates);

	if(Opt->DataType == DISCRETE)
	{
		fprintf(Str, "Base frequency (PI's):           ");
		switch(Opt->PiTypes)
		{
			case PI_NONE:
				fprintf(Str, "None\n");
				break;
			case PI_EMP:
				fprintf(Str, "Empirical: ");
				PrintEmpPis(Str, Opt);
				fprintf(Str, "\n");
				break;
			case PI_UNI:
				fprintf(Str, "Uniform\n");
				break;
		}
		
		fprintf(Str, "Pis used in ancestral state estimation:      ");
		if(Opt->UsePisInAncStates == TRUE)
			fprintf(Str, "Yes\n");
		else
			fprintf(Str, "No\n");
		
		fprintf(Str, "Character Symbols:               ");

		if(Opt->Model == M_MULTISTATE)
		{
			NOS = Opt->Trees->NoStates;
			if(Opt->UseCovarion == TRUE)
				NOS = NOS / 2;
			for(Index=0;Index<NOS-1;Index++)
				fprintf(Str, "%c,", Opt->Trees->SymbolList[Index]);
			fprintf(Str, "%c\n", Opt->Trees->SymbolList[Index]);
		}
		else
			fprintf(Str, "00,01,10,11\n");
	}

	if(Opt->DataType == CONTINUOUS)
	{
		fprintf(Str, "Test for trait correlation:      ");
		if(Opt->TestCorrel == TRUE)
			fprintf(Str, "True\n");
		else
			fprintf(Str, "False\n");

		fprintf(Str, "Kappa:                           ");
		PrintConPar(Str, Opt->UseKappa, Opt->EstKappa, Opt->FixKappa);

		fprintf(Str, "Delta:                           ");
		PrintConPar(Str, Opt->UseDelta, Opt->EstDelta, Opt->FixDelta);

		fprintf(Str, "Lambda:                          ");
		PrintConPar(Str, Opt->UseLambda, Opt->EstLambda, Opt->FixLambda);

		fprintf(Str, "OU:                              ");
		PrintConPar(Str, Opt->UseOU, Opt->EstOU, Opt->FixOU);
		
		if(Opt->Analsis == ANALMCMC)
			PrintRJLocalTrans(Str, Opt);

		fprintf(Str, "RJ Landscape Rate Group:         ");
		if(Opt->UseRJLandscapeRateGroup == TRUE)
			fprintf(Str, "True\n");
		

		if(Opt->AlphaZero == TRUE)
			fprintf(Str, "Alpha through zero:              True\n");

		if(Opt->NodeData == TRUE)
			fprintf(Str, "Model for Node Data:             True\n");

		if(Opt->NodeBLData == TRUE)
			fprintf(Str, "Model for Node BLS Data:         True\n");
	
	}
	else
	{
		if(Opt->UseKappa == TRUE)
		{
			fprintf(Str, "Kappa:                           ");
			PrintConPar(Str, Opt->UseKappa, Opt->EstKappa, Opt->FixKappa);
		}

		if(Opt->UseGamma == TRUE)
		{
			fprintf(Str, "Gamma:                           ");
			PrintConPar(Str, Opt->UseGamma, Opt->EstGamma, Opt->FixGamma);
			fprintf(Str, "Gamma Categories:                %d\n", Opt->GammaCats);
		}
		
		fprintf(Str, "Using a covarion model:          ");
		if(Opt->UseCovarion == TRUE)
			fprintf(Str, "True\n");
		else
			fprintf(Str, "False\n");

		if(Opt->UseRModel == TRUE)
		{
			fprintf(Str, "Using R Model:\n");
			if(Opt->RModelP != -1)
				fprintf(Str, "R Model Rates fixed to:          %f", Opt->RModelP);
		}
	}
	
	if(Opt->DataType == DISCRETE)
	{
		fprintf(Str, "Normalise Q Matrix:              ");
		if(Opt->NormQMat == TRUE)
			fprintf(Str, "True\n");
		else
			fprintf(Str, "False\n");
	}

	PrintTags(Str, Opt);

	if(Opt->NoLocalTransforms > 0)
		PrintLocalTransforms(Str, Opt->LocalTransforms, Opt->NoLocalTransforms);

	if(Opt->SaveModels == TRUE)
		fprintf(Str, "Save Model:                      %s\n", Opt->SaveModelsFN);

	if(Opt->LoadModels == TRUE)
		fprintf(Str, "Load Model:                      %s\n", Opt->LoadModelsFN);

	PrintEstData(Str, Opt);

	if(Opt->ScaleTrees != -1)
		fprintf(Str, "Scale Tree:                      %f\n", Opt->ScaleTrees);

	if(Opt->AnalyticalP == TRUE)
		fprintf(Str, "Analytical P:                    True\n");
	
	if(Opt->Model == M_DISC_HET)
	{
		fprintf(Str, "Tree 1 Partitions :			   \t");
	//	PrintTreePart(Str, Opt->Trees, 0);
	}

	PrintOptRes(Str, Opt);
	if(Opt->Analsis == ANALMCMC)
		PrintPriorOpt(Str, Opt);

	if(Opt->Model != M_CONTINUOUS_RR && Opt->Model != M_CONTINUOUS_DIR)
		PrintRecNodes(Str, Opt);

	if(Opt->Trees->NoOfRemovedTaxa != 0)
	{
		fprintf(Str, "Removed taxa:\n");
		for(Index=0;Index<Opt->Trees->NoOfRemovedTaxa;Index++)
			fprintf(Str, "          %s\n", Opt->Trees->RemovedTaxa[Index]);
	}
	
	if(Opt->Model == M_DISC_HET)
		PrintHetMap(Str, Opt, Opt->Trees);
	
	if(Opt->DistData != NULL)
		PrintDistData(Str, Opt->DistData);

	if(Opt->NoCShed > 0)
	{
		fprintf(Str, "Custom Schedule:\n");
		PrintCustomSchedule(Str, Opt->NoCShed, Opt->CShedList);
	}
	


	PrintTimeSlices(Str, Opt->TimeSlices);
	PrintPatterns(Str, Opt->NoPatterns, Opt->PatternList);

	PrintTreesInfo(Str, Opt->Trees, Opt->DataType);
	fflush(Str);
}

void	FreeOptions(OPTIONS *Opt, int NoSites)
{
	int		Index;
	
	for(Index=0;Index<Opt->NoOfRates;Index++)
		free(Opt->RateName[Index]);
	free(Opt->RateName);

	for(Index=0;Index<Opt->DefNoRates;Index++)
		free(Opt->DefRateNames[Index]);
	free(Opt->DefRateNames);

	if(Opt->AllPriors != NULL)
	{
		for(Index=0;Index<Opt->NoAllPriors;Index++)
			FreePrior(Opt->AllPriors[Index]);
		free(Opt->AllPriors);
	}

	if(Opt->EstDataSites != NULL)
		free(Opt->EstDataSites);
	
	free(Opt->DataFN);
	free(Opt->TreeFN);
	free(Opt->BaseOutputFN);
	fclose(Opt->LogFile);

	if(Opt->MLAlg != NULL)
		free(Opt->MLAlg);

	if(Opt->LogFileRead != NULL)
		fclose(Opt->LogFileRead);

	if(Opt->LogFileBuffer != NULL)
		free(Opt->LogFileBuffer);

	if(Opt->PassedOut != NULL)
		free(Opt->PassedOut);

	free(Opt->ResTypes);
	free(Opt->ResNo);
	free(Opt->ResConst);
	 
	if(Opt->SaveInitialTrees != NULL)
		free(Opt->SaveInitialTrees);
	
	FreeRecNodes(Opt, NoSites);

	if(Opt->SaveModelsFN != NULL)
		free(Opt->SaveModelsFN);

	if(Opt->LoadModelsFN != NULL)
		free(Opt->LoadModelsFN);

	if(Opt->Stones != NULL)
		FreeStones(Opt->Stones);

	for(Index=0;Index<Opt->NoLocalTransforms;Index++)
		FreeLocalTransforms(Opt->LocalTransforms[Index]);

	if(Opt->LocalTransforms != NULL)
		free(Opt->LocalTransforms);
	
	for(Index=0;Index<Opt->NoTags;Index++)
		FreeTag(Opt->TagList[Index]);
	if(Opt->TagList != NULL)
		free(Opt->TagList);

	if(Opt->DistData != NULL)
		FreeDistData(Opt->DistData);

	if(Opt->NoCShed > 0)
	{
		for(Index=0;Index<Opt->NoCShed;Index++)
			FreeCustomSchedule(Opt->CShedList[Index]);

		free(Opt->CShedList);
	}

	FreeTimeSlices(Opt->TimeSlices);

	for(Index=0;Index<Opt->NoPatterns;Index++)
		FreePattern(Opt->PatternList[Index]);
	
	if(Opt->PatternList != NULL)
		free(Opt->PatternList);


	if(Opt->RateScalars != NULL)
		free(Opt->RateScalars);

	free(Opt);
}

char*	CreatRateName(char N1, char N2)
{
	char	*Ret;

	Ret = (char*)SMalloc(sizeof(char) * 4);

	sprintf(Ret, "q%c%c", N1, N2);
	
	return Ret;
}

char**	ModelARateName(OPTIONS* Opt)
{
	char**	Ret;
	char*	Buffer;
	int		Index;

	Opt->NoOfRates = Opt->Trees->NoSites;

	Ret = (char**)SMalloc(sizeof(char*)*Opt->NoOfRates);
	Buffer = (char*)SMalloc(sizeof(char) * BUFFERSIZE);

	for(Index=0;Index<Opt->NoOfRates;Index++)
	{
		sprintf(Buffer, "alpha-%d", Index+1);
		Ret[Index] = StrMake(Buffer);
	}

	free(Buffer);
	return Ret;
}

char**	ModelBRateName(OPTIONS* Opt)
{
	char**	Ret;
	char*	Buffer;
	int		No, Index;

	Ret = (char**)SMalloc(sizeof(char*)*Opt->NoOfRates);
	Buffer = (char*)SMalloc(sizeof(char) * BUFFERSIZE);

	for(Index=0;Index<Opt->NoOfRates;Index++)
	{
		if(Index<Opt->Trees->NoSites)
		{
			No = Index+1;
			sprintf(Buffer, "Alpha-%d", No);
		}
		else
		{
			No = (Index - Opt->Trees->NoSites) + 1;
			sprintf(Buffer, "Beta-%d", No);
		}

		Ret[Index] = StrMake(Buffer);
	}

	free(Buffer);
	return Ret;
}

char**	RetModelRateName(OPTIONS* Opt)
{
	char**	Ret;
	char*	Buffer;
	int		Index;
	
	Ret = (char**)SMalloc(sizeof(char*) * Opt->NoOfRates);
	Buffer = (char*)SMalloc(sizeof(char) * BUFFERSIZE);

	sprintf(Buffer, "Alpha");
	Ret[0] = StrMake(Buffer);
	
	for(Index=1;Index<Opt->NoOfRates;Index++)
	{
		sprintf(Buffer, "Beta-%d", Index+1);
		Ret[Index] = StrMake(Buffer);
	}
	 
	free(Buffer);
	return Ret;
}

char**	ContrastRateNames(OPTIONS *Opt)
{
	char	**Ret;
	char	*Buffer;
	int		Index, NOS, i;
	
	NOS = Opt->Trees->NoSites;
	 
	Ret = (char**)SMalloc(sizeof(char**) * Opt->NoOfRates);
	Buffer = (char*)SMalloc(sizeof(char*) * BUFFERSIZE);
	
	i = 0;
	for(Index=0;Index<Opt->Trees->NoSites;Index++)
	{
		sprintf(Buffer, "Alpha-%d", Index+1);
		Ret[i++] = StrMake(Buffer);
	}

	free(Buffer);
	return Ret;
}

char**	ContrastFullRateNames(OPTIONS *Opt)
{
	char	**Ret;
	char	*Buffer;
	int		Index, NOS, i;
	
	NOS = Opt->Trees->NoSites;
	 
	Ret = (char**)SMalloc(sizeof(char**) * Opt->NoOfRates);
	Buffer = (char*)SMalloc(sizeof(char*) * BUFFERSIZE);
	
	i = 0;
	for(Index=0;Index<NOS;Index++)
	{
		sprintf(Buffer, "Alpha-%d", Index+1);
		Ret[i++] = StrMake(Buffer);
	}

	for(Index=0;Index<NOS;Index++)
	{
		sprintf(Buffer, "Sigma-%d", Index+1);
		Ret[i++] = StrMake(Buffer);
	}

	free(Buffer);
	return Ret;
}

char**	ContrastRegRateNames(OPTIONS *Opt)
{
	char	**Ret;
	char	*Buffer;
	int		Index, Pos;
		 
	Ret = (char**)SMalloc(sizeof(char**) * Opt->NoOfRates);
	Buffer = (char*)SMalloc(sizeof(char*) * BUFFERSIZE);
	
	Pos = 0;
	for(Index=1;Index<Opt->Trees->NoSites;Index++)
	{
		sprintf(Buffer, "Beta-%d", Index);
		Ret[Pos++] = StrMake(Buffer);
	}

	free(Buffer);
	return Ret;
}

char**	FatTailRateNames(OPTIONS *Opt)
{
	char	**Ret;
	char	*Buffer;
	int		Index, Pos;
	 
	Ret = (char**)SMalloc(sizeof(char**) * Opt->NoOfRates);
	Buffer = (char*)SMalloc(sizeof(char*) * BUFFERSIZE);
		
	Pos = 0;
	for(Index=0;Index<Opt->Trees->NoSites;Index++)
	{
		sprintf(Buffer, "Sig2-%d", Index+1);
		Ret[Pos++] = StrMake(Buffer);
	}

	free(Buffer);
	return Ret;
}

char**	GeoRateNames(OPTIONS *Opt)
{
	char	**Ret;
	char	*Buffer;

	Ret = (char**)SMalloc(sizeof(char**) * Opt->NoOfRates);
	Buffer = (char*)SMalloc(sizeof(char*) * BUFFERSIZE);

	sprintf(Buffer, "Alpha");
	Ret[0] = StrMake(Buffer);

	sprintf(Buffer, "Scale");
	Ret[1] = StrMake(Buffer);

	free(Buffer);
	return Ret;
}

char**	CreatContinusRateName(OPTIONS* Opt)
{

	Opt->NoOfRates  = FindNoConRates(Opt);

	switch(Opt->Model)
	{
		case M_CONTINUOUS_RR:
			return ModelARateName(Opt);

		case M_CONTINUOUS_DIR:
			return ModelBRateName(Opt);

		case M_CONTINUOUS_REG:
			return RetModelRateName(Opt);

		case M_CONTRAST:
			return ContrastFullRateNames(Opt);

		case M_CONTRAST_CORREL:
			return ContrastRateNames(Opt);

		case M_CONTRAST_REG:
			return ContrastRegRateNames(Opt);

		case M_FATTAIL:
			return FatTailRateNames(Opt);

		case M_GEO:
			return GeoRateNames(Opt);
		// Keep CLang happy
		default:
			break;
	}

	return NULL;
}

void	SetOptRateNamesFixed(OPTIONS *Opt, int NoRates, char *Rates[])
{
	int Index;

	Opt->NoOfRates = NoRates;
	Opt->RateName = (char**)SMalloc(sizeof(char*) * NoRates);

	for(Index=0;Index<NoRates;Index++)
		Opt->RateName[Index] = StrMake(Rates[Index]);
}

double*	AllocDefRateScalars(int NoPram)
{
	double *Ret;
	int Index;

	Ret = (double*)SMalloc(sizeof(double) * NoPram);

	for(Index=0;Index<NoPram;Index++)
		Ret[Index] = 1.0;

	return Ret;
}

void	SetOptRates(OPTIONS* Opt, int NOS, char *SymbolList)
{
	int		Index, Inner, Outter;

	if(Opt->DataType == CONTINUOUS)
	{
		Opt->RateName	= CreatContinusRateName(Opt);
		return;
	}

	if(Opt->Model == M_DISC_INDEP)
		SetOptRateNamesFixed(Opt, 4, INDEPPRAMS);

	if(Opt->Model == M_DISC_DEP)
		SetOptRateNamesFixed(Opt, 8, DEPPRAMS);

	if(Opt->Model == M_DISC_CV)
		SetOptRateNamesFixed(Opt, 14, DEPCVPRAMS);
	
	if(Opt->Model == M_DISC_HET)
		SetOptRateNamesFixed(Opt, 12, DEPHETROPRAMS);
	
	if(Opt->Model == M_MULTISTATE)
	{
		Opt->NoOfRates	= (NOS * NOS) - NOS;
		Opt->RateName	= (char**)SMalloc(sizeof(char*)*Opt->NoOfRates);

		for(Outter=0,Index=0;Outter<NOS;Outter++)
		{
			for(Inner=0;Inner<NOS;Inner++)
			{
				if(Outter != Inner)
				{
					Opt->RateName[Index] = CreatRateName(SymbolList[Outter], SymbolList[Inner]);
					Index++;
				}
			}
		}
	}

	Opt->RateScalars = AllocDefRateScalars(Opt->NoOfRates);
}

void		AllocRestictions(OPTIONS *Opt)
{
	int	Index;

	if(Opt->ResTypes != NULL)
		free(Opt->ResTypes);
	
	if(Opt->ResNo != NULL)
		free(Opt->ResNo);
 
	if(Opt->ResConst != NULL)
		free(Opt->ResConst);
	
	Opt->ResTypes		= (RESTYPES*)SMalloc(sizeof(RESTYPES) * Opt->NoOfRates);
	Opt->ResNo			= (int*)SMalloc(sizeof(int) * Opt->NoOfRates);
	Opt->ResConst		= (double*)SMalloc(sizeof(double) * Opt->NoOfRates);

	for(Index=0;Index<Opt->NoOfRates;Index++)
	{
		Opt->ResTypes[Index] = RESNONE;
		Opt->ResNo[Index] = -1;
		Opt->ResConst[Index] = -1;
	}
}

void	SetFatTailPrior(OPTIONS *Opt)
{
	int Index, Pos;
	PRIOR *Prior;


	Pos = 0;
	
	for(Index=0;Index<Opt->Trees->NoSites;Index++)
	{
//		Prior = CreateUniformPrior(Opt->RateName[Pos], 0.2, 2.0);
//		AddPriorToOpt(Opt, Prior);

		Prior = CreateUniformPrior(Opt->RateName[Pos], 0.0, 100.0);
		AddPriorToOpt(Opt, Prior);
	}
}

void	GetGeoPriors(OPTIONS *Opt)
{
	PRIOR *Prior;

	Prior = CreateUniformPrior(Opt->RateName[0], 0.2, 2.0);
	AddPriorToOpt(Opt, Prior);

	Prior = CreateUniformPrior(Opt->RateName[1], 0.0, 1000000.0);
	AddPriorToOpt(Opt, Prior);
}

void	SetAnsStatesEst(OPTIONS *Opt, TREES *Trees)
{
	int TIndex, SIndex;
	TAXA *Taxa;

	if(EstData(Trees) == FALSE)
		return;

	for(TIndex=0;TIndex<Trees->NoTaxa;TIndex++)
	{
		Taxa = Trees->Taxa[TIndex];
		for(SIndex=0;SIndex<Trees->NoSites;SIndex++)
			if(Taxa->EstDataP[SIndex] == TRUE)
				AddConAnsStatePrior(Opt, SIndex+1);
	}
}

void	AllocRatePriors(OPTIONS *Opt, TREES *Trees)
{
	int		Index;
	PRIOR	*Prior;
	
	if(Opt->ModelType == MT_CONTINUOUS)
		SetAnsStatesEst(Opt, Trees);

	if(Opt->Model == M_FATTAIL)
	{
		SetFatTailPrior(Opt);
		return;
	}

	if(Opt->Model == M_GEO)
	{
		GetGeoPriors(Opt);
		return;
	}

	for(Index=0;Index<Opt->NoOfRates;Index++)
	{
		if(Opt->ModelType == DISCRETE)
			Prior = CreateUniformPrior(Opt->RateName[Index], 0, 100);
		else
			Prior = CreateUniformPrior(Opt->RateName[Index], -100, 100);
				
		if(Opt->Model == M_CONTRAST && Index >= Opt->Trees->NoSites)
			Prior->DistVals[0] = 0;
		
		AddPriorToOpt(Opt, Prior);
	}
/*
	if(Opt->Model == M_CONTRAST_REG)
	{
		Prior = CreateUniformPrior("Var", -100, 100);
		AddPriorToOpt(Opt, Prior);
	}
*/
}

MODEL_TYPE	GetModelType(MODEL Model)
{
	switch(Model)
	{
		case	M_MULTISTATE:		return MT_DISCRETE; break;
		case	M_DISC_INDEP:		return MT_DISCRETE; break;
		case	M_DISC_DEP:			return MT_DISCRETE; break;
		case	M_CONTINUOUS_RR:	return MT_CONTINUOUS; break;
		case	M_CONTINUOUS_DIR:	return MT_CONTINUOUS; break;
		case	M_CONTINUOUS_REG:	return MT_CONTINUOUS; break;
		case	M_CONTRAST_CORREL:	return MT_CONTRAST; break;
		case	M_CONTRAST_REG:		return MT_CONTRAST; break;
		case	M_CONTRAST:			return MT_CONTRAST; break;
		case	M_DISC_CV:			return MT_DISCRETE; break;
		case	M_DISC_HET:			return MT_DISCRETE; break;
		case	M_FATTAIL:			return MT_FATTAIL; break;
		case	M_GEO:				return MT_FATTAIL; break;
	}

	printf("Unkown model type (%s::%d)\n", __FILE__, __LINE__);
	exit(1);

	return MT_DISCRETE;
}

void		SetDefRates(OPTIONS *Opt)
{
	int Index;

	Opt->DefNoRates = Opt->NoOfRates;
	Opt->DefRateNames = (char**)SMalloc(sizeof(char*) * Opt->DefNoRates);

	for(Index=0;Index<Opt->DefNoRates;Index++)
		Opt->DefRateNames[Index] = StrMake(Opt->RateName[Index]);
}

OPTIONS*	CreatOptions(MODEL Model, ANALSIS Analsis, int NOS, char *TreeFN, char *DataFN, char *SymbolList, TREES* Trees)
{
	OPTIONS *Ret;
	int		Index;
	char	*Buffer;

	Ret = (OPTIONS*)SMalloc(sizeof(OPTIONS));

	Buffer = (char*)SMalloc(sizeof(char) * BUFFERSIZE);

	Ret->Trees		= Trees;
	Ret->Model		= Model;
	Ret->Analsis	= Analsis;

	Ret->ModelType	= GetModelType(Model);

	Ret->TestCorrel	= FALSE;
	Ret->UseCovarion= FALSE;

	Ret->NodeData	= FALSE;
	Ret->NodeBLData = FALSE;
	Ret->AlphaZero	= FALSE;
	Ret->HPDev		= 1;

	Ret->OutTrees	= NULL;
	Ret->VarRatesLog= NULL;

	Ret->LogFatTail	= NULL;
	Ret->LogLandscapeGroups	 = NULL;
	
	Ret->UseRModel	= FALSE;
	Ret->RModelP	= -1;
	Ret->EstData	= FALSE;


	Ret->NoEstDataSite	=	0;
	Ret->EstDataSites	=	NULL;
	Ret->NoEstChanges	=	5;

	Ret->NOSPerSite		=	FALSE;

	if(Ret->ModelType == MT_DISCRETE)
		Ret->DataType = DISCRETE;

	if(	Ret->ModelType == MT_CONTINUOUS	|| 
		Ret->ModelType == MT_CONTRAST	||
		Ret->ModelType	== MT_FATTAIL)
	{
		Ret->TestCorrel = TRUE;
		Ret->DataType	= CONTINUOUS;
	}

	Ret->ResConst	= NULL;
	Ret->ResTypes	= NULL;
	Ret->ResNo		= NULL;

	Ret->DefNoRates	= -1;
	Ret->DefRateNames= NULL;

	Ret->RateScalars = NULL;

	SetOptRates(Ret, NOS, SymbolList);

	SetDefRates(Ret);
	
	AllocRestictions(Ret);
			
	Ret->TreeFN = StrMake(TreeFN);
	Ret->DataFN = StrMake(DataFN);

	Ret->BaseOutputFN = StrMake(DataFN);

	Ret->LogFile		= NULL;
	Ret->LogFileRead	= NULL;
	Ret->LogFileBuffer	= NULL;
	Ret->PassedOut		= NULL;		
	Ret->UseSchedule	= FALSE;

	Ret->MCMCMLStart	= FALSE; 


	Ret->MLTries		=	-1;
	Ret->MLTol			=	-1;
	Ret->MLMaxEVals		=	-1;
	Ret->MLAlg			=	NULL;



	Ret->NoAllPriors	= 0;
	Ret->AllPriors		= NULL;
	
	Ret->PriorCats		= -1;

	Ret->RateMin		= RATE_MIN;
	Ret->RateMax		= RATE_MAX;

	if(Ret->Analsis == ANALML)
	{
		Ret->Itters		=	-1;
		Ret->Sample		=	-1;
		Ret->BurnIn		=	-1;

		Ret->MLTries		=	10;
		Ret->MLTol			=	0.000001;
		Ret->MLMaxEVals		=	20000;
		Ret->MLAlg			=	StrMake("BOBYQA");
	}
	
	if(Ret->Analsis == ANALMCMC)
	{
		
		Ret->Itters		=	5050000;
		Ret->BurnIn		=	50000;
		
		Ret->Itters		=	1010000;
		Ret->BurnIn		=	10000;
		Ret->Sample		=	1000;
		
		Ret->PriorCats	=	100;
		
		AllocRatePriors(Ret, Trees);

		Ret->UseSchedule	= TRUE;
	}


	Ret->RecNodeList	=	NULL;
	Ret->NoOfRecNodes	=	0;
	Ret->Summary		=	FALSE;
	Ret->PiTypes		=	PI_NONE;

	Ret->UseKappa		=	FALSE;
	Ret->UseDelta		=	FALSE;
	Ret->UseLambda		=	FALSE;
	Ret->UseGamma		=	FALSE;
	Ret->UseOU			=	FALSE;

	Ret->EstKappa		=	FALSE;
	Ret->EstDelta		=	FALSE;
	Ret->EstLambda		=	FALSE;
	Ret->EstGamma		=	FALSE;
	Ret->EstOU			=	FALSE;


	Ret->FixKappa		=	-1;
	Ret->FixDelta		=	-1;
	Ret->FixLambda		=	-1;
	Ret->FixGamma		=	-1;
	Ret->FixOU			=	-1;

	Ret->InvertV		=	FALSE;

	Ret->UseRJMCMC		=	FALSE;
	Ret->CapRJRatesNo	=	-1;

	Ret->FindCF			=	FALSE;
	Ret->CFRate			=	NULL;

//	Ret->Headers		=	TRUE;

	Ret->AnalyticalP	=	FALSE;


	Ret->Seed			=	GetSeed();

	Ret->MakeUM			=	FALSE;

	Ret->UseEqualTrees	=	FALSE;
	Ret->ETreeBI		=	-1;

	Ret->SaveInitialTrees=	NULL;
	Ret->SaveTrees		=	FALSE;

	Ret->Precision		=	sizeof(double) * 8;

#ifdef BIG_LH
	Ret->Precision		=	256;
#endif

	Ret->Cores			=	GetMaxThreads();

	Ret->SaveModels		=	FALSE;
	Ret->SaveModelsFN	=	NULL;

	Ret->LoadModels		=	FALSE;
	Ret->LoadModelsFN	=	NULL;

	Ret->Stones			=	NULL;
	Ret->RJDummy		=	FALSE;
	Ret->RJDummyLog		=	NULL;

	Ret->RJDummyBetaDev =	0.1;

	Ret->ScaleTrees		=	-1;

	Ret->FatTailNormal	=	FALSE;

	for(Index=0;Index<NO_RJ_LOCAL_SCALAR;Index++)
		Ret->UseRJLocalScalar[Index] = FALSE;

	Ret->NoTags	= 0;
	Ret->TagList = NULL;

	Ret->NoLocalTransforms = 0;
	Ret->LocalTransforms = NULL;

	Ret->UseDistData = FALSE;
	Ret->DistData = NULL;

	Ret->NoLh = FALSE;

	Ret->NoCShed = 0;
	Ret->CShedList = NULL;
	
	Ret->TimeSlices = CreateTimeSlices();


	Ret->NoPatterns = 0;
	Ret->PatternList = NULL;

	Ret->MinTransTaxaNo = MIN_NO_TAXA_RJ_LOCAL_TRANS;

	Ret->NormQMat = FALSE;
	
	Ret->NoSliceSampleSteps = 100;

	Ret->UsePisInAncStates = TRUE;
	Ret->RJZero			   = TRUE;
	
	Ret->UseRJLandscapeRateGroup	 = FALSE;
	Ret->NoLandscapeRateGroup = -1;

	free(Buffer);

	return Ret; 
}

void	PrintModelChoic(TREES *Trees)
{
	printf("Please select the model of evolution to use.\n");
	printf("1)	MultiState\n");
	if(Trees->NoSites == 2 && Trees->NoStates == 2)
	{
		printf("2)\tDiscrete: Independent\n");
		printf("3)\tDiscrete: Dependant\n");
	}

	if(Trees->ValidCData == TRUE)
	{
		printf("4)\tContinuous: Random Walk (Model A)\n");
		printf("5)\tContinuous: Directional (Model B)\n");

		if(Trees->NoSites > 1)
			printf("6)\tContinuous: Regression\n");

		if(EstData(Trees) == FALSE)
		{
			printf("7)\tIndependent Contrast\n");

			printf("8)\tIndependent Contrast: Correlation\n");

			if(Trees->NoSites > 1)
				printf("9)\tIndependent Contrast: Regression\n");
		}
	}

	if(Trees->NoSites == 2 && Trees->NoStates == 2)
	{
		printf("10)\tDiscrete: Covarion\n");
		if(Trees->NoTrees == 1)
			printf("11)	Discrete: Heterogeneous \n");
	}

	if(Trees->ValidCData == TRUE)
		printf("12)\tFat Tail\n");
	
	if(Trees->ValidCData == TRUE && Trees->NoSites == 2)
		printf("13)\tGeo\n");
}

int		GetModelInt()
{
	char	*Buffer;
	int		Ret;

	Ret = -1;

	Buffer = (char*)SMalloc(sizeof(char) * BUFFERSIZE);

	if(fgets(Buffer, BUFFERSIZE, stdin) == NULL)
	{
		printf("Fatal error Reading model choice\n");
		exit(0);
	}
	
	ReplaceChar('\n', '\0', Buffer);

	if(IsValidInt(Buffer) == FALSE)
		printf("%s is not a valid model choice\n", Buffer);
	else
		Ret = atoi(Buffer);
		
	free(Buffer);

	return Ret;
}

int		ValidModelChoice(TREES *Trees, MODEL Model)
{
	if(Model == M_MULTISTATE)
		return TRUE;

	if(Model == M_DISC_INDEP || Model == M_DISC_DEP || Model == M_DISC_CV || Model == M_DISC_HET)
	{
		if(Trees->NoSites != 2 || Trees->NoStates != 2)
		{
			printf("Discrete analisis requiers two two state characters\n");
			printf("There are %d states and %d sites in the current data set.\n", Trees->NoStates, Trees->NoSites);
			return FALSE;
		}

		if(Model == M_DISC_HET && Trees->NoTrees != 1)
		{
			printf("Discrete: Heterogeneous requires a single tree.\n");
			return FALSE;
		}

		return TRUE;
	}

	if(Trees->ValidCData == FALSE)
	{
		printf("Model requires continues data.\n");
		return FALSE;
	}

	if(Model == M_CONTINUOUS_REG || Model == M_CONTRAST_REG)
	{
		if(Trees->NoSites < 2)
		{
			printf("Regression, requires two or more sites.\n");
			return FALSE;
		}
	}

	if(EstData(Trees) == TRUE)
	{
		if(GetModelType(Model) == MT_CONTRAST)
		{
			printf("Currently estimating unknown data values is not supported in independent contrast methods, please use the continues models.\n");
			return FALSE;
		}
	}
	
	return TRUE;
}

MODEL	IntToModel(int No)
{
	if(No == 1)
		return M_MULTISTATE;

	if(No == 2)
		return M_DISC_INDEP;

	if(No == 3)
		return M_DISC_DEP;

	if(No == 4)
		return M_CONTINUOUS_RR;

	if(No == 5)
		return M_CONTINUOUS_DIR;

	if(No == 6)
		return M_CONTINUOUS_REG;

	if(No == 7)
		return M_CONTRAST;

	if(No == 8)
		return M_CONTRAST_CORREL;
	
	if(No == 9)
		return M_CONTRAST_REG;

	if(No == 10)
		return M_DISC_CV;

	if(No == 11)
		return M_DISC_HET;

	if(No == 12)
		return M_FATTAIL;

	if(No == 13)
		return M_GEO;

	printf("Unknown model\n");
	exit(0);
}

void	GetModelChoice(TREES *Trees, MODEL *Model)
{
	int		No;
	MODEL	M;

	
	No = GetModelInt();
	if(No == -1)
		return;	

	M = IntToModel(No);

	*Model = M;

}

MODEL	GetModel(TREES *Trees)
{
	MODEL	Model;

	PrintModelChoic(Trees);
		
	GetModelChoice(Trees, &Model);
	
	return Model;
}


ANALSIS	GetAnalsis(TREES *Trees)
{
	char	*Buffer;
	int		No;
	
	Buffer =  (char*)SMalloc(sizeof(char) * BUFFERSIZE);

	printf("Please select the analysis method to use.\n");
	printf("1)	Maximum Likelihood.\n");
	printf("2)	MCMC\n");
	
	fgets(&Buffer[0], BUFFERSIZE, stdin);
	Buffer[BUFFERSIZE-1] = '\0';

	No = atoi(Buffer);
	free(Buffer);

	if(No != 1 && No != 2)
	{
		printf("Invalid analysis choice\n");
		exit(1);
	}
	
	if(No == 1)
		return ANALML;

	return ANALMCMC;
}

COMMANDS	StringToCommand(char *Str)
{
	int	CIndex=0;
	int	SIndex=0;

	if(Str[0] == '#')
		return CCOMMENT;
	do
	{
		if(StrICmp(Str, COMMANDSTRINGS[SIndex]) == 0 || StrICmp(Str, COMMANDSTRINGS[SIndex+1]) == 0)
			return (COMMANDS)CIndex;
		
		SIndex+=2;
		CIndex++;
	} while(strcmp(COMMANDSTRINGS[SIndex], "") != 0);

	return CUNKNOWN;
}

int		StrToRate(OPTIONS* Opt, char* Str)
{
	int	Index;
	
	for(Index=0;Index<Opt->NoOfRates;Index++)
	{
		if(strcmp(Str, Opt->RateName[Index]) == 0)
			return Index;
	}

	return -1;
}

int		RestrictToConst(OPTIONS *Opt, int Tokes, char *argv[], double Const)
{
	int	Index;
	int	RateNo;
	
	for(Index=1;Index<Tokes-1;Index++)
	{
		RateNo = StrToRate(Opt, argv[Index]);
		if(RateNo==-1)
		{
			printf("Rate paramtier: %s is unknown\n", argv[Index]);
			return FALSE;
		}

		Opt->ResTypes[RateNo] = RESCONST;
		Opt->ResConst[RateNo] = Const;
	}

	return TRUE;
}

int		RecRes(OPTIONS *Opt, int RateNo, int Verb)
{
	int	Current;

	if(Verb == TRUE)
		printf("%s -> ", Opt->RateName[RateNo]);


	Current = Opt->ResNo[RateNo];

	for(;;)
	{
		if(Verb == TRUE)
			printf("%s -> ", Opt->RateName[Current]);

		if((Opt->ResTypes[Current] == RESNONE) || (Opt->ResTypes[Current] == RESCONST))
			return FALSE;
		
		if(Current == RateNo)
		{
			if(Verb == TRUE)
				printf("Infinite\n");
			return TRUE;
		}

		Current = Opt->ResNo[Current];
	}

	return TRUE;
}

int		RestrictToRate(OPTIONS *Opt, int Tokes, char *argv[], char* Rate)
{
	int	ToNo;
	int	RateNo;
	int	Index;

	ToNo = StrToRate(Opt, Rate);
	if(ToNo == -1)
	{
		printf("Could not convert %s to a valid rate paramiter\n", Rate);
		return FALSE;
	}

	for(Index=1;Index<Tokes-1;Index++)
	{
		RateNo = StrToRate(Opt, argv[Index]);
		if(RateNo == -1)
		{
			printf("Could not restrict %s to a valid rate paramtier\n", argv[Index]);
			return FALSE;
		}

		Opt->ResTypes[RateNo] = RESRATE;
		Opt->ResNo[RateNo] = ToNo;

		if(RecRes(Opt, RateNo, FALSE) == TRUE)
		{
			printf("Restcting %s to %s cause a recusrive restriction\n", Opt->RateName[RateNo], Opt->RateName[ToNo]);

			RecRes(Opt, RateNo, TRUE);
			return FALSE;
		}
	}
	
	return TRUE;
}

void		RestrictAll(OPTIONS *Opt, char *To)
{
	double	Const;
	int		RateTo;
	int		Index;

	Const = atof(To);
	if((Const != 0) || (strcmp(To, "0")==0))
	{
		for(Index=0;Index<Opt->NoOfRates;Index++)
		{
			Opt->ResTypes[Index] = RESCONST;
			Opt->ResConst[Index] = Const;
		}
		return;
	}
	else
	{
		RateTo = StrToRate(Opt, To);
		if(RateTo == -1)
		{
			printf("Could not conver %s to a valid rate paramiter\n", To);
			return;
		}

		for(Index=0;Index<Opt->NoOfRates;Index++)
		{
			if(Index != RateTo)
			{
				Opt->ResTypes[Index]	= RESRATE;
				Opt->ResNo[Index]		= RateTo;
			}
			else
			{
				Opt->ResTypes[Index]	= RESNONE;
				Opt->ResConst[Index]	= -1;
				Opt->ResNo[Index]		= -1;
			}
		}
	}
}


void		Restrict(OPTIONS *Opt, int Tokes, char *argv[])
{
	RESTYPES	*BackRes;
	double		*BackConst;
	int			*BackNo;
	double		Const;
	int			Safe;
	
	BackRes		= (RESTYPES*)SMalloc(sizeof(RESTYPES) * Opt->NoOfRates);
	BackConst	= (double*)SMalloc(sizeof(double) * Opt->NoOfRates);
	BackNo		= (int*)SMalloc(sizeof(int) * Opt->NoOfRates);
	
	memcpy(BackRes, Opt->ResTypes, sizeof(RESTYPES) * Opt->NoOfRates); 
	memcpy(BackConst, Opt->ResConst, sizeof(double) * Opt->NoOfRates);
	memcpy(BackNo, Opt->ResNo, sizeof(int) * Opt->NoOfRates);

	Const = atof(argv[Tokes-1]);

	if((Const != 0) || (strcmp(argv[Tokes-1], "0")==0))
	{
		Safe = RestrictToConst(Opt, Tokes, argv, Const);
	}
	else
	{
		Safe = RestrictToRate(Opt, Tokes, argv, argv[Tokes-1]);
	}

	if(Safe == TRUE)
	{
		free(BackRes);
		free(BackConst);
		free(BackNo);
	}
	else
	{
		free(Opt->ResConst);
		free(Opt->ResNo);
		free(Opt->ResTypes);

		Opt->ResTypes		= BackRes;
		Opt->ResNo			= BackNo;
		Opt->ResConst		= BackConst;
	}
}

void	UnRestict(OPTIONS *Opt, char* Rate)
{
	int RateNo;
	
	RateNo = StrToRate(Opt, Rate);
	if(RateNo == -1)
	{
		printf("Could not conver %s to a valid rate paramtier\n", Rate);
		return;
	}
	else
	{
		Opt->ResTypes[RateNo]	= RESNONE;
		Opt->ResConst[RateNo]	= -1;
		Opt->ResNo[RateNo]		= -1;
	}
}

void	UnRestictAll(OPTIONS *Opt)
{
	int	Index;

	for(Index=0;Index<Opt->NoOfRates;Index++)
	{
		Opt->ResTypes[Index]	= RESNONE;
		Opt->ResConst[Index]	= -1;
		Opt->ResNo[Index]		= -1;
	}
}

void	SetPrior(OPTIONS *Opt, char *Name, int Tokes, char **argv)
{
	PRIOR	*CPrior, *NPrior;
	
	CPrior = GetPriorFromName(Name, Opt->AllPriors, Opt->NoAllPriors);

	if(CPrior == NULL)
	{
		printf("Cannot find prior %s, priors are case sensitive.\n", Name);
		exit(0);
	}

	NPrior = CreatePriorFromStr(Name, Tokes, argv);

	ReplacePrior(Opt, NPrior);

	FreePrior(CPrior);
}


void	SetPriorCmd(OPTIONS *Opt, int Tokes, char **argv)
{
	char	*Name;

	if(Tokes < 4)
	{
		printf("The prior command must take a parameter name, distribution and distribution parameters.");
		exit(0);
	}

	Name = argv[1];

	SetPrior(Opt, Name, Tokes-2, &argv[2]);
}

void	SetAllRatePriors(OPTIONS *Opt, int Tokes, char **argv)
{
	int Index;

	if(Tokes < 3)
	{
		printf("The set all prior command takes a prior.\n");
		exit(0);
	}

	for(Index=0;Index<Opt->NoOfRates;Index++)
		SetPrior(Opt, Opt->RateName[Index], Tokes-1, &argv[1]);
}

void	SetHyperPrior(OPTIONS *Opt, char *Name, int Tokes, char **argv)
{
	PRIOR	*CPrior, *NPrior;

	CPrior = GetPriorFromName(Name, Opt->AllPriors, Opt->NoAllPriors);

	if(CPrior == NULL)
	{

		printf("Cannot find prior %s, priors are case sensitive.\n", Name);
		exit(0);
	}

	NPrior = CreateHyperPriorFromStr(Name, Tokes, argv);

	ReplacePrior(Opt, NPrior);

	FreePrior(CPrior);
}


void	SetHyperPriorCmd(OPTIONS *Opt, int Tokes, char **argv)
{
	char *Name;

	if(Tokes < 4)
	{
		printf("The Hyper prior command must take a parameter name, distribution and min / max valus for distribution parameters.");
		exit(0);
	}

	Name = argv[1];

	SetHyperPrior(Opt, Name, Tokes-2, &argv[2]);
}

void	SetHyperPriorAllCmd(OPTIONS *Opt, int Tokes, char **argv)
{
	int	Index;

	for(Index=0;Index<Opt->NoOfRates;Index++)
		SetHyperPrior(Opt, Opt->RateName[Index], Tokes-1, &argv[1]);
}

void	PrintUnderNode(NODE N)
{
	int	Index;

	if(N->Tip == TRUE)
		printf("%s\t", N->Taxa->Name);
	else
	{
		for(Index=0;Index<N->NoNodes;Index++)
			PrintUnderNode(N->NodeList[Index]);
	}
}


RECNODE*	OptFindRecNode(OPTIONS *Opt, char* Name)
{
	RECNODE *Ret;
	int Index;

	for(Index=0;Index<Opt->NoOfRecNodes;Index++)
	{
		Ret = Opt->RecNodeList[Index];

		if(strcmp(Name, Ret->Name)==0)
			return Ret;
	}

	return NULL;
}

/*
X 	Likilhood values unchanged
-	Likilhood set to zero


Symbol	0,0	0,1	1,0	1,1
0		X	-	-	-
1		-	X	-	-
2		-	-	X	-
3		-	-	-	X
				
10		X	X	-	-
11		X	-	X	-
12		X	-	-	X
13		-	X	X	-
14		-	X	-	X
15		-	-	X	X
				
20		X	X	X	-
21		X	X	-	X
22		X	-	X	X
23		-	X	X	X
*/

void	ValidDescFossileState(int FNo)
{
	if((FNo >= 0) && (FNo <= 3))
		return;

	if((FNo >= 10) && (FNo <= 15))
		return;

	if((FNo >= 20) && (FNo <= 23))
		return ;

	printf("Unknown discrete fossilisation sate.\nKnown values are\n");

	printf("X	Likilhood values unchanged\n");
	printf("-	Likilhood set to zero\n\n");
	printf("Symbol     0,0   0,1   1,0   1,1\n");
	printf("0          X     -     -     -\n");
	printf("1          -     X     -     -\n");
	printf("2          -     -     X     -\n");
	printf("3          -     -     -     X\n");
	printf("\n");
	printf("10         X     X     -     -\n");
	printf("11         X     -     X     -\n");
	printf("12         X     -     -     X\n");
	printf("13         -     X     X     -\n");
	printf("14         -     X     -     X\n");
	printf("15         -     -     X     X\n");
	printf("\n");
	printf("20         X     X     X     -\n");
	printf("21         X     X     -     X\n");
	printf("22         X     -     X     X\n");
	printf("23         -     X     X     X\n");


	exit(1);
}

/* Will have to add support for extra states (10 to 23) */
int*	DesFossilSate(char* State, OPTIONS *Opt)
{
	int FState;
	int	*Ret;

	if(IsValidInt(State) == FALSE)
	{
		printf("Could not convert %s to a valid Discrete state", State);
		exit(1);
	}

	FState = atoi(State);
	ValidDescFossileState(FState);

	Ret = (int*)SMalloc(sizeof(int));
	Ret[0] = FState;

	return Ret;
}

int		MSStateToNo(char State, char *SList, int NoS)
{
	int Index;

	for(Index=0;Index<NoS;Index++)
		if(State == SList[Index])
			return Index;

	printf("Cannot convert %c to a valid multi-state value.\n", State);
	exit(1);
	return -1;
}

int*	CrateMSFossilStateList(char *List, OPTIONS *Opt, int *No)
{
	int *Ret, Index;
	
	*No = (int)strlen(List);

	Ret = (int*)SMalloc(sizeof(int) * *No);

	for(Index=0;Index<*No;Index++)
		Ret[Index] = MSStateToNo(List[Index], Opt->Trees->SymbolList, Opt->Trees->NoStates);

	return Ret;
}

int*	FossilState(char *States, OPTIONS *Opt, int *No)
{
	if(Opt->Model == M_MULTISTATE)
		return CrateMSFossilStateList(States, Opt, No);

	*No = 1;
	return DesFossilSate(States, Opt);
}

int		GetTaxaNoFormName(char* Name, TREES* Trees, int *No)
{
	int	Index;
	
	for(Index=0;Index<Trees->NoTaxa;Index++)
	{
		if(strcmp(Name, Trees->Taxa[Index]->Name) == 0)
		{
			*No = Trees->Taxa[Index]->No;
			return TRUE;
		}
	}

	*No = -1;
	return FALSE;
}

int		ValidTaxaList(char** List, int Start, int No, OPTIONS *Opt)
{
	int		Index;
	int		OK;
	int		TaxaNo;
		
	for(Index=Start;Index<No;Index++)
	{
		OK = FALSE;
		/* Check to see if the taxa is a number */
		if(IsValidInt(List[Index]) == TRUE)
		{
			OK = TRUE;
			TaxaNo = atoi(List[Index]);

			if(GetTaxaFromID(TaxaNo, Opt->Trees->Taxa, Opt->Trees->NoTaxa) == NULL)
			{
				printf("Error: Could not convert %s to a valid taxa number.\n", List[Index]);
				exit(0);
				return FALSE;
			}
		}
		else
		{
			if(GetTaxaNoFormName(List[Index], Opt->Trees, &TaxaNo) == FALSE)
			{
				printf("Error: Could not convert %s to a valid taxa name.\n", List[Index]);
				exit(0);
				return FALSE;
			}
		}
	}
	
	return TRUE;
}
/*
TAXA*	GetTaxaFromName(char *ID, TREES* Trees)
{
	int	Index;

	for(Index=0;Index<Trees->NoTaxa;Index++)
		if(strcmp(ID, Trees->Taxa[Index]->Name) == 0)
			return Trees->Taxa[Index];

	return NULL;
}*/		  

char**	SetConFState(OPTIONS *Opt, NODETYPE NodeType, char **argv)
{
	int		Index;
	char	**Ret;

	Ret = (char**)SMalloc(sizeof(char*) * Opt->Trees->NoSites);

	for(Index=0;Index<Opt->Trees->NoSites;Index++)
	{
		if(NodeType != FOSSIL)
			Ret[Index] = StrMake(ESTDATAPOINT); 
		else
			Ret[Index] = StrMake(argv[Index]);
	}

	return Ret;
}
/*
void	AddRecNode(OPTIONS *Opt, NODETYPE NodeType, int Tokes, char *argv[])
{
	RECNODE		*RNode;
	int			Index, NoTaxa;
	int			*FStates, NoFStates;
	char**		ConFState;
	
	
	ConFState = NULL;
	FStates = NULL;

	if(OptFindRecNode(Opt, argv[1]) != NULL)
	{
		printf("Node name %s is allready is use please chose another\n", argv[1]);
		exit(0);
		return;
	}

	Index=2;

	if(NodeType == FOSSIL)
	{
		if(Opt->DataType == DISCRETE)
		{
			FStates = FossilState(argv[2], Opt, &NoFStates);
			Index++;
		}
		else
			Index += Opt->Trees->NoSites;
	}

	if(Opt->DataType == CONTINUOUS)		
	{
		ConFState = SetConFState(Opt, NodeType, &argv[Index]);
		if(NodeType == FOSSIL)
			Index += Opt->Trees->NoSites;
	}
	
	if(ValidTaxaList(argv, Index, Tokes, Opt) == FALSE)
		return;

	RNode = (RECNODE*)SMalloc(sizeof(RECNODE));

	RNode->Part				= NULL;
	RNode->FossilStates		= NULL;
	RNode->NoFossilStates	= -1;
	RNode->ConData			= NULL;
	if(Opt->DataType == CONTINUOUS)		
		RNode->ConData	= ConFState;

	RNode->Name = StrMake(argv[1]);
		
	RNode->NodeType		= NodeType;

	if(RNode->NodeType == FOSSIL)
	{
		RNode->FossilStates = FStates;
		RNode->NoFossilStates = NoFStates;
		NoTaxa = Tokes - 3;
	}
	else
		NoTaxa= Tokes - 2;

	RNode->Taxa = (TAXA**)SMalloc(sizeof(TAXA*)*NoTaxa);

	for(Index=0;Index<NoTaxa;Index++)
	{
		if(NodeType == FOSSIL)
			RNode->Taxa[Index] = GetTaxaFromNameNo(argv[Index+3], Opt->Trees);
		else
			RNode->Taxa[Index] = GetTaxaFromNameNo(argv[Index+2], Opt->Trees);
	}

	RNode->Part = CreatPart(NoTaxa);
	
	RNode->TreeNodes = (NODE*)SMalloc(sizeof(NODE)*Opt->Trees->NoTrees);

	SetRecNodes(RNode, Opt->Trees);

	Opt->RecNodeList = (RECNODE**)AddToList(&Opt->NoOfRecNodes, (void**)Opt->RecNodeList, RNode);
}
*/

RECNODE*	CreateRecNode(void)
{
	RECNODE *Ret;

	Ret = (RECNODE*)SMalloc(sizeof(RECNODE));

	Ret->FossilStates	= NULL;
	Ret->NoFossilStates	= -1;
	Ret->ConData		= NULL;
	Ret->Name			= NULL;
	Ret->Hits			= 0;
	Ret->Tag			= NULL;

	return Ret;
}

void	DumpLineError(int Tokes, char **argv)
{
	int Index;
	
	printf("\n\nError with line:\n");
	for(Index=0;Index<Tokes;Index++)
		printf("%s ", argv[Index]);
	printf("\n");
}

void	AddRecNodeCheck(OPTIONS *Opt, NODETYPE NodeType, int Tokes, char **argv)
{
	if(NodeType == MRCA || NodeType == NODEREC)
	{
		if(Tokes != 3)
		{
			DumpLineError(Tokes, argv);
			printf("Reconstructing an internal node requires a unique name and tag.\n");
			exit(1);
		}
	}

	if(NodeType == FOSSIL)
	{
		if(Tokes < 4)
		{
			DumpLineError(Tokes, argv);
			printf("The fossil command requires a unique name, a tag and the states to fossilise.\n");
			exit(1);
		}
	}
	
	if(OptFindRecNode(Opt, argv[1]) != NULL)
	{
		printf("Node name %s is allready is use please chose another\n", argv[1]);
		exit(1);
	}

	if(GetTagFromName(Opt, argv[2]) == NULL)
	{
		printf("Could not find tag %s.", argv[2]);
		exit(1);
	}
}

void	AddConAnsStatePrior(OPTIONS *Opt, int SiteNo)
{
	char *Buffer;
	PRIOR *Prior;

	Buffer = (char*)SMalloc(sizeof(char) * 128);

 	if(Opt->Model == M_CONTINUOUS_REG && SiteNo == 1)
		sprintf(Buffer, "AncState-Dep");
	else
		sprintf(Buffer, "AncState-%d", SiteNo);

	RemovePriorFormOpt(Buffer, Opt);
	Prior = CreateUniformPrior(Buffer, -100, 100);
	AddPriorToOpt(Opt, Prior);
	
	free(Buffer);
}

void	SetConAnsStatesPriors(OPTIONS *Opt, TREES *Trees)
{
	int Index;
		
	for(Index=0;Index<Trees->NoSites;Index++)
		AddConAnsStatePrior(Opt, Index+1);
}

void	AddRecNode(OPTIONS *Opt, NODETYPE NodeType, int Tokes, char **argv)
{
	RECNODE	*RNode;
	
	AddRecNodeCheck(Opt, NodeType, Tokes, argv);
	
	RNode = CreateRecNode();
	
	RNode->Name = StrMake(argv[1]);
	RNode->Tag	= GetTagFromName(Opt, argv[2]);
	RNode->NodeType	 = NodeType;
	
	if(NodeType == FOSSIL && Opt->DataType == DISCRETE)
		RNode->FossilStates = FossilState(argv[3], Opt, &RNode->NoFossilStates);

	if(Opt->DataType == CONTINUOUS)	
	{
		RNode->ConData = SetConFState(Opt, NodeType, &argv[3]);
		SetConAnsStatesPriors(Opt, Opt->Trees);
	}
		
	Opt->RecNodeList = (RECNODE**)AddToList(&Opt->NoOfRecNodes, (void**)Opt->RecNodeList, RNode);
}

void	SetEvenRoot(TREES *Trees)
{
	int		TIndex;
	int		NIndex;
	double	t;
	NODE	Root;

	for(TIndex=0;TIndex<Trees->NoTrees;TIndex++)
	{
		Root = Trees->Tree[TIndex]->Root;
		t = 0;

		for(NIndex=0;NIndex<Root->NoNodes;NIndex++)
			t += Root->NodeList[NIndex]->Length;

		t = t / (double)Root->NoNodes;
		for(NIndex=0;NIndex<Root->NoNodes;NIndex++)
			Root->NodeList[NIndex]->Length = t;
	}
}

void	SetLogFile(OPTIONS *Opt, int Tokes, char **Passed)
{
	if(Tokes != 2)
	{
		printf("LogFile take the base name to call output files, the default is the data file name.\n");
		exit(1);
	}
	
	free(Opt->BaseOutputFN);

	Opt->BaseOutputFN = StrMake(Passed[1]);
}

void	PreSet(OPTIONS *Opt, int Tokes, char **Passed)
{
	char *PSet;

	if(Tokes != 2)
	{
		printf("The PreSet command take a preset\n");
		exit(0);
	}
	
	PSet = Passed[1];

	MakeLower(PSet);

	if(strcmp(PSet , "m1p")==0)
	{
		RestrictAll(Opt, Opt->RateName[0]);
		Opt->AnalyticalP = TRUE;
		return;
	}

	printf("Unknown preset.\n");
	exit(0);
}

void	GetBasePis(OPTIONS *Opt, char* Type)
{
	MakeLower(Type);

	if(strcmp(Type, "emp") == 0 || strcmp(Type, "empirical") == 0)
	{
		Opt->PiTypes = PI_EMP;
		return;
	}

	if(strcmp(Type, "uni") == 0 || strcmp(Type, "uniform") == 0)
	{
		Opt->PiTypes = PI_UNI;
		return;
	}

	if(strcmp(Type, "none") == 0)
	{
		Opt->PiTypes = PI_NONE;
		return;
	}

	printf("The option %s, is unknown. Valid options are empirical, uni and none\n", Type);
	exit(0);
}

int		CmdVailWithDataType(OPTIONS *Opt, COMMANDS	Command)
{
#ifdef BTOCL
	if(Command == CCOVARION)
	{
		printf("Covarion is not currently supported with OpenCL\n");
		return FALSE;
	}
#endif

	if(Opt->DataType == CONTINUOUS)
	{
		if(Opt->Model != M_FATTAIL && Command == CNOSLICESAMPLESTEPS)
			return FALSE;

		if(Command == CVARRATES)
		{
			if(Opt->ModelType == MT_CONTINUOUS)
				return FALSE;

			if(Opt->Analsis == ANALMCMC)
				return TRUE;

			return FALSE;
		}

		if(Command == CFOSSIL)
		{
			if(Opt->ModelType == MT_CONTRAST)
				return FALSE;
		}

		if(Command == CDISTDATA)
		{
			if(Command == MT_CONTINUOUS)
				return FALSE;
		}
		
		if((Opt->Model != M_CONTRAST_REG) && (Command == CRJDUMMY))
			return FALSE;

		if(Opt->Model == M_CONTRAST_CORREL)
		{
			if(Command == CNODE)
				return TRUE;
		}

		if( (Command == CNODE)		||
			(Command == CADDTAXA)	||
			(Command == CDELTAXA)	||
			(Command == CCOVARION)	||
			(Command == CRES)		||
			(Command == CRESALL)	||
			(Command == CGAMMA)		||
			(Command == CRMODEL)	||
			(Command == CREVJUMP)	||
			(Command == CHYPERPRIOR)||
			(Command == CHPRJ)		||
			(Command == CHPALL)		||
			(Command == CPIS)		||
			(Command == CPRECISION) ||
			(Command == CNOSPERSITE)|| 
			(Command == CSYMMETRICAL) ||
			(Command == CADDPATTERN) ||
			Command == CRATESCALARS
			)
		{
			printf("Command %s (%s) is not valid with the current model\n", COMMANDSTRINGS[Command*2], COMMANDSTRINGS[(Command*2)+1]);
			return FALSE;
		}
	}
	else
	{
		if(
			Command == CDELTA		|| 
			Command == CLAMBDA		||
			Command == COU			||
			Command == CALPHAZERO	||
			Command == CNODEBLDATA	||
			Command == CNODEDATA	||
			Command == CRJDUMMY		||
			Command == CRJLOCALTRANSFORM ||
			Command == CDISTDATA	||
			Command == CRJLANDSCAPERATEGROUP
			)
		{
			printf("Command %s (%s) is not valid with Discrete data\n", COMMANDSTRINGS[Command*2], COMMANDSTRINGS[(Command*2)+1]);
			return FALSE;
		}
	}

	if(Opt->Analsis == ANALMCMC)
	{
		if(	Command == CCI ||
			Command == CMLTOL ||
			Command == CMLTRIES ||
			Command == CMLEVAL ||
			Command == CMLALG ||
			Command == CSETMINMAXRATE
			)
		{
			printf("Command %s (%s) is not valid with the MCMC model\n", COMMANDSTRINGS[Command*2], COMMANDSTRINGS[(Command*2)+1]);
			return FALSE;
		}
	}

	if(Opt->Analsis == ANALML)
	{
		if(
			Command ==	CITTERS		||
			Command ==	CBURNIN		||
			Command ==	CSAMPLE		||
			Command ==	CHYPERPRIOR	||
			Command ==	CHPRJ		||
			Command == CPRIOR		||
			Command ==	CHPALL		||
			Command ==	CREVJUMP	||
			Command == CMCMCMLSTART	||
			Command == CCAPRJRATES	||
			Command == CSAVEMODELS	||
			Command == CLOADMODELS	||
			Command == CSHEDULE		||
			Command == CRJDUMMY		||
			Command == CVARRATES	||
			Command == CDISTDATA	||
			Command == CNOLH		||
			Command == CSTONES		||
			Command == CCSCHED		||
			Command == CRJLANDSCAPERATEGROUP
			)
		{
			printf("Command %s (%s) is not valid with the ML model\n", COMMANDSTRINGS[Command*2], COMMANDSTRINGS[(Command*2)+1]);
			return FALSE;
		}
	}

	return TRUE;
}


int		SetConVar(int Tokes, char** Passed, int *InUse, int *Est, double *Const)
{
	double	Temp;

	if(Tokes == 1)
	{
		*Const	= -1;
		if((*InUse) == TRUE)
		{
			*InUse	= FALSE;
			*Est	= FALSE;
			
		}
		else
		{
			*InUse	= TRUE;
			*Est	= TRUE;
		}
		return TRUE;
	}

	if(IsValidDouble(Passed[1])==TRUE)
	{
		Temp = atof(Passed[1]);
		*Est	= FALSE;
		*InUse	= TRUE;
		*Const = Temp;
	}
	
	return TRUE;
}

void	SetKappa(OPTIONS *Opt, int Tokes, char **Passed)
{
	double Val;

	if(Opt->Analsis == ANALMCMC)
		RemovePriorFormOpt("Kappa", Opt);

	if(Tokes == 2)
	{
		if(IsValidDouble(Passed[1]) == FALSE)
		{
			printf("Cannot convert %s to a valid kappa value", Passed[1]);
			exit(1);
		}

		Val = atof(Passed[1]);
		if(Val < MIN_KAPPA)
		{
			printf("%f is lower than the minium kappa %f", Val, MIN_KAPPA);
			exit(1);
		}

		Opt->EstKappa = FALSE;
		Opt->UseKappa = TRUE;
		Opt->FixKappa = Val;
		return;
	}

	if(Opt->UseKappa == TRUE)
	{
		Opt->UseKappa = FALSE;
		Opt->EstKappa = FALSE;
		return;
	}

	Opt->UseKappa = TRUE;
	Opt->EstKappa = TRUE;
	
	SetLocalTransformPrior(Opt, VR_KAPPA);
}

void	SetLambda(OPTIONS *Opt, int Tokes, char **Passed)
{
	double Val;

	if(Opt->Analsis == ANALMCMC)
		RemovePriorFormOpt("Lambda", Opt);

	if(Tokes == 2)
	{
		if(IsValidDouble(Passed[1]) == FALSE)
		{
			printf("Cannot convert %s to a valid lambda value", Passed[1]);
			exit(1);
		}

		Val = atof(Passed[1]);
		if(Val < MIN_LAMBDA)
		{
			printf("%f is lower than the minium lambda %12.12f", Val, MIN_LAMBDA);
			exit(1);
		}

		Opt->EstLambda = FALSE;
		Opt->UseLambda = TRUE;
		Opt->FixLambda = Val;
		return;
	}

	if(Opt->UseLambda == TRUE)
	{
		Opt->UseLambda = FALSE;
		Opt->EstLambda = FALSE;
		return;
	}

	Opt->UseLambda = TRUE;
	Opt->EstLambda = TRUE;

	SetLocalTransformPrior(Opt, VR_LAMBDA);
}

void	SetDelta(OPTIONS *Opt, int Tokes, char **Passed)
{
	double Val;

	if(Opt->Analsis == ANALMCMC)
		RemovePriorFormOpt("Detla", Opt);

	if(Tokes == 2)
	{
		if(IsValidDouble(Passed[1]) == FALSE)
		{
			printf("Cannot convert %s to a valid detla value", Passed[1]);
			exit(1);
		}

		Val = atof(Passed[1]);
		if(Val < MIN_DELTA)
		{
			printf("%f is lower than the minium detla %f", Val, MIN_DELTA);
			exit(1);
		}

		Opt->EstDelta = FALSE;
		Opt->UseDelta = TRUE;
		Opt->FixDelta = Val;
		return;
	}

	if(Opt->UseDelta == TRUE)
	{
		Opt->UseDelta = FALSE;
		Opt->EstData = FALSE;
		return;
	}

	Opt->UseDelta = TRUE;
	Opt->EstDelta = TRUE;

	SetLocalTransformPrior(Opt, VR_DELTA);
}

void	SetOU(OPTIONS *Opt, int Tokes, char **Passed)
{
	double Val;

	if(Opt->Analsis == ANALMCMC)
		RemovePriorFormOpt("OU", Opt);

	if(Tokes == 2)
	{
		if(IsValidDouble(Passed[1]) == FALSE)
		{
			printf("Cannot convert %s to a valid OU value", Passed[1]);
			exit(1);
		}

		Val = atof(Passed[1]);
		if(Val < MIN_OU)
		{
			printf("%f is lower than the minium OU %f", Val, MIN_OU);
			exit(1);
		}

		Opt->EstOU = FALSE;
		Opt->UseOU = TRUE;
		Opt->FixOU = Val;
		return;
	}

	if(Opt->UseOU == TRUE)
	{
		Opt->UseOU = FALSE;
		Opt->EstOU = FALSE;
		return;
	}

	Opt->UseOU = TRUE;
	Opt->EstOU = TRUE;

	SetLocalTransformPrior(Opt, VR_OU);
}

void	ExcludeTaxa(OPTIONS *Opt, int Tokes, char **Passed)
{
	int		Index;
	char	*Name;
	TAXA	*Taxa;

	for(Index=0;Index<Tokes;Index++)
	{
		Taxa = GetTaxaFromName(Passed[Index], Opt->Trees->Taxa, Opt->Trees->NoTaxa);

		if(Taxa == NULL)
		{
			printf("Exclude Taxa, invalid taxa name %s\n", Passed[Index]);
			exit(0);
		}

		CheckDelTaxa(Opt, Opt->Trees, Passed[Index]);

		Name = Taxa->Name;

		RemoveTaxa(Opt->Trees, Name);
	}

	SetParts(Opt->Trees);
}

void	RemoveRatePriors(OPTIONS *Opt)
{
	int Index;

	for(Index=0;Index<Opt->NoOfRates;Index++)
		RemovePriorFormOpt(Opt->RateName[Index], Opt);
}

void	SetRJMCMC(OPTIONS *Opt, int Tokes, char** Passed)
{
	PRIOR		*Prior;

	RemoveRatePriors(Opt);
	RemovePriorFormOpt("RJRates", Opt);

	Opt->UseRJMCMC = TRUE;

	Prior = CreatePriorFromStr("RJRates", Tokes, Passed);
	AddPriorToOpt(Opt, Prior);

	

	Opt->UseRJMCMC = TRUE;
}

void	SetRJMCMCHP(OPTIONS *Opt, int Tokes, char** Passed)
{
	PRIOR		*Prior;

	RemoveRatePriors(Opt);

	RemovePriorFormOpt("RJRates", Opt);

	Opt->UseRJMCMC = TRUE;

	Prior = CreateHyperPriorFromStr("RJRates", Tokes, Passed);
	AddPriorToOpt(Opt, Prior);

	Opt->UseRJMCMC = TRUE;
}



void	SetGamma(OPTIONS *Opt, char** Passed, int Tokes)
{
	int		GammaCats;
	double	Value;
	PRIOR	*Prior;

	if(Opt->Analsis == ANALMCMC)
		RemovePriorFormOpt("Gamma", Opt);

	if(Tokes == 1)
	{
		if(Opt->UseGamma == TRUE)
		{
			Opt->UseGamma	= FALSE;
			Opt->FixGamma	= -1;
			Opt->EstGamma	= FALSE;
			return;
		}
		else
		{
			printf("The Gamma command take 0, 1 or 2 parameters.\n 0 to turn Gamma off\n1 to estermate Gamma\n 2 to fix it to a constant\n");
			exit(0);
		}
	}

	if(IsValidInt(Passed[1]) == FALSE)
	{
		printf("Could not convert %s to a valid number of categories to divide the gamma disruption up into.\n", Passed[1]);
		exit(0);
	}

	GammaCats = atoi(Passed[1]);

	if(GammaCats < 2 || GammaCats > 8)
	{
		printf("The number of gamma catergoires must be grater than 1 and less than 8\n");
		exit(0);
	}

	if(Tokes == 2)
	{
		Opt->UseGamma = TRUE;
		Opt->EstGamma = TRUE;
		Opt->FixGamma = -1;
	
		Opt->GammaCats = GammaCats;

		if(Opt->Analsis == ANALMCMC)
		{
			Prior	= CreateUniformPrior("Gamma", MIN_GAMMA, MAX_GAMMA);
			AddPriorToOpt(Opt, Prior);
		}
		
		return;
	}

	if(Tokes == 3)
	{
		if(IsValidDouble(Passed[2]) == FALSE)
		{
			printf("Could not convert %s to a valid gamma shapre parmiter\n", Passed[2]);
			exit(0);
		}
		Value = atof(Passed[2]);

		if((Value < MIN_GAMMA) || (Value > MAX_GAMMA))
		{
			printf("Gamma shape parmiter must be grater than %f and less than %f\n", (double)MIN_GAMMA, (double)MAX_GAMMA);
			return;
		}


		Opt->UseGamma = TRUE;
		Opt->FixGamma = Value;
		Opt->EstGamma = FALSE;

		Opt->GammaCats= GammaCats;
		
		return;
	}

	printf("The Gamma command take 0, 1 or 2 parameters.\n 0 to turn Gamma off\n1 to estermate Gamma\n 2 to fix it to a constant\n");
	exit(0);
}

void	SetCI(OPTIONS *Opt, char *Rate)
{
	int	Index;

	if(Opt->FindCF == TRUE)
	{
		Opt->FindCF = FALSE;
		Opt->CFRate = NULL;
		return;
	}

	MakeLower(Rate);

	for(Index=0;Index<Opt->NoOfRates;Index++)
	{
		
	}
}


void	FreeRecNode(RECNODE *R, int NoSites)
{
	int Index;

	if(R->ConData != NULL)
	{
		for(Index=0;Index<NoSites;Index++)
			free(R->ConData[Index]);
		free(R->ConData);
	}

	free(R->Name);

	if(R->FossilStates != NULL)
		free(R->FossilStates);

	free(R);
}

void	FreeRecNodes(OPTIONS *Opt, int NoSites)
{
	int	Index;
	RECNODE	*R;

	if(Opt->RecNodeList == NULL)
		return;

	for(Index=0;Index<Opt->NoOfRecNodes;Index++)
	{
		R = Opt->RecNodeList[Index];
		FreeRecNode(R, NoSites);
	}
	
	free(Opt->RecNodeList);
	Opt->NoOfRecNodes = 0;
	Opt->RecNodeList = NULL;
}


void	SetNOSPerSiteOpt(OPTIONS *Opt)
{
	if(Opt->NOSPerSite == FALSE)
	{
		Opt->NOSPerSite = TRUE;
		Opt->AnalyticalP = TRUE;
		RestrictAll(Opt, Opt->RateName[0]);
	}
	else
	{
		Opt->NOSPerSite = FALSE;
		Opt->AnalyticalP = FALSE;
	}
}

void	OptSetSeed(OPTIONS *Opt, char	*CSeed)
{
	if(IsValidInt(CSeed) == FALSE)
	{
		printf("%s is not a valid Unsigned integer.\n", CSeed);
		return;
	}

	sscanf(CSeed, "%lu", &Opt->Seed);

//	Opt->Seed = atoi(CSeed);
}


void	SetEqualTrees(OPTIONS *Opt, int Tokes, char **Passed)
{
	int NoTreeBI;
	char *TreeBI;

	if(Tokes == 1)
	{
		Opt->UseEqualTrees = FALSE;
		return;
	}
	
	if(Tokes != 2)
	{
		printf("Equal trees takes no parameters to toggle off or 1 parameter, tree specific burn-in.\n");
		return;
	}

	TreeBI = Passed[1];

	if(IsValidInt(TreeBI) == FALSE)
	{
		printf("Cound not convert %s to a valid tree specific burn-in number.\n", TreeBI);
		return;
	}

	NoTreeBI = atoi(TreeBI);
	if(NoTreeBI < 0)
	{
		printf("Tree specific burn-in must be greater than 0\n");
		return;
	}

	Opt->UseEqualTrees = TRUE;
	Opt->ETreeBI = NoTreeBI;
}

void	SetPrecision(OPTIONS *Opt, char *Token)
{
	int Pre;

	if(IsValidInt(Token) == FALSE)
	{
		printf("%s is not a valid precision. Precision must be an integer >= 64.\n", Token);
		return;
	}

	Pre = atoi(Token);
	if(Pre <= sizeof(double) * 8)
	{
		printf("%s is not a valid precision. Precision must be an integer >= 64.\n", Token);
		return;
	}

	Opt->Precision = Pre;
}

void	SetCores(OPTIONS *Opt, int Tokes, char** Passed)
{
	int Cores;


#ifndef CLIK_P
#ifndef OPENMP_THR
	printf("Cores is not valid with this build. please use the threaded build.");
	return;
#endif
#endif

	if(Tokes != 2)
	{
		printf("Cores takes the number of cores to use.\n");
		return;
	}

	if(IsValidInt(Passed[1]) == FALSE)
	{
		printf("Could not covert %s to a valid number of cores\n", Passed[0]);
		return;
	}

	Cores = atoi(Passed[1]);
	if(Cores < 1)
	{
		printf("the number of course must be >= 0\n");
		return;
	}

	Opt->Cores = Cores;
}

void	ResRateNo(OPTIONS *Opt, int From, int To)
{
	Opt->ResTypes[From] = RESRATE;
	Opt->ResNo[From] = To;
}

void	SetMSSymmetrical(OPTIONS *Opt)
{
	int		x, y, NOS, From, To;
	char	FromS[4], ToS[4];
		
	NOS = Opt->Trees->NoStates;

	for(x=0;x<NOS;x++)
	{
		for(y=0;y<x;y++)
		{
			if(x != y)
			{
				sprintf(&FromS[0], "q%c%c", Opt->Trees->SymbolList[x], Opt->Trees->SymbolList[y]);
				sprintf(&ToS[0], "q%c%c", Opt->Trees->SymbolList[y], Opt->Trees->SymbolList[x]);
				
				From	= StrToRate(Opt, &FromS[0]);
				To		= StrToRate(Opt, &ToS[0]);
				
				ResRateNo(Opt, From, To);
			}
		}
	}	
}

void	SetSymmetrical(OPTIONS *Opt)
{
	UnRestictAll(Opt);

	if(Opt->Model == M_MULTISTATE)
	{
		SetMSSymmetrical(Opt);
	}

	if(Opt->Model == M_DISC_INDEP)
	{
		/* Beta 1 = Alpha 1 */
		ResRateNo(Opt, 1, 0);

		/* Beta 2 = Alpha 2 */
		ResRateNo(Opt, 3, 2);
	}

	if(Opt->Model == M_DISC_DEP)
	{
		/* q21 = q12 */
		ResRateNo(Opt, 2, 0);

		/* q31 = q13 */
		ResRateNo(Opt, 4, 1);
		
		/* q42 = q24 */ 
		ResRateNo(Opt, 6, 3);
		
		/* q43 = q34 */
		ResRateNo(Opt, 7, 5);
	}
}

void	SetMCMCMLStart(OPTIONS *Opt)
{
	if(Opt->MCMCMLStart == TRUE)
		Opt->MCMCMLStart = FALSE;
	else
		Opt->MCMCMLStart = TRUE;	
}

void	SetTestCorrel(OPTIONS *Opt)
{
	if(Opt->TestCorrel == FALSE)
		Opt->TestCorrel = TRUE;
	else
		Opt->TestCorrel = FALSE;
}

void	CapRJRatesNo(OPTIONS *Opt, int Tokes, char **Passed)
{
	int Cap;

	if(Tokes > 2)
	{
		printf("CapRJRates takes the maximum number of RJ Rates to alow.\n");
		return;
	}

	if(Tokes == 1)
	{
		Opt->CapRJRatesNo = -1;
		return;
	}

	if(IsValidInt(Passed[1]) == FALSE)
	{
		printf("%s not a valid cap for RJ Rates.\n", Passed[1]);
		return;
	}

	Cap = atoi(Passed[1]);
	if(Cap < 1)
	{
		printf("%s not a valid cap for RJ Rates.\n", Passed[1]);
		return;
	}

	Opt->CapRJRatesNo = Cap;
}

void	SetSaveModels(OPTIONS *Opt, int Tokes, char **Passed)
{
	if(Tokes > 2)
	{
		printf("SaveModels command take zero paramiters to turn it off, or a file name to save the models to.\n");
		return;
	}

	if(Opt->SaveModelsFN != NULL)
		free(Opt->SaveModelsFN);

	if(Tokes == 1)
	{
		Opt->SaveModels = FALSE;
		return;
	}

	Opt->SaveModels = TRUE;
	Opt->SaveModelsFN = StrMake(Passed[1]);
}

void	SetLoadModels(OPTIONS *Opt, int Tokes, char **Passed)
{
	if(Tokes > 2)
	{
		printf("LoadModels command take zero paramiters to turn it off, or a file name to save the models to.\n");
		return;
	}

	if(Opt->LoadModelsFN != NULL)
		free(Opt->LoadModelsFN);

	if(Tokes == 1)
	{
		Opt->LoadModels = FALSE;
		return;
	}

	Opt->LoadModels = TRUE;
	Opt->LoadModelsFN = StrMake(Passed[1]);
}


void	LineAddErr(TREES *Trees, char *Line)
{
	char *Buffer;
	char **Passed;
	int		Tokes, ID;
	double	Err;

	Buffer = StrMake(Line);
	Passed = (char**)SMalloc(sizeof(char*) * strlen(Line));

	Tokes = MakeArgv(Buffer, Passed, (int)strlen(Line));

	if(Tokes == 0)
		return;

	if(Tokes != 2)
	{
		printf("AddErr: %s is not a valid line.\n", Line);
		printf("AddErr: Each line should contain a  taxa name and standard error.\n");
		return;
	}

	if(GetTaxaNoFormName(Passed[0], Trees, &ID) == FALSE)
	{
		printf("AddErr: %s is an invalid taxa name\n", Passed[0]);
		return;
	}

	if(IsValidDouble(Passed[1]) == FALSE)
	{
		printf("AddErr: Could not convert %s to a valid error\n", Passed[1]);
		return;
	}

	Err = atof(Passed[1]);

	AddTaxaErr(Trees, ID, Err);

	free(Passed);
	free(Buffer);	
}

void	LoadAddErr(OPTIONS *Opt, int Tokes, char **argv)
{
	TEXTFILE *TF;
	int	Index;
	char *FName;

	if(Tokes != 2)
	{
		printf("AddErr requires one parameter, a file with taxa names and error.\n");
		printf("File names cannot contain spaces.\n");
		exit(0);
	}

	FName = argv[1];

	TF = LoadTextFile(FName, FALSE);

	for(Index=0;Index<TF->NoOfLines;Index++)
		LineAddErr(Opt->Trees, TF->Data[Index]);
	
	FreeTextFile(TF);
}

void	SetSteppingstone(OPTIONS *Opt, char **Passed, int Tokes)
{
	int K, Sample;
	double	Alpha, Beta;

	if(Tokes == 1)
	{
		if(Opt->Stones != NULL)
			FreeStones(Opt->Stones);
		Opt->Stones = NULL;
	}

	Beta = 1.0;
	Alpha= 0.4;

	if((Tokes != 3) && (Tokes != 5))
	{
		printf("Stones takes the number of stones and length to sample each stone.\n");
		printf("Or\nStones takes the number of stones, length to sample each stone, alpha and beta values of each stone.\n");
		return;
	}
	
	if(IsValidInt(Passed[1]) == FALSE)
	{
		printf("Stones: could not convert %s to a valid number of stones.\n", Passed[1]);
		return;
	}

	K = atoi(Passed[1]);
	if(K < 1)
	{
		printf("Stones: could not convert %s to a valid number of stones.\n", Passed[1]);
		return;
	}

	if(IsValidInt(Passed[2]) == FALSE)
	{
		printf("Stones: could not convert %s to a valid number of itterations per stone.\n", Passed[2]);
		return;
	}

	Sample = atoi(Passed[2]);
	if(Sample < 1)
	{
		printf("Stones: could not convert %s to a valid number of iterations per stone \n", Passed[2]);
		return;
	}	

	if(Tokes == 5)
	{
		if(IsValidDouble(Passed[3]) == FALSE)
		{
			printf("Stones: could not convert %s to a valid alpha.\n", Passed[3]);
			return;
		}

		Alpha = atof(Passed[3]);
		if(Alpha < 0)
		{
			printf("Stones: could not convert %s to a valid alpha.\n", Passed[3]);
			return;
		}		
		
		if(IsValidDouble(Passed[4]) == FALSE)
		{
			printf("Stones: could not convert %s to a valid beta.\n", Passed[4]);
			return;
		}

		Beta = atof(Passed[4]);
		if(Beta < 0)
		{
			printf("Stones: could not convert %s to a valid beta.\n", Passed[4]);
			return;
		}		
	}

	if(Opt->Stones != NULL)
		FreeStones(Opt->Stones);

	Opt->Stones = CratesStones(K,  Sample,  Alpha,  Beta);
}

void	SetRJDummy(OPTIONS *Opt, char **Passed, int Tokes)
{
	if(Opt->Trees->NoTrees != 1)
	{
		printf("RJ Dummy coding only works on a singel tree.\n");
		return;
	}

	if(Opt->Trees->NoSites > 2)
	{
		printf("RJ Dummy coding does not currently work with multiple regressions.\n");
		return;
	}

	if(Opt->RJDummy == TRUE)
		Opt->RJDummy = FALSE;
	else
		Opt->RJDummy = TRUE;
}

void	SetScaleTree(OPTIONS *Opt, char **Passed, int Tokes)
{
	double S;

	if(Tokes == 1)
	{
		Opt->ScaleTrees = FindTreeNormalise(Opt->Trees);
		return;
	}

	if(Tokes > 2)
	{
		printf("Scale Tree take a scale to multiple the branch lengths by.\n");
		return;
	}

	if(IsValidDouble(Passed[1]) == FALSE)
	{
		printf("Could not convert %s to a valid scalar.\n", Passed[1]);
		return;
	}

	S = atof(Passed[1]);
	if(S <= 0.0)
	{
		printf("Scalar has to be greater than zero.\n");
		return;
	}

	Opt->ScaleTrees = S;
}


int		ValidRJLocalScalarModel(OPTIONS *Opt, char **Passed, int Tokes)
{
	int Err;

	if(Tokes != 2)
	{
		printf("RJ Local Scalar take a scalar names (kappa, lambda, delta, OU).\n");
		return FALSE;
	}

	if(Opt->ModelType != MT_CONTRAST || Opt->Analsis == ANALML)
	{
		printf("RJ Local Scalar is only valid with MCMC and a contrast model.\n");
		return FALSE;
	}

	NameToRJLocalType(Passed[1], &Err);

	if(Err == TRUE)
	{
		printf("invalid transform name, valid scalars are (kappa, lambda, delta, OU).\n");
		return FALSE;
	}

	return TRUE;
}

void	SetLocalTransformPrior(OPTIONS *Opt, TRANSFORM_TYPE	Type)
{
	PRIOR *Prior;

	if(Opt->Analsis == ANALML)
		return;

	if(Type == VR_KAPPA)
	{
		RemovePriorFormOpt("Kappa", Opt);
		Prior = CreateUniformPrior("Kappa", MIN_KAPPA, MAX_KAPPA);
		AddPriorToOpt(Opt, Prior);
	}
	
	if(Type == VR_LAMBDA)
	{
		RemovePriorFormOpt("Lambda", Opt);
		Prior = CreateUniformPrior("Lambda", MIN_LAMBDA, MAX_LAMBDA);
		AddPriorToOpt(Opt, Prior);
	}

	if(Type == VR_DELTA)
	{
		RemovePriorFormOpt("Delta", Opt);
		Prior = CreateUniformPrior("Delta", MIN_DELTA, MAX_DELTA);
		AddPriorToOpt(Opt, Prior);
	}

	if(Type == VR_OU)
	{
		RemovePriorFormOpt("OU", Opt);
		Prior = CreateExpPrior("OU", 1.0);
		AddPriorToOpt(Opt, Prior);
	}

	if(Type == VR_BL)
	{
		RemovePriorFormOpt("VRBL", Opt);
		Prior = CreateSGammaPrior("VRBL", VARRATES_ALPHA, VARRATES_BETA);
		AddPriorToOpt(Opt, Prior);
	}

	if(Type == VR_NODE)
	{
		RemovePriorFormOpt("VRNode", Opt);
		Prior = CreateSGammaPrior("VRNode", VARRATES_ALPHA, VARRATES_BETA);
		AddPriorToOpt(Opt, Prior);
	}

	if(Type == VR_LS_BL)
	{
		RemovePriorFormOpt("VR_LS_BL", Opt);
//		Prior = CreateUniformPrior("VR_LS_BL", -5, 5);
		Prior = CreateNormalPrior("VR_LS_BL", 0, 2.0);
		AddPriorToOpt(Opt, Prior);
	}
}

void	SetRJLocalTransform(OPTIONS *Opt, char **Passed, int Tokes)
{
	TRANSFORM_TYPE	Type;
	int Err;

	if(ValidRJLocalScalarModel(Opt, Passed, Tokes) == FALSE)
		exit(1);

	Type = NameToRJLocalType(Passed[1], &Err);

	Opt->UseRJLocalScalar[Type]	= TRUE;

	SetLocalTransformPrior(Opt, Type);
	Opt->SaveTrees = TRUE;
}

void	SetFatTailNormal(OPTIONS *Opt)
{
	if(Opt->ModelType != MT_FATTAIL)
	{
		printf("Fat Tail Normal can only be used with Fat Tail models.\n");
		return;
	}

	if(Opt->FatTailNormal == FALSE)
		Opt->FatTailNormal = TRUE;
	else
		Opt->FatTailNormal = FALSE;
}



double	ValidLocalRateScalar(char *Str)
{
	double Ret;

	if(IsValidDouble(Str) == FALSE)
	{
		printf("Cannot convert %s to a valid scalar.\n", Str);
		exit(0);
	}

	Ret = atof(Str);
	if(Ret < 0)
	{
		printf("Scalars must be > 0.\n");
		exit(0);
	}

	return Ret;
}

TAG**	GetTagListFromNames(OPTIONS *Opt, char **NList, int Tokes, int *NoTags)
{
	TAG **Ret;
	TAG *Tag;
	int Index;

	Ret = (TAG**)SMalloc(sizeof(TAG*) * Tokes);

	*NoTags = 0;

	for(Index=0;Index<Tokes;Index++)
	{
		Tag = GetTagFromName(Opt, NList[Index]);
		if(Tag != NULL)
		{
			Ret[*NoTags] = Tag;
			(*NoTags)++;
		}
		else
			return Ret;
	}

	printf("LocalTransform takes a name, a list of tags, A transform type (node, bl, kappa, lambda, delta, OU) and an optional fixed scalar");
	exit(0);

	return NULL;
}

void	AddLocalTransform(OPTIONS *Opt, int Tokes, char **Passed)
{
	TRANSFORM_TYPE		Type;
	double				Scale;
	int					Est, Pos, NoTags;
	char				*Name;
	TAG					**Tags;
	LOCAL_TRANSFORM		*UVR;

	if(Tokes < 4)
	{
		printf("LocalTransform takes a name, a list of tags, A transform type (node, branch, kappa, lambda, delta, OU) and an optional fixed scalar");
		exit(1);
	}

	Name = Passed[1];

	Scale = -1;
	Est = TRUE;

	NoTags = Tokes;

	Tags = GetTagListFromNames(Opt, &Passed[2], Tokes-2, &NoTags);

	Pos = 2 + NoTags;

	Type = StrToVarRatesType(Passed[Pos++]);

	if(Tokes != Pos)
	{
		Scale = ValidLocalRateScalar(Passed[Pos]);
		Est = FALSE;
	}

	
	UVR = CreateLocalTransforms(Name, Tags, NoTags, Type, Est, Scale);
	Opt->LocalTransforms = (LOCAL_TRANSFORM**)AddToList(&Opt->NoLocalTransforms, (void**)Opt->LocalTransforms, UVR);

	SetLocalTransformPrior(Opt, Type);

	free(Tags);

	Opt->SaveTrees = TRUE;
}

void	SetDistData(OPTIONS *Opt, int Tokes, char **Passed)
{
	if(Tokes != 2)
	{
		printf("Data Dist takes a data distribution file.\n");
		exit(0);
	}

	if(Opt->DistData != NULL)
		FreeDistData(Opt->DistData);

	Opt->DistData = LoadDistData(Opt, Opt->Trees, Passed[1]);
	Opt->UseDistData = TRUE;
}

int		CompCShed(const void *CS1, const void *CS2)
{
	CUSTOM_SCHEDULE **S1, **S2;
	long long Diff;

	S1 = (CUSTOM_SCHEDULE**)CS1;
	S2 = (CUSTOM_SCHEDULE**)CS2;

	Diff  = (*S1)->Iteration - (*S2)->Iteration;

	if(Diff == 0)
		return 0;

	if(Diff < 0)
		return -1;

	return 1;
}

void	AddTimeSlicePriors(TIME_SLICE *TS, OPTIONS *Opt)
{
	PRIOR *Prior;

	if(TS->FixedTime == FALSE)
	{
		RemovePriorFormOpt("TimeSlice-Time", Opt);
		Prior = CreateUniformPrior("TimeSlice-Time", 0.0, 1.0);
		AddPriorToOpt(Opt, Prior);
	}

	if(TS->FixedScale == FALSE)
	{
		RemovePriorFormOpt("TimeSlice-Scale", Opt);
		Prior = CreateSGammaPrior("TimeSlice-Scale", VARRATES_ALPHA, VARRATES_BETA);
		AddPriorToOpt(Opt, Prior);
	}

}


void	OptAddTimeSlice(OPTIONS *Opt, int Tokes, char **Passed)
{
	TIME_SLICE *TS;

	double Time, Scale;
	char *Name;

	Time = Scale = -1;

	if(Tokes != 2 && Tokes != 3 && Tokes != 4)
	{ 
		printf("AddTimeSlice requires a name and an optional time point and scale value.");
		exit(1);
	}

	Name = Passed[1];

	if(GetTimeSlice(Opt->TimeSlices, Name)  != NULL)
	{
		printf("Time slices %s is allready defined.\n", Name);
		exit(1);
	}

	if(Tokes >= 3)
	{
		if(IsValidDouble(Passed[2]) == FALSE)
		{
			printf("Cannot convert %s to a valid time point, must be between 0-1.\n", Passed[2]);
			exit(1);
		}
		Time = atof(Passed[2]);

		if(Time < 0 || Time > 1.0)
		{
			printf("Cannot convert %s to a valid time point, must be between 0-1.\n", Passed[2]);
			exit(1);
		}
	}

	if(Tokes == 4)
	{
		if(IsValidDouble(Passed[3]) == FALSE)
		{
			printf("Cannot convert %s to a valid scalar, must be a number >0.\n", Passed[2]);
			exit(1);
		}

		Scale = atof(Passed[3]);
		if(Scale < 0)
		{
			printf("Cannot convert %s to a valid scalar, must be a number >0.\n", Passed[2]);
			exit(1);
		}
	}
	
	TS = AddTimeSlice(Opt->TimeSlices, Name, Time, Scale);
	Opt->SaveTrees = TRUE;

	if(Opt->Analsis == ANALMCMC)
		AddTimeSlicePriors(TS, Opt);
}

void	RemoveRateNamePriors(OPTIONS *Opt)
{
	int Index;

	for(Index=0;Index<Opt->NoOfRates;Index++)
		RemovePriorFormOpt(Opt->RateName[Index], Opt);
}

void	AddRateNamePriors(OPTIONS *Opt)
{
	int Index;
	PRIOR *Prior;

	for(Index=0;Index<Opt->NoOfRates;Index++)
	{
		Prior = CreateUniformPrior(Opt->RateName[Index], 0, 100);
		AddPriorToOpt(Opt, Prior);
	}
}

void	OptAddPattern(OPTIONS *Opt, int Tokes, char **Passed)
{
	int Index;
	TAG *Tag;

	if(Tokes < 3)
	{
		printf("AddPattern takes a pattern name and one or more tags.\n");
		exit(0);
	}
	
	for(Index=2;Index<Tokes;Index++)
	{
		Tag = GetTagFromName(Opt, Passed[Index]);
		if(Tag == NULL)
		{
			printf("%s is not a valid tag.\n", Passed[Index]);
			exit(1);
		}
	}

	AddPattern(Opt, Passed[1], Tokes-2, &Passed[2]);

	if(Opt->Analsis == ANALMCMC)
		RemoveRateNamePriors(Opt);

	SetPatternRateNames(Opt);
	AllocRestictions(Opt);

	if(Opt->Analsis == ANALMCMC)
		AddRateNamePriors(Opt);
}

void SetOptCustomSchedule(OPTIONS *Opt, int Tokes, char **Passed)
{
	int Index;
	CUSTOM_SCHEDULE *CShed;
	double Freq;

	if(!(Tokes == 2 || Tokes == NO_SCHEDULE_OPT + 2))
	{
		printf("Custom Schedule requires in iteration number and a vector of %d frequencies to schedule each operation, or just an iteration number to set the default schedule.\n", NO_SCHEDULE_OPT);
		printf("Operators are :\n");
		for(Index=0;Index<NO_SCHEDULE_OPT;Index++)
			printf("\t%d\t%s\n", Index, SHEDOP[Index]);
		
		exit(0);
	}

	CShed = AllocCustomSchedule();

	sscanf(Passed[1], "%lld", &CShed->Iteration);
	if(CShed->Iteration < 0)
	{
		printf("Iteration number must be >0.\n");
		exit(1);
	}

	if(Tokes == 2)
		CShed->Default = TRUE;
	else
	{
		for(Index=0;Index<NO_SCHEDULE_OPT;Index++)
		{
			if(IsValidDouble(Passed[Index+2]) == FALSE)
			{
				printf("Cannot covert %s to a valid frequncy.\n", Passed[Index+2]);
				exit(0);
			}

			Freq = atof(Passed[Index+2]);

			if(Freq < 0)
			{
				printf("Frequncyies my be > 0, value is %s.\n", Passed[Index+2]);
				exit(0);
			}

			CShed->Frequencies[Index] = Freq;
		}

		NormaliseVector(CShed->Frequencies, NO_SCHEDULE_OPT);
	}

	Opt->CShedList =(CUSTOM_SCHEDULE**) AddToList(&Opt->NoCShed, (void**)Opt->CShedList, (void*)CShed);
	qsort(Opt->CShedList, Opt->NoCShed, sizeof(CUSTOM_SCHEDULE*), CompCShed);
}

void	SetNoLh(OPTIONS *Opt)
{
	if(Opt->NoLh == TRUE)
		Opt->NoLh = FALSE;
	else
		Opt->NoLh = TRUE;
}

void	SetSaveTrees(OPTIONS *Opt)
{
	if(Opt->SaveTrees == TRUE)
		Opt->SaveTrees = FALSE;
	else
		Opt->SaveTrees = TRUE;
}

void	SetBurnIn(OPTIONS *Opt, int Tokes, char **Passed)
{
	long long TBurnIn;

	if(Tokes != 2)
	{
		printf("The Burn In command take the number of iterations before sampling.\n");
		return;
	}
	
	if(IsValidInt(Passed[1]) == FALSE)
	{
		printf("%s is not a valid number of iterations for Burn In.\n", Passed[1]);
		return;
	}

	sscanf(Passed[1], "%lld", &TBurnIn);
	if(TBurnIn < 0)
	{
		printf("Burn In period must be greater than 0.\n");
		return;
	}

	Opt->BurnIn = TBurnIn;
}


void	SetItters(OPTIONS *Opt, int Tokes, char **Passed)
{
	long long TItter;

	if(Tokes != 2)
	{
		printf("The iterations command take the number of iterations to run the chain for.\n");
		return;
	}
	
	if(IsValidInt(Passed[1]) == FALSE)
	{
		printf("%s is not a valid number of iterations.\n", Passed[1]);
		return;
	}

	sscanf(Passed[1], "%lld", &TItter);
	if(TItter < -1) 
	{
		printf("Number of itterations must be 0 or greater, or -1 to run an infinite chain.\n");
		return;
	}

	Opt->Itters = TItter;
}

void	SetVarRatesOpt(OPTIONS *Opt)
{
	if(	Opt->UseRJLocalScalar[VR_BL]	== TRUE ||
		Opt->UseRJLocalScalar[VR_NODE]	== TRUE)
	{
		printf("Node or branch RJ allready set.\n");
		exit(1);
	}

	SetLocalTransformPrior(Opt, VR_BL);
	SetLocalTransformPrior(Opt, VR_NODE);

	Opt->UseRJLocalScalar[VR_BL]	= TRUE;
	Opt->UseRJLocalScalar[VR_NODE]	= TRUE;
	
	Opt->SaveTrees = TRUE;
}


void	SetLandscape(OPTIONS *Opt, int Tokes, char **Passed)
{
	if(Tokes != 1)
	{
		printf("Landscape Does not take any parameters.\n");
		exit(1);
	}

	if(	Opt->UseRJLocalScalar[VR_NODE]	== TRUE ||
		Opt->UseRJLocalScalar[VR_LS_BL]	== TRUE)
	{
		printf("Node or Landscape BL allready set.\n");
		exit(1);
	}

	SetLocalTransformPrior(Opt, VR_LS_BL);
	SetLocalTransformPrior(Opt, VR_NODE);

	Opt->UseRJLocalScalar[VR_LS_BL]	= TRUE;
	Opt->UseRJLocalScalar[VR_NODE]	= TRUE;

	Opt->SaveTrees = TRUE;
}

void	SetRateScalars(OPTIONS *Opt, int Tokes, char **Passed)
{
	int Index;
	double Val;

	if(Tokes == 1)
	{
		for(Index=0;Index<Opt->NoOfRates;Index++)
			Opt->RateScalars[Index] = 1.0;
		return;
	}

	if(Tokes - 1 != Opt->NoOfRates)
	{
		printf("RateScalars takes takes %d rate scalars.\n", Opt->NoOfRates);
		exit(1);
	}

	for(Index=0;Index<Opt->NoOfRates;Index++)
	{
		if(IsValidDouble(Passed[Index+1]) == FALSE)
		{
			printf("Cannot convert %s to a valid rate.\n", Passed[Index+1]);
			exit(1);
		}

		Val = atof(Passed[Index+1]);
		if(Val < 0.0)
		{
			printf("Rates sclars (%s) must be > 0.\n", Passed[Index+1]);
			exit(1);
		}

		Opt->RateScalars[Index] = Val;
	}
}

void	SetRJLandscapeRateGroup(OPTIONS *Opt, int Tokes, char **Passed)
{
	PRIOR *Prior;

	RemovePriorFormOpt("RJ_Landscape_Rate_Group", Opt);

	if(Opt->UseRJLandscapeRateGroup == FALSE)
	{
		Opt->UseRJLandscapeRateGroup = TRUE;
		Prior = CreateNormalPrior("RJ_Landscape_Rate_Group", 0, 2.0);
		AddPriorToOpt(Opt, Prior);

		Opt->NoLandscapeRateGroup = atoi(Passed[1]);
	}
	else
		Opt->UseRJLandscapeRateGroup = FALSE;
}

void	SaveInitialTrees(OPTIONS *Opt, int Tokes, char **Passed)
{
	if(Opt->SaveInitialTrees != NULL)
	{
		free(Opt->SaveInitialTrees);
		Opt->SaveInitialTrees = NULL;
	}

	if(Tokes == 2)
		Opt->SaveInitialTrees = StrMake(Passed[1]);
	else
	{
		printf("Save Initial Trees requies a file name.\n");
		exit(1);
	}
}

void	SetMLTol(OPTIONS *Opt, int Tokes, char **Passed)
{

	if(Tokes != 2)
	{
		printf("MLTol Takes a likelihood tolerance used as a terminate criteria.\n");
		exit(1);
	}

	if(IsValidDouble(Passed[1]) == FALSE)
	{
		printf("Cannot convert %s to a valid tolerance.\n", Passed[1]);
		exit(1);
	}

	Opt->MLTol = atof(Passed[1]);

	if(Opt->MLTol < 0)
	{
		printf("MLTol must be >0.\n");
		exit(1);
	}
}

void	SetMLMaxEval(OPTIONS *Opt, int Tokes, char **Passed)
{
	if(Tokes != 2)
	{
		printf("MLMaxEVal takes the maximum number of likelihood evaluations to try.");
		exit(1);
	}

	if(IsValidInt(Passed[1]) == FALSE)
	{
		printf("Cannot convert maximum number of evaluations to a valid number.\n");
		exit(1);
	}

	Opt->MLMaxEVals = atoi(Passed[1]);

	if(Opt->MLMaxEVals <= -1 || Opt->MLMaxEVals == 0)
	{
		printf("maximum number of evaluations must be greater than zero, or -1 for no limit.\n");
		exit(1);
	}
}

void	SetMLAlg(OPTIONS *Opt, int Tokes, char **Passed)
{
	if(Tokes != 2)
	{
		printf("MLAlg takes an algorithm name, valid names are ");
		PrintAlgNames();
		exit(1);
	}

	if(ValidMLAlgName(Passed[1]) == FALSE)
	{
		printf("Invalid algorithm name %s\n", Passed[1]);
		PrintAlgNames();
		exit(1);
	}

	if(Opt->MLAlg != NULL)
		free(Opt->MLAlg);

	Opt->MLAlg = StrMake(Passed[1]);
}

void	SetMinTaxaNoTrans(OPTIONS *Opt, int Tokes, char **Passed)
{
	int No;

	if(Tokes != 2)
	{
		printf("SetMinTransTaxaNo takes a minimum number of taxa to transfom a node using RJ kappa, lambda, delta and OU.\n");
		exit(1);
	}

	if(IsValidInt(Passed[1]) == FALSE)
	{
		printf("Cannot Convert %s to a valid integer.\n", Passed[1]);
		exit(1);
	}

	No = atoi(Passed[1]);

	if(No < 1)
	{
		printf("Number of taxa must be greater than 1.\n");
		exit(1);
	}

	Opt->MinTransTaxaNo = No;
}

void SetMinMaxRate(OPTIONS *Opt, int Tokes, char **Passed)
{
	double Min, Max;

	if(Tokes != 3)
	{
		printf("SetMinMaxRate takes the a minimum and maximum rate.\n");
		exit(0);
	}

	if(IsValidDouble(Passed[1]) == FALSE || IsValidDouble(Passed[2]) == FALSE)
	{
		printf("Cannot convert values to rates.\n");
		exit(0);
	}

	Min = atof(Passed[1]);
	Max = atof(Passed[2]);

	if(Min >= Max)
	{
		printf("Min %f is not valid", Min);
		exit(0);
	}

	if(Min < RATE_MIN)
		Min = RATE_MIN;

	Opt->RateMin = Min;
	Opt->RateMax = Max;
}

void	SetNormQMatrix(OPTIONS *Opt, int Tokes, char **Passed)
{
	PRIOR *Prior;

	if(Opt->Model != M_MULTISTATE)
	{
		printf("Only multistate models can be normalised.\n");
		exit(0);
	}

	if(Tokes != 1)
	{
		printf("Normalise Q Matrix does not take any parameters.\n");
		exit(0);
	}

	if(Opt->NormQMat == TRUE)
	{
		Opt->NormQMat = FALSE;
		RemovePriorFormOpt("GlobalRate", Opt);
	}
	else
	{
		Prior = CreateUniformPrior("GlobalRate", 0, 100);
		AddPriorToOpt(Opt, Prior);
		Opt->NormQMat = TRUE;

	}
}

void	SetNoSliceSampleSteps(OPTIONS *Opt, int Tokes, char **Passed)
{
	int NoSteps;

	if(Tokes != 2)
	{
		printf("NoSliceSampleSteps take a number of steps.\n");
		exit(1);
	}

	if(IsValidInt(Passed[1]) == FALSE)
	{
		printf("Cannot convert %s to a valid number of steps", Passed[1]);
		exit(1);
	}

	NoSteps = atoi(Passed[1]);

	if(NoSteps < 2)
	{
		printf("Cannot convert %s to a valid number of steps", Passed[1]);
		exit(1);
	}

	Opt->NoSliceSampleSteps = NoSteps;
}

void	SetPisAncStates(OPTIONS *Opt, int Tokes, char **Passed)
{
	if(Tokes != 1)
	{
		printf("PisAncStates Does not take any parameters.\n");
		exit(1);
	}
	
	if(Opt->UsePisInAncStates == TRUE)
		Opt->UsePisInAncStates = FALSE;
	else
		Opt->UsePisInAncStates = TRUE;
}

void	SetRJZero(OPTIONS *Opt, int Tokes, char **Passed)
{
	if(Tokes != 1)
	{
		printf("RJZero Does not take any parameters.\n");
		exit(1);
	}

	if(Opt->RJZero == TRUE)
		Opt->RJZero = FALSE;
	else
		Opt->RJZero = TRUE;
}



int		PassLine(OPTIONS *Opt, char *Buffer, char **Passed)
{
	int			Tokes;
	COMMANDS	Command;
	int			Index;
	int			Temp;
	
	ReplaceChar(';', ' ', Buffer);
	ReplaceChar('=', ' ', Buffer);
	ReplaceChar(',', ' ', Buffer);
	RemoveChar('\n', Buffer);

	
	Tokes = MakeArgv(Buffer, Passed, BUFFERSIZE);

	if(Tokes == BUFFERSIZE)
	{
		printf("%s - Command line is too long, please contract developers for a solution.\n", Buffer);
		exit(0);
	}

	if(Tokes >= 1)
		MakeLower(Passed[0]);

	if(Tokes <= 0)
		return FALSE;

	Command = StringToCommand(Passed[0]);

	if(Command == CUNKNOWN)
	{
		printf("Unknown command: %s\n",Passed[0]);
		exit(1);
	}

	if(CmdVailWithDataType(Opt,Command) == FALSE)
	{
		printf("Command is not valid with data / model.\n");
		exit(1);
	}

	if(Command == CRUN)
		return TRUE;

	if(Command == CRES)
	{
		if(Tokes >= 3)
			Restrict(Opt, Tokes, Passed);
		else
			printf("The Restrict command takes two parimiters, a rate to restict and a constant or rate to restict it to.\n");
	}

	if(Command == CUNRES)
	{
		if(Tokes == 2)
			UnRestict(Opt, Passed[1]);
		else
			printf("The unresict command takes one paramtier, the name of the rate to unresict\n");
	}

	if(Command == CRESALL)
	{
		if(Tokes == 2)
			RestrictAll(Opt, Passed[1]);
		else
			printf("The RestrictAll command takes one parimiter eather a rate or a constant\n");
	}

	if(Command == CUNRESALL)
		UnRestictAll(Opt);
	
	if(Command == CPRIOR)
		SetPriorCmd(Opt, Tokes, Passed);

	if(Command == CPRIORALL)
		SetAllRatePriors(Opt, Tokes, Passed);
	
	if(Command == CHYPERPRIOR)
		SetHyperPriorCmd(Opt, Tokes, Passed);

	if(Command == CHPALL)
		SetHyperPriorAllCmd(Opt, Tokes, Passed);
	
	if(Command == CITTERS)
		SetItters(Opt, Tokes, Passed);
		
	if(Command == CSAMPLE)
	{
		if(Tokes == 2)
		{
			Temp = atoi(Passed[1]);
			if(Temp <= 0)
				printf("Could not convert %s to a valid samile frequncy", Passed[1]);
			else
				Opt->Sample= Temp;
		}
		else
			printf("Sample requires a number that spesifies the sample frequncy.\n");

	}
		
	
	if(Command == CPRIORCAT)
	{
		if(Tokes == 2)
		{
			Temp = atoi(Passed[1]);
			if(Temp == 0)
				printf("Could not convert %s to a valid number of descrete prior catagiores\n", Passed[1]);
			else
				Opt->PriorCats = Temp;
		}
		else
			printf("Itters requires a number that spesifies the number of descrete prior catagiores.\n");
	}

	if(Command == CMLTRIES)
	{
		if(Tokes == 2)
		{
			Temp = atoi(Passed[1]);
			if(Temp == 0)
				printf("Could not convert %s to a valid number of run of the optermising function (per tree).\n", Passed[1]);
			else
				Opt->MLTries = Temp;
		}
		else
			printf("Itters requires a number that spesifies the number of times to run the optermising function (per tree).\n");
	}
	
	if(Command == CMLTOL)
		SetMLTol(Opt, Tokes, Passed);
	
	if(Command == CMLEVAL)
		SetMLMaxEval(Opt, Tokes, Passed);

	if(Command == CMLALG)
		SetMLAlg(Opt, Tokes, Passed);

	if(Command == CINFO)
		PrintOptions(stdout, Opt);
		
	if(Command == CHELP)
	{
		Index=0;
		do
		{
			printf("%s\t%s\n", COMMANDSTRINGS[Index], COMMANDSTRINGS[Index+1]);
			Index+=2;
		}while(COMMANDSTRINGS[Index][0] != '\0');
	}

	if(Command == CNODE)
		AddRecNode(Opt, NODEREC, Tokes, Passed);

	if(Command == CMRCA)
		AddRecNode(Opt, MRCA, Tokes, Passed);


	if(Command == CADDTAXA)
	{
		printf("The AddTaxa command is no longer supported.\n");
		exit(1);
	}
	
	if(Command == CDELTAXA)
	{
		printf("The DelTaxa command is no longer supported.\n");
		exit(1);
	}

	if(Command == CEVENROOT)
		SetEvenRoot(Opt->Trees);

	if(Command == CLOGFILE)
		SetLogFile(Opt, Tokes, Passed); 


	if(Command == CPRESET)
		PreSet(Opt, Tokes, Passed);

	if(Command == CSUMMARY)
	{
		if(Opt->Summary == FALSE)
			Opt->Summary = TRUE;
		else
			Opt->Summary = FALSE;
	}

	if(Command == CBURNIN)
		SetBurnIn(Opt, Tokes, Passed);
	

	if(Command == CPIS)
	{
		if(Tokes == 2)
		{
			GetBasePis(Opt, Passed[1]);
		}
		else
			printf("PiTypes requres a type for the base fequncys, uni, none, est, emp\n");
	}

	if(Command == CKAPPA)
		SetKappa(Opt, Tokes, Passed);
	
	if(Command == CDELTA)
		SetDelta(Opt, Tokes, Passed);

	if(Command == CLAMBDA)
		SetLambda(Opt, Tokes, Passed);

	if(Command == COU)
		SetOU(Opt, Tokes, Passed);

	if(Command == CEXTTAXA)
	{
		if(Tokes > 1)
		{
			FreeTreeParts(Opt->Trees);
			FreeRecNodes(Opt, Opt->Trees->NoSites);
			ExcludeTaxa(Opt, Tokes-1, &Passed[1]);
			SetParts(Opt->Trees);
		}
		else
		{
			printf("The exclude taxa (et) command requires one or more taxa names or numbers\n");
		}
	}

	if(Command == CSAVEINITIALTREES)
		SaveInitialTrees(Opt, Tokes, Passed);

	if(Command == CTESTCORREL)
	{
		SetTestCorrel(Opt);
	}

	if(Command == CCOVARION)
	{
		if(Opt->UseCovarion == TRUE)
			Opt->UseCovarion = FALSE;
		else
			Opt->UseCovarion = TRUE;
	}

	if(Command == CREVJUMP)
		SetRJMCMC(Opt, Tokes-1, &Passed[1]);
	

	if(Command == CHPRJ)
		SetRJMCMCHP(Opt, Tokes-1, &Passed[1]);
	
	if(Command == CEXIT)
	{
		exit(0);
	}

	if(Command == CFOSSIL)
		AddRecNode(Opt, FOSSIL, Tokes, Passed);


	if(Command == CNODEDATA)
	{
		if(Opt->NodeData == FALSE)
		{
			Opt->NodeData = TRUE;
			Opt->NodeBLData = FALSE;
		}
		else
			Opt->NodeData = FALSE;
	}

	if(Command == CALPHAZERO)
	{
		if(Opt->AlphaZero == FALSE)
			Opt->AlphaZero = TRUE;
		else
			Opt->AlphaZero = FALSE;
	}



	if(Command == CNODEBLDATA)
	{
		if(Opt->NodeBLData == FALSE)
		{
			Opt->NodeBLData = TRUE;
			Opt->NodeData	= FALSE;
		}
		else
			Opt->NodeBLData	= FALSE;
	}

	if(Command == CGAMMA)
		SetGamma(Opt, Passed, Tokes);

	if(Command == CCI)
	{
		if(Tokes == 2)
			SetCI(Opt, Passed[1]);
		else
		{
			if(Opt->FindCF == TRUE)
			{
				Opt->FindCF = FALSE;
				Opt->CFRate = NULL;
			}
			else
				printf("Confidence intervals requires a paramiter\n");
		}
	}


	if(Command == CRMODEL)
	{
		if(Opt->UseRModel == TRUE)
		{
			Opt->UseRModel = FALSE;
			Opt->RModelP = -1;
		}
		else
		{
			if(Tokes == 1)
			{
				Opt->UseRModel = TRUE;
				Opt->RModelP = -1;
			}

			if(Tokes == 2)
			{
				if(IsValidDouble(Passed[1]) == FALSE)
					printf("Could not convert %s to a valid R Model rate\n", Passed[1]);
				else
				{
					Opt->UseRModel = TRUE;
					Opt->RModelP = atof(Passed[1]);
				}
			}
		}
	}



	if(Command == CSETSEED)
	{
		if(Tokes == 2)
			OptSetSeed(Opt, Passed[1]);
		else
			printf("SetSeed take an unsinged intger.\n");
		
	}

	if(Command == CMAKEUM)
		MakeUM(Opt->Trees);

	if(Command == CVARRATES)
		SetVarRatesOpt(Opt);
	
	if(Command == CEQUALTREES)
		SetEqualTrees(Opt, Tokes, Passed);
	
	if(Command == CPRECISION)
	{
#ifndef BIG_LH
		printf("Precision is only valid with the Big Lh build of BayesTraits.\n");
		exit(1);
		return FALSE;
#endif
		SetPrecision(Opt, Passed[1]);
	}

	if(Command == CCORES)
		SetCores(Opt, Tokes, Passed);

	if(Command == CSYMMETRICAL)
		SetSymmetrical(Opt);

	if(Command == CMCMCMLSTART)
		SetMCMCMLStart(Opt);

	if(Command == CCAPRJRATES)
		CapRJRatesNo(Opt, Tokes ,Passed);

	if(Command == CSAVEMODELS)
		SetSaveModels(Opt, Tokes, Passed);

	if(Command == CLOADMODELS)
		SetLoadModels(Opt, Tokes, Passed);

	if(Command == CADDERR)
		LoadAddErr(Opt, Tokes, Passed);
	
	if(Command == CSTONES)
		SetSteppingstone(Opt, Passed, Tokes);
	

	if(Command == CSHEDULE)
	{
		if(Opt->UseSchedule == FALSE)
			Opt->UseSchedule = TRUE;
		else
			Opt->UseSchedule = FALSE;
	}

	if(Command == CRJDUMMY)
		SetRJDummy(Opt, Passed, Tokes);
		
	if(Command == CSCALETREES)
		SetScaleTree(Opt, Passed, Tokes);
	
	if(Command == CRJLOCALTRANSFORM)
		SetRJLocalTransform(Opt, Passed, Tokes);
	
	if(Command == CFATTAILNORMAL)
		SetFatTailNormal(Opt);

	if(Command == CADDTAG)
		AddTag(Opt, Tokes, Passed);

	if(Command == CLOCALTRANSFORM)
		AddLocalTransform(Opt, Tokes, Passed);

	if(Command == CDISTDATA)
		SetDistData(Opt, Tokes, Passed);

	if(Command == CNOLH)
		SetNoLh(Opt);

	if(Command == CSAVETREES)
		SetSaveTrees(Opt);

	if(Command == CCSCHED)
		SetOptCustomSchedule(Opt, Tokes, Passed);

	if(Command == CADDTIMESLICE)
		OptAddTimeSlice(Opt, Tokes, Passed);

	if(Command == CADDPATTERN)
		OptAddPattern(Opt, Tokes, Passed);

	if(Command == CSETMINTAXATRANS)
		SetMinTaxaNoTrans(Opt, Tokes, Passed);

	if(Command == CSETMINMAXRATE)
		SetMinMaxRate(Opt, Tokes, Passed);

	if(Command == CNORMQMAT)
		SetNormQMatrix(Opt, Tokes, Passed);

	if(Command == CNOSLICESAMPLESTEPS)
		SetNoSliceSampleSteps(Opt, Tokes, Passed);

	if(Command == CPISANCSTATES)
		SetPisAncStates(Opt, Tokes, Passed);

	if(Command == CRJZERO)
		SetRJZero(Opt, Tokes, Passed);

	if(Command == CLANDSCAPE)
		SetLandscape(Opt, Tokes, Passed);
	
	if(Command == CRATESCALARS)
		SetRateScalars(Opt, Tokes, Passed);

	if(Command == CRJLANDSCAPERATEGROUP)
		SetRJLandscapeRateGroup(Opt, Tokes, Passed);

	return FALSE;
}

void	GetOptionsArry(OPTIONS *Opt, int Size, char** OptStr)
{
	int		Index;
	char	**Passed;

	Passed = (char**)SMalloc(sizeof(char*)*BUFFERSIZE);

	for(Index=0;Index<Size;Index++)
		PassLine(Opt, OptStr[Index], Passed);

	free(Passed);
}

void	GetOptions(OPTIONS *Opt)
{
	char	*Buffer;
	char	**Passed;

	Passed = (char**)SMalloc(sizeof(char*) * BUFFERSIZE);
	Buffer = (char*)SMalloc(sizeof(char) * BUFFERSIZE);

	do
	{
		fgets(Buffer, BUFFERSIZE, stdin); 
	} while(PassLine(Opt,Buffer, Passed) == FALSE);
	
	free(Buffer);
	free(Passed);
}

void	CheckOptions(OPTIONS *Opt)
{
	int NoFreeP;
	printf("\n");

	if(Opt->LoadModels == TRUE)
	{
		if(Opt->UseRJMCMC == TRUE)
		{
			printf("RJ MCMC and the use of a model file are mutuality exclusive.\n");
			exit(0);
		}
	}

	NoFreeP = FindNoOfRates(Opt);

	if(Opt->UseRJMCMC == FALSE && Opt->DataType == DISCRETE)
	{
		if(FindNoOfRates(Opt) > 25)
		{
			printf("Too many free parameter to estimate (%d), try reducing the model or using RJ MCMC\n", NoFreeP);
			printf("If you believe you data can support this number of free parameter please contact the developers to have this limitation removed.\n");
			exit(0);
		}	
	}

	if(Opt->Stones != NULL && Opt->Itters == -1)
	{
		printf("Stepping stone sampler is not valid with an infinite number of iterations.\n");
		exit(0);
	}
} 