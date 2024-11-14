/*
*  BayesTriats 4.0
*
*  copyright 2022
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
#include <math.h>
#include <string.h>

#include "TypeDef.h"
#include "Initialise.h"
#include "TypeDef.h"
#include "Trees.h"
#include "Data.h"
#include "Options.h"
#include "Rates.h"
#include "Likelihood.h"
#include "RandLib.h"
#include "Priors.h"
#include "MCMC.h"
#include "Praxis.h"
#include "ML.h"
#include "GenLib.h"
#include "Continuous.h"
#include "Initialise.h"
#include "VarRates.h"
#include "BigLh.h"
#include "PTrees.h"
#include "Threaded.h"
#include "Contrasts.h"
#include "FatTail.h"
#include "Fossil.h"
#include "Pattern.h"
#include "RestrictionMap.h"
#include "Part.h"

#ifdef BTOCL
#include "btocl_discrete.h"
#endif

OPTIONS*	SetUpOptions(TREES* Trees, char	*TreeFN, char *DataFN)
{
	OPTIONS		*Opt;
	MODEL		Model;
	ANALSIS		Analsis;


	Model	= GetModel(Trees);
	Analsis = GetAnalsis(Trees);

	if(GetModelType(Model) == MT_FATTAIL && Analsis == ANALYSIS_ML)
	{
		printf("Fat Tail models require MCMC analysis.\n");
		exit(1);
	}

	CheckDataWithModel(DataFN, Trees, Model);

	PreProcessDataWithModel(Trees, Model);

	Opt = CreatOptions(Model, Analsis, Trees->NoStates, TreeFN, DataFN, Trees->SymbolList, Trees);
	
	return Opt;
}

void	CheckFabricGlobalRates(OPTIONS *Opt, TREES *Trees)
{
	if(!(Opt->Model == M_CONTRAST || Opt->Model == M_CONTRAST_CORREL))
		return;

	if(Trees->NoSites == 1)
		return;

	if(DataModifiedOptions(Opt) == FALSE)
		return;

	printf("The model cannot be use with multiple sites and either the directional betas or a global trend.\n");
	exit(1);
}

void	CheckMissingDataAndNonParametricMethods(OPTIONS *Opt, TREES *Trees)
{
	if(UseNonParametricMethods(Opt) == TRUE && Trees->NoOfRemovedTaxa > 0)
	{
		printf("\n\nNon-parametric methods, Variable rates, fabric model ect cannot be used with missing data.\nAs it effects the post processing. Delete the missing data from the tree.\n");

		printf("To quickly delete taxa with missing data, run a basic ML Brownian model and use the save trees (SaveTrees) command output the tree with missing data remove. For example\n\n");
		printf("7\n1\nSaveTrees\nRun\n\n");

		printf("this will produce and outtput file with the extension .Output.trees\n");
		exit(1);
	}
}

void CheckRegOpt(OPTIONS *Opt, TREES *Trees)
{
	int TaxaIndex,SiteIndex;
	TAXA	*Taxa;


	if(Opt->Model != M_CONTINUOUS_REG)
		return;

	if(Opt->NoOfRecNodes > 0)
	{
		printf("Ancestral states cannot be estimated for GLM regression models, are it requires the estimation of independent data. To estimate the dependent value, add a dummy node to the tree and include a independent variables in the data file, using ? for the dependent values.\n");
		exit(1);
	}

	for(TaxaIndex=0;TaxaIndex<Trees->NoTaxa;TaxaIndex++)
	{
		Taxa = Trees->Taxa[TaxaIndex];

		for(SiteIndex=1;SiteIndex<Trees->NoSites;SiteIndex++)
		{
			if(Taxa->EstDataP[SiteIndex] == TRUE)
			{
				printf("Independent sites %d (for taxa %s) cannot be estimated under a GLM regression. If a value is needed, try estimating under a non-regression GLM. Dependent sites can be estimated.\n", SiteIndex+1, Taxa->Name);
				exit(1);
			}
		}
	}
}

void	CheckOptions(OPTIONS *Opt, TREES *Trees)
{
	int NoFreeP;

	CheckMissingDataAndNonParametricMethods(Opt, Trees);
	
	if(Opt->LoadModels == TRUE)
	{
		if(Opt->UseRJMCMC == TRUE)
		{
			printf("RJ MCMC and the use of a model file are mutuality exclusive.\n");
			exit(1);
		}
	}

	NoFreeP = FindNoOfRates(Opt);

	if(Opt->UseRJMCMC == FALSE && Opt->DataType == DISCRETE)
	{
		if(FindNoOfRates(Opt) > 25)
		{
			printf("Too many free parameter to estimate (%d), try reducing the model or using RJ MCMC\n", NoFreeP);
			printf("If you believe you data can support this number of free parameter please contact the developers to have this limitation removed.\n");
			exit(1);
		}	
	}

	if(Opt->StoneOptions != NULL && Opt->Itters == -1)
	{
		printf("Stepping stone sampler is not valid with an infinite number of iterations.\n");
		exit(1);
	}

	if(Opt->RJLockModel == TRUE && Opt->StoneOptions == NULL)
	{
		printf("RJLockModel is only valid with Stepping Stones.\n");
		exit(1);
	}


	CheckFabricGlobalRates(Opt, Trees);

	CheckResMaps(Opt->RestrictionMaps, Opt->NoRestrictionMaps);


	CheckRegOpt(Opt, Trees);
} 


void	CheckUndefinedPrior(OPTIONS *Opt)
{
	int Index;
	PRIOR *Prior;

	for(Index=0;Index<Opt->NoAllPriors;Index++)
	{
		Prior = Opt->AllPriors[Index];
		if(Prior->Dist == PDIST_UNDEFINED)
		{
			if(strcmp("FabricBeta", Prior->Name) == 0)
			{
				printf("The priors on the betas x branch (FabricBeta) are undefined.\nThis prior is specific to the data under analysis, currently there is no generic prior that can be used.\nPlease see the manual section \"Some guidelines for developing prior distributions for directional effects\" and the paper \"General statistical model shows that macroevolutionary patterns and processes are consistent with Darwinian gradualism\" for more information.\n");
				exit(1);
			}

			printf("prior on %s is undefined.\n", Prior->Name);
			exit(1);
		}
	}
}

void	SetLockRJBranch(OPTIONS* Opt, TREES* Trees)
{
	TREE *Tree;
	NODE Node;
	TAG *Tag;
	int TagIndex, TIndex, NIndex;


	if(Opt->NoLockedRJBL == 0)
		return;

	for(TIndex=0;TIndex<Trees->NoTrees;TIndex++)
	{
		Tree = Trees->Tree[TIndex];
		for(NIndex=0;NIndex<Tree->NoNodes;NIndex++)
			Tree->NodeList[NIndex]->RJLockNode = FALSE;
	}

	for(TagIndex=0;TagIndex<Opt->NoLockedRJBL;TagIndex++)
	{
		Tag = Opt->LockedRJBL[TagIndex];

		for(TIndex=0;TIndex<Trees->NoTrees;TIndex++)
		{
			Tree = Trees->Tree[TIndex];
			Node = PartGetMRCA(Tree, Tag->Part); 
			Node->RJLockNode = TRUE;
		}
	}
}


void	PreProcess(OPTIONS *Opt, TREES* Trees)
{
	int		Index, ID;


	// Some other checks. 
	CheckOptions(Opt, Trees);

	SetTreesDistToRoot(Trees);
	
	CheckSingleDescendent(Trees);

	SetPatternNo(Opt, Trees);

	if(Opt->ScaleTrees != -1.0)
		ScaleUserTrees(Trees, Opt->ScaleTrees);

	if(Opt->Model == M_CONTRAST_REG)
	{
		if(Opt->TestCorrel == FALSE)
			SetDataRegTC(Opt, Trees);
	}

	if(Opt->NormQMat == TRUE && Opt->NoPatterns > 0)
	{
		printf("Normalisation and multiple patters are not supported together.\n");
		exit(0);
	}
	
	SetNoOfThreads(Opt->Cores);

//	Test for stones here. 
	
	for(Index=0;Index<Trees->NoTrees;Index++)
	{
		ID = 0;
		SetNodeIDs(Trees->Tree[Index]);
	}

		
	Trees->UseCovarion	= Opt->UseCovarion;

	SetPTrees(Opt, Trees);

	SetTreesInternalNodes(Trees);


	if(Opt->ModelType == MT_CONTINUOUS)
		InitContinus(Opt, Trees);

	if(Opt->ModelType == MT_CONTRAST)
		InitContrastAll(Opt, Trees);

	if(Opt->ModelType == MT_FATTAIL)
		InitFatTailTrees(Opt, Trees);

	if(Opt->ModelType == MT_DISCRETE)
	{
//		NormaliseTrees(Trees->NormConst, Trees);
		
		if(Opt->UseCovarion == TRUE)
			Trees->NoStates = Trees->NoStates * 2;

		if(Opt->Model == M_DISC_CV)
			Trees->NoStates = Trees->NoStates * 2;

		if(Opt->UseKappa == TRUE && Opt->FixKappa != -1)
		{
			for(Index=0;Index<Trees->NoTrees;Index++)
				TreeBLToPower(Trees, Trees->Tree[Index], Opt->FixKappa);

			Opt->FixKappa = -1;
			Opt->UseKappa = FALSE;
		}

		AllocPartial(Opt, Trees, Opt->UseGamma);
		AllocLHInfo(Trees, Opt);

		SetFossils(Trees, Opt);
		
		SetNOSPerSite(Opt, Trees);

		InitTreeBigLh(Opt, Trees);

#ifdef BTOCL
		btocl_AllocPMatrixInfo(Trees);
		//btocl_AllocLhInfo(Trees);
#endif
	}

		
	if(FindNoEstDataPoint(Trees) > 0)
		Opt->EstData = TRUE;
	else
		Opt->EstData = FALSE;
		
	if(Opt->SaveInitialTrees != NULL)
		SaveTrees(Opt->SaveInitialTrees, Trees);
	
	SaveUserBrachLengths(Trees);
	
	MapResMaps(Opt, Trees, Opt->RestrictionMaps, Opt->NoRestrictionMaps);

	CheckUndefinedPrior(Opt);

	SetLockRJBranch(Opt, Trees);

}



void Finalise(OPTIONS *Opt, TREES* Trees)
{
	if(Opt->OutTrees != NULL)
	{
		fprintf(Opt->OutTrees, "end;");
		fclose(Opt->OutTrees);
	}
}