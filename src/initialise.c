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
#include <math.h>

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
#include "QuadDouble.h"
#include "Contrasts.h"
#include "FatTail.h"
#include "Fossil.h"
#include "Pattern.h"

#ifdef BTOCL
#include "btocl_discrete.h"
#endif

OPTIONS*	SetUpOptions(TREES* Trees, char	*TreeFN, char *DataFN)
{
	OPTIONS		*Opt;
	MODEL		Model;
	ANALSIS		Analsis;


	Model	= GetModel(Trees);
	
	if(GetModelType(Model) == MT_FATTAIL)
		Analsis = ANALMCMC;
	else
		Analsis = GetAnalsis(Trees);

	CheckDataWithModel(DataFN, Trees, Model);

	PreProcessDataWithModel(Trees, Model);

	Opt = CreatOptions(Model, Analsis, Trees->NoStates, TreeFN, DataFN, Trees->SymbolList, Trees);
	
	return Opt;
}

void	PreProcess(OPTIONS *Opt, TREES* Trees)
{
	int		Index, ID;
	
	CheckSingleDescendent(Trees);

	SetPatternNo(Opt, Trees);

	if(Opt->ScaleTrees != -1.0)
		ScaleUserTrees(Trees, Opt->ScaleTrees);

	if(Opt->Model == M_CONTRAST_REG)
	{
		if(Opt->TestCorrel == FALSE)
			SetDataRegTC(Opt);
	}

	if(Opt->NormQMat == TRUE && Opt->NoPatterns > 0)
	{
		printf("Normalisation and multiple patters are not supported together.\n");
		exit(0);
	}
	
	SetNoOfThreads(Opt->Cores);

	if(Opt->Stones != NULL)
		Opt->Stones->ItStart = Opt->Itters + 1;
	
	for(Index=0;Index<Trees->NoTrees;Index++)
	{
		ID = 0;
		SetNodeIDs(Trees->Tree[Index]);
	}

	Opt->LogFile		= OpenWriteWithExt(Opt->BaseOutputFN, OUTPUT_EXT_LOG);

	#ifdef JNIRUN
		Opt->LogFileRead = OpenRead(Opt->LogFN);
		Opt->LogFileBuffer = (char*)SMalloc(sizeof(char) * LOGFILEBUFFERSIZE);
		Opt->PassedOut = (char**)SMalloc(sizeof(char*) * LOGFILEBUFFERSIZE);
	#endif

	Trees->UseCovarion	= Opt->UseCovarion;

	SetPTrees(Opt, Trees);


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

		if(Opt->Model == M_DESCCV)
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
		
		SetNOSPerSite(Opt);

		InitTreeBigLh(Opt, Trees);

#ifdef QUAD_DOUBLE
		InitQuadDoubleLh(Opt, Trees);
#endif

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
		SaveTrees(Opt->SaveInitialTrees, Opt->Trees);

	if(Opt->SaveTrees == TRUE)
		InitialiseOutputTrees(Opt, Trees);
	
	SaveUserBrachLengths(Trees);
}



void Finalise(OPTIONS *Opt, TREES* Trees)
{
	if(Opt->SaveTrees == TRUE)
	{
		fprintf(Opt->OutTrees, "end;");
		fclose(Opt->OutTrees);
	}
}