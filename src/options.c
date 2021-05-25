#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#include "typedef.h"
#include "genlib.h"
#include "options.h"
#include "data.h"
#include "trees.h"
#include "priors.h"
#include "treenode.h"
#include "rand.h"

#define	RATEOUTPUTLEN	33
#define	RIGHTINDENT		4

void	FreeRecNodes(OPTIONS *Opt);

char*	FormatRateName(char* RateName)
{
	char* Ret;
	int	Index;

	Ret = (char*)malloc(sizeof(char) * RATEOUTPUTLEN+1);
	if(Ret == NULL)
		MallocErr();

	for(Index=0;Index<RATEOUTPUTLEN;Index++)
		Ret[Index] = ' ';
	Ret[RATEOUTPUTLEN] = '\0';

	sprintf(&Ret[RIGHTINDENT], "%s", RateName);
	Ret[RIGHTINDENT+strlen(RateName)] = ' ';
	
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
		FRateName = FormatRateName(Opt->RateName[Index]);
		fprintf(Str, "%s", FRateName);

		free(FRateName);

		if(Opt->UseRJMCMC == TRUE)
		{
			fprintf(Str, "RJ MCMC\n");
		}
		else
		{
			switch(Opt->ResTypes[Index])
			{
				case RESNONE:
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
}

void	PrintPriorVals(PRIORS	*P, FILE* Str)
{
	int		VIndex;

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

int		DepPrior(OPTIONS* Opt, PRIORS	*P)
{
	char*	Buffer;

	if(Opt->Model != CONTINUOUSREG)
		return FALSE;

	Buffer = (char*)malloc(sizeof(char) * BUFFERSIZE);
	if(Buffer == NULL)
		MallocErr();

	sprintf(Buffer, "Beta-%d", Opt->DependantSite + 1);
	
	if(strcmp(Buffer, P->RateName) == 0)
	{
		free(Buffer);
		return TRUE;
	}

	free(Buffer);
	return FALSE;
}

void	PrintPriorOpt(FILE* Str, OPTIONS *Opt)
{
	int		Index;
	int		VIndex;
	PRIORS	*P;
	char	*FRateName;

	fprintf(Str, "Prior Information:\n");
	fprintf(Str, "    Prior Categories:            %d\n", Opt->PriorCats);

	if(Opt->UseRJMCMC == TRUE)
	{
		fprintf(Str, "    RJ Prior                     %s ", DISTNAMES[(int)Opt->RJPrior->Dist]);
		if(Opt->RJPrior->HP == FALSE)
		{
			for(VIndex=0;VIndex<DISTPRAMS[Opt->RJPrior->Dist];VIndex++)
				fprintf(Str, "%2.2f ", Opt->RJPrior->DistVals[VIndex]);
		}
		else
		{
			for(VIndex=0;VIndex<DISTPRAMS[Opt->RJPrior->Dist];VIndex++)
				fprintf(Str, "(%2.2f,%2.2f) ", Opt->RJPrior->HP[VIndex*2], Opt->RJPrior->HP[(VIndex*2)+1]);
		}
		fprintf(Str, "\n");
	}

	for(Index=0;Index<Opt->NoOfRates;Index++)
	{
		FRateName = FormatRateName(Opt->RateName[Index]);
		fprintf(Str, "%s", FRateName);
		free(FRateName);

		P = Opt->Priors[Index];

		if((Opt->ResTypes[Index] == RESNONE) && (Opt->UseRJMCMC == FALSE) && (DepPrior(Opt, P) == FALSE))
		{
			PrintPriorVals(P, Str);
		}
		else
		{
			fprintf(Str, "N\\A\n");
		}
	}

	if(Opt->EstGamma == TRUE)
	{
		fprintf(Str, "        Gamma                    ");
		P = Opt->PriorGamma;

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

}

double	FindAveNodeDepth(RECNODE RNode, OPTIONS *Opt)
{
	double Ret=0;

	Ret = Opt->Trees->NoOfTaxa * Opt->Trees->NoOfTrees;
	Ret = Ret - RNode->Hits;
	Ret = Ret / Opt->Trees->NoOfTrees;
	Ret = Opt->Trees->NoOfTaxa - Ret;

	return Ret;
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

	if(Opt->DataType == CONTINUOUS)
	{
		fprintf(Str, "Data Deviation:                  %f\n", Opt->EstDataDev);
		fprintf(Str, "Estimating values for taxa and Sites\n");
	}
	else
		fprintf(Str, "Estimating values for taxa and Sites\n");

	for(TIndex=0;TIndex<Trees->NoOfTaxa;TIndex++)
	{
		Taxa = &Trees->Taxa[TIndex];
		if(Taxa->EstData == TRUE)
		{
			fprintf(Str, "\t");
			PrintFixSize(Taxa->Name, 20, Str);
			for(SIndex=0;SIndex<Trees->NoOfSites;SIndex++)
			{
				if(Taxa->EstDataP[SIndex] == TRUE)
					fprintf(Str, "%d ", SIndex+1);
			}
			fprintf(Str, "\n");
		}
	}
}

void	PrintOptions(FILE* Str, OPTIONS *Opt)
{
	RECNODE	RNode;
	int		Index;
	int		NOS;

	fprintf(Str, "Options:\n");

	fprintf(Str, "Model:                           %s\n", MODELNAMES[Opt->Model]);
	
	fprintf(Str, "Tree File Name:                  %s\n", Opt->TreeFN);
	fprintf(Str, "Data File Name:                  %s\n", Opt->DataFN);
	fprintf(Str, "Log File Name:                   %s\n", Opt->LogFN);

	if(Opt->Headers == FALSE)
	fprintf(Str, "Output Headers                   FALSE\n");

	fprintf(Str, "Summary:                         ");
	if(Opt->Summary == FALSE)
		fprintf(Str, "False\n");
	else
		fprintf(Str, "True\n");
	
	fprintf(Str, "Seed                             %lu\n", Opt->Seed);

	if(Opt->MakeUM == TRUE)
		fprintf(Str, "Make UM                      True\n");

	if(Opt->Analsis == ANALML)
	{
		fprintf(Str, "Analsis Type:                    Maximum Likelihood\n" );
		fprintf(Str, "ML attempt per tree:             %d\n", Opt->MLTries);
	}

	if(Opt->Analsis == ANALMCMC)
	{
		fprintf(Str, "Analysis Type:                   MCMC\n" );
		fprintf(Str, "Sample Period:                   %d\n", Opt->Sample);
		fprintf(Str, "Iterations:                      %d\n", Opt->Itters);
		fprintf(Str, "Burn in:                         %d\n", Opt->BurnIn);
		fprintf(Str, "Rate Dev:                        %f\n", Opt->RateDev);
		if(Opt->UseSchedule	== TRUE)
			fprintf(Str, "Schedule File:                   %s\n", Opt->ScheduleFile);
		if(Opt->DataType == CONTINUOUS)
		{
			for(Index=0;Index<Opt->NoOfRates;Index++)
			{
				fprintf(Str, "    ");
				PrintFixSize(Opt->RateName[Index], 29, Str);
				fprintf(Str, "%f\n", Opt->RateDevList[Index]);
			}
		}

		if(Opt->UseModelFile == TRUE)
		fprintf(Str, "Models taken form:               %s\n", Opt->ModelFile);

		fprintf(Str, "Solo Tree Move:                  ");
		if(Opt->SoloTreeMove == TRUE)
			fprintf(Str, "True\n");
		else
			fprintf(Str, "False\n");
	}

	if(Opt->NOSPerSite == TRUE)
		fprintf(Str, "Fit no of states per Site:       Yes\n");
	fprintf(Str, "No of Rates:                     %d\n", Opt->NoOfRates);

	if(Opt->DataType == DISCRETE)
	{
		fprintf(Str, "Base frequency (PI's)            ");
		switch(Opt->PiTypes)
		{
			case PINONE:
				fprintf(Str, "None\n");
				break;
			case PIEST:
				fprintf(Str, "Estimate\n");
				break;
			case PIEMP:
				fprintf(Str, "Empirical\n");
				break;
			case PIUNI:
				fprintf(Str, "Uniform\n");
				break;
		}

		fprintf(Str, "Character Symbols                ");

		if(Opt->Model == MULTISTATE)
		{
			NOS = Opt->Trees->NoOfStates;
			if(Opt->UseCovarion == TRUE)
				NOS = NOS / 2;
			for(Index=0;Index<NOS	-1;Index++)
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

		if(Opt->Model == CONTINUOUSREG)
			fprintf(Str, "Dependant Site:                  %d\n", Opt->DependantSite+1);

		fprintf(Str, "Kappa                            ");
		PrintConPar(Str, Opt->UseKappa, Opt->EstKappa, Opt->FixKappa);

		fprintf(Str, "Delta                            ");
		PrintConPar(Str, Opt->UseDelta, Opt->EstDelta, Opt->FixDelta);

		fprintf(Str, "Lambda                           ");
		PrintConPar(Str, Opt->UseLambda, Opt->EstLambda, Opt->FixLambda);

		if(Opt->AlphaZero == TRUE)
			fprintf(Str, "Alpha through zero:              True\n");

		if(Opt->NodeData == TRUE)
			fprintf(Str, "Model for Node Data:             True\n");

		if(Opt->NodeBLData == TRUE)
			fprintf(Str, "Model for Node BLS Data:         True\n");

		if(Opt->UseVarData == TRUE)
			fprintf(Str, "Varable Data form file:          %s\n", Opt->VarDataFile);
	
		if(Opt->UsePhyloPlasty == TRUE)
			fprintf(Str, "Using PhyloPlasty:               True\n");
	
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
			fprintf(Str, "Using R Model\n");
			if(Opt->RModelP != -1)
				fprintf(Str, "R Model Rates fixed to:          %f", Opt->RModelP);
				
		}
	}

	PrintEstData(Str, Opt);

	if(Opt->AnalyticalP == TRUE)
		fprintf(Str, "Analytical P:                    True\n");
	
	PrintOptRes(Str, Opt);
	if(Opt->Analsis == ANALMCMC)
		PrintPriorOpt(Str, Opt);

	RNode = Opt->RecNode;
	while((RNode != NULL) && (Opt->DataType != CONTINUOUS))
	{
		if(RNode->NodeType == MRCA)
			fprintf(Str, "NRCA:         %s                    %f\n", RNode->Name, FindAveNodeDepth(RNode, Opt));
			
		
		if(RNode->NodeType == NODEREC)
			fprintf(Str, "Node:         %s                    %f\n", RNode->Name, ((double)RNode->Hits / Opt->Trees->NoOfTrees)*100);
		

		if(RNode->NodeType == FOSSIL)
			fprintf(Str, "Fossil:       %s                    %f (%d)\n", RNode->Name, FindAveNodeDepth(RNode, Opt), RNode->FossilState);
		
		for(Index=0;Index<RNode->NoOfTaxa;Index++)
			fprintf(Str, "             %d\t%s\n", RNode->Taxa[Index]->No, RNode->Taxa[Index]->Name);

		RNode = RNode->Next;
	}

	if(Opt->Trees->NoOfRemovedTaxa != 0)
	{
		fprintf(Str, "Removed taxa:\n");
		for(Index=0;Index<Opt->Trees->NoOfRemovedTaxa;Index++)
			fprintf(Str, "          %s\n", Opt->Trees->RemovedTaxa[Index]);
	}


	PrintTrees(Str, Opt->Trees, Opt->DataType);
	fflush(stdout);
}

void	FreeOptions(OPTIONS *Opt)
{
	int		Index;
	
	if((Opt->Model == MULTISTATE) || (Opt->DataType == CONTINUOUS))
	{
		for(Index=0;Index<Opt->NoOfRates;Index++)
			free(Opt->RateName[Index]);
		free(Opt->RateName);
	}

	if(Opt->Analsis == ANALMCMC)
	{
		for(Index=0;Index<Opt->NoOfRates;Index++)
			FreePrior(Opt->Priors[Index]);
		
		free(Opt->Priors);

		FreePrior(Opt->RJPrior);
	}

	if(Opt->ModelFile != NULL)
		free(Opt->ModelFile);

	if(Opt->EstDataSites != NULL)
		free(Opt->EstDataSites);

	if(Opt->RateDevList != NULL)
		free(Opt->RateDevList);

	if(Opt->ScheduleFile != NULL)
		free(Opt->ScheduleFile);

	free(Opt->DataFN);
	free(Opt->TreeFN);
	free(Opt->LogFN);
	fclose(Opt->LogFile);

	if(Opt->LogFileRead != NULL)
		fclose(Opt->LogFileRead);

	if(Opt->LogFileBuffer != NULL)
		free(Opt->LogFileBuffer);

	if(Opt->PassedOut != NULL)
		free(Opt->PassedOut);

	free(Opt->ResTypes);
	free(Opt->ResNo);
	free(Opt->ResConst);

	FreeRecNodes(Opt);

	free(Opt);
}

char*	CreatRateName(char N1, char N2)
{
	char	Buffer[128];
	char	*Ret;

	sprintf(&Buffer[0], "q%c%c", N1, N2);
	Ret = (char*)malloc(sizeof(char)*strlen(&Buffer[0]) + 1);
	strcpy(Ret, &Buffer[0]);
	
	return Ret;
}

char**	ModelARateName(OPTIONS* Opt)
{
	char**	Ret;
	char*	Buffer;
	int		Index;

	Opt->NoOfRates = Opt->Trees->NoOfSites;

	Ret = (char**)malloc(sizeof(char*)*Opt->NoOfRates);
	if(Ret == NULL)
		MallocErr();

	Buffer = (char*)malloc(sizeof(char) * BUFFERSIZE);
	if(Buffer == NULL)
		MallocErr();

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
	int		Index;

	Opt->NoOfRates = Opt->Trees->NoOfSites * 2;

	Ret = (char**)malloc(sizeof(char*)*Opt->NoOfRates);
	if(Ret == NULL)
		MallocErr();

	Buffer = (char*)malloc(sizeof(char) * BUFFERSIZE);
	if(Buffer == NULL)
		MallocErr();

	for(Index=0;Index<Opt->NoOfRates;Index++)
	{
		if(Index<Opt->Trees->NoOfSites)
			sprintf(Buffer, "Alpha-%d", Index+1);
		else
			sprintf(Buffer, "Beta-%d", Index+1);

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

	Opt->NoOfRates = Opt->Trees->NoOfSites + 1;

	Ret = (char**)malloc(sizeof(char*)*Opt->NoOfRates);
	if(Ret == NULL)
		MallocErr();

	Buffer = (char*)malloc(sizeof(char) * BUFFERSIZE);
	if(Buffer == NULL)
		MallocErr();

	sprintf(Buffer, "Alpha");
	Ret[0] = StrMake(Buffer);
	
	for(Index=1;Index<Opt->NoOfRates;Index++)
	{
		sprintf(Buffer, "Beta-%d", Index);
		Ret[Index] = StrMake(Buffer);
	}

	free(Buffer);
	return Ret;
}

char**	ContrastRateNames(OPTIONS *Opt)
{
	char	**Ret;
	char	*Buffer;
	int		Index;
	int		i;
	
	Opt->NoOfRates = Opt->Trees->NoOfSites * 2;

	Ret = (char**)malloc(sizeof(char**) * Opt->NoOfRates);
	Buffer = (char*)malloc(sizeof(char*) * BUFFERSIZE);
	if((Ret == NULL) || (Buffer == NULL))
		MallocErr();
	
	i = 0;
	for(Index=0;Index<Opt->Trees->NoOfSites;Index++)
	{
		sprintf(Buffer, "Alpha-%d", Index+1);
		Ret[i++] = StrMake(Buffer);

		sprintf(Buffer, "Sigma-%d", Index+1);
		Ret[i++] = StrMake(Buffer);
	}

	free(Buffer);
	return Ret;
}

char**	CreatContinusRateName(OPTIONS* Opt)
{
	switch(Opt->Model)
	{
		case CONTINUOUSRR:
			return ModelARateName(Opt);
		case CONTINUOUSDIR:
			return ModelBRateName(Opt);
		case CONTINUOUSREG:
			return RetModelRateName(Opt);
		case CONTRASTM:
			return ContrastRateNames(Opt);
	}

	return NULL;
}

void	SetOptRates(OPTIONS* Opt, int NOS, char *SymbolList)
{
	int		Inner;
	int		Outter;
	int		Index;
	
	if(Opt->DataType == CONTINUOUS)
	{
		Opt->RateName	= CreatContinusRateName(Opt);
		return;
	}

	if(Opt->Model == DESCINDEP)
	{
		Opt->NoOfRates	= 4;
		Opt->RateName	= INDEPPRAMS;
		return;
	}

	if(Opt->Model == DESCDEP)
	{
		Opt->NoOfRates	= 8;
		Opt->RateName	= DEPPRAMS;
	}

	if(Opt->Model == MULTISTATE)
	{
		Opt->NoOfRates	= (NOS * NOS) - NOS;
		Opt->RateName	= (char**)malloc(sizeof(char*)*Opt->NoOfRates);
		if(Opt->RateName == NULL)
			MallocErr();

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
}

void		AllocRestictions(OPTIONS *Opt)
{
	int	Index;

	Opt->ResTypes		= (RESTYPES*)malloc(sizeof(RESTYPES) * Opt->NoOfRates);
	Opt->ResNo			= (int*)malloc(sizeof(int) * Opt->NoOfRates);
	Opt->ResConst		= (double*)malloc(sizeof(double) * Opt->NoOfRates);

	if( (Opt->ResTypes == NULL) ||
		(Opt->ResNo == NULL) ||
		(Opt->ResConst == NULL))
		MallocErr();

	for(Index=0;Index<Opt->NoOfRates;Index++)
	{
		Opt->ResTypes[Index] = RESNONE;
		Opt->ResNo[Index] = -1;
		Opt->ResConst[Index] = -1;
	}
}

PRIORS*	CreatDefPrior(OPTIONS *Opt)
{
	PRIORS *P;

	P = (PRIORS*)malloc(sizeof(PRIORS));
	if(P == NULL)
		MallocErr();

	P->Dist = UNIFORM;
	P->DistVals = (double*)malloc(sizeof(double) * 2);
	if(P->DistVals == NULL)
		MallocErr();

	if(Opt->DataType == DISCRETE)
		P->DistVals[0] = 0;
	else
		P->DistVals[0] = -100;

	P->DistVals[1] = 100;

	P->RateNo = -1;

	P->HP		= NULL;
	P->UseHP	= FALSE;
	P->OffSet	= 0;
	P->RateName	= NULL;

	return P;
}

void	AllocPrios(OPTIONS *Opt)
{
	int		Index;

	Opt->Priors		= (PRIORS**)malloc(sizeof(PRIORS*) * Opt->NoOfRates);
	if(Opt->Priors == NULL)
		MallocErr();

	for(Index=0;Index<Opt->NoOfRates;Index++)
	{
		Opt->Priors[Index] = CreatDefPrior(Opt);
		Opt->Priors[Index]->RateName = StrMake(Opt->RateName[Index]);
	}

	Opt->RJPrior = CreatDefPrior(Opt);
}

void	SetAllRateDevs(OPTIONS *Opt, double Dev)
{
	int Index;

	for(Index=0;Index<Opt->NoOfRates;Index++)
		Opt->RateDevList[Index] = Dev;
	Opt->RateDev = Dev;
}

OPTIONS*	CreatOptions(MODEL Model, ANALSIS Analsis, int NOS, char *TreeFN, char *DataFN, char *SymbolList, TREES* Trees)
{
	OPTIONS *Ret=NULL;
	char	Buffer[1024];

	if((Model == DESCDEP) || (Model == DESCINDEP))
		SquashDep(Trees);

	if((Model == CONTINUOUSRR) || (Model == CONTINUOUSDIR) || (Model == CONTINUOUSREG))
		RemoveConMissingData(Trees);


	Ret = (OPTIONS*)malloc(sizeof(OPTIONS));
	if(Ret == NULL)
		MallocErr();

	Ret->Trees		= Trees;
	Ret->Model		= Model;
	Ret->Analsis	= Analsis;
	Ret->TestCorrel	= FALSE;
	Ret->UseCovarion= FALSE;

	Ret->NodeData	= FALSE;
	Ret->NodeBLData = FALSE;
	Ret->AlphaZero	= FALSE;
	Ret->HPDev		= 1;
	Ret->DependantSite= -1;

	Ret->ModelFile	= NULL;
	Ret->UseModelFile= FALSE;

	Ret->UseRModel = FALSE;
	Ret->RModelP	= -1;
	Ret->EstDataDev	= 0.2;

	Ret->NoEstDataSite	=	0;
	Ret->EstDataSites	=	NULL;
	Ret->NoEstChanges	=	5;

	Ret->NOSPerSite		=	FALSE;
	Ret->RateDevList	=	NULL;


	if((Ret->Model == MULTISTATE) ||
		(Ret->Model == DESCINDEP) ||
		(Ret->Model == DESCDEP))
		Ret->DataType = DISCRETE;

	if( (Ret->Model == CONTINUOUSRR) || 
		(Ret->Model == CONTINUOUSDIR) || 
		(Ret->Model == CONTINUOUSREG) ||
		(Ret->Model == CONTRASTM))
	{
		Ret->TestCorrel = TRUE;
		Ret->DataType	= CONTINUOUS;

		if(Ret->Model == CONTINUOUSREG)
			Ret->DependantSite = 0;
	}

	SetOptRates(Ret, NOS, SymbolList);

	
	AllocRestictions(Ret);

	Ret->TreeFN = (char*)malloc(sizeof(char)*strlen(TreeFN)+1);
	Ret->DataFN = (char*)malloc(sizeof(char)*strlen(DataFN)+1);

	if((Ret->TreeFN == NULL) || (Ret->DataFN == NULL))
		MallocErr();

	strcpy(Ret->TreeFN, TreeFN);
	strcpy(Ret->DataFN, DataFN);

	sprintf(&Buffer[0], "%s.%s", DataFN, LOGFILEEXT);

	Ret->LogFN = (char*)malloc(sizeof(char)*strlen(Buffer)+1);
	if(Ret->LogFN== NULL)
		MallocErr();
	strcpy(Ret->LogFN, &Buffer[0]);

	Ret->LogFile		= NULL;
	Ret->LogFileRead	= NULL;
	Ret->LogFileBuffer	= NULL;
	Ret->PassedOut		= NULL;		

	if(Ret->Analsis == ANALML)
	{
		Ret->MLTries	=	10;
		Ret->Itters		=	-1;
		Ret->Sample		=	-1;
		Ret->Priors		=	NULL;
		Ret->PriorCats	=	-1;
		Ret->RateDev	=	-1;
		Ret->BurnIn		=	-1;
		Ret->RJPrior		=	NULL;
	}
	
	if(Ret->Analsis == ANALMCMC)
	{
		Ret->MLTries	=	-1;
		Ret->Itters		=	5050000;
		Ret->BurnIn		=	50000;

		Ret->Itters		=	1050000;
		Ret->BurnIn		=	50000;
		Ret->Sample		=	100;

		Ret->PriorCats	=	100;
		Ret->RateDev	=	2;
		Ret->EstDataDev	=	.05;
		AllocPrios(Ret);
	}

	if(Ret->DataType == CONTINUOUS)
	{
		Ret->RateDevList = (double*)malloc(sizeof(double)  * Ret->NoOfRates);
		if(Ret->RateDevList == NULL)
			MallocErr();
		SetAllRateDevs(Ret, Ret->RateDev);
	}

	Ret->RecNode		=	NULL;
	Ret->RecNodeList	=	NULL;
	Ret->NoOfRecNodes	=	0;
	Ret->Summary		=	FALSE;
	Ret->PiTypes		=	PINONE;

	Ret->UseKappa		=	FALSE;
	Ret->UseDelta		=	FALSE;
	Ret->UseLambda		=	FALSE;
	Ret->UseGamma		=	FALSE;

	Ret->EstKappa		=	FALSE;
	Ret->EstDelta		=	FALSE;
	Ret->EstLambda		=	FALSE;
	Ret->EstGamma		=	FALSE;

	Ret->FixKappa		=	-1;
	Ret->FixDelta		=	-1;
	Ret->FixLambda		=	-1;
	Ret->FixGamma		=	-1;

	Ret->InvertV		=	FALSE;

	Ret->PriorGamma		=	NULL;
	Ret->PriorKappa		=	NULL;
	Ret->PriorLambda	=	NULL;
	Ret->PriorDelta		=	NULL;

	Ret->UseRJMCMC		=	FALSE;

	Ret->FindCF			=	FALSE;
	Ret->CFRate			=	NULL;

	Ret->Headers		=	TRUE;

	Ret->VarData		=	NULL;
	Ret->UseVarData		=	FALSE;
	Ret->VarDataFile	=	NULL;

	Ret->AnalyticalP	=	FALSE;

	Ret->UseSchedule	=	FALSE;
	Ret->ScheduleFile	=	NULL;

	Ret->SoloTreeMove	=	FALSE;
	Ret->Seed			=	GetSeed();

	Ret->MakeUM			=	FALSE;

	Ret->UsePhyloPlasty	=	FALSE; 
	return Ret; 
}

MODEL	GetModel(TREES *Trees)
{
	char	Buffer[1024];
	int		Comment;

	Comment = FALSE;
	for(;;)
	{
		if(Comment == FALSE)
		{
			printf("Please Select the model of evolution to use.\n");
			printf("1)	MultiState.\n");
			printf("2)	Discrete: Independent\n");
			printf("3)	Discrete: Dependant\n");

			if(Trees->ValidCData == TRUE)
			{
				printf("4)	Continuous: Random Walk (Model A)\n");
				printf("5)	Continuous: Directional (Model B)\n");

				if(Trees->NoOfSites > 1)
					printf("6)	Continuous: Regression\n");

				printf("7)	Independent contrast\n");
			}
		}

		if(fgets(&Buffer[0], 1024, stdin) != NULL)
		{
			if(Buffer[0] != '#')
			{
				Comment = FALSE;
				switch (atoi(&Buffer[0]))
				{
					case 1:	
						return 	MULTISTATE;
					break;

					case 2:
						if((Trees->NoOfSites != 2) || (Trees->NoOfStates > 2))
						{
							printf("Discrete analisis requiers two two state characters\n");
							printf("There are %d states and %d sites in the current data set.\n", Trees->NoOfStates, Trees->NoOfSites);
							break;
						}
						return DESCINDEP;
					break;

					case 3:
						if((Trees->NoOfSites != 2) || (Trees->NoOfStates > 2))
						{
							printf("Descete analisis requiers two two stat characters\n");
							printf("There are %d states and %d sites in the current data set.\n", Trees->NoOfStates, Trees->NoOfSites);
							break;
						}

						return DESCDEP;
					break;

					case 4:
						if(Trees->ValidCData == TRUE)
							return CONTINUOUSRR;
					break;

					case 5:
						if(Trees->ValidCData == TRUE)
							return CONTINUOUSDIR;
					break;

					case 6:
						if(Trees->ValidCData == TRUE)
						{
							if(Trees->NoOfSites > 1)
								return CONTINUOUSREG;
							else
								printf("Continuous: Regression model requires more than one site\n");
						}
					break;


					case 7:
						if(Trees->ValidCData == TRUE)
							return CONTRASTM;
					break;


					default:
						if(Trees->ValidCData == FALSE)
							printf("%s is not a valid choice please enter a number between 1 and 3 followed by a return.\n", Buffer);
						else
							printf("%s is not a valid choice please enter a number between 1 and 5 followed by a return.\n", Buffer);
				}
			}
			else
				Comment = TRUE;
		}
	}

	return -1;
}


ANALSIS	GetAnalsis(TREES *Trees)
{
	char	Buffer[1024];
	int		Comment;

	Comment = FALSE;

	for(;;)
	{
		if(Comment == FALSE)
		{
			printf("Please Select the analysis method to use.\n");
			printf("1)	Maximum Likelihood.\n");
			printf("2)	MCMC\n");
		}
		
		fgets(&Buffer[0], 1024, stdin);

		if(Buffer[0] != '#')
		{
			Comment = FALSE;
			switch (atoi(&Buffer[0]))
			{
				case 1:
					return ANALML;
				break;

				case 2:
					return ANALMCMC;
				break;
			}

			printf("Unknown choice: %s\n", Buffer);
		}
		else
			Comment = TRUE;
	}
}

COMMANDS	StringToCommand(char *Str)
{
	int	CIndex=0;
	int	SIndex=0;

	if(Str[0] == '#')
		return CCOMMENT;

	do
	{
		if((strcmp(Str, COMMANDSTRINGS[SIndex]) == 0) || (strcmp(Str, COMMANDSTRINGS[SIndex+1]) == 0))
			return CIndex;
		
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
	
	printf("Restring to const\n");

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
	RESTYPES	*BackRes=NULL;
	double		*BackConst=NULL;
	int			*BackNo=NULL;
	double		Const;
	int			Safe;

	BackRes		= (RESTYPES*)malloc(sizeof(RESTYPES) * Opt->NoOfRates);
	BackConst	= (double*)malloc(sizeof(double) * Opt->NoOfRates);
	BackNo		= (int*)malloc(sizeof(int) * Opt->NoOfRates);
	
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

PRIORDIST	StrToPriorDist(char* Str)
{
	int			Index=0;
	
	do
	{
		if(strcmp(Str, DISTNAMES[Index])==0)
			return BETA+Index;
		Index++;
	} while(DISTNAMES[Index][0] != '\0');

	return -1;
}

int	SetPriorNo(OPTIONS *Opt, int RateNo, int Tokes, char *argv[])
{
	PRIORDIST	Dist;
	double		Num;
	int			Index;
	PRIORS		*P;

	Dist = StrToPriorDist(argv[0]);

	if(Dist == -1)
	{
		printf("Could not conver %s to a valid distriubntion\n", argv[0]);
		printf("Valid distiubtions are 	beta, gamma, uniform, exp\n");
		return FALSE;
	}

	if(Tokes!=DISTPRAMS[Dist]+1)
	{
		printf("Prior %s equires %d parmeters\n", DISTNAMES[Dist], DISTPRAMS[Dist]);
		return FALSE;
	}

	for(Index=1;Index<Tokes;Index++)
	{
		Num = atof(argv[Index]);
		if((Num == 0) && (strcmp(argv[Index], "0")!=0))
		{
			printf("Could not conver Prior paramiter %s to a valid number\n", argv[Index]);
			return FALSE;
		}
	}

	P = Opt->Priors[RateNo];

	free(P->DistVals);

	P->DistVals = (double*)malloc(sizeof(double) * DISTPRAMS[Dist]);
	if(P->DistVals == NULL)
		MallocErr();

	P->Dist = Dist;

	for(Index=0;Index<DISTPRAMS[Dist];Index++)
		P->DistVals[Index] = atof(argv[1+Index]);

	if(P->UseHP == TRUE)
	{
		free(P->HP);
		P->HP = NULL;
		P->UseHP = FALSE;
	}

	return TRUE;
}

void	SetPrior(OPTIONS *Opt, int Tokes, char *argv[])
{
	int			Rate;

	if(Opt->Analsis != ANALMCMC)
	{
		printf("Priors can only be set of MCMC analis\n");
		return;
	}

	Rate = StrToRate(Opt, argv[1]);
	if(Rate == -1)
	{
		printf("Could not convert %s to a valid rate paramiter\n", argv[1]);
		return;
	}

	SetPriorNo(Opt, Rate, Tokes-2, &argv[2]);
}

void	SetAllPriors(OPTIONS *Opt, int Tokes, char *argv[])
{
	int			Index;

	if(Opt->Analsis != ANALMCMC)
	{
		printf("Priors can only be set of MCMC analis\n");
		return;
	}

	for(Index=0;Index<Opt->NoOfRates;Index++)
		SetPriorNo(Opt, Index, Tokes-1, &argv[1]);

}

void	PrintUnderNode(NODE N)
{
	if(N->Tip == TRUE)
		printf("%s\t", N->Taxa->Name);
	else
	{
		PrintUnderNode(N->Left);
		PrintUnderNode(N->Right);
	}
}

void	FreeRecNode(RECNODE RNode)
{
		
	free(RNode->Name);
	free(RNode->Taxa);
	free(RNode->TreeNodes);
	free(RNode);
}

RECNODE	OptFindRecNode(OPTIONS *Opt, char* Name)
{
	RECNODE Ret;

	Ret = Opt->RecNode;
	while(Ret!=NULL)
	{
		if(strcmp(Name, Ret->Name)==0)
			return Ret;
		Ret = Ret->Next;
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

int	ValidDescFossileState(int FNo)
{
	if((FNo >= 0) && (FNo <= 3))
		return TRUE;

	if((FNo >= 10) && (FNo <= 15))
		return TRUE;

	if((FNo >= 20) && (FNo <= 23))
		return TRUE;

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


	return FALSE;
}

/* Will have to add support for extra states (10 to 23) */
int	DesFossilSate(char* State, OPTIONS *Opt)
{
	int	Ret;

	if(IsValidInt(State) == FALSE)
	{
		printf("Could not convert %s to a valid Discrete state", State);
		return -1;
	}

	Ret = atoi(State);
	if(ValidDescFossileState(Ret) == FALSE)
		return -1;

	return Ret;
}

int	FossilState(char *State, OPTIONS *Opt)
{
	int	Index;

	if(Opt->Model == MULTISTATE)
	{
		for(Index=0;Index<Opt->Trees->NoOfStates;Index++)
		{
			if(State[0] == Opt->Trees->SymbolList[Index])
				return Index;
		}

		printf("State %s is invalid\n", State);
		return -1;
	}

	return DesFossilSate(State, Opt);
}

int		GetTaxaNoFormName(char* Name, TREES* Trees, int *No)
{
	int	Index;
	
	for(Index=0;Index<Trees->NoOfTaxa;Index++)
	{
		if(strcmp(Name, Trees->Taxa[Index].Name) == 0)
		{
			*No = Trees->Taxa[Index].No;
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

			if(GetTaxaFromID(TaxaNo, Opt->Trees->Taxa, Opt->Trees->NoOfTaxa) == NULL)
			{
				printf("Could not convert %s to a valid taxa number\n", List[Index]);
				return FALSE;
			}
		}
		else
		{
			if(GetTaxaNoFormName(List[Index], Opt->Trees, &TaxaNo) == FALSE)
			{
				printf("Could not convert %s to a valid taxa name\n", List[Index]);
				return FALSE;
			}
		}
	}
	
	return TRUE;
}

TAXA*	GetTaxaFromNameNo(char *ID, TREES* Trees)
{
	int	Index;

	if(IsValidInt(ID) == TRUE)
		return GetTaxaFromID(atoi(ID), Trees->Taxa, Trees->NoOfTaxa);

	for(Index=0;Index<Trees->NoOfTaxa;Index++)
		if(strcmp(ID, Trees->Taxa[Index].Name) == 0)
			return &Trees->Taxa[Index];

	return NULL;
}					  

int	CheckConFState(char *State, int No)
{
	if(strcmp(State, ESTDATAPOINT) == 0)
		return TRUE;

	if(IsValidDouble(State) == FALSE)
	{
		
	}
}

char**	SetConFState(OPTIONS *Opt, NODETYPE NodeType, int Tokes, char *argv[])
{
	int		Index;
	char	**Ret;

	if(Opt->DataType == DISCRETE)
		return NULL;

	Ret = (char**)malloc(sizeof(char*) * Opt->Trees->NoOfSites);
	if(Ret == NULL)
		MallocErr();

	if(NodeType != FOSSIL)
	{
		for(Index=0;Index<Opt->Trees->NoOfSites;Index++)
			Ret[Index] = StrMake(ESTDATAPOINT); 
	}

	for(Index=0;Index<Opt->Trees->NoOfSites;Index++)
	{
		
	}
}

void	AddRecNode(OPTIONS *Opt, NODETYPE NodeType, int Tokes, char *argv[])
{
	RECNODE		RNode;
	int			Index;
	int			FState;
	char**		ConFState;
	

	FState = -1;
	
	ConFState = NULL;

	if(OptFindRecNode(Opt, argv[1]) != NULL)
	{
		printf("Node name %s is allready is use please chose another or delete the node\n", argv[1]);
		return;
	}

	Index=2;

	if(NodeType == FOSSIL)
	{
		if(Opt->DataType == DISCRETE)
		{
			FState = FossilState(argv[2], Opt);
			if(FState == -1)
				return;
				
			Index++;
		}
		else
		{
			Index += Opt->Trees->NoOfSites;
		}
	}

	
	if(ValidTaxaList(argv, Index, Tokes, Opt) == FALSE)
		return;

	RNode = (RECNODE)malloc(sizeof(struct RNODE));
	if(RNode == NULL)
		MallocErr();

	RNode->TaxaID	= NULL;
	RNode->ConData	= NULL;

	RNode->Name = (char*)malloc(sizeof(char)*strlen(argv[1])+1);
	if(RNode->Name == NULL)
		MallocErr();
	strcpy(RNode->Name, argv[1]);

	RNode->NodeType		= NodeType;
	RNode->FossilState	= FState;

	if(RNode->NodeType == FOSSIL)
		RNode->NodeType = MRCA;

	if(NodeType == FOSSIL)
		RNode->NoOfTaxa = Tokes - 3;
	else
		RNode->NoOfTaxa = Tokes - 2;

	RNode->Taxa = (TAXA**)malloc(sizeof(TAXA*)*RNode->NoOfTaxa);
	if(RNode->Taxa == NULL)
		MallocErr();

	for(Index=0;Index<RNode->NoOfTaxa;Index++)
	{
		if(NodeType == FOSSIL)
			RNode->Taxa[Index] = GetTaxaFromNameNo(argv[Index+3], Opt->Trees);
		else
			RNode->Taxa[Index] = GetTaxaFromNameNo(argv[Index+2], Opt->Trees);
	}

	RNode->Next = Opt->RecNode;
	Opt->RecNode = RNode;

	RNode->TreeNodes = (NODE*)malloc(sizeof(NODE)*Opt->Trees->NoOfTrees);
	if(RNode->TreeNodes == NULL)
		MallocErr();

	SetRecNodes(RNode, Opt);
	
	Opt->NoOfRecNodes++;

	RNode->NodeType = NodeType;
}

void	DelRecNode(OPTIONS *Opt, char* NodeName)
{
	int		Found;
	RECNODE	RNode=NULL;
	RECNODE Last;

	Found = FALSE;

	RNode = Opt->RecNode;

	while(RNode != NULL)
	{
		if(strcmp(RNode->Name, NodeName)==0)
			Found = TRUE;

		RNode = RNode->Next;
	}

	if(Found == FALSE)
	{
		printf("Could not find node %s\n", NodeName);
		return;
	}

	RNode = Opt->RecNode;
	if(strcmp(RNode->Name, NodeName)==0)
		Opt->RecNode = Opt->RecNode->Next;
	else
	{
		Last = RNode;
		RNode = RNode->Next;
		Found = FALSE;
		while(Found == FALSE)
		{
			if(strcmp(RNode->Name, NodeName) == 0)
			{
				Last->Next = RNode->Next;
				Found = TRUE;
			}
			else
			{	
				Last = RNode;
				RNode = RNode->Next;
			}
		}
	}

	FreeRecNode(RNode);
	Opt->NoOfRecNodes--;
}


void	AddToRecNode(OPTIONS *Opt, int Tokes, char *argv[])
{
	RECNODE RNode=NULL;
	int		TIndex;
	int		AddIndex;
	TAXA	**NewTaxaList;
	
	RNode = OptFindRecNode(Opt, argv[1]);
	if(RNode == NULL)
	{
		printf("Node %s could not be found\n", argv[1]);
		return;
	}

	for(AddIndex=2;AddIndex<Tokes;AddIndex++)
	{
		if((atoi(argv[AddIndex]) == 0) && (strcmp(argv[AddIndex], "0") != 0))
		{
			printf("Could not convert %s to a valid taxa number\n", argv[AddIndex]);
			return;
		}
		if(GetTaxaFromID(atoi(argv[AddIndex]), Opt->Trees->Taxa, Opt->Trees->NoOfTaxa) == NULL)
		{
			printf("Could not convert %s to a valid taxa number\n", argv[AddIndex]);
			return;
		}

		for(TIndex=0;TIndex<RNode->NoOfTaxa;TIndex++)
		{
			if(atof(argv[AddIndex]) == RNode->Taxa[TIndex]->No)
			{
				printf("Could not added taxa %s as it is allready in the list\n", argv[AddIndex]);
				return;
			}
		}
	}

	NewTaxaList = (TAXA**)malloc(sizeof(TAXA*) * ((Tokes-2) + (RNode->NoOfTaxa)));
	if(NewTaxaList == NULL)
		MallocErr();

	TIndex=0;
	for(AddIndex=0;AddIndex<RNode->NoOfTaxa;AddIndex++,TIndex++)
		NewTaxaList[TIndex] = GetTaxaFromID(RNode->Taxa[TIndex]->No, Opt->Trees->Taxa, Opt->Trees->NoOfTaxa);

	for(AddIndex=2;AddIndex<Tokes;TIndex++,AddIndex++)
		NewTaxaList[TIndex] = GetTaxaFromID(atoi(argv[AddIndex]), Opt->Trees->Taxa, Opt->Trees->NoOfTaxa);

	RNode->NoOfTaxa += Tokes-2;
	free(RNode->Taxa);
	RNode->Taxa = NewTaxaList;

	SetRecNodes(RNode, Opt);

	fflush(stdout);
}

void	DelToRecNode(OPTIONS *Opt, int Tokes, char *argv[])
{
	RECNODE RNode=NULL;
	int		TIndex;
	int		ArgIndex;
	int		No;
	int		Valid;
	TAXA	**NewTaxaList=NULL;

	RNode = OptFindRecNode(Opt, argv[1]);

	if(RNode == NULL)
	{
		printf("Could not find node %s\n", argv[1]);
		return;
	}

	if(RNode->NoOfTaxa - (Tokes -2) < 2)
	{
		printf("Each node must have two or more taxa in it. Use the DelNode command to remove the Node\n");
		return;
	}

	for(ArgIndex=2;ArgIndex<Tokes;ArgIndex++)
	{
		No = atoi(argv[ArgIndex]);

		if((No == 0) && (strcmp(argv[ArgIndex], "0")!=0))
		{
			printf("Could not conver %s to a valid taxa number\n", argv[ArgIndex]);
			return;
		}


		Valid = FALSE;
		for(TIndex=0;TIndex<RNode->NoOfTaxa;TIndex++)
			if(RNode->Taxa[TIndex]->No == No)
				Valid = TRUE;

		if(Valid == FALSE)
		{
			printf("Taxa no %d is not a member of the Node\n", No);
			return;
		}
	}

	NewTaxaList = (TAXA**)malloc(sizeof(TAXA*)*(RNode->NoOfTaxa - (Tokes - 2)));
	if(NewTaxaList == NULL)
		MallocErr();

	No = 0;
	for(TIndex=0;TIndex<RNode->NoOfTaxa;TIndex++)
	{	Valid = TRUE;
		for(ArgIndex=2;ArgIndex<Tokes;ArgIndex++)
			if(atoi(argv[ArgIndex]) == RNode->Taxa[TIndex]->No)
				Valid = FALSE;
		if(Valid == TRUE)
			NewTaxaList[No++] = RNode->Taxa[TIndex];
	}

	free(RNode->Taxa);
	RNode->Taxa = NewTaxaList;

	RNode->NoOfTaxa = RNode->NoOfTaxa - (Tokes - 2);

	SetRecNodes(RNode, Opt);
}

void	SetEvenRoot(TREES *Trees)
{
	int		TIndex;
	double	t;
	NODE	Root;

	for(TIndex=0;TIndex<Trees->NoOfTrees;TIndex++)
	{
		Root = Trees->Tree[TIndex].Root;
		t = 0;

		t = (Root->Right->Length + Root->Left->Length) / 2;
		Root->Left->Length = t;
		Root->Right->Length = t;
	}
}

void	LogFile(OPTIONS *Opt, char *LogFN)
{
	FILE*	TempLogFile=NULL;

	TempLogFile = fopen(LogFN, "w");
	if(TempLogFile == NULL)
	{
		printf("Could not open file %s for writting\n", LogFN);
		return;
	}
	fclose(TempLogFile);

	free(Opt->LogFN);

	Opt->LogFN = StrMake(LogFN);
}

void	PreSet(OPTIONS *Opt, char* Set)
{
	MakeLower(Set);

	if(strcmp(Set, "lang1ml")==0)
	{
		RestrictAll(Opt, Opt->RateName[0]);
		Opt->AnalyticalP = TRUE;
	}

	if(strcmp(Set, "lang1mcmc")==0)
	{
		RestrictAll(Opt, Opt->RateName[0]);
		Opt->AnalyticalP = TRUE;
	}
}

void	GetBasePis(OPTIONS *Opt, char* Type)
{
	MakeLower(Type);

	if(strcmp(Type, "est")==0)
	{
		Opt->PiTypes = PIEST;
		return;
	}

	if(strcmp(Type, "emp")==0)
	{
		Opt->PiTypes = PIEMP;
		return;
	}

	if(strcmp(Type, "uni")==0)
	{
		Opt->PiTypes = PIUNI;
		return;
	}

	if(strcmp(Type, "none")==0)
	{
		Opt->PiTypes = PINONE;
		return;
	}

	printf("The option %s, is unknown. Valid options are est, emp, uni and none\n");
}

int		CmdVailWithDataType(OPTIONS *Opt, COMMANDS	Command)
{
	if(Opt->DataType == CONTINUOUS)
	{
		if( (Command == CNODE)		||
			(Command == CDELNODE)	||
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
			(Command == CFOSSIL)	||
			(Command == CPRIORALL)	||
			(Command == CPRIOR)		||
			(Command == CPIS)		||
			(Command == CNOSPERSITE)
			)
		{
			printf("Command %s (%s) is not valid with the current model\n", COMMANDSTRINGS[Command*2], COMMANDSTRINGS[(Command*2)+1]);
			return FALSE;
		}
	}
	else
	{
		if(
			(Command == CDELTA)		|| 
			(Command == CLAMBDA)	||
			(Command == CALPHAZERO)	||
			(Command == CNODEBLDATA)||
			(Command == CNODEDATA)  ||
			(Command == CDEPSITE)   ||
			(Command == CDATADEV)	||
			(Command == CPHYLOPLASTY) 
			)
		{
			printf("Command %s (%s) is not valid with Discrete data\n", COMMANDSTRINGS[Command*2], COMMANDSTRINGS[(Command*2)+1]);
			return FALSE;
		}
	}

	if(Opt->Analsis == ANALMCMC)
	{
		if(Command == CCI)
		{
			printf("Confidence intervals can only be established under ML.\n");
			return FALSE;
		}
	}

	if(Opt->Analsis == ANALML)
	{
		if(
			(Command ==	CHYPERPRIOR)||
			(Command ==	CHPRJ)		||
			(Command ==	CHPALL)		||
			(Command ==	CREVJUMP)   ||
			(Command == CMODELFILE) ||
			(Command == CSCHEDULE)	||
			(Command == CSTREEMOVE) ||
			(Command == CPHYLOPLASTY)
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

void	ExcludeTaxa(OPTIONS *Opt, int	Tokes, char **Passed)
{
	int		Index;
	int		No;
	char	*Name;
	TAXA	*Taxa;

	for(Index=0;Index<Tokes;Index++)
	{
		
		if(IsValidInt(Passed[Index]))
		{
			No	 = atoi(Passed[Index]);
			Taxa = GetTaxaFromID(No, Opt->Trees->Taxa, Opt->Trees->NoOfTaxa);
		}
		else
			Taxa = GetTaxaFromName(Passed[Index], Opt->Trees->Taxa, Opt->Trees->NoOfTaxa);

		if(Taxa == NULL)
		{
			printf("Paramiter %s cannot be converted to a valid taxa number of name\n", Passed[Index]);
			return;
		}

		Name = Taxa->Name;

		if(RemoveTaxa(Opt, Opt->Trees, Name) == FALSE)
			return;
	}
}

void	PrintTaxaInfo(OPTIONS *Opt)
{
	int		TIndex;
	TREES	*Trees;
	char	Buffer[64];

	Trees = Opt->Trees;

	printf("Taxa Info\n");

	for(TIndex=0;TIndex<Trees->NoOfTaxa;TIndex++)
	{
		printf("    ");
		sprintf(&Buffer[0], "%d", Trees->Taxa[TIndex].No);
		PrintFixSize(&Buffer[0], 5, stdout);
		printf("%s", Trees->Taxa[TIndex].Name);
		printf("\n");
	}
}

void	SetRJMCMC(OPTIONS *Opt, int Tokes, char** Passed)
{
	int			Index;
	PRIORDIST	Dist;

	if(Opt->UseRJMCMC == TRUE)
	{
		Opt->UseRJMCMC = FALSE;
		return;
	}

	if(Tokes < 1)
	{
		printf("To turn RJ MCMC on prior distrusions has to be spesified\n");
		return;
	}

	Dist = StrToPriorDist(Passed[0]);

	if(Dist == -1)
	{
		printf("Could not convter %s to a valid distrubtion\n", Passed[0]);
		return;
	}


	if(Tokes - 1 != DISTPRAMS[Dist])
	{
		printf("Prior %s take %d parmeters\n", DISTNAMES[Dist], DISTPRAMS[Dist]);
		return;
	}

	for(Index=0;Index<DISTPRAMS[Dist];Index++)
	{
		if(IsValidDouble(Passed[Index + 1]) == FALSE)
		{
			printf("Could not conver %s to a valid prior paramiter\n", Passed[Index + 1]);
			return;
		}
	}

	free(Opt->RJPrior->DistVals);
	Opt->RJPrior->DistVals = (double*)malloc(sizeof(double) * DISTPRAMS[Dist]);
	if(Opt->RJPrior->DistVals == NULL)
		MallocErr();
	
	for(Index=0;Index<DISTPRAMS[Dist];Index++)
		Opt->RJPrior->DistVals[Index] = atof(Passed[Index+1]);

	Opt->RJPrior->Dist = Dist;
	Opt->UseRJMCMC = TRUE;
}

PRIORS*	NameToPrior(OPTIONS *Opt, char* Name)
{
	int	Index=0;

	for(Index=0;Index<Opt->NoOfRates;Index++)
		if(strcmp(Opt->RateName[Index], Name)==0)
			return Opt->Priors[Index];

	return NULL;
}

void	SetHyperPrior(OPTIONS *Opt, char* PName, char** Passed, int Tokes, PRIORS*	P)
{
	PRIORS*		Prior;
	PRIORDIST	Dist;
	int			Index;
	int			PIndex;
	double		Low, High;

	if(P == NULL)
		Prior = NameToPrior(Opt, PName);
	else
		Prior = P;

	if(Prior == NULL)
	{
		printf("Could not convert %s to a valid rate name\n", PName);
		return;
	}

	Dist = StrToPriorDist(Passed[0]);

	if(Dist == -1)
	{
		printf("Could not convert %s to a valid distibution type\n", Passed[0]);
		return;
	}

	if(Tokes - 1 != DISTPRAMS[Dist] * 2)
	{
		printf("The hyper prior command require an upper and lower bound for each paramtier in the distribution\n");
		return;
	}

	for(Index=0;Index<DISTPRAMS[Dist] * 2;Index++)
	{
		if(IsValidDouble(Passed[Index+1]) == FALSE)
		{
			printf("Could not convert %s to a valid double\n", Passed[Index+1]);
		}
	}
	
	PIndex=1;
	for(Index=0;Index<DISTPRAMS[Dist];Index++)
	{
		Low		= atof(Passed[PIndex]);
		High	= atof(Passed[PIndex+1]);
		
		if(High <= Low)
		{
			printf("%f must be grater than %f\n", High, Low);
			return;
		}

		PIndex+=2;
	}


	free(Prior->DistVals);

	Prior->UseHP = TRUE;
	if(Prior->HP != NULL)
		free(Prior->HP);
	Prior->HP = (double*)malloc(sizeof(double) * (DISTPRAMS[Dist] * 2));
	Prior->DistVals = (double*)malloc(sizeof(double) * DISTPRAMS[Dist]);

	if((Prior->HP == NULL) || (Prior->DistVals == NULL))
		MallocErr();

	for(Index=0;Index<DISTPRAMS[Dist] * 2;Index++)
		Prior->HP[Index] = atof(Passed[Index+1]);

	for(Index=0;Index<DISTPRAMS[Dist];Index++)
	{
		Prior->DistVals[Index] = Prior->HP[(Index * 2) + 1] - Prior->HP[(Index * 2)];
		Prior->DistVals[Index] = Prior->DistVals[Index] / 2;
		Prior->DistVals[Index] += Prior->HP[(Index * 2)];
	}

	Prior->Dist		= Dist;
	Prior->UseHP	= TRUE;
}

void	SetHPAll(OPTIONS *Opt, char** Passed, int NoOfTokes)
{
	int	Index;

	for(Index=0;Index<Opt->NoOfRates;Index++)
	{
		SetHyperPrior(Opt, Opt->RateName[Index], Passed, NoOfTokes, NULL);
	}
}

void	SetGamma(OPTIONS *Opt, char** Passed, int Tokes)
{
	int		GammaCats;
	double	Value;

	if(Tokes == 1)
	{
		if(Opt->UseGamma == TRUE)
		{
			Opt->UseGamma	= FALSE;
			Opt->FixGamma	= -1;
			Opt->EstGamma	= FALSE;
			if(Opt->PriorGamma != NULL)
			{
				FreePrior(Opt->PriorGamma);
				Opt->PriorGamma = NULL;
			}
			return;
		}
		else
		{
			printf("The Gamma command take 0, 1 or 2 parameters.\n 0 to turn Gamma off\n1 to estermate Gamma\n 2 to fix it to a constant\n");
			return;
		}
	}

	if(IsValidInt(Passed[1]) == FALSE)
	{
		printf("Could not convert %s to a valid number of categories to divide the gamma disruption up into.\n", Passed[1]);
		return;
	}

	GammaCats = atoi(Passed[1]);

	if((GammaCats < 2) || (GammaCats > 8))
	{
		printf("The number of gamma catergoires must be grater than 1 and less than 8\n");
		return;
	}

	if(Tokes == 2)
	{
		if(Opt->PriorGamma != NULL)
		{
			FreePrior(Opt->PriorGamma);
			Opt->PriorGamma = NULL;
		}

		Opt->UseGamma = TRUE;
		Opt->EstGamma = TRUE;
		Opt->FixGamma = -1;
	
		Opt->GammaCats = GammaCats;

		if(Opt->Analsis == ANALMCMC)
		{
			Opt->PriorGamma = CreatDefPrior(Opt);
			Opt->PriorGamma->RateNo = -1;
		}
		
		return;
	}

	if(Tokes == 3)
	{
		if(IsValidDouble(Passed[2]) == FALSE)
		{
			printf("Could not convert %s to a valid gamma shapre parmiter\n", Passed[2]);
			return;
		}
		Value = atof(Passed[2]);

		if((Value < GAMMAMIN) || (Value > GAMMAMAX))
		{
			printf("Gamma shape parmiter must be grater than %f and less than %f\n", GAMMAMIN, GAMMAMAX);
			return;
		}

		if(Opt->PriorGamma != NULL)
		{
			FreePrior(Opt->PriorGamma);
			Opt->PriorGamma = NULL;
		}

		Opt->UseGamma = TRUE;
		Opt->FixGamma = Value;
		Opt->EstGamma = FALSE;

		Opt->GammaCats= GammaCats;
		
		return;
	}

	printf("The Gamma command take 0, 1 or 2 parameters.\n 0 to turn Gamma off\n1 to estermate Gamma\n 2 to fix it to a constant\n");
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

void	SetModelFile(OPTIONS *Opt, int Tokes, char **Passed)
{
	if(Tokes == 1)
	{
		if(Opt->ModelFile != NULL)
		{
			free(Opt->ModelFile);
			Opt->ModelFile = NULL;
		}

		Opt->UseModelFile = FALSE;
		return;
	}

	if(Tokes == 2)
	{
		if(Opt->ModelFile != NULL)
			free(Opt->ModelFile);

		Opt->ModelFile = StrMake(Passed[1]);
		Opt->UseModelFile = TRUE;
	}
	else
	{
		printf("ModelFile (MF) command take 0 or 1 parameters, 0 to turn the use of a model file to off. 1 to specify the model file name.\n");
		return;
	}
}

void	SetVarDataFile(OPTIONS *Opt, int Tokes, char** Passed)
{
	if(Tokes == 1)
	{
		if(Opt->VarDataFile != NULL)
		{
			free(Opt->VarDataFile);
			Opt->VarDataFile = NULL;
		}

		Opt->UseVarData = FALSE;
		return;
	}

	if(Tokes == 2)
	{
		if(Opt->VarDataFile != NULL)
			free(Opt->VarDataFile);

		Opt->VarDataFile = StrMake(Passed[1]);
		Opt->UseVarData = TRUE;
	}
	else
	{
		printf("VarData (VD) command take 0 or 1 parameters, 0 to turn the use of a data file to off. 1 to specify the data file name.\n");
		return;
	}
}


void	FlattenRecNode(OPTIONS *Opt)
{
	int		Index;
	RECNODE	RNode;

	if(Opt->NoOfRecNodes == 0)
		return;
	
	Opt->RecNodeList = (RECNODE*)malloc(sizeof(struct RNODE**) * Opt->NoOfRecNodes);
	if(Opt->RecNodeList == NULL)
		MallocErr();

	for(Index=0,RNode=Opt->RecNode;Index<Opt->NoOfRecNodes;Index++, RNode = RNode->Next)
		Opt->RecNodeList[Index] = RNode;
}

void	FreeRecNodes(OPTIONS *Opt)
{
	int	Index;
	RECNODE	R;

	if(Opt->RecNodeList != NULL)
		free(Opt->RecNodeList);

	FlattenRecNode(Opt);

	for(Index=0;Index<Opt->NoOfRecNodes;Index++)
	{
		R = Opt->RecNodeList[Index];


		free(R->TreeNodes);
		free(R->Name);
		free(R->TaxaID);
		free(R->Taxa);		

		free(R);
	}
	
	free(Opt->RecNodeList);
	Opt->NoOfRecNodes = 0;
	Opt->RecNode = NULL;
}

void	SetDataDev(OPTIONS *Opt, int Tokes, char ** Passed)
{
	double	NewVal;

	if(Opt->Analsis == ANALML)
	{
		printf("Missing Data can only be estimated under MCMC.\n");
		return;
	}
/*
	if(EstData(Opt->Trees) == FALSE)
	{
		printf("No Data is being estimated, DataDev cannot be set.\n");
		return;
	}
*/
	if(IsValidDouble(Passed[1]) == FALSE)
	{
		printf("Could not convert %s to a valid DataDev\n", Passed[1]);
		return;
	}

	NewVal = atof(Passed[1]);
	if(NewVal <= 0)
	{
		printf("DataDev must be gratern then 0\n");
		return;
	}

	Opt->EstDataDev = NewVal;
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

void	SetScheduleUse(OPTIONS *Opt)
{
	char *Temp;

	if(Opt->UseSchedule == TRUE)
	{
		free(Opt->ScheduleFile);
		Opt->UseSchedule = FALSE;
		return;
	}
	
	Temp = (char*)malloc(sizeof(char) * BUFFERSIZE);
	if(Temp == NULL)
		MallocErr();

	sprintf(Temp, "%s.Schedule.txt", Opt->DataFN);
	Opt->ScheduleFile = StrMake(Temp);
	free(Temp);

	Opt->UseSchedule = TRUE;
}

void	SetConRateDev(OPTIONS *Opt, char *PName, char *PRateDev)
{
	int		RPos;
	double	RD;

	RPos = StrToRate(Opt, PName);
	if(RPos == -1)
	{
		printf("Could not find parameter %s\n", PName);
		return;
	}

	if(IsValidDouble(PRateDev) == FALSE)
	{
		printf("Could not conver %s to a valid doble\n", PRateDev);
		return;
	}

	RD = atof(PRateDev);
	if(RD < 0)
	{
		printf("%s has to be a value grater then 0.\n", RPos);
		return;
	}

	Opt->RateDevList[RPos] = RD;
}

void	OptSetSeed(OPTIONS *Opt, char	*CSeed)
{
	if(IsValidInt(CSeed) == FALSE)
	{
		printf("%s is not a valid Unsigned integer.\n", CSeed);
		return;
	}

	IntSetSeed(atoi(CSeed));
	Opt->Seed = GetSeed();
}

int	FossilNoPramOK(OPTIONS *Opt, int Tokes)
{
	return TRUE;
}


int		PassLine(OPTIONS *Opt, char *Buffer)
{
	char		*Passed[1024];
	int			Tokes;
	COMMANDS	Command;
	int			Index;
	int			Temp;
	double		TempDouble;

	ReplaceChar(';', ' ', Buffer);
	ReplaceChar('=', ' ', Buffer);
	ReplaceChar(',', ' ', Buffer);
	RemoveChar('\n', Buffer);

	Tokes = MakeArgv(Buffer, Passed, 1024);

	if(Tokes >= 1)
		MakeLower(Passed[0]);

	if(Tokes <= 0)
		return FALSE;

	Command = StringToCommand(Passed[0]);

	if(Command == CUNKNOWN)
		printf("Unknown command: %s\n",Passed[0]);

	if(CmdVailWithDataType(Opt,Command)==FALSE)
		return FALSE;

	if(Command == CRUN)
		return TRUE;

	if(Command == CRES)
	{
		if(Tokes >= 3)
			Restrict(Opt, Tokes, Passed);
		else
			printf("The Restrict command takes two parimiters a paramiter to restict and a constant or rate to restict it to.\n");

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
	{
		UnRestictAll(Opt);
	}

	if(Command == CPRIOR)
	{
		if(Tokes >= 4)
			SetPrior(Opt, Tokes, Passed);
		else
		{
			printf("Prior set the prior values, requires a rate parmeters, adistribution type and a number of parmeters\n");
			printf("E.G., Prior q01 Beta 6 24.5\n");
		}
	}	

	if(Command == CITTERS)
	{
		if(Tokes == 2)
		{
			Temp = atoi(Passed[1]);
			if(Temp == 0)
				printf("Could not convert %s to a valid number of Itterashions. Use -1 for infinite\n", Passed[1]);
			else
				Opt->Itters = Temp;
		}
		else
			printf("Itters requires a number that spesifies the number of itters to run the chain for, use -1 of infineite\n");
	}
	
	
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
		
	if(Command == CINFO)
	{
		PrintOptions(stdout, Opt);
	}
		
	if(Command == CPRIORALL)
	{
		if(Tokes >= 3)
			SetAllPriors(Opt, Tokes, Passed);
		else
			printf("Set all prionrs take a prior distruntion and a number of paramter\n");
	}
	
	if(Command == CHELP)
	{
		Index=0;
		do
		{
			printf("%s\t%s\n", COMMANDSTRINGS[Index], COMMANDSTRINGS[Index+1]);
			Index+=2;
		}while(COMMANDSTRINGS[Index][0] != '\0');
	}

	if((Command == CNODE) || (Command == CMRCA))
	{
		if(Tokes >= 4)
		{
			if(Command == CNODE)
				AddRecNode(Opt, NODEREC, Tokes, Passed);

			if(Command == CMRCA)
				AddRecNode(Opt, MRCA, Tokes, Passed);
		}
		else
		{
			printf("Command must supplie a name and two or more taxa numbers\n");
		}
	}

	if(Command == CDELNODE)
	{
		if(Tokes == 2)
			DelRecNode(Opt, Passed[1]);
		else
			printf("The del node command remove a node, it take the name of the node\n");

	}

	if(Command == CADDTAXA)
	{
		if(Tokes >= 3)
			AddToRecNode(Opt, Tokes, Passed);
		else
			printf("The AddNode command takes at least two parmeters a Node Name and taxa number/s\n");
	}


	if(Command == CDELTAXA)
	{
		if(Tokes >= 3)
			DelToRecNode(Opt, Tokes, Passed);
		else
			printf("The DelTaxa command requies 2 or more parmeters and Node Name and a list of taxa numbers to remove from that node\n");
	}

	if(Command == CEVENROOT)
	{
		SetEvenRoot(Opt->Trees);
	}

	if(Command == CLOGFILE)
	{
		if(Tokes == 2)
			LogFile(Opt, Passed[1]); 
		else
			printf("The LogFile command requires a file name to use a log file\n");
	}

	if(Command == CRATEDEV)
	{
		if(Tokes == 2)
		{
			TempDouble = atof(Passed[1]);
			if(TempDouble <= 0)
				printf("Could not convert %s to a valid Rate deveation\n");
			else
				Opt->RateDev  = TempDouble;

			if(Opt->DataType == CONTINUOUS)
				SetAllRateDevs(Opt, Opt->RateDev);
		}
		else if ((Tokes == 3) && (Opt->DataType == CONTINUOUS))
			SetConRateDev(Opt, Passed[1], Passed[2]);
		else
			printf("The RateDev command requires a floating point number\n");
	}

	if(Command == CPRESET)
	{
		if(Tokes == 2)
		{
			PreSet(Opt, Passed[1]);
		}
		else
			printf("The PreSet command take a preset\n");
	}

	if(Command == CSUMMARY)
	{
		if(Opt->Summary == FALSE)
			Opt->Summary = TRUE;
		else
			Opt->Summary = FALSE;
	}

	if(Command == CBURNIN)
	{
		if((Tokes == 2) && (Opt->Analsis == ANALMCMC))
		{
			if(IsValidInt(Passed[1]) == TRUE)
			{
				Temp = atoi(Passed[1]);
				Opt->BurnIn = Temp;
			}
			else
				printf("The value %s cannot be converted to a valid ingteger\n", Passed[1]);
		}
		else
			printf("The Burn In command take an interger, the number of itteraions needed befor convergence. Only allpicable for MCMC\n");
	}

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
	{
		if(SetConVar(Tokes, Passed, &Opt->UseKappa, &Opt->EstKappa, &Opt->FixKappa ) == FALSE)
			printf("The Kappa command take 0 or 1 parameters, 0 to toggle kappa (on / off) and 1 to fix it to a constant\n");
				
	}

	if(Command == CDELTA)
	{
		if(SetConVar(Tokes, Passed, &Opt->UseDelta, &Opt->EstDelta, &Opt->FixDelta ) == FALSE)
			printf("The Delta command take 0 or 1 parameters, 0 to toggle Delta (on / off) and 1 to fix it to a constant\n");
				
	}

	if(Command == CLAMBDA)
	{
		if(SetConVar(Tokes, Passed, &Opt->UseLambda, &Opt->EstLambda, &Opt->FixLambda ) == FALSE)
			printf("The Lambda command take 0 or 1 parameters, 0 to toggle Lambda (on / off) and 1 to fix it to a constant\n");
	}

	if(Command == CEXTTAXA)
	{
		if(Tokes > 1)
		{
			FreePartitions(Opt->Trees);
			FreeRecNodes(Opt);
			ExcludeTaxa(Opt, Tokes-1, &Passed[1]);
			SetPartitions(Opt->Trees);
		}
		else
		{
			printf("The exclude taxa (et) command requires one or more taxa names or numbers\n");
		}
	}

	if(Command == CTAXAINFO)
	{
		PrintTaxaInfo(Opt);
	}

	if(Command == CSAVETREES)
	{
		if(Tokes == 2)
		{
			PrintTree(Passed[1], Opt->Trees, Opt);
		}
		else
			printf("Save trees requies a file name.\n");
	}

	if(Command == CTESTCORREL)
	{
		if(Opt->TestCorrel == FALSE)
			Opt->TestCorrel = TRUE;
		else
			Opt->TestCorrel = FALSE;
	}

	if(Command == CSURFACE)
	{
	}

	if(Command == CCOVARION)
	{
		if(Opt->UseCovarion == TRUE)
			Opt->UseCovarion = FALSE;
		else
			Opt->UseCovarion = TRUE;
	}

	if(Command == CREVJUMP)
	{
		SetRJMCMC(Opt, Tokes-1, &Passed[1]);
	}

	if(Command == CEXIT)
	{
		exit(0);
	}

	if(Command == CFOSSIL)
	{
		if(FossilNoPramOK(Opt, Tokes) == TRUE)
			AddRecNode(Opt, FOSSIL, Tokes, Passed);
		else
		{
			if(Opt->DataType == DISCRETE)
			{
				printf("Error: The fossil command take a node name, a state to fossilise in and a list of taxa that define the node of interest.\n");
			}
			else
			{
				printf("Error: The fossil command take a node name, a state to fossilise in for each site and a list of taxa that define the node of interests.\n");
			}
		}
	}


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

	if(Command == CHYPERPRIOR)
	{
		if(Tokes > 4)
		{
			SetHyperPrior(Opt, Passed[1], &Passed[2], Tokes-2, NULL);
		}
		else
			printf("HyperPrior requires a rate, a distrubion and a set of upper and lower values for the paramtiers\n");
	}

	if(Command == CHPRJ)
	{
		if(Tokes > 3)
		{
			SetHyperPrior(Opt, NULL, &Passed[1], Tokes-1, Opt->RJPrior);
			Opt->UseRJMCMC = TRUE;
		}
		else
		{
			if(Opt->UseRJMCMC == TRUE)
				Opt->UseRJMCMC = FALSE;
			else
				printf("HPRevJump rquires a distrubion and a set of upper and lower values for the paramtiers\n");
		}
	}

	if(Command == CHPALL)
	{
		if(Tokes > 3)
		{
			SetHPAll(Opt, &Passed[1], Tokes-1);	
		}
		else
		{
			printf("HyperPriorAll rquires a distrubion and a set of upper and lower values for the paramtiers\n");
		}
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
	{
		SetGamma(Opt, Passed, Tokes);
	}

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

	if(Command == CDEPSITE)
	{
		if(Tokes == 2)
		{
			if(IsValidInt(Passed[1]) == TRUE)
			{
				Temp = atoi(Passed[1]);
				if((Temp < 1) || (Temp > Opt->Trees->NoOfSites))
					printf("The designated dependent site must be > 0 and <= %d\n", Opt->Trees->NoOfSites);
				else
					Opt->DependantSite = Temp-1;
			}
			else
				printf("DepSite requires an integers to designated as the dependent value.\n");
		}
		else
			printf("DepSite requires a site number to designated as the dependent value.\n");
	}

	if(Command == CHEADERS)
	{
		if(Opt->Headers == TRUE)
			Opt->Headers = FALSE;
		else
			Opt->Headers = TRUE;

	}

	if(Command == CMODELFILE)
		SetModelFile(Opt, Tokes, Passed);

	if(Command == CVARDATA)
		SetVarDataFile(Opt, Tokes, Passed);

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

	if(Command == CDATADEV)
		SetDataDev(Opt, Tokes, Passed);

	if(Command == CNOSPERSITE)
		SetNOSPerSiteOpt(Opt);

	if(Command == CSCHEDULE)
		SetScheduleUse(Opt); 

	if(Command == CSTREEMOVE)
	{
		if(Opt->SoloTreeMove == TRUE)
			Opt->SoloTreeMove = FALSE;
		else
			Opt->SoloTreeMove = TRUE;
	}

	if(Command == CSETSEED)
	{
		if(Tokes == 2)
			OptSetSeed(Opt, Passed[1]);
		else
			printf("SetSeed take an unsinged intger.\n");
		
	}

	if(Command == CMAKEUM)
	{

		MakeUM(Opt->Trees);
	/*	if(Opt->MakeUM == TRUE)
			Opt->MakeUM = FALSE;
		else
			Opt->MakeUM = TRUE;
	*/
	}

	if(Command == CPHYLOPLASTY) 
	{
		if(Opt->Trees->NoOfTrees > 1)
		{
			printf("PhyloPlasty can only be used on a single tree.\n");
			return FALSE;
		}

		if(Opt->UsePhyloPlasty == FALSE)
			Opt->UsePhyloPlasty = TRUE;
		else
			Opt->UsePhyloPlasty = FALSE;
	}

	return FALSE;
}

void	GetOptionsArry(OPTIONS *Opt, int Size, char** OptStr)
{
	int	Index;

	for(Index=0;Index<Size;Index++)
		PassLine(Opt, OptStr[Index]);
}

void	GetOptions(OPTIONS *Opt)
{
	char	*Buffer;

	Buffer = (char*)malloc(sizeof(char) * BUFFERSIZE);
	if(Buffer == NULL)
		MallocErr();

	do
	{
		printf("> ");
		fgets(Buffer, BUFFERSIZE, stdin); 
	} while(PassLine(Opt,Buffer) == FALSE);
	
	free(Buffer);
}