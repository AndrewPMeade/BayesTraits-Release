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
#include "RandLib.h"
#include "threaded.h"
#include "part.h"
#include "rates.h"

#define	RATEOUTPUTLEN	33
#define	RIGHTINDENT		4

void	FreeRecNodes(OPTIONS *Opt, int NoSites);

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

		if((Opt->ResTypes[Index] == RESNONE) && (Opt->UseRJMCMC == FALSE))
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
		Taxa = Trees->Taxa[TIndex];
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
	
	fprintf(Str, "Precision:                       %d bits\n", Opt->Precision);
	fprintf(Str, "Cores:                           %d\n", Opt->Cores);

	if(Opt->Analsis == ANALMCMC)
	{
		fprintf(Str, "Analysis Type:                   MCMC\n" );
		fprintf(Str, "Sample Period:                   %d\n", Opt->Sample);
		fprintf(Str, "Iterations:                      %d\n", Opt->Itters);
		fprintf(Str, "Burn in:                         %d\n", Opt->BurnIn);

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
		}

		if(Opt->UseSchedule	== TRUE)
			fprintf(Str, "Schedule File:                   %s.Schedule.txt\n", Opt->LogFN);

		if(Opt->AutoTuneRD == TRUE)
			fprintf(Str, "Rate Dev:                        AutoTune\n");
		else
		{
			fprintf(Str, "Rate Dev:                        %f\n", Opt->RateDev);

			if(Opt->DataType == CONTINUOUS)
			{
				for(Index=0;Index<Opt->NoOfRates;Index++)
				{
					fprintf(Str, "    ");
					PrintFixSize(Opt->RateName[Index], 29, Str);
					fprintf(Str, "%f\n", Opt->RateDevList[Index]);
				}
			}
		}
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
			case PIEMP:
				fprintf(Str, "Empirical\n");
				break;
			case PIUNI:
				fprintf(Str, "Uniform\n");
				break;
		}

		fprintf(Str, "Character Symbols:               ");

		if(Opt->Model == M_MULTISTATE)
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

		
//		fprintf(Str, "Normalisation Constant:          %20.20f\n", Opt->Trees->NormConst);
	}

	if(Opt->DataType == CONTINUOUS)
	{
		fprintf(Str, "Test for trait correlation:      ");
		if(Opt->TestCorrel == TRUE)
			fprintf(Str, "True\n");
		else
			fprintf(Str, "False\n");

		fprintf(Str, "Kappa                            ");
		PrintConPar(Str, Opt->UseKappa, Opt->EstKappa, Opt->FixKappa);

		fprintf(Str, "Delta                            ");
		PrintConPar(Str, Opt->UseDelta, Opt->EstDelta, Opt->FixDelta);

		fprintf(Str, "Lambda                           ");
		PrintConPar(Str, Opt->UseLambda, Opt->EstLambda, Opt->FixLambda);

		fprintf(Str, "OU                               ");
		PrintConPar(Str, Opt->UseOU, Opt->EstOU, Opt->FixOU);
				
		if(Opt->AlphaZero == TRUE)
			fprintf(Str, "Alpha through zero:              True\n");

		if(Opt->NodeData == TRUE)
			fprintf(Str, "Model for Node Data:             True\n");

		if(Opt->NodeBLData == TRUE)
			fprintf(Str, "Model for Node BLS Data:         True\n");

		if(Opt->UseVarData == TRUE)
			fprintf(Str, "Varable Data form file:          %s\n", Opt->VarDataFile);
	
		if(Opt->UseVarRates == TRUE)
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

	if(Opt->SaveModels == TRUE)
		fprintf(Str, "Save Model:                      %s\n", Opt->SaveModelsFN);

	if(Opt->LoadModels == TRUE)
		fprintf(Str, "Load Model:                      %s\n", Opt->LoadModelsFN);

	PrintEstData(Str, Opt);

	if(Opt->AnalyticalP == TRUE)
		fprintf(Str, "Analytical P:                    True\n");
	
	if(Opt->Model == M_DESCHET)
	{
		fprintf(Str, "Tree 1 Partitions :			   \t");
	//	PrintTreePart(Str, Opt->Trees, 0);
	}

	PrintOptRes(Str, Opt);
	if(Opt->Analsis == ANALMCMC)
		PrintPriorOpt(Str, Opt);

	RNode = Opt->RecNode;
	while((RNode != NULL) && (Opt->Model != M_CONTINUOUSRR) && (Opt->Model != M_CONTINUOUSDIR))
	{
		if(RNode->NodeType == MRCA)
			fprintf(Str, "MRCA:         %s                    %f\n", RNode->Name, FindAveNodeDepth(RNode, Opt));
			
		
		if(RNode->NodeType == NODEREC)
			fprintf(Str, "Node:         %s                    %f\n", RNode->Name, ((double)RNode->Hits / Opt->Trees->NoOfTrees)*100);
		

		if(RNode->NodeType == FOSSIL)
			fprintf(Str, "Fossil:       %s                    %f (%d)\n", RNode->Name, FindAveNodeDepth(RNode, Opt), RNode->FossilState);
		
		for(Index=0;Index<RNode->Part->NoTaxa;Index++)
			fprintf(Str, "             %d\t%s\n", RNode->Taxa[Index]->No, RNode->Taxa[Index]->Name);

		RNode = RNode->Next;
	}

	if(Opt->Trees->NoOfRemovedTaxa != 0)
	{
		fprintf(Str, "Removed taxa:\n");
		for(Index=0;Index<Opt->Trees->NoOfRemovedTaxa;Index++)
			fprintf(Str, "          %s\n", Opt->Trees->RemovedTaxa[Index]);
	}


	PrintTreesInfo(Str, Opt->Trees, Opt->DataType);
	fflush(stdout);
}

void	FreeOptions(OPTIONS *Opt, int NoSites)
{
	int		Index;
	
	if((Opt->Model == M_MULTISTATE) || (Opt->DataType == CONTINUOUS))
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


	if(Opt->EstDataSites != NULL)
		free(Opt->EstDataSites);

	if(Opt->RateDevList != NULL)
		free(Opt->RateDevList);

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

	if(Opt->SaveTrees != NULL)
		free(Opt->SaveTrees);
	
	FreeRecNodes(Opt, NoSites);

	if(Opt->SaveModelsFN != NULL)
		free(Opt->SaveModelsFN);

	if(Opt->LoadModelsFN != NULL)
		free(Opt->LoadModelsFN);

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
	int		No, Index;


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
		{
			No = Index+1;
			sprintf(Buffer, "Alpha-%d", No);
		}
		else
		{
			No = (Index - Opt->Trees->NoOfSites) + 1;
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

	Opt->NoOfRates = Opt->Trees->NoOfSites;

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
	

	Opt->NoOfRates = Opt->Trees->NoOfSites;

	NOS = Opt->Trees->NoOfSites;
	 
	Ret = (char**)malloc(sizeof(char**) * Opt->NoOfRates);
	Buffer = (char*)malloc(sizeof(char*) * BUFFERSIZE);
	if((Ret == NULL) || (Buffer == NULL))
		MallocErr();
	
	i = 0;
	for(Index=0;Index<Opt->Trees->NoOfSites;Index++)
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
	
	Opt->NoOfRates = Opt->Trees->NoOfSites * 2;

	NOS = Opt->Trees->NoOfSites;
	 
	Ret = (char**)malloc(sizeof(char**) * Opt->NoOfRates);
	Buffer = (char*)malloc(sizeof(char*) * BUFFERSIZE);
	if((Ret == NULL) || (Buffer == NULL))
		MallocErr();
	
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
	
	Opt->NoOfRates = Opt->Trees->NoOfSites;
	 
	Ret = (char**)malloc(sizeof(char**) * Opt->NoOfRates);
	Buffer = (char*)malloc(sizeof(char*) * BUFFERSIZE);
	if((Ret == NULL) || (Buffer == NULL))
		MallocErr();
	
	Pos = 0;
	sprintf(Buffer, "Alpha");
	Ret[Pos++] = StrMake(Buffer);

	for(Index=1;Index<Opt->Trees->NoOfSites;Index++)
	{
		sprintf(Buffer, "Beta-%d", Index);
		Ret[Pos++] = StrMake(Buffer);
	}

	free(Buffer);
	return Ret;
}

char**	CreatContinusRateName(OPTIONS* Opt)
{
	switch(Opt->Model)
	{
		case M_CONTINUOUSRR:
			return ModelARateName(Opt);

		case M_CONTINUOUSDIR:
			return ModelBRateName(Opt);

		case M_CONTINUOUSREG:
			return RetModelRateName(Opt);

		case M_CONTRAST_STD:
			return ContrastRateNames(Opt);

		case M_CONTRAST_REG:
			return ContrastRegRateNames(Opt);

		case M_CONTRAST_FULL:
			return ContrastFullRateNames(Opt);
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

	if(Opt->Model == M_DESCINDEP)
	{
		Opt->NoOfRates	= 4;
		Opt->RateName	= INDEPPRAMS;
		return;
	}

	if(Opt->Model == M_DESCDEP)
	{
		Opt->NoOfRates	= 8;
		Opt->RateName	= DEPPRAMS;
	}

	if(Opt->Model == M_DESCCV)
	{
		Opt->NoOfRates = 14;
		Opt->RateName = DEPCVPRAMS;
	}

	if(Opt->Model == M_DESCHET)
	{
		Opt->NoOfRates = 12;
		Opt->RateName = DEPHETROPRAMS;
	}

	if(Opt->Model == M_MULTISTATE)
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
	{
		P->DistVals[0] = 0;
		P->DistVals[1] = 100;
	}
	else
	{
		P->DistVals[0] = -100;
		P->DistVals[1] = 100;
	}
	
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

MODEL_TYPE	GetModelType(MODEL Model)
{
	switch(Model)
	{
		case	M_MULTISTATE:		return MT_DISCRETE; break;
		case	M_DESCINDEP:		return MT_DISCRETE; break;
		case	M_DESCDEP:			return MT_DISCRETE; break;
		case	M_CONTINUOUSRR:		return MT_CONTINUOUS; break;
		case	M_CONTINUOUSDIR:	return MT_CONTINUOUS; break;
		case	M_CONTINUOUSREG:	return MT_CONTINUOUS; break;
		case	M_CONTRAST_STD:		return MT_CONTRAST; break;
		case	M_CONTRAST_REG:		return MT_CONTRAST; break;
		case	M_CONTRAST_FULL:	return MT_CONTRAST; break;
		case	M_DESCCV:			return MT_DISCRETE; break;
		case	M_DESCHET:			return MT_DISCRETE; break;
	}

	printf("Unkown model type (%s::%d)\n", __FILE__, __LINE__);
	exit(1);

	return MT_DISCRETE;
}

OPTIONS*	CreatOptions(MODEL Model, ANALSIS Analsis, int NOS, char *TreeFN, char *DataFN, char *SymbolList, TREES* Trees)
{
	OPTIONS *Ret;
	char	*Buffer;

		
	if((Model == M_DESCDEP) || (Model == M_DESCINDEP) || (Model == M_DESCCV) || (Model == M_DESCHET))
		SquashDep(Trees);

	if((GetModelType(Model) == MT_CONTINUOUS) || (GetModelType(Model) == MT_CONTRAST))
		RemoveConMissingData(Trees);

	Ret = (OPTIONS*)malloc(sizeof(OPTIONS));
	if(Ret == NULL)
		MallocErr();

	Buffer = (char*)malloc(sizeof(char) * BUFFERSIZE);
	if(Buffer == NULL)
		MallocErr();

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

	Ret->PPTree		= NULL;
	Ret->PPLog		= NULL;

	
	Ret->UseRModel	= FALSE;
	Ret->RModelP	= -1;
	Ret->EstDataDev	= 0.2;
	
	Ret->VarRatesScaleDev = PPSCALEDEV;
	Ret->AutoTuneVarRates = FALSE;

	Ret->NoEstDataSite	=	0;
	Ret->EstDataSites	=	NULL;
	Ret->NoEstChanges	=	5;

	Ret->NOSPerSite		=	FALSE;
	Ret->RateDevList	=	NULL;


	if(Ret->ModelType == MT_DISCRETE);
		Ret->DataType = DISCRETE;

	if((Ret->ModelType == MT_CONTINUOUS) || (Ret->ModelType == MT_CONTRAST))
	{
		Ret->TestCorrel = TRUE;
		Ret->DataType	= CONTINUOUS;
	}

	SetOptRates(Ret, NOS, SymbolList);
	
	AllocRestictions(Ret);
		
	Ret->TreeFN = StrMake(TreeFN);
	Ret->DataFN = StrMake(DataFN);

	sprintf(Buffer, "%s.%s", DataFN, LOGFILEEXT);
	Ret->LogFN = StrMake(Buffer);

	Ret->LogFile		= NULL;
	Ret->LogFileRead	= NULL;
	Ret->LogFileBuffer	= NULL;
	Ret->PassedOut		= NULL;		
	Ret->UseSchedule	= FALSE;

	Ret->MLTries		= 10;
	Ret->MCMCMLStart	= FALSE; 
	Ret->AutoTuneRD		= FALSE;
	Ret->AutoTuneDD		= FALSE;
	Ret->RateDevPerParm	= TRUE;

	if(Ret->Analsis == ANALML)
	{
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
		
		Ret->Itters		=	5050000;
		Ret->BurnIn		=	50000;

		Ret->Itters		=	1010000;
		Ret->BurnIn		=	10000;
		Ret->Sample		=	1000;
		/*
		Ret->BurnIn		=	10000;
		Ret->Sample		=	100;
		Ret->Itters		=	1000;
		Ret->Itters		=	101000;
		*/

		Ret->PriorCats	=	100;
		Ret->RateDev	=	1;
		Ret->AutoTuneRD	=	TRUE;
		Ret->EstDataDev	=	.05;
		AllocPrios(Ret);

		Ret->UseSchedule	= TRUE;
	}

	if(Ret->DataType == CONTINUOUS)
	{
		Ret->RateDevList = (double*)malloc(sizeof(double)  * Ret->NoOfRates);
		if(Ret->RateDevList == NULL)
			MallocErr();
		SetAllRateDevs(Ret, Ret->RateDev);
	}
	else
	{
		Ret->RateDevList = (double*)malloc(sizeof(double));
		if(Ret->RateDevList == NULL)
			MallocErr();
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
	Ret->UseOU			=	FALSE;

	Ret->EstKappa		=	FALSE;
	Ret->EstDelta		=	FALSE;
	Ret->EstLambda		=	FALSE;
	Ret->EstGamma		=	FALSE;
	Ret->EstOU			=	FALSE;
//	Ret->EstOU			=	TRUE;

	Ret->FixKappa		=	-1;
	Ret->FixDelta		=	-1;
	Ret->FixLambda		=	-1;
	Ret->FixGamma		=	-1;
	Ret->FixOU			=	-1;

	Ret->RateDevKappa	=	1.0;
	Ret->RateDevLambda	=	1.0;
	Ret->RateDevDelta	=	1.0;
	Ret->RateDevOU		=	1.0;

	Ret->InvertV		=	FALSE;

	Ret->PriorGamma		=	NULL;
	Ret->PriorKappa		=	NULL;
	Ret->PriorLambda	=	NULL;
	Ret->PriorDelta		=	NULL;

	Ret->UseRJMCMC		=	FALSE;
	Ret->CapRJRatesNo	=	-1;

	Ret->FindCF			=	FALSE;
	Ret->CFRate			=	NULL;

	Ret->Headers		=	TRUE;

	Ret->VarData		=	NULL;
	Ret->UseVarData		=	FALSE;
	Ret->VarDataFile	=	NULL;
	

	Ret->AnalyticalP	=	FALSE;


	Ret->Seed			=	GetSeed();

	Ret->MakeUM			=	FALSE;

	Ret->UseVarRates	=	FALSE; 

	Ret->UseEqualTrees	=	FALSE;
	Ret->ETreeBI		=	-1;

	Ret->SaveTrees		=	NULL;


	Ret->Precision		=	sizeof(double) * 8;
#ifdef BIG_LH
	Ret->Precision		=	256;
#endif

	Ret->Cores			=	GetMaxThreads();

	Ret->SaveModels		=	FALSE;
	Ret->SaveModelsFN	=	NULL;

	Ret->LoadModels		=	FALSE;
	Ret->LoadModelsFN	=	NULL;


	free(Buffer);
	return Ret; 
}

void	PrintModelChoic(TREES *Trees)
{
	printf("Please Select the model of evolution to use.\n");
	printf("1)	MultiState\n");
	if((Trees->NoOfSites == 2) && (Trees->NoOfStates == 2))
	{
		printf("2)	Discrete: Independent\n");
		printf("3)	Discrete: Dependant\n");
	}

	if(Trees->ValidCData == TRUE)
	{
		printf("4)	Continuous: Random Walk (Model A)\n");
		printf("5)	Continuous: Directional (Model B)\n");

		if(Trees->NoOfSites > 1)
			printf("6)	Continuous: Regression\n");

		printf("7)	Independent Contrast\n");

		if(Trees->NoOfSites > 1)
			printf("8)	Independent Contrast: Regression\n");

		printf("9)	Independent Contrast: Full\n");
	}

#ifndef PUBLIC_BUILD
	printf("8)	Discrete: Covarion\n");
	printf("9)	Discrete: Heterogeneous \n");
#endif
}

int		GetModelInt()
{
	char	*Buffer;
	int		Ret;

	Ret = -1;

	Buffer = (char*)malloc(sizeof(char) * BUFFERSIZE);
	if(Buffer == NULL)
		MallocErr();

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

int		ValidModelChoice(TREES *Trees, int ModelNo)
{

	int MaxModelNo;

	
	MaxModelNo = 10;

#ifdef PUBLIC_BUILD
	MaxModelNo = 9;
#endif

	if((ModelNo < 1) || (ModelNo > MaxModelNo))
	{
		printf("Model must be 1-8.\n");
		return FALSE;
	}

	if(ModelNo == 1)
		return TRUE;

	if((ModelNo == 2) || (ModelNo == 3) || (ModelNo == 10) || (ModelNo == 11))
	{
		if((Trees->NoOfSites != 2) || (Trees->NoOfStates != 2))
		{
			printf("Discrete analisis requiers two two state characters\n");
			printf("There are %d states and %d sites in the current data set.\n", Trees->NoOfStates, Trees->NoOfSites);
			return FALSE;
		}
		return TRUE;
	}

	if(Trees->ValidCData == FALSE)
	{
		printf("Model %d requires continues data.\n", ModelNo);
		return FALSE;
	}

	if((ModelNo == 6) && (Trees->NoOfSites < 2))
	{
		printf("Continuous Regression, requires two or more sites.\n");
		return FALSE;
	}

	if((ModelNo == 8) && (Trees->NoOfSites < 2))
	{
		printf("Regression, requires two or more sites.\n");
		return FALSE;
	}

	return TRUE;
}

MODEL	IntToModel(int No)
{
	if(No == 1)
		return M_MULTISTATE;

	if(No == 2)
		return M_DESCINDEP;

	if(No == 3)
		return M_DESCDEP;

	if(No == 4)
		return M_CONTINUOUSRR;

	if(No == 5)
		return M_CONTINUOUSDIR;

	if(No == 6)
		return M_CONTINUOUSREG;

	if(No == 7)
		return M_CONTRAST_STD;

	if(No == 8)
		return M_CONTRAST_REG;

	if(No == 9)
		return M_CONTRAST_FULL;

	if(No == 10)
		return M_DESCCV;

	if(No == 11)
		return M_DESCHET;

	printf("Unkown model\n");
	exit(0);
}


//./Seq/MamBigTrim.trees ./Seq/MamBigTrim.txt < BigMamIn.txt > sout.txt
void	GetModelChoice(TREES *Trees, MODEL *Model, int *Valid)
{
	int		No;

	*Valid = FALSE;

	No = GetModelInt();
	if(No == -1)
		return;	

	if(ValidModelChoice(Trees, No) == FALSE)
		return;

	*Valid = TRUE;
	*Model = IntToModel(No);
}

MODEL	GetModel(TREES *Trees)
{
	MODEL	Model;
	int		Valid;

	do
	{
		PrintModelChoic(Trees);
		
		GetModelChoice(Trees, &Model, &Valid);
	} while(Valid == FALSE);

	
	return Model;
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
	
	BackRes		= (RESTYPES*)malloc(sizeof(RESTYPES) * Opt->NoOfRates);
	BackConst	= (double*)malloc(sizeof(double) * Opt->NoOfRates);
	BackNo		= (int*)malloc(sizeof(int) * Opt->NoOfRates);

	if((BackRes == NULL) || (BackConst == NULL) || (BackNo == NULL))
		MallocErr();
	
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
	int			Index;

	MakeLower(Str);
	
	Index = 0;
	do
	{
		if(strcmp(Str, DISTNAMES[Index])==0)
			return (PRIORDIST)(BETA+Index);
		Index++;
	} while(DISTNAMES[Index][0] != '\0');

	return (PRIORDIST)-1;
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

	if((Opt->DataType == CONTINUOUS) && (Dist != UNIFORM))
	{
		printf("Only uniform priors can be used with continuous data.\n");
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
	int	Index;

	if(N->Tip == TRUE)
		printf("%s\t", N->Taxa->Name);
	else
	{
		for(Index=0;Index<N->NoNodes;Index++)
			PrintUnderNode(N->NodeList[Index]);
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

	if(Opt->Model == M_MULTISTATE)
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
		if(strcmp(ID, Trees->Taxa[Index]->Name) == 0)
			return Trees->Taxa[Index];

	return NULL;
}					  

char**	SetConFState(OPTIONS *Opt, NODETYPE NodeType, char *argv[])
{
	int		Index;
	char	**Ret;

	if(Opt->DataType == DISCRETE) 
		return NULL;

	Ret = (char**)malloc(sizeof(char*) * Opt->Trees->NoOfSites);
	if(Ret == NULL)
		MallocErr();


	for(Index=0;Index<Opt->Trees->NoOfSites;Index++)
	{
		if(NodeType != FOSSIL)
			Ret[Index] = StrMake(ESTDATAPOINT); 
		else
			Ret[Index] = StrMake(argv[Index]);
	}

	return Ret;
}

void	AddRecNode(OPTIONS *Opt, NODETYPE NodeType, int Tokes, char *argv[])
{
	RECNODE		RNode;
	int			Index, FState, NoTaxa;
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

	if(Opt->DataType == CONTINUOUS)		
	{
		ConFState = SetConFState(Opt, NodeType, &argv[Index]);
		if(NodeType == FOSSIL)
			Index += Opt->Trees->NoOfSites;
	}
	
	if(ValidTaxaList(argv, Index, Tokes, Opt) == FALSE)
		return;

	RNode = (RECNODE)malloc(sizeof(struct RNODE));
	if(RNode == NULL)
		MallocErr();

	RNode->Part		= NULL;
	RNode->ConData	= NULL;
	if(Opt->DataType == CONTINUOUS)		
		RNode->ConData	= ConFState;

	RNode->Name = (char*)malloc(sizeof(char)*(strlen(argv[1])+1));
	if(RNode->Name == NULL)
		MallocErr();
	strcpy(RNode->Name, argv[1]);

	RNode->NodeType		= NodeType;
	RNode->FossilState	= FState;

	if(RNode->NodeType == FOSSIL)
		RNode->NodeType = MRCA;

	if(NodeType == FOSSIL)
		NoTaxa = Tokes - 3;
	else
		NoTaxa= Tokes - 2;

	RNode->Taxa = (TAXA**)malloc(sizeof(TAXA*)*NoTaxa);
	if(RNode->Taxa == NULL)
		MallocErr();

	for(Index=0;Index<NoTaxa;Index++)
	{
		if(NodeType == FOSSIL)
			RNode->Taxa[Index] = GetTaxaFromNameNo(argv[Index+3], Opt->Trees);
		else
			RNode->Taxa[Index] = GetTaxaFromNameNo(argv[Index+2], Opt->Trees);
	}

	RNode->Part = CreatPart(NoTaxa);
	
	RNode->Next = Opt->RecNode;
	Opt->RecNode = RNode;

	RNode->TreeNodes = (NODE*)malloc(sizeof(NODE)*Opt->Trees->NoOfTrees);
	if(RNode->TreeNodes == NULL)
		MallocErr();

	SetRecNodes(RNode, Opt->Trees);
	
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


void	SetEvenRoot(TREES *Trees)
{
	int		TIndex;
	int		NIndex;
	double	t;
	NODE	Root;

	for(TIndex=0;TIndex<Trees->NoOfTrees;TIndex++)
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

void	LogFile(OPTIONS *Opt, char *LogFN)
{
	FILE*	TempLogFile;

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

	if(strcmp(Set, "m1p")==0)
	{
		RestrictAll(Opt, Opt->RateName[0]);
		Opt->AnalyticalP = TRUE;
	}
}

void	GetBasePis(OPTIONS *Opt, char* Type)
{
	MakeLower(Type);

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

	printf("The option %s, is unknown. Valid options are est, emp, uni and none\n", Type);
}

int		CmdVailWithDataType(OPTIONS *Opt, COMMANDS	Command)
{
	if(Opt->DataType == CONTINUOUS)
	{
		if(Opt->Model == M_CONTRAST_STD)
		{
			if(Command == CNODE)
				return TRUE;
		}

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
//			(Command == CPRIORALL)	||
//			(Command == CPRIOR)		||
			(Command == CPIS)		||
			(Command == CPRECISION) ||
			(Command == CNOSPERSITE)|| 
			(Command == CSYMMETRICAL)
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
			(Command == COU)		||
			(Command == CALPHAZERO)	||
			(Command == CNODEBLDATA)||
			(Command == CNODEDATA)  ||
			(Command == CDEPSITE)   ||
			(Command == CDATADEV)	||
			(Command == CPHYLOPLASTY) ||
			(Command == CVARRATES)
	
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
			(Command ==	CITTERS)		||
			(Command ==	CBURNIN)		||
			(Command ==	CSAMPLE)		||
			(Command ==	CHYPERPRIOR)	||
			(Command ==	CRATEDEV)		||
			(Command ==	CHPRJ)			||
			(Command ==	CHPALL)			||
			(Command ==	CREVJUMP)		||
			(Command == CMCMCMLSTART)	||
			(Command == CPHYLOPLASTY)	||
			(Command == CCAPRJRATES)	||
			(Command == CVARRATES)		||
			(Command == CSAVEMODELS)	||
			(Command == CLOADMODELS)
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
		sprintf(&Buffer[0], "%d", Trees->Taxa[TIndex]->No);
		PrintFixSize(&Buffer[0], 5, stdout);
		printf("%s", Trees->Taxa[TIndex]->Name);
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
	int	Index;


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
			printf("Gamma shape parmiter must be grater than %f and less than %f\n", (double)GAMMAMIN, (double)GAMMAMAX);
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

void	FreeRecNodes(OPTIONS *Opt, int NoSites)
{
	int	Index, CIndex;
	RECNODE	R;

	if(Opt->RecNodeList != NULL)
		free(Opt->RecNodeList);

	FlattenRecNode(Opt);

	for(Index=0;Index<Opt->NoOfRecNodes;Index++)
	{
		R = Opt->RecNodeList[Index];

		if(R->ConData != NULL)
		{
			for(CIndex=0;CIndex<NoSites;CIndex++)
				free(R->ConData[CIndex]);
			free(R->ConData);
		}

		free(R->TreeNodes);
		free(R->Name);
		FreePart(R->Part);
		
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

	if(Tokes == 0)
	{
		Opt->AutoTuneDD = TRUE;
		return;
	}
	
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
	Opt->AutoTuneDD = FALSE;
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
		printf("Could not convert %s to a valid double\n", PRateDev);
		return;
	}

	RD = atof(PRateDev);
	if(RD < 0)
	{
		printf("%s has to be a value grater then 0.\n", PRateDev);
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

	Opt->Seed = atoi(CSeed);
}

int	FossilNoPramOK(OPTIONS *Opt, char **Passed, int Tokes)
{
	if(Tokes < 5)
		return FALSE;

	return TRUE;
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

//#if !(CLIK_P || THREADED)
#ifndef CLIK_P
#ifndef THREADED
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
		
	NOS = Opt->Trees->NoOfStates;

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

	if(Opt->Model == M_DESCINDEP)
	{
		/* Beta 1 = Alpha 1 */
		ResRateNo(Opt, 1, 0);

		/* Beta 2 = Alpha 2 */
		ResRateNo(Opt, 3, 2);
	}

	if(Opt->Model == M_DESCDEP)
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

void	SetRateDev(OPTIONS *Opt, int Tokes, char **Passed)
{
	double Temp;

	if(Tokes == 0)
	{
		Opt->AutoTuneRD = TRUE;
		return;
	}

//	MakeLower(Passed[0]);
		
	if(Tokes == 1)
	{
		Temp = atof(Passed[0]);
	
		if(Temp <= 0)
			printf("Could not convert %s to a valid Rate deveation\n", Passed[1]);
		else
			Opt->RateDev  = Temp;

		if(Opt->DataType == CONTINUOUS)
			SetAllRateDevs(Opt, Opt->RateDev);

		Opt->AutoTuneRD = FALSE;
	}
	else if ((Tokes == 2) && (Opt->DataType == CONTINUOUS))
	{
		SetConRateDev(Opt, Passed[0], Passed[1]);
		Opt->AutoTuneRD = FALSE;
	}
	else
		printf("The RateDev command requires a floating point number\n");
}

void	LineAddErr(TREES *Trees, char *Line)
{
	char *Buffer;
	char **Passed;
	int		Tokes, ID;
	double	Err;

	Buffer = StrMake(Line);
	Passed = (char**)malloc(sizeof(char*) * strlen(Line));
	if(Passed == NULL)
		MallocErr();

	Tokes = MakeArgv(Buffer, Passed, strlen(Line));

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

void	LoadAddErr(OPTIONS *Opt, char *FName)
{
	TEXTFILE *TF;
	int	Index;

	TF = LoadTextFile(FName, FALSE);

	for(Index=0;Index<TF->NoOfLines;Index++)
		LineAddErr(Opt->Trees, TF->Data[Index]);
	
	FreeTextFile(TF);
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
//	Tokes = MakeArgv(Buffer, Passed, 1024);

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
		printf("Unknown command: %s\n",Passed[0]);

	if(CmdVailWithDataType(Opt,Command) == FALSE)
		return FALSE;

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
		printf("The AddTaxa command is no longer supported.\n");
	
	/*	if(Tokes >= 3)
			AddToRecNode(Opt, Tokes, Passed);
		else
			printf("The AddNode command takes at least two parmeters a Node Name and taxa number/s\n");*/
	}


	if(Command == CDELTAXA)
	{
		printf("The DelTaxa command is no longer supported.\n");
/*		if(Tokes >= 3)
			DelToRecNode(Opt, Tokes, Passed);
		else
			printf("The DelTaxa command requies 2 or more parmeters and Node Name and a list of taxa numbers to remove from that node\n");*/
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
		SetRateDev(Opt, Tokes-1, &Passed[1]);
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

	if(Command == COU)
	{
		if(SetConVar(Tokes, Passed, &Opt->UseOU, &Opt->EstOU, &Opt->FixOU) == FALSE)
			printf("The OU command take 0 or 1 parameters, 0 to toggle OU (on / off) and 1 to fix it to a constant\n");
	}

	if(Command == CEXTTAXA)
	{
		if(Tokes > 1)
		{
			FreeParts(Opt->Trees);
			FreeRecNodes(Opt, Opt->Trees->NoOfSites);
			ExcludeTaxa(Opt, Tokes-1, &Passed[1]);
			SetParts(Opt->Trees);
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
			Opt->SaveTrees = StrMake(Passed[1]);
		//	PrintTree(Passed[1], Opt->Trees, Opt);
		}
		else
		{
			if(Opt->SaveTrees != NULL)
			{
				free(Opt->SaveTrees);
				Opt->SaveTrees = NULL;
			}
			else
				printf("Save trees requies a file name.\n");
		}
	}

	if(Command == CTESTCORREL)
	{
		SetTestCorrel(Opt);
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
		if(FossilNoPramOK(Opt, Passed, Tokes) == TRUE)
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

	}

	if(Command == CHEADERS)
	{
		if(Opt->Headers == TRUE)
			Opt->Headers = FALSE;
		else
			Opt->Headers = TRUE;

	}

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

	if((Command == CPHYLOPLASTY) || (Command == CVARRATES))
	{
		if(Opt->Trees->NoOfTrees > 1)
		{
			printf("VarRates can only be used on a single tree.\n");
			
		}
		else
		{
			if(Opt->UseVarRates == FALSE)
			{
				Opt->UseVarRates = TRUE;
				Opt->AutoTuneVarRates = TRUE;
			}
			else
			{
				Opt->UseVarRates = FALSE;
				Opt->AutoTuneVarRates = FALSE;
			}
		}
	}

	if(Command == CEQUALTREES)
	{
		SetEqualTrees(Opt, Tokes, Passed);
	}

	if(Command == CPRECISION)
	{
#ifndef BIG_LH
		printf("Precision is only valid with the Big Lh build of BayesTraits.\n");
		return FALSE;
#endif
		SetPrecision(Opt, Passed[1]);
	}

	if(Command == CCORES)
	{
		SetCores(Opt, Tokes, Passed);
	
		return FALSE;
	}

	if(Command == CSYMMETRICAL)
	{
		SetSymmetrical(Opt);
		return FALSE;
	}

	if(Command == CMCMCMLSTART)
	{
		SetMCMCMLStart(Opt);
		return FALSE;
	}

	if(Command == CCAPRJRATES)
	{
		CapRJRatesNo(Opt, Tokes ,Passed);
		return FALSE;
	}

	if(Command == CSAVEMODELS)
	{
		SetSaveModels(Opt, Tokes, Passed);
		return FALSE;
	}

	if(Command == CLOADMODELS)
	{
		SetLoadModels(Opt, Tokes, Passed);
		return FALSE;
	}

	if(Command == CADDERR)
	{
		if(Tokes != 2)
		{
			printf("AddErr requires one parameter, a file with taxa names and error.\n");
			printf("File names cannot contain spaces.\n");
			return FALSE;
		}
		
		LoadAddErr(Opt, Passed[1]);
	}

	return FALSE;
}

void	GetOptionsArry(OPTIONS *Opt, int Size, char** OptStr)
{
	int		Index;
	char	**Passed;

	Passed = (char**)malloc(sizeof(char*)*BUFFERSIZE);
	if(Passed == NULL)
		MallocErr();

	for(Index=0;Index<Size;Index++)
		PassLine(Opt, OptStr[Index], Passed);

	free(Passed);
}

void	GetOptions(OPTIONS *Opt)
{
	char	*Buffer;
	char	**Passed;

	Passed = (char**)malloc(sizeof(char*) * BUFFERSIZE);
	Buffer = (char*)malloc(sizeof(char) * BUFFERSIZE);
	if((Buffer == NULL) || (Passed == NULL))
		MallocErr();

	do
	{
		printf("> ");
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

	if((Opt->UseRJMCMC == FALSE) && (Opt->DataType == DISCRETE))
	{
		if(FindNoOfRates(Opt) > 25)
		{
			printf("To many free parameter to estimate (%d), try reducing the model or using RJ MCMC\n", NoFreeP);
			printf("If you believe you data can support this number of free parameter please contact the developers to have this limitation removed.\n");
			exit(0);
		}	
	}
}