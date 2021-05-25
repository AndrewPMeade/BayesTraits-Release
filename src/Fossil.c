#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "typedef.h"
#include "Fossil.h"
#include "genlib.h"


void	MakeMSFossilMask(RECNODE *RNode, int *Mask, int NOS)
{
	int Index;

	for(Index=0;Index<RNode->NoFossilStates;Index++)
		Mask[RNode->FossilStates[Index]] = 1;
}

void	MakeDiscretePattern(int *Mask, int s00, int s01, int s10, int s11)
{
	if(s00 == 1)
		Mask[0] = 1;

	if(s01 == 1)
		Mask[1] = 1;

	if(s10 == 1)
		Mask[1] = 1;

	if(s11 == 1)
		Mask[3] = 1;
}

void	MakeDiscreteFossilMask(int *Mask, int FState)
{
	switch(FState)
	{
	case 0:
		MakeDiscretePattern(Mask, 1, 0, 0, 0);
		break;

	case 1:
		MakeDiscretePattern(Mask, 0, 1, 0, 0);
		break;

	case 2:
		MakeDiscretePattern(Mask, 0, 0, 1, 0);
		break;

	case 3:
		MakeDiscretePattern(Mask, 0, 0, 0, 1);
		break;

	case 10:
		MakeDiscretePattern(Mask, 1, 1, 0, 0);
		break;

	case 11:
		MakeDiscretePattern(Mask, 1, 0, 1, 0);
		break;

	case 12:
		MakeDiscretePattern(Mask, 1, 0, 0, 1);
		break;

	case 13:
		MakeDiscretePattern(Mask, 0, 1, 1, 0);
		break;

	case 14:
		MakeDiscretePattern(Mask, 0, 1, 0, 1);
		break;

	case 15:
		MakeDiscretePattern(Mask, 0, 0, 1, 1);
		break;

	case 20:
		MakeDiscretePattern(Mask, 1, 1, 1, 0);
		break;

	case 21:
		MakeDiscretePattern(Mask, 1, 1, 0, 1);
		break;

	case 22:
		MakeDiscretePattern(Mask, 1, 0, 1, 1);
		break;

	case 23:
		MakeDiscretePattern(Mask, 0, 1, 1, 1);
		break;
	}
}

int*	MakeFossilMask(RECNODE *RNode, int NOS, MODEL M)
{
	int *Ret, Index;

	Ret = (int*)SMalloc(sizeof(int) * NOS);

	for(Index=0;Index<NOS;Index++)
		Ret[Index] = 0;

	if(M == M_MULTISTATE)
	{
		MakeMSFossilMask(RNode, Ret, NOS);
		return Ret;
	}

	if(M == M_DESCINDEP || M == M_DESCDEP)
	{
		MakeDiscreteFossilMask(Ret, RNode->FossilStates[0]);
		return Ret;
	}

	printf("CV not supported with Fossil, Please contact support if you need this feature\n");
	exit(1);

	return NULL;
}


void	SetFossils(TREES *Trees, OPTIONS *Opt)
{
	RECNODE	*RNode;
	int		Index, TIndex;
	NODE	N;

	if(Opt->UseCovarion == TRUE)
	{
		printf("CV not supported with Fossil, Please contact support if you need this feature\n");
		exit(1);
	}


	for(Index=0;Index<Opt->NoOfRecNodes;Index++)
	{
		RNode = Opt->RecNodeList[Index];
		if(RNode->NodeType == FOSSIL)
		{
			for(TIndex=0;TIndex<Trees->NoOfTrees;TIndex++)
			{
				N = RNode->TreeNodes[TIndex];
							
				N->FossilMask = MakeFossilMask(RNode, Trees->NoOfStates, Opt->Model);
			}
		}

	}
}


void	FossilLh(NODE N, OPTIONS *Opt, TREES *Trees, int SiteNo)
{
	int	Index, NOS;
	int *Mask;
	
	NOS = Trees->NoOfStates;

	Mask = N->FossilMask;

	for(Index=0;Index<NOS;Index++)
	{
#ifdef BIG_LH
		FossilLhBig(N, Trees, Mask,  SiteNo);
		return;
#endif

#ifdef QUAD_DOUBLE
		N->BigPartial[SiteNo][Index] = N->BigPartial[SiteNo][Index] * (double)Mask[Index];
#else
		N->Partial[SiteNo][Index] = N->Partial[SiteNo][Index] * (double)Mask[Index];
#endif
	}

}