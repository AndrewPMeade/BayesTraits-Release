#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "GenLib.h"
#include "TypeDef.h"
#include "Matrix.h"
#include "AncStateV.h"
#include "Trees.h"

#include "RandDists.h"

void	SetNodeID(TREE *Tree)
{
	int NIndex, NID;
	NODE N;

	NID = 0;

	for(NIndex=0;NIndex<Tree->NoNodes;NIndex++)
	{
		N = Tree->NodeList[NIndex];
		if(N->Ans != NULL)
		{
			N->ID = NID++;
			printf("%f\t",N->DistToRoot);

			if(N->Tip == TRUE)
				printf("%s\n", N->Taxa->Name);
			else
				printf("x%d\n", N->ID);
		}
		else
			N->ID = -1;

	}
}

NODE	GetNodeFromID(TREE *Tree, int ID)
{
	int NIndex;

	for(NIndex=0;NIndex<Tree->NoNodes;NIndex++)
		if(Tree->NodeList[NIndex]->ID == ID)
			return Tree->NodeList[NIndex];

	return NULL;
}


double		NodeDistToRoot(NODE N)
{
	return N->DistToRoot;
}

void		RetSetNodeVisit(TREE *Tree)
{
	int Index;

	for(Index=0;Index<Tree->NoNodes;Index++)
		Tree->NodeList[Index]->Visited = FALSE;
}

void		TraceNode(NODE N)
{
	N->Visited = TRUE;
	if(N->Ans == NULL)
		return;
	TraceNode(N->Ans);
}


NODE		GetCommonNode(NODE N)
{
	do
	{
		if(N->Visited == TRUE)
			return N;
		N = N->Ans;
	}while(N);

	return NULL;
}

double		NodeToNodeDist(TREE *Tree, int ID_x, int ID_y)
{
	NODE X, Y, C;

	X = GetNodeFromID(Tree, ID_x);
	if(ID_x == ID_y)
		return NodeDistToRoot(X);

	Y = GetNodeFromID(Tree, ID_y);

	RetSetNodeVisit(Tree);
	
	TraceNode(X);
	C = GetCommonNode(Y);

	return NodeDistToRoot(C); 
}


void	CaclAncVMatrix(double **V, int N, TREE *Tree)
{
	int x,y;


	for(x=0;x<N;x++)
	{
		for(y=0;y<N;y++)
		{
			V[x][y] = NodeToNodeDist(Tree, x,y);
		//	printf("%f\t",V[x][y]);
		}
		//printf("\n");
	}
}

void	PrintXVect(TREE *Tree)
{
	int Index;
	NODE N;

	printf("X = {");
	for(Index=1;Index<Tree->NoNodes;Index++)
	{
		N = Tree->NodeList[Index];

		if(N->Tip == TRUE)
			printf("%f", N->Taxa->ConData[0]);
		else
			printf("x%d", N->ID);

		if(Index != Tree->NoNodes-1)
			printf(",");
	}

	printf("};\n");
}

void PrintEstValues(TREE *Tree)
{
	int Index;
	NODE N;

	
	for(Index=1;Index<Tree->NoNodes;Index++)
	{
		N = Tree->NodeList[Index];

		if(N->Tip == FALSE)
		{
			printf("x%d=.;", N->ID);
			
		}
	}
	printf("\n");
}

int*	GetEstNodeID(TREE *Tree, int *NoEstNodes)
{
	int Index, NoEst;
	int *Ret;
	NODE Node;

	NoEst = 0;
	Ret = (int*)SMalloc(sizeof(int) * Tree->NoNodes);

	for(Index=1;Index<Tree->NoNodes;Index++)
	{
		Node = Tree->NodeList[Index];

		if(Node->Tip == FALSE)
			Ret[NoEst++] = Node->ID;
	}

	(*NoEstNodes) = NoEst;
	return Ret;
}

void	PrintEstMax(int NoEstNodes, int *EstID)
{
	int Index;

	printf("FindMaximum[Lh, {");

	for(Index=0;Index<NoEstNodes-1;Index++)
		printf("x%d,", EstID[Index]);
	printf("x%d", EstID[Index]);

	printf("}]\n");
}


void	MultiVarNorma(MATRIX *V)
{
	MATRIX *Res;
	
	Res = MultivariateNormal(100, V);

	PrintMatrix(Res, "Res=", stdout);
}



void MakeAncVMatrix(OPTIONS *Opt, TREES *Trees)
{
	MATRIX *VMat;
	int N;
	TREE *Tree;
	int NoEstNodes;
	int *EstID;

	Tree = Trees->Tree[0];

	SetTreesDistToRoot(Trees);

	N = Tree->NoNodes - 1;

	VMat = AllocMatrix(N,N);
	

	SetNodeID(Tree);

	CaclAncVMatrix(VMat->me, N, Tree);

	PrintMathematicaMatrix(VMat, "V=", stdout);
		
	PrintEstValues(Tree);
	PrintXVect(Tree);
	
	printf("InvV = Inverse[V];\n");
	printf("Sig2 = 2.0;u=0;\n");
	printf("Xu = X - u;\n");

	printf("Lh = Log[(1/Sqrt[2.0*Pi*Sig2])*E^((-Xu.InvV.Xu) / (2. Sig2))]\n");
	
	EstID = GetEstNodeID(Tree, &NoEstNodes);
	PrintEstMax(NoEstNodes, EstID);



	MultiVarNorma(VMat);

	FreeMatrix(VMat);
	

	exit(0);
}