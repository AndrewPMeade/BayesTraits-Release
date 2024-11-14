
#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <math.h>

#include "IntraNode.h"
//#include "GenLib.h"
#include "Trees.h"
#include "TypeDef.h"
#include "FatTail.h"
#include "Geo.h"
#include "RestrictionMap.h"

void RecGetNoIntraNode(NODE Node, double IntraNodeDist, double Height, double Point, int *NoNodes)
{
	int Index;
	
	while(Point > Height && Point < Height + Node->Length)
	{
		(*NoNodes)++;
		Point += IntraNodeDist;
	}
	
	for(Index=0;Index<Node->NoNodes;Index++)
		RecGetNoIntraNode(Node->NodeList[Index], IntraNodeDist, Height + Node->Length, Point, NoNodes);
}

int GetNoIntraNode(NODE Root, double IntraNodeDist)
{
	int Index, NoNodes;
	double Dist;

	Dist = 0.0;
	NoNodes = 0;

	for(Index=0;Index<Root->NoNodes;Index++)
		RecGetNoIntraNode(Root->NodeList[Index], IntraNodeDist, 0.0, IntraNodeDist, &NoNodes);
		
	return NoNodes;
}

void RecLinkIntraNodesToTree(NODE Node, INTRA_NODE *IntraNodePtr, int *NodeOffSet, double IntraNodeDist, double Height, double Point)
{
	int Index, NoNodes;

	NoNodes = 1;
	Node->IntraNodes =  &IntraNodePtr[*NodeOffSet];

	while(Point > Height && Point < Height + Node->Length)
	{
		NoNodes++;
		Point += IntraNodeDist;
	}

	Node->NoIntraNodes = NoNodes;
	(*NodeOffSet) += NoNodes;

	for(Index=0;Index<Node->NoNodes;Index++)
		RecLinkIntraNodesToTree(Node->NodeList[Index], IntraNodePtr, NodeOffSet, IntraNodeDist, Height + Node->Length, Point);
}

void LinkIntraNodesToTree(TREE *Tree, double IntraNodeDist)
{
	NODE Root;
	int Index, NodeOffSet;
	double Dist;

	Dist = 0.0;
	
	Root = Tree->Root;
	
	Root->IntraNodes = &Tree->IntraNodes[0];
	Root->NoIntraNodes = 1;
	NodeOffSet = 1;

	for(Index=0;Index<Root->NoNodes;Index++)
		RecLinkIntraNodesToTree(Root->NodeList[Index], Tree->IntraNodes, &NodeOffSet, IntraNodeDist, 0.0, IntraNodeDist);
}

void SetTreeInterNodes(TREE *Tree, RATES *Rates)
{
	memcpy(Tree->IntraNodes, Rates->IntraNodes, sizeof(INTRA_NODE) * Rates->NoIntraNodes);
}

void GetTreeInterNodes(TREE *Tree, RATES *Rates)
{
	memcpy(Rates->IntraNodes, Tree->IntraNodes,  sizeof(INTRA_NODE) * Rates->NoIntraNodes);
}

void RecSetInterNodesHeights(NODE Node, double IntraNodeDist, double Height, double Point)
{
	int Index;
	
	Index=0;
	
	while(Point > Height && Point < Height + Node->Length)
	{
		Node->IntraNodes[Index].Height = Point;
		Point += IntraNodeDist;
		Index++;
	}
	
	Node->IntraNodes[Index].Height = Height + Node->Length;

	for(Index=0;Index<Node->NoNodes;Index++)
		RecSetInterNodesHeights(Node->NodeList[Index], IntraNodeDist, Height + Node->Length, Point);
}

double GetMaxIntraNodeHeight(TREE *Tree)
{
	int NIndex;
	double Max;
	NODE Node;

	Max = 0;

	for(NIndex=0;NIndex<Tree->NoNodes;NIndex++)
	{
		Node = Tree->NodeList[NIndex];
		if(Node->IntraNodes[Node->NoIntraNodes-1].Height > Max)
			Max = Node->IntraNodes[Node->NoIntraNodes-1].Height;

	}

	return Max;
}

void RecSetHeightMaxDist(NODE Node, double MaxHeight)
{
	int Index;

	for(Index=0;Index<Node->NoIntraNodes;Index++)
		Node->IntraNodes[Index].Height = MaxHeight - Node->IntraNodes[Index].Height;

	for(Index=0;Index<Node->NoNodes;Index++)
		RecSetHeightMaxDist(Node->NodeList[Index], MaxHeight);
}

void SetInterNodesHeights(TREE *Tree, double IntraNodeDist)
{
	int Index;
	NODE Root;
	double Max;

	Root = Tree->Root;

	Root->IntraNodes[0].Height = 0;

	for(Index=0;Index<Root->NoNodes;Index++)
		RecSetInterNodesHeights(Root->NodeList[Index],  IntraNodeDist, 0.0, IntraNodeDist);

	Max = GetMaxIntraNodeHeight(Tree);

	RecSetHeightMaxDist(Tree->Root, Max);

//	Node->Height = MaxDistToRoot - Node->DistToRoot;


}

void SetInterNodesNodeLengths(NODE Node)
{
	int Index;
	double LastHight;

	LastHight = 0;
	
	if(Node->Ans != NULL)
		LastHight = Node->Ans->IntraNodes[Node->Ans->NoIntraNodes-1].Height;
	else
		LastHight = Node->IntraNodes[0].Height;

	for(Index=0;Index<Node->NoIntraNodes;Index++)
	{
		Node->IntraNodes[Index].Length = LastHight - Node->IntraNodes[Index].Height;
		Node->IntraNodes[Index].InputLength = Node->IntraNodes[Index].Length;

		LastHight = Node->IntraNodes[Index].Height;	
	}
}

void SetInterNodesLengths(TREE *Tree)
{
	int Index;
	
	for(Index=0;Index<Tree->NoNodes;Index++)
		SetInterNodesNodeLengths(Tree->NodeList[Index]);
}

void PrintInterNode(NODE Node)
{
	int Index;
	INTRA_NODE *INode;
	double Long, Lat;


	printf("%d\t", Node->NoIntraNodes);

	for(Index=0;Index<Node->NoIntraNodes;Index++)
	{
		INode = &Node->IntraNodes[Index];

		XYZToLongLat(INode->X, INode->Y, INode->Z, &Long, &Lat);


		printf("%f\t%f\t%f\t%f\t%f\t%f\t%f\t", INode->Height, INode->Length, INode->X, INode->Y, INode->Z, Long, Lat);
	}

	RecPRintNodeTaxa(Node, '\t');
	printf("\n");

	for(Index=0;Index<Node->NoNodes;Index++)
		PrintInterNode(Node->NodeList[Index]);
}

void SetIntraNodeResMaps(OPTIONS *Opt, TREES *Trees, RATES *Rates)
{
	int Index;
	INTRA_NODE *INode;
	
	if(Opt->NoRestrictionMaps == 0)
		return;

	for(Index=0;Index<Rates->NoIntraNodes;Index++)
	{
		INode = &Rates->IntraNodes[Index];
		INode->ResMap = GetMapFromHeight(INode->Height, Opt->RestrictionMaps, Opt->NoRestrictionMaps);
	}

}

void CreateIntraNode(OPTIONS *Opt, TREES *Trees, RATES *Rates)
{
	TREE *Tree;
	int Index;

	Tree = Trees->Tree[0];

	Rates->NoIntraNodes = GetNoIntraNode(Tree->Root, Opt->IntraNodeDist);
	Rates->NoIntraNodes += Tree->NoNodes;

	Rates->IntraNodes = (INTRA_NODE*)SMalloc(sizeof(INTRA_NODE) * Rates->NoIntraNodes);
	for(Index=0;Index<Rates->NoIntraNodes;Index++)
		Rates->IntraNodes[Index].ResMap = NULL;

	if(Tree->IntraNodes == NULL)
		Tree->IntraNodes = (INTRA_NODE*)SMalloc(sizeof(INTRA_NODE) * Rates->NoIntraNodes);
	
	SetTreeInterNodes(Tree, Rates);
	LinkIntraNodesToTree(Tree, Opt->IntraNodeDist);

	SetInterNodesHeights(Tree, Opt->IntraNodeDist);

	SetInterNodesLengths(Tree);

	GetTreeInterNodes(Tree, Rates);

	SetIntraNodeResMaps(Opt, Trees, Rates);
}

void CopyIntraNode(RATES *A, RATES *B)
{
	memcpy(A->IntraNodes, B->IntraNodes, sizeof(INTRA_NODE) * A->NoIntraNodes);
}

void GetNodeCoordinates(NODE Node)
{
	INTRA_NODE *INode;

	INode = &Node->IntraNodes[Node->NoIntraNodes-1];

	NodeToXYZ(Node, &INode->X, &INode->Y, &INode->Z);
}

/* 
	Gets the node Coordinates, form the Anc State and places it into the INTRA_NODE.
*/
void GetAllNodeCoordinates(TREE *Tree)
{
	int Index;

	for(Index=0;Index<Tree->NoNodes;Index++)
		GetNodeCoordinates(Tree->NodeList[Index]);
}

void SetNodeCoordinates(NODE Node)
{
	INTRA_NODE *INode;

	INode = &Node->IntraNodes[Node->NoIntraNodes-1];

	XYZToNode(Node, INode->X, INode->Y, INode->Z);


//	printf("%f\t%f\t%f\n", INode->X, INode->Y, INode->Z);


}
/* 
	Set the node Coordinates, form the the INTRA_NODE and places it in the Anc State.
*/
void SetAllNodeCoordinates(TREE *Tree)
{
	int Index;

	for(Index=0;Index<Tree->NoNodes;Index++)
		SetNodeCoordinates(Tree->NodeList[Index]);
}

void GetRandIntLXYZ(gsl_rng *RNG, double X, double Y, double Z, double *RX, double *RY, double *RZ, double Radius, NODE Node)
{
	int Tried;
	double Long, Lat;
	double RLong, RLat;


	XYZToLongLat(X, Y, Z, &Long, &Lat);

	Tried = 0;
	do
	{
		if(Tried < INIT_RAND_TRIES)
			GetRandLongLat(RNG, Long, Lat, &RLong, &RLat, Radius);
		else
			GetGloablRandLongLat(RNG, &RLong, &RLat);

		Tried++;
	}while(ValidGeoNodeLongLat(Node, RLong, RLat) == FALSE);

	
	LongLatToXYZ(RLong, RLat, RX, RY, RZ);
}


void SetInitIntaNode(NODE Node)
{
	int Index;
	INTRA_NODE *Last, *Current;
	
	for(Index=Node->NoIntraNodes-2;Index>=0;Index--)
	{
		Last = &Node->IntraNodes[Index+1];
		Current = &Node->IntraNodes[Index];

		GetRandIntLXYZ(Node->RNG, Last->X, Last->Y, Last->Z, &Current->X, &Current->Y, &Current->Z, INIT_RAND_DIST, Node);
	}
}

void DumpNodeXYZ(TREE *Tree)
{
	NODE Node;
	INTRA_NODE *INode;
	double X,Y,Z;
	int Index;

	printf("\n\n\n");
	for(Index=0;Index<Tree->NoNodes;Index++)
	{
		Node = Tree->NodeList[Index];
		INode = &Node->IntraNodes[Node->NoIntraNodes-1];
		NodeToXYZ(Node, &X, &Y, &Z);

		printf("%d\t%f\t%f\t%f\t%f\t%f\t%f\n", Index, X, Y, Z, INode->X, INode->Y,INode->Z);
	}

	printf("\n\n\n");
}

void SetInitIntraNodeLocations(OPTIONS *Opt, TREES *Trees, RATES *Rates)
{
	TREE *Tree;
	FATTAILRATES *FTR;
	int Index;

	Tree = Trees->Tree[0];

	FTR = Rates->FatTailRates;

	CheckRestictedMaps(Trees, FTR);

	FatTailSetAnsSates(Tree, Trees->NoSites, FTR);
	SetTreeInterNodes(Tree, Rates);
	GetAllNodeCoordinates(Tree);

	for(Index=0;Index<Tree->NoNodes;Index++)
		SetInitIntaNode(Tree->NodeList[Index]);


	SetAllNodeCoordinates(Tree);
	GetTreeInterNodes(Tree, Rates);
	FatTailGetAnsSates(Tree, Trees->NoSites, FTR);

	CheckIntraNodeRestictedMaps(Trees, Rates);
}

double CalcIntraNodePairLh(double Scale, INTRA_NODE *A, INTRA_NODE *B, double Length)
{
	double Lh;

	Lh = StableDistTPDF(Scale, A->X - B->X, Length); 
	Lh += StableDistTPDF(Scale, A->Y - B->Y, Length); 
	Lh += StableDistTPDF(Scale, A->Z - B->Z, Length); 

	return Lh;
}

double CalcAllIntraNodeNodeLh(double Scale, NODE Node)
{
	double Lh;
	int Index;
	INTRA_NODE *A, *B;

	Lh = 0;
	if(Node->Tip == FALSE)
	{
		A = &Node->IntraNodes[Node->NoIntraNodes-1];
		for(Index=0;Index<Node->NoNodes;Index++)
		{
			B = &Node->NodeList[Index]->IntraNodes[0];
			Lh += CalcIntraNodePairLh(Scale, A, B, B->Length);
		}
	}

	for(Index=Node->NoIntraNodes-1;Index>=1;Index--)
	{
		A = &Node->IntraNodes[Index-1];
		B = &Node->IntraNodes[Index];

		Lh += CalcIntraNodePairLh(Scale, A, B, B->Length);
	}

	return Lh;
}

double CalcIntraNodeLh(TREES *Trees, RATES *Rates)
{
	TREE *Tree;
	double Lh;
	int Index;

	Tree = Trees->Tree[0];

	SetTreeInterNodes(Tree, Rates);
	GetAllNodeCoordinates(Tree);

	Lh = 0;
	for(Index=0;Index<Tree->NoNodes;Index++)
		Lh += CalcAllIntraNodeNodeLh(Rates->Rates[0], Tree->NodeList[Index]);

	return Lh;
}

double	CalcIntraNodeNodeLh(NODE Node, INTRA_NODE *A, double Scale)
{
	double Ret;
	int Index;
	INTRA_NODE *B;

	Ret = 0;
	B = NULL;
	if(Node->NoIntraNodes != 1)
	{
		B = &Node->IntraNodes[Node->NoIntraNodes-2];
	}
	else
	{
		if(Node->Ans != NULL)
			B = &Node->Ans->IntraNodes[Node->Ans->NoIntraNodes-1];
	}


	if(B != NULL)
		Ret += CalcIntraNodePairLh(Scale, A, B, A->Length);
	
	for(Index=0;Index<Node->NoNodes;Index++)
	{
		B = &Node->NodeList[Index]->IntraNodes[0];
		Ret += CalcIntraNodePairLh(Scale, A, B, B->Length);
	}

	return Ret;
}

void	SetIntraNode(INTRA_NODE *A, double X, double Y, double Z)
{
	A->X = X;
	A->Y = Y;
	A->Z = Z;
}

void	ChangeIntraNodeNode(NODE Node, INTRA_NODE *A, double Scale)
{
	double RX, RY, RZ;
	double TX, TY, TZ;
	double CLh, NLh;
	int Tried, Changed;

	Changed = FALSE;
	Tried = FALSE;
	
	CLh = CalcIntraNodeNodeLh(Node, A, Scale);

	TX = A->X; 	TY = A->Y; 	TZ = A->Z;

	do
	{
		SetIntraNode(A, TX, TY, TZ);

		GetRandIntLXYZ(Node->RNG, A->X, A->Y, A->Z, &RX, &RY, &RZ, CHANGE_RAND_DIST, Node);
	
		SetIntraNode(A, RX, RY, RZ);

		NLh = CalcIntraNodeNodeLh(Node, A, Scale);

		if(log(gsl_rng_uniform_pos(Node->RNG)) < (NLh - CLh))
			Changed = TRUE;
		else
			Tried++;
	
	}while(Changed == FALSE && Tried < TRIES_PER_NODE_CHANGE);

	if(Changed == FALSE)
		SetIntraNode(A, TX, TY, TZ);

}

double CalcIntraNodeIntaLh(NODE Node, INTRA_NODE *A, INTRA_NODE *B, INTRA_NODE *C, double Scale)
{
	return  CalcIntraNodePairLh(Scale, C, B, C->Length) + CalcIntraNodePairLh(Scale, A, B, B->Length);
}

// B is the node that is changeing, A is the rootward node and C is the tipward one. 
void	ChangeIntraNodeIntra(NODE Node, INTRA_NODE *A, INTRA_NODE *B, INTRA_NODE *C, double Scale)
{
	double RX, RY, RZ;
	double TX, TY, TZ;
	double CLh, NLh;
	int Tried, Changed;

	Changed = FALSE;
	Tried = 0;
	
	CLh = CalcIntraNodeIntaLh(Node, A, B, C, Scale);

	TX = B->X; 	TY = B->Y; 	TZ = B->Z;

	do
	{
		SetIntraNode(B, TX, TY, TZ);

		GetRandIntLXYZ(Node->RNG, B->X, B->Y, B->Z, &RX, &RY, &RZ, CHANGE_RAND_DIST, Node);
	
		SetIntraNode(B, RX, RY, RZ);

		NLh = CalcIntraNodeIntaLh(Node, A, B, C, Scale);

		if(log(gsl_rng_uniform_pos(Node->RNG)) < (NLh - CLh))
			Changed = TRUE;
		else
			Tried++;
	
	}while(Changed == FALSE && Tried < TRIES_PER_NODE_CHANGE);

	if(Changed == FALSE)
		SetIntraNode(B, TX, TY, TZ);

}

void	RecChangeAllIntraNodeLh(NODE Node, double Scale)
{
	int Index;
	INTRA_NODE *A, *B, *C;

	for(Index=0;Index<Node->NoNodes;Index++)
		RecChangeAllIntraNodeLh(Node->NodeList[Index], Scale);

	if(Node->Tip == FALSE)
		ChangeIntraNodeNode(Node, &Node->IntraNodes[Node->NoIntraNodes-1], Scale);

	for(Index=Node->NoIntraNodes-2;Index>=0;Index--)
	{
		if(Index != 0)
			A = &Node->IntraNodes[Index-1];
		else
			A = &Node->Ans->IntraNodes[Node->Ans->NoIntraNodes-1]; 
		B = &Node->IntraNodes[Index];
		C = &Node->IntraNodes[Index+1];

		ChangeIntraNodeIntra(Node, A, B, C, Scale);
	}
}

void ChangeAllIntraNodeLh(OPTIONS *Opt, TREES *Trees, RATES *Rates)
{
	TREE *Tree;
	FATTAILRATES *FTR;

	Tree = Trees->Tree[Rates->TreeNo];
	FTR = Rates->FatTailRates;

	
	SetTreeInterNodes(Tree, Rates);
	GetAllNodeCoordinates(Tree);

//	CheckIntraNodeRestictedMaps(Trees, Rates);
		

	RecChangeAllIntraNodeLh(Tree->Root, Rates->Rates[0]);

	SetAllNodeCoordinates(Tree);
	GetTreeInterNodes(Tree, Rates);
}


void	OutputIntraNodeHeader(FILE *Str, TREES *Trees, NODE Node, INTRA_NODE *INode, int NodeIndex, int INodeIndex)
{
	int Index;
	TAXA *Taxa;


	fprintf(Str, "Node-%05d-%05d\t", NodeIndex, INodeIndex);
	fprintf(Str, "%f\t", Node->Length);
	fprintf(Str, "%f\t%f\t", INode->Length, INode->Height);

	if(INode->ResMap == NULL)
		fprintf(Str, "None\t");
	else
		fprintf(Str, "%s\t", INode->ResMap->FileName);

	for(Index=0;Index<Node->Part->NoTaxa;Index++)
	{
		Taxa = Trees->Taxa[Node->Part->Taxa[Index]];
		fprintf(Str, "%s\t", Taxa->Name);
	}

	fprintf(Str, "\n");
}

void	InitIntraNodeFile(OPTIONS *Opt, TREES *Trees)
{
	TREE *Tree;
	int Index, IIndex;
	NODE Node;
	INTRA_NODE *INode;

	Opt->LogIntraNode = OpenWriteWithExt(Opt->BaseOutputFN, OUTPUT_EXT_INTRA_NODE);
	Tree = Trees->Tree[0];

	fprintf(Opt->LogIntraNode, "Node\tBranch Length\tSegment Length\tHeight\tRestriction Map\tTaxa\n");
	for(Index=0;Index<Tree->NoNodes;Index++)
	{
		Node = Tree->NodeList[Index];
		for(IIndex=0;IIndex<Node->NoIntraNodes;IIndex++)
		{
			INode = &Node->IntraNodes[IIndex];
			OutputIntraNodeHeader(Opt->LogIntraNode, Trees, Node, INode, Index, IIndex);
		}
	}
	

	fprintf(Opt->LogIntraNode, "Itter\tLh\tScale\t");
	for(Index=0;Index<Tree->NoNodes;Index++)
	{
		Node = Tree->NodeList[Index];
		for(IIndex=0;IIndex<Node->NoIntraNodes;IIndex++)
		{
			fprintf(Opt->LogIntraNode, "Node-%05d-%05d - Seg Length\t", Index, IIndex);
			fprintf(Opt->LogIntraNode, "Node-%05d-%05d - Long\t", Index, IIndex);
			fprintf(Opt->LogIntraNode, "Node-%05d-%05d - Lat\t", Index, IIndex);
		}	
	}

	fprintf(Opt->LogIntraNode, "\n");
	fflush(Opt->LogIntraNode);
}


void	OutputIntraNode(size_t Itter, OPTIONS *Opt, TREES *Trees, RATES *Rates, FILE *Str)
{
	TREE *Tree;
	int Index, IIndex;
	NODE Node;
	double Long, Lat;
	INTRA_NODE *IntraNode;

//	Str = Opt->LogIntraNode;
	
	Tree = Trees->Tree[Rates->TreeNo];


	FatTailSetAnsSates(Tree, Trees->NoSites, Rates->FatTailRates);
	SetTreeInterNodes(Tree, Rates);
	GetAllNodeCoordinates(Tree);
	
	fprintf(Str, "%zu\t%f\t%f\t", Itter, Rates->Lh, Rates->Rates[0]);	


	for(Index=0;Index<Tree->NoNodes;Index++)
	{
		Node = Tree->NodeList[Index];
		for(IIndex=0;IIndex<Node->NoIntraNodes;IIndex++)
		{
			IntraNode = &Node->IntraNodes[IIndex];

			XYZToLongLat(IntraNode->X, IntraNode->Y, IntraNode->Z, &Long, &Lat);
			fprintf(Str, "%f\t%f\t%f\t", IntraNode->Length, Long, Lat);
		}
	}

	fprintf(Str, "\n");
	fflush(Str);
}

void	SetIntraNodeTransformTree(OPTIONS *Opt, TREES *Trees, RATES *Rates)
{
	NODE Node;
	INTRA_NODE *INode;
	int Index, IIndex;
	TREE *Tree;
	double Ratio;

	Tree = Trees->Tree[Rates->TreeNo];

	SetTreeInterNodes(Tree, Rates);

	for(Index=0;Index<Tree->NoNodes;Index++)
	{
		Node = Tree->NodeList[Index];
		Ratio = Node->Length / Node->UserLength;

		for(IIndex=0;IIndex<Node->NoIntraNodes;IIndex++)
		{
			INode = &Node->IntraNodes[IIndex];

			INode->Length = Ratio * INode->InputLength;
		}
	}

	GetTreeInterNodes(Tree, Rates);
}
/*
void	CheckIntraNodeNodeRestictedMapsLongLat(RESTRICTION_MAP *ResMap, double Long, double Lat)
{
	int Valid;

	Valid = ValidResPoint(ResMap, Long, Lat);
	if(Valid == FALSE)
	{
		printf("Invalid Node:\t%f\t%f\n", Long, Lat);
		exit(1);
	}

}

int	CheckIntraNodeNodeRestictedMaps(NODE Node, INTRA_NODE *INode)
{
	int Valid;
	double Long, Lat;

	XYZToLongLat(INode->X, INode->Y, INode->Z, &Long, &Lat);

	Valid = ValidResPoint(INode->ResMap, Long, Lat);

	return Valid;
}
*/
void	CheckIntraNodeRestictedMaps(TREES *Trees, RATES *Rates)
{
	TREE *Tree;
	int Index, Valid, IIndex, MaxINode;
	NODE Node;
	INTRA_NODE *INode;

	double Long, Lat;

	Tree = Trees->Tree[0];

	FatTailSetAnsSates(Tree, Trees->NoSites, Rates->FatTailRates);
	SetTreeInterNodes(Tree, Rates);

	for(Index=0;Index<Tree->NoNodes;Index++)
	{
		Node = Tree->NodeList[Index];
		MaxINode = Node->NoIntraNodes;

		if(Node->Tip == TRUE)
			MaxINode--;

		for(IIndex=0;IIndex<MaxINode;IIndex++)
		{
			INode = &Node->IntraNodes[IIndex];
					

			XYZToLongLat(INode->X, INode->Y, INode->Z, &Long, &Lat);

			printf("Coude will not work with flipped nodes, its not clear if Intra nodes will go forword.\n");
			exit(1);
//			Valid = ValidResPoint(INode->ResMap, Long, Lat);

			if(Valid == FALSE)
			{
				printf("%d\t%d\t%d\t%f\t%f\n", Index, Valid, Node->Tip, Long, Lat);
				printf("Invalid point.\n");
				exit(1);
			}
		}
	}
} 

