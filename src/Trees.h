#if !defined LOADTREES
#define LOADTREES

#include "TypeDef.h"

TREES*	LoadTrees(char* FileName);

void	ReSetTaxaID(TREES *Trees);

void	FreeTrees(TREES* Trees, OPTIONS *Opt);
void	AllocPartial(OPTIONS *Opt, TREES* Trees, int Gamma);
TAXA*	GetTaxaFromID(int ID, TAXA **Taxa, int NoTaxa);
TAXA*	GetTaxaFromName(char *Name, TAXA **Taxa, int NoTaxa);

void	CTaxaBelow(NODE N, int *No);
void	PrintTreesInfo(FILE*	Str, TREES *Trees, DATATYPE DataType);

/* void	PrintNodeRec(FILE *Str, NODE Node, int NOS, int NoOfSites, RATES* Rates); */
double	GetStateProbPct(int State, int NoStates, double *Part);

int		SymbolToPos(char Symbol, char *List);
int		SiteHadUnKnownState(char *StatList);
void	RemoveTaxa(TREES *Trees, char *TName);
void	CheckDelTaxa(OPTIONS *Opt, TREES *Trees, char *TName);

void	SaveTrees(char	*FileName, TREES* Trees);

void	SetFossils(TREES *Trees, OPTIONS *Opt);

void	SetMinBL(TREES *Trees);
void	SetNOSPerSite(OPTIONS *Opt);

void	AddNewRecNode(TREES* Trees, RECNODE *RecNode);
void	SetNodeTipData(OPTIONS *Opt, NODE N, TREE* Tree, TREES *Trees);

void	MakeUM(TREES* Trees);

void	SetNodeIDs(TREE* Tree);

void	ListOddPPTaxa(TREES *Trees);

double	FindTreeNormalise(TREES *Trees);
void	NormaliseTrees(double NormC, TREES *Trees);

void	SetVisitedNode(NODE N, int Val);
void	SetVisitedTree(TREE *Tree, int Val);

double*	GetPhyChanges(TREES *Trees, TREE *Tree, double RateDev, RANDSTATES *RS);

void	AddTaxaErr(TREES *Trees, int TaxaID, double Err);

int		TaxaIndexToNo(TREES *Trees, int Index);
int		TaxaNoToIndex(TREES *Trees, int ID);

void	SaveUserBrachLengths(TREES *Trees);
void	SetUserBranchLength(TREE *Tree);
void	SetTreesDistToRoot(TREES *Trees);
void	SetTreeDistToRoot(TREE *Tree);
void	RecSetDistToRoot(NODE N);

void	RecScaleSubTree(NODE N, double Scale);
double	SumNodeBL(NODE N);
void	ScaleSubTree(NODE N, double Scale);
void	ScaleTrees(TREES *Trees, double Scale);
void	ScaleUserTrees(TREES *Trees, double Scale);

NODE	GetTreeTaxaNode(TREE *Tree, int TaxaNo);

void	InitialiseOutputTrees(OPTIONS *Opt, TREES *Trees);
void	OutputTree(OPTIONS *Opt, TREES *Trees, RATES *Rates, long long No, FILE *Out);

void	PrintTreeBL(TREE *Tree);
void	PrintTreeNode(TREES *Trees, TREE *Tree);
void	RecPrintNode(NODE N);



void	CheckSingleDescendent(TREES *Trees);

int		GetNoInternalNodes(TREE *Tree);

#endif
