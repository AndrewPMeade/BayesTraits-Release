#if !defined (LOADTREES)
#define LOADTREES

#include "typedef.h"

TREES*	LoadTrees(char* FileName);

void	FreeTrees(TREES* Trees, OPTIONS *Opt);
void	AllocPartial(TREES* Trees, int Gamma);
TAXA*	GetTaxaFromID(int ID, TAXA *Taxa, int NoOfTaxa);
TAXA*	GetTaxaFromName(char *Name, TAXA *Taxa, int NoOfTaxa);

void	PrintTrees(FILE*	Str, TREES *Trees, DATATYPE DataType);

/* void	PrintNodeRec(FILE *Str, NODE Node, int NOS, int NoOfSites, RATES* Rates); */
double	GetStateProbPct(int State, int NoOfStates, double *Part);

int		SymbolToPos(char Symbol, char *List);
int		SiteHadUnKnownState(char *StatList);
void	AllocConVarCoVar(TREES* Trees);
int		RemoveTaxa(OPTIONS *Opt, TREES *Trees, char *TName);
void	PrintTree(char	*FileName, TREES* Trees, OPTIONS *Opt);
void	SetFossiles(TREES *Trees, OPTIONS *Opt);

char*	TipIDToTaxaName(int TID, TAXA* Taxa);
void	SetMinBL(TREES *Trees);
void	SetNOSPerSite(OPTIONS *Opt);

void	AddNewRecNode(TREES* Trees, RECNODE RecNode);
void	SetNodeTipData(NODE N, TREE* Tree, TREES *Trees);

void	MakeUM(TREES* Trees);

#endif
