#ifndef CONTINUOUSHEADDER
#define CONTINUOUSHEADDER

void	InitContinus(OPTIONS *Opt, TREES* Trees);
void	InitContinusTree(OPTIONS *Opt, TREES* Trees, int TreeNo);
NODE	TaxaToNode(TREES* Trees, TREE *Tree, TAXA *Taxa);
double	MLFindAlphaMean(TREES* Trees, TREE *Tree, int Site);
double	LHRandWalk(OPTIONS *Opt, TREES* Trees, RATES* Rates);
void	FreeConVar(CONVAR* ConVar, int NoTaxa);
void	TreeBLToPower(TREES *Trees, TREE *Tree, double Power);
int		FindInvV(TREES *Trees, TREE* Tree);

void	CheckZeroTaxaBL(TREES *Trees);

void	CalcZ(TREES* Trees, TREE *Tree, OPTIONS *Opt);
/*	double	MLFindAlphaReg(TREES* Trees, TREE *Tree);	*/
double	MLFindAlphaReg(TREES* Trees, TREE *Tree, double *Data);

MATRIX*	FindRegVar(TREES *Trees, RATES* Rates, int AlphaZero);

void	FreeTempConVars(TEMPCONVAR* TempCon);

void	SetNodesVPos(TREES *Trees, TREE *Tree);

void	MapPartID(TREES* Trees, PART *Part, int *Map);
int		GetMapID(TREES *Trees, int ID);

void	VToTree(MATRIX *V, TREE *Tree);
void	TreeToV(TREES* Trees, TREE *Tree, MATRIX *V);
#endif
