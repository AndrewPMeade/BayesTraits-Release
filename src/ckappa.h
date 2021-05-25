#if !defined (CKAPPAHEADDER)
#define CKAPPAHEADDER

void	InitCKappa(OPTIONS	*Opt, TREES	*Trees);

void	InitCKappaTree(TREES *Trees, TREE* Tree);
void	FreeCKappaTree(CONVAR *ConVar, int NoOfTaxa);

void	KappaVarCoVar(TREES	*Trees, TREE *Tree);

void	MakeKappaV(TREES *Trees, TREE	*Tree, double Kappa);



#endif

