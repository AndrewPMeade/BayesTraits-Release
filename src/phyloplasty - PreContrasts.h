#if !defined (PHYLOPLASTYHEADDER)
#define PHYLOPLASTYHEADDER

#include "typedef.h"

PPCOVARV*		InitPhyloPlastyConVar(TREES	*Trees, TREE *Tree);
void			MapPhyloPlastyToV(TREES	*Trees, TREE *Tree);

void			PhyloPlasyMove(RATES *Rates);
PHYLOPLASTY*	CreatPhyloPlasty(OPTIONS *Opt, RATES *Rates);

void			SetPhyloPlastyV(TREES* Trees, RATES* Rates);
void			CopyPhyloPlasty(PHYLOPLASTY *A, PHYLOPLASTY *B);

void			BlankPP(PHYLOPLASTY *PP);
void			AddPPNode(RATES *Rates, NODE N, double Scale);

void			PhyloPlasyMove(RATES *Rates);

#endif



