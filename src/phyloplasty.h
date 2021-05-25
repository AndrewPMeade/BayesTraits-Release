#if !defined (PHYLOPLASTYHEADDER)
#define PHYLOPLASTYHEADDER

#include "typedef.h"

PPCOVARV*		InitPhyloPlastyConVar(TREES	*Trees, TREE *Tree);
void			MapPhyloPlastyToV(TREES	*Trees, TREE *Tree);

void			PhyloPlasyMove(RATES *Rates);
PHYLOPLASTY*	CreatPhyloPlasty(OPTIONS *Opt, RATES *Rates);

void			SetPhyloPlastyV(TREES* Trees, RATES* Rates);
#endif



