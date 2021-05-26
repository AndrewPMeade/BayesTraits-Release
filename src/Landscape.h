#ifndef LANDSCAPE_H
#define LANDSCAPE_H

#include "TypeDef.h"

LANDSCAPE*	CreateLandScape(TREES *Tree);
void		FreeLandScape(LANDSCAPE *Landscape);

int			UseLandScapeModel(OPTIONS* Opt);

void		ResetTreeLandscape(TREE *Tree);
void		MapLandscape(OPTIONS *Opt, TREES *Trees, RATES *Rates);

void		MLLandscape(OPTIONS *Opt, TREES *Trees, RATES *Rates);

#endif 
