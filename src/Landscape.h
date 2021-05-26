#ifndef LANDSCAPE_H
#define LANDSCAPE_H

#include "TypeDef.h"

int			UseLandscapeBeta(OPTIONS* Opt);

void		ResetTreeLandscape(TREE *Tree);
void		MapLandscape(OPTIONS *Opt, TREES *Trees, RATES *Rates);

void		MLLandscape(OPTIONS *Opt, TREES *Trees, RATES *Rates);

#endif 
