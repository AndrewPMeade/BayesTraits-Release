#ifndef LANDSCAPE_H
#define LANDSCAPE_H

#include "TypeDef.h"

LANDSCAPE*	CreateLandscape(OPTIONS *Opt);
void		FreeLandscape(LANDSCAPE *Landscape);

void		BlankLandscape(LANDSCAPE *Land);
void		CopyLandscape(LANDSCAPE *A, LANDSCAPE *B, int NoTrees);

void		MapLandscape(OPTIONS *Opt, TREES *Trees, RATES *Rates);

void		AddLandscapeFromPart(LANDSCAPE *Land, PART *Part, TREES *Trees, double Beta);

#endif 
