#ifndef CONTRASTS_TRANS_H
#define CONTRASTS_TRANS_H

#include "typedef.h"


void	TransformTree(OPTIONS *Opt, TREES *Trees, RATES *Rates, int Norm);


void	TransformTreeDelta(NODE N, double Delta, int Norm);
void	TransformTreeKappa(NODE N, double Kappa, int Norm);
void	TransformTreeOU(TREES *Trees, NODE N, double OU, int Norm);
void	TransformTreeLambda(NODE N, double Lambda, int Norm);

int		NeedToReSetBL(OPTIONS *Opt, RATES *Rates);

double	GetTransformDefValue(TRANSFORM_TYPE TranType);

#endif
