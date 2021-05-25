#ifndef CONTRASTS_TRANS_H
#define CONTRASTS_TRANS_H

#include "typedef.h"

void	TransformContrastTree(OPTIONS *Opt, TREES *Trees, RATES *Rates);
void	TransformContrastTreeFixed(OPTIONS *Opt, TREES *Trees);

void	TransContNodeDelta(NODE N, double Delta, int Norm);
void	TransContNodeKappa(NODE N, double Kappa, int Norm);
void	TransContNodeOU(NODE N, double OU, int Norm);
void	TransContNodeLambda(NODE N, double Lambda, int Norm);

int		NeedToReSetBL(OPTIONS *Opt);

#endif
