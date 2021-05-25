#ifndef QUAD_DOUBLE_H
#define QUAD_DOUBLE_H

void	InitQuadDoubleLh(OPTIONS *Opt, TREES *Trees);
void	FreeQuadLh(OPTIONS *Opt, TREES *Trees);

void	NodeLhQuadDouble(NODE N, OPTIONS *Opt, TREES *Trees, int SiteNo);

void	FossilDepLhQuadDobule(NODE N, int SiteNo, int s00, int s01, int s10, int s11);


double	CombineQuadDoubleLh(RATES* Rates, TREES *Trees, OPTIONS *Opt, int SiteNo, int NOS);
void	SetQuadDoubleNodeRec(NODE N, int NOS, int NoSites, RATES *Rates, OPTIONS *Opt);
#endif
