#ifndef INTRA_NODE_H
#define INTRA_NODE_H

#include "GenLib.h"
#include "TypeDef.h"

#define INIT_RAND_DIST 200
#define INIT_RAND_TRIES 100

#define CHANGE_RAND_DIST 2000
#define TRIES_PER_NODE_CHANGE 10



void CreateIntraNode(OPTIONS *Opt, TREES *Trees, RATES *Rates);
void CopyIntraNode(RATES *A, RATES *B);

void SetInitIntraNodeLocations(OPTIONS *Opt, TREES *Trees, RATES *Rates);

double CalcIntraNodeLh(TREES *Trees, RATES *Rates);

void ChangeAllIntraNodeLh(OPTIONS *Opt, TREES *Trees, RATES *Rates);

void	InitIntraNodeFile(OPTIONS *Opt, TREES *Trees);
void	OutputIntraNode(size_t Itter, OPTIONS *Opt, TREES *Trees, RATES *Rates, FILE *Str);

void	SetIntraNodeTransformTree(OPTIONS *Opt, TREES *Trees, RATES *Rates);

//void	CheckIntraNodeNodeRestictedMapsLongLat(RESTRICTION_MAP *ResMap, double Long, double Lat);
//int	CheckIntraNodeNodeRestictedMaps(NODE Node, INTRA_NODE *INode);
void	CheckIntraNodeRestictedMaps(TREES *Trees, RATES *Rates);

#endif
