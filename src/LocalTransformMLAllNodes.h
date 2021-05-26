#ifndef		LOCAL_TRANSFORMS_ML_ALL_NODES_H
#define		LOCAL_TRANSFORMS_ML_ALL_NODES_H

#include "TypeDef.h"

#define MAX_LT_NAME_SIZE 128
#define ML_CUT_POINT	1.96

typedef struct
{
	TAG **TagList;
	int	*ValidNodes;
	int NoNodes;

} LT_ALL_NODES;

void	LocalTransformMLAllNodes(OPTIONS *Opt, TREES *Trees, RATES *Rates);

#endif