#ifndef BTLAPACK_INTERFACE_H
#define BTLAPACK_INTERFACE_H

#ifdef BTLAPACK

#include "btlapack.h"

void	btlapack_FindInvV(TREES *Trees, TREE* Tree);
void	btlapack_InitConTree(TREES *Trees, TREE* Tree);
void	btlapack_FreeConTree(TREES *Trees, TREE* Tree);
#endif

#endif
