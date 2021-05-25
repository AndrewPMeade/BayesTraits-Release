#ifndef BTOCL_CONTINUOUS_H
#define BTOCL_CONTINUOUS_H

#ifdef BTOCL_CON
#include "btocl_runtime.h"
#include "btocl_runtime_kernels.h"
#include "btocl_lin.h"

void	btocl_FindInvV(TREES *Trees, TREE* Tree);
void    btocl_VectByKroneckerMult(TREE* Tree);
void	btocl_AllocConVar(CONVAR* ConVar, TREES *Trees);
void	btocl_FreeConVar(CONVAR* ConVar);
#endif

#endif
