#ifndef INVMKL_H
#define INVMKL_H

#include "typedef.h"
#ifdef USE_MLK

#define MKL_INT int
#include <mkl.h>

int			InvMLK(TREES *Trees, TREE *Tree);

#endif
#endif