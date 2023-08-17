#ifndef INITIALISE_H
#define INITIALISE_H


#include "TypeDef.h"

OPTIONS*	SetUpOptions(TREES* Trees, char	*TreeFN, char *DataFN);
void		PreProcess(OPTIONS *Opt, TREES* Trees);
void		Finalise(OPTIONS *Opt, TREES* Trees); 

#endif
