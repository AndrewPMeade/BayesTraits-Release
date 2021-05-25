#ifndef TAG_H
#define TAG_H

#include "typedef.h"

TAG*	GetTagFromName(OPTIONS *Opt, char *Name);
void	AddTag(OPTIONS *Opt, int Tokes, char **Passed);

#endif