#ifndef TAG_H
#define TAG_H

#include "TypeDef.h"
TAG*	GetTagFromName(OPTIONS *Opt, char *Name);

void	AddTag(OPTIONS *Opt, int Tokes, char **Passed);
TAG*	CreateTag(TREES *Trees, char *Name, int NoTaxa, char **TaxaNames);

void	FreeTag(TAG *Tag);

void	PrintTags(FILE *Out, OPTIONS *Opt);
void	PrintTag(FILE *Out, TAG *Tag);



#endif