#ifndef OPTIONS_H
#define OPTIONS_H

#include "TypeDef.h"

OPTIONS*	CreatOptions(MODEL Model, ANALSIS Analsis, int NOS, char *TreeFN, char *DataFN, char *SymbolList, TREES* Trees);
void		FreeOptions(OPTIONS *Opt, int NoSites);

MODEL		GetModel(TREES *Trees);
ANALSIS		GetAnalsis(TREES *Trees);
MODEL_TYPE	GetModelType(MODEL Model);
void		GetOptions(OPTIONS *Opt, TREES *Trees);

int			ValidModelChoice(TREES *Trees, MODEL Model);

void		PrintOptions(FILE* Str, OPTIONS *Opt, TREES *Trees);
void		GetOptionsArry(OPTIONS *Opt, int Size, char** OptStr, TREES *Trees);

MODEL		IntToModel(int No);

int			DataModifiedOptions(OPTIONS *Opt);

#endif
