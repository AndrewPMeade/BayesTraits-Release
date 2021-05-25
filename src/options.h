#if !defined (OPTIONSHEADDER)
#define OPTIONSHEADDER

#include "typedef.h"

void		FlattenRecNode(OPTIONS *Opt);

OPTIONS*	CreatOptions(MODEL Model, ANALSIS Analsis, int NOS, char *TreeFN, char *DataFN, char *SymbolList, TREES* Trees);
MODEL		GetModel(TREES *Trees);
ANALSIS		GetAnalsis(TREES *Trees);
MODEL_TYPE	GetModelType(MODEL Model);
void		GetOptions(OPTIONS *Opt);

int			ValidModelChoice(TREES *Trees, int ModelNo);

void		PrintOptions(FILE* Str, OPTIONS *Opt);
void		FreeOptions(OPTIONS *Opt, int NoSites);
void		GetOptionsArry(OPTIONS *Opt, int Size, char** OptStr);

void		CheckOptions(OPTIONS *Opt);

MODEL		IntToModel(int No);

#endif
