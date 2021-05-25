#if !defined (OPTIONSHEADDER)
#define OPTIONSHEADDER

#include "typedef.h"

void	FlattenRecNode(OPTIONS *Opt);

OPTIONS*	CreatOptions(MODEL Model, ANALSIS Analsis, int NOS, char *TreeFN, char *DataFN, char *SymbolList, TREES* Trees);
MODEL		GetModel(TREES *Trees);
ANALSIS		GetAnalsis(TREES *Trees);
void		GetOptions(OPTIONS *Opt);

void		PrintOptions(FILE* Str, OPTIONS *Opt);
void		FreeOptions(OPTIONS *Opt);
void		GetOptionsArry(OPTIONS *Opt, int Size, char** OptStr);

#endif



