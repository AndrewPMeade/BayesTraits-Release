#ifndef FATTAIL_HEADDER
#define FATTAIL_HEADDER

#include "TypeDef.h"

//#define FAT_TAIL_ML_PARAM
#ifdef FAT_TAIL_ML_PARAM
	#define FAT_TAIL_ML_SIG2	1.874589821111	
	#define FAT_TAIL_ML_ROOT	0.073277330319
#endif

void			InitFatTailTrees(OPTIONS *Opt, TREES *Trees);

void			CheckFatTailBL(TREES *Trees);

void			MapRatesToFatTailRate(RATES *Rates, FATTAILRATES *FatTailRates);
void			MapFatTailRateToRates(RATES *Rates, FATTAILRATES *FatTailRates);

void			FatTailSetAnsSates(TREE *Tree, int NoSites, FATTAILRATES *FTR);
void 			FatTailGetAnsSates(TREE *Tree, int NoSites, FATTAILRATES *FTR);

FATTAILRATES*	CreateFatTailRates(OPTIONS *Opt, TREES *Trees);
void			FreeFatTailRates(FATTAILRATES* FTR, int NoSites, int NoCores);

void			MutateFatTailRates(OPTIONS *Opt, TREES* Trees, RATES* Rates, SCHEDULE*	Shed);

void			FreeFatTailTree(FATTAILTREE *FatTailTree);

void			CopyFatTailRates(TREES *Trees, FATTAILRATES *A, FATTAILRATES *B);

double			CalcTreeStableLh(OPTIONS *Opt, TREES *Trees, RATES *Rates);

void			SSAnsStatesFatTail(OPTIONS *Opt, TREES *Trees, RATES *Rates);
void			SSAllAnsStatesFatTail(OPTIONS *Opt, TREES *Trees, RATES *Rates);

void			InitFattailFile(OPTIONS *Opt, TREES *Trees);
void			OutputFatTail(long long Itter, OPTIONS *Opt, TREES *Trees, RATES *Rates);

void			FatTailTest(int argc, char **argv);

void			LoadRatesFromStr(char *Str, RATES *Rates, OPTIONS *Opt, TREES *Trees);

#endif