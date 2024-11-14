#ifndef MODELFILE_H
#define MODELFILE_H

#include "TypeDef.h"

FILE*			InitSaveModelFile(char *FName, OPTIONS *Opt, TREES *Trees, RATES *Rates);

void			SaveModelFile(FILE *MFile, OPTIONS *Opt, TREES *Trees, RATES *Rates);
MODELFILE*		LoadModelFile(char *FileName, OPTIONS *Opt, TREES *Trees, RATES *Rates);
void			FreeModelFile(MODELFILE *MF);


void		TestModelFile(OPTIONS *Opt, TREES *Trees, RATES *Rates);
void		ChangeModelFile(RATES *Rates, gsl_rng *RNG);
void		MapModelFile(OPTIONS *Opt, RATES *Rates, TREES *Trees);

#endif