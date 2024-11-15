#ifndef RJDUMMY_H
#define RJDUMMY_H

#include "TypeDef.h"

int			GetMaxDummy(OPTIONS *Opt, TREES *Trees);

RJDUMMY*	CreatRJDummyCode(OPTIONS *Opt, TREES *Trees);
void		FreeRJDummyCode(RJDUMMY *RJDummy);
void		RJDummyCopy(RATES *A, RATES *B);

void		MapDummyCodes(TREES *Trees, RATES *Rates);

void		RJDummyMove(size_t Itters, OPTIONS *Opt, TREES *Trees, RATES *Rates);
void		RJDummyMoveNode(OPTIONS *Opt, TREES *Trees, RATES *Rates);

void		InitRJDummyFile(OPTIONS *Opt, TREES *Trees);
void		PrintRJDummy(size_t Itter, OPTIONS *Opt, TREES *Trees, RATES *Rates);

void		BuildDummyCodeBeta(RJDUMMY *RJDummy);

void		RJDummyChange(OPTIONS *Opt, TREES *Trees, RATES *Rates);

void		TestDummyCodeSig(OPTIONS *Opt, TREES *Trees, RATES *Rates);

#endif