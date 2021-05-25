#ifndef PART_H
#define PART_H

#include "TypeDef.h"


void	SetParts(TREES *Trees);
void	FreeParts(TREES *Trees);

void	FreePart(PART *Part);
PART*	CreatPart(int NoTaxa);

void	GetPartDiff(PART *Ans, PART *Cur, PART *Diff);

void	PrintPart(FILE *Str, TREES *Trees, PART *Part);
void	PrintPartTaxaOnly(FILE *Str, TREES *Trees, PART *Part);

PART*	CreatePart(TREES *Trees, int NoTaxa, char **TaxaList);

NODE	PartGetMRCA(TREE *Tree, PART *Part);

int		TaxaInPart(int TaxaIndex, PART *P);
int		PartEqual(PART *A, PART *B);

int		PartSubSet(PART *A, PART *B);


#endif