#ifndef PART_H
#define PART_H

#include "typedef.h"


void	SetParts(TREES *Trees);
void	FreeParts(TREES *Trees);

void	FreePart(PART *Part);
PART*	CreatPart(int NoTaxa);

void	SetRecNodes(RECNODE RNode, TREES *Trees);
NODE	FindNode(RECNODE RNode, TREE *Tree, int *Depth);

void	GetPartDiff(PART *Ans, PART *Cur, PART *Diff);

void	PrintPart(FILE *Str, TREES *Trees, PART *Part);

PART*	CreatePart(TREES *Trees, int NoTaxa, char **TaxaList);

NODE	PartGetMRCA(TREE *Tree, PART *Part);


#endif