#ifndef STONES_H
#define STONES_H

#include "TypeDef.h"

void	PrintStones(FILE *Str, STONES *Stones);
void	OutputStoneHeadder(FILE *Out, STONES *Stones);

int		StonesStarted(STONES *Stones, long long Itter);

double	GetStoneHeat(STONES *Stones, long long Itter, double Heat);

void	FreeStones(STONES *Stones);
STONES*	CratesStones(int NoS, int Sample, double Alpha, double Beta);

int		ChangeSample(STONES *Stones, int Itters);
void	StoneItter(STONES *Stones, long long Itter, double Lh, FILE *Out);
int		StoneExit(STONES *Stones, long long Itters);

#endif