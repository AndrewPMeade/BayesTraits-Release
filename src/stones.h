#ifndef STONES_H
#define STONES_H

#include "typedef.h"
void	PrintStones(FILE *Str, STONES *Stones);
void	OutputStoneHeadder(FILE *Out, STONES *Stones);

int		StonesStarted(STONES *Stones, int Itter);

double	GetStoneHeat(STONES *Stones, int Itter, double Heat);

void	FreeStones(STONES *Stones);
STONES*	CratesStones(int Start, int K, int Sample, double Alpha, double Beta);

int		ChangeSample(STONES *Stones, int Itters);
void	StoneItter(STONES *Stones, int Itter, double Lh, FILE *Out);
int		StoneExit(STONES *Stones, int Itters);


#endif