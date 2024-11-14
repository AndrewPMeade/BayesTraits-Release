#ifndef STONES_H
#define STONES_H

#include "TypeDef.h"

STONES_OPTIONS*	CrateStonesOptions(size_t NoStones, size_t ItPerStone, double Alpha, double Beta);
void	FreeStonesOptions(STONES_OPTIONS* StoneOpt);

void	PrintStones(FILE *Str, STONES *Stones);
void	OutputStoneHeadder(FILE *Out, STONES *Stones);

int		StonesStarted(STONES *Stones, size_t Itter);

double	GetStoneHeat(STONES *Stones, size_t Itter, double Heat);

void	FreeStones(STONES *Stones);
STONES*	CratesStones(size_t NoS, size_t Sample, double Alpha, double Beta);

int		ChangeSample(STONES *Stones, size_t Itters);
void	StoneItter(STONES *Stones, size_t Itter, double Lh, FILE *Out);
int		StoneExit(STONES *Stones, size_t Itters);

#endif