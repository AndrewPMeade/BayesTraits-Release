#ifndef PATTERN_H
#define PATTERN_H

#include "typedef.h"
#include "genlib.h"

void	AddPattern(OPTIONS *Opt, char *Name, int NoTags, char **TagNameList);
void	FreePattern(PATTERN *Pattern);

void	PrintPatterns(FILE *Str, int NoPatterns, PATTERN **PList);

void	SetPatternRateNames(OPTIONS *Opt);

void	SetPatternNo(OPTIONS *Opt, TREES *Trees);


#endif