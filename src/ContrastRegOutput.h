#ifndef CON_REG_OUTPUT_H
#define CON_REG_OUTPUT_H

#include "TypeDef.h"
#include "GenLib.h"

typedef struct
{
	double	*Predicated;
	double	*Residuals;

	double	**Contrasts;

	int		NoSites;
	int		NoContrasts;

	double	MSE;
	double	*SSE;
	double	*R2;
	double	*SSRes;
	double	*SSTotal;

	double *SE;

	double	GR2;

} REG_CON_POST;

void	OutputConReg(FILE *Str, OPTIONS *Opt, TREES *Trees, RATES *Rates);

#endif
