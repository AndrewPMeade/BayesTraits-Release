#ifndef VAR_RATES_ML_LOCAL_TRANSFOMR_H
#define VAR_RATES_ML_LOCAL_TRANSFOMR_H

#include "TypeDef.h"



//	After reading in a var rates check point and with a local transfom this function find ML values for each local transfom in tern. 
void	VarRatesMLLocalTransfom(OPTIONS *Opt, TREES *Trees, RATES *Rates);

#endif