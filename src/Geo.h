#ifndef GEO_HEADDER
#define GEO_HEADDER

#include "typedef.h"

// Data should be Long / Lat

void	ValidGeoData(TREES *Trees);
void	PreProcessGeoData(TREES *Trees);

void	LongLatToXYZ(double Long, double Lat, double *X, double *Y, double *Z);
void	XYZToLongLat(double X, double Y, double Z, double *Long, double *Lat);

void	NodeToLongLat(NODE N, double *Long, double *Lat);

void	CorrectIntGeoNodes(TREE *Tree);

void	GeoTest(OPTIONS *Opt, TREES *Trees, RATES *Rates);

void	GeoUpDateAllAnsStates(OPTIONS *Opt, TREES *Trees, RATES *Rates);

void	LoadGeoData(char *Str, OPTIONS *Opt, TREES *Trees, RATES *CRates);
#endif