#ifndef GEO_HEADDER
#define GEO_HEADDER

#include "TypeDef.h"

// Data should be Long / Lat

void	ValidGeoData(TREES *Trees);
void	PreProcessGeoData(TREES *Trees);

void	LongLatToXYZ(double Long, double Lat, double *X, double *Y, double *Z);
void	XYZToLongLat(double X, double Y, double Z, double *Long, double *Lat);

void	NodeToLongLat(NODE N, double *Long, double *Lat);

void	CorrectIntGeoNodes(TREE *Tree);

void	GeoUpDateAllAnsStates(OPTIONS *Opt, TREES *Trees, RATES *Rates);
void	GeoUpDateAnsStates (OPTIONS *Opt, TREES *Trees, RATES *Rates);

void	LoadGeoData(OPTIONS *Opt, TREES *Trees, RATES *CRates, char *Str);


void	NodeToXYZ(NODE N, double *X, double *Y, double *Z);
void	XYZToNode(NODE N, double X, double Y, double Z);
void	NodeToLongLat(NODE N, double *Long, double *Lat);
void	LongLatToNode(NODE N, double Long, double Lat);

int		ValidGeoNodeLongLat(NODE Node, double Long, double Lat);


void	GetRandLongLat(gsl_rng* RNG, double Long, double Lat, double* RLong, double* RLat, double Radius);
void	GetGloablRandLongLat(gsl_rng *RNG, double *RLong, double *RLat);


void	CheckRestictedMaps(TREES *Trees, FATTAILRATES *FTR);

void	SimGeoData(OPTIONS *Opt, TREES *Trees, RATES *Rates);

#endif