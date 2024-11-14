#ifndef GEO_MAP_H
#define GEO_MAP_H

#include "TypeDef.h"

#define GLOBAL_AGE_VAL -1


RESTRICTION_MAP*	LoadResMap(char *FileName, double AgeMin, double AgeMax);
void		FreeResMap(RESTRICTION_MAP* GeoMap);
void		PrintResMap(RESTRICTION_MAP* GeoMap);

RESTRICTION_MAP* GetMapFromHeight(double Height, RESTRICTION_MAP** Maps, int NoMaps);


RESTRICTION_POINT*	GetClosestPoint(RESTRICTION_MAP* GeoMap, double Long, double Lat);
RESTRICTION_POINT*	GetClosestPointBox(RESTRICTION_MAP* GeoMap, double Long, double Lat);

double	GeoDist(double Long1, double Lat1, double Long2, double Lat2);
double	GeoPointDist(RESTRICTION_POINT* Point, double Long, double Lat);
//int		ValidResPoint(RESTRICTION_MAP* Map, double Long, double Lat);
int		ValidNodeResPoint(NODE_RES_MAP *NodeResMap, double Long, double Lat);

NODE_RES_MAP*	CrateNodeResMap(RESTRICTION_MAP *ResMap);
void			FreeNodeResMap(NODE_RES_MAP* NodeResMap);

void	CheckResMaps(RESTRICTION_MAP** Maps, int NoMaps);
void	MapResMaps(OPTIONS *Opt, TREES *Trees, RESTRICTION_MAP** Maps, int NoMaps);

void	CreateGeoGrid(RESTRICTION_MAP* GeoMap, size_t NoLong, size_t NoLat);
//void	BuildBinResMaps( int argc, char **argv);

#endif