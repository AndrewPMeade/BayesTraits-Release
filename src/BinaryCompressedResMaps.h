#ifndef BINCOMPRESMAP_H
#define BINCOMPRESMAP_H


typedef struct
{
	char *FileName;
	double AgeMin, AgeMax;
} BIN_MAP_INFO;

typedef struct
{
	int NoMaps;
	BIN_MAP_INFO** BinMaps;

	char* FileName;

	int NoLatBoxes;
	int NoLongBoxes;

} MAP_GROUP;

void	BuildBinResMaps(char *MapFileName);
RESTRICTION_MAP**	LoadBinCompResMaps(char *FName, int *Size);


#endif