#ifndef		LOCAL_TRANSFORMS_H
#define		LOCAL_TRANSFORMS_H

#include "typedef.h"

void				ApplyLocalTransforms(RATES *Rates, TREES *Trees, OPTIONS *Opt, int Norm);

LOCAL_TRANSFORM*	CreateLocalTransforms(char *Name, TAG **TagList, int NoTags, TRANSFORM_TYPE Type, int Est, double Scale);
void				FreeLocalTransforms(LOCAL_TRANSFORM* LTrans);

void				CopyLocalTransforms(LOCAL_TRANSFORM* ATrans, LOCAL_TRANSFORM* BTrans);
LOCAL_TRANSFORM*	CloneLocalTransform(LOCAL_TRANSFORM* LTrans);

void				PrintLocalTransform(FILE *Str, LOCAL_TRANSFORM* Trans);
void				PrintLocalTransforms(FILE *Str, LOCAL_TRANSFORM** List, int NoTrans);

int					EstLocalTransforms(LOCAL_TRANSFORM** List, int NoTrans);
int					NoEstLocalTransform(LOCAL_TRANSFORM** List, int NoTrans);

void				ChangeLocalTransform(OPTIONS *Opt, TREES *Trees, RATES *Rates, SCHEDULE *Shed);

int					GetNoTransformType(TRANSFORM_TYPE TType, RATES *Rates);

double				CaclLocalTransformsPrior(RATES *Rates);

double				ChangeLocalScale(RANDSTATES	*RS, double Scale, double Dev);

#endif