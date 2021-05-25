#if !defined (BAYESTRAITSJNIH)
#define BAYESTRAITSJNIH

#include "typedef.h"

int		mainJNI(JNIEnv *Env, jobject Obj, int Size, char** RunP);
void	SetProgress(JNIEnv *Env, jobject Obj, int Progress);
void	CheckStop(JNIEnv *Env, jobject Obj, TREES *Trees);

void	ProcessHeaders(JNIEnv *Env, jobject Obj, OPTIONS *Opt);
void	AddResults(JNIEnv *Env, jobject Obj, OPTIONS *Opt);
#endif

