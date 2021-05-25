#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

//#define DllExport __declspec(dllexport)
//#define DllImport __declspec(dllimport)


#include "JBayesTraits_JAnalysis.h"

#include "typedef.h"
#include "trees.h"
#include "data.h"
#include "options.h"
#include "rates.h"
#include "likelihood.h"
#include "rand.h"
#include "priors.h"
#include "mcmc.h"
#include "praxis.h"
#include "ml.h"
#include "genlib.h"
#include "continuous.h"
#include "initialise.h"
#include "rand.h"

#include "./MathLib/dcalc.h"
#include "./MathLib/mconf.h"

#include <jni.h>

#define	MAXJANALYSIS 128
			 

int		SetOnce = FALSE;
/*
char	*DeadAnalysis[MAXJANALYSIS];

void	InitBTJNI()
{
	int i;
	if(SetOnce == TRUE)
		return;
	SetOnce = TRUE;

	for(i=0;i<MAXJANALYSIS;i++)
		DeadAnalysis[i] = NULL;
}

*/

MODEL	GetModelFromNo(int ModelNo)
{
	switch(ModelNo)
	{
		case 0:	return 	MULTISTATE;
		case 1: return	DESCINDEP;
		case 2: return 	DESCDEP;
		case 3: return	CONTINUOUSRR;
		case 4: return	CONTINUOUSDIR;
		case 5: return 	CONTINUOUSREG;
	}

	exit(0);
}

ANALSIS	AnalsisFromNo(int AnalisisNo)
{
	if(AnalisisNo == 0)
		return ANALML;

	if(AnalisisNo == 1)
		return ANALMCMC;

	exit(0);
}


void	SetProgress(JNIEnv *Env, jobject Obj, int Progress)
{
	jclass		classID;
	jmethodID	mID;

	classID = (*Env)->GetObjectClass(Env, Obj);
	mID		= (*Env)->GetMethodID(Env, classID, "SetProgress", "(I)V");
	(*Env)->CallVoidMethod(Env, Obj, mID, Progress);
}

void	CheckStop(JNIEnv *Env, jobject Obj, TREES *Trees)
{
	jclass		classID;
	jfieldID	fID;
	jboolean	Stop;

	classID = (*Env)->GetObjectClass(Env, Obj);
	fID		= (*Env)->GetFieldID(Env, classID, "Stop", "Z");
	Stop	= (*Env)->GetBooleanField(Env, Obj, fID);

	if(Stop == JNI_TRUE)
		Trees->JStop = TRUE;
}

int mainJNI(JNIEnv *Env, jobject Obj, int Size, char** RunP)
{
	TREES*		Trees=NULL;
	OPTIONS*	Opt=NULL; 
	char*		TreeFN;
	char*		DataFN;
	MODEL		Model;
	ANALSIS		Analysis;

	
	TreeFN	= RunP[0];
	DataFN	= RunP[1];
	Model	= GetModelFromNo(atoi(RunP[2]));
	Analysis= AnalsisFromNo(atoi(RunP[3]));
	

	if(SetOnce == FALSE)
	{
		SetSeed();
		SetOnce = TRUE;
	}

	Trees  = LoadTrees(TreeFN);
	LoadData(DataFN, Trees);
	
	Opt = CreatOptions(Model, Analysis, Trees->NoOfStates, TreeFN, DataFN, Trees->SymbolList, Trees);
	GetOptionsArry(Opt, Size-4, &RunP[4]);

	PreProcess(Opt, Trees);

	if(Opt->Analsis == ANALMCMC)
		MCMC(Opt, Trees, Env, Obj);

	if(Opt->Analsis == ANALML)
		FindML(Opt, Trees, Env, Obj);

	FreeTrees(Trees, Opt);
	FreeOptions(Opt);

	return 0;
}


/*

#ifdef _MANAGED
#pragma managed(push, off)
#endif

*/

JNIEXPORT void JNICALL Java_JBayesTraits_JAnalysis_RunBayesTraits
/* JNIEXPORT void JNICALL Java_JBayesTraits_JModelRun_RunBayesTraits */
  (JNIEnv *Env, jobject Obj, jobjectArray Arr)
{
	int		Index;
	int		ArrLen;
	char**	Commands;
	jstring strElem;

	ArrLen = (*Env)->GetArrayLength(Env, Arr);

	Commands = (char**)malloc(sizeof(char*) * ArrLen);
	if(Commands == NULL)
		MallocErr();
	
	for(Index=0;Index<ArrLen;Index++) 
	{
		strElem = (jstring)(*Env)->GetObjectArrayElement(Env, Arr, Index);
		if(strElem != NULL) 
		{
			const char *strTemp = (*Env)->GetStringUTFChars(Env, strElem, NULL);
			Commands[Index] = StrMake(strTemp);
			(*Env)->ReleaseStringChars(Env, strElem, (jchar*)strTemp);
			(*Env)->DeleteLocalRef(Env, strElem);
		}
	}

//	for(Index=0;;Index++)Index=0;
	mainJNI(Env, Obj, ArrLen, Commands);

	for(Index=0;Index<ArrLen;Index++)
		free(Commands[Index]);
	free(Commands);
}

/*
BOOL APIENTRY DllMain( HMODULE hModule,
                       DWORD  ul_reason_for_call,
                       LPVOID lpReserved
					 )
{
    return TRUE;
}

#ifdef _MANAGED
#pragma managed(pop)
#endif
*/
