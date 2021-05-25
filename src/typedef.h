#if !defined (TYPEDEFS)
#define TYPEDEFS

#pragma warning(disable : 4996)
#include <stdio.h>
#include "matrix.h"
/*
#define	JNIRUN
*/
#ifdef	 JNIRUN
	#include "jni.h"
/*	Build 
	gcc-4.2 -dynamiclib -lm -O3 -w -o libBayesTraits.jnilib *.c ./MathLib/*.c -static-libgcc

	trying
	gcc-4.2 -dynamiclib -lm -O3 -w -o libBayesTraits.jnilib *.c ./MathLib/*.c -nodefaultlibs -static-libgcc -lsystem.B -nostatfiles

	find out what dynamic libs are in play. 
	otool -L libBayesTraits.jnilib

	its normaly, will have to sort this out.  
	/usr/lib/libSystem.B.dylib
*/
#endif

/* Cost of Uniform prior */
#define PPUNICOST		2.302585093

#define PPJACOBIAN	1

/* Max uniform vlaue */
#define PPMAXSCALE	10
//#define	PPSCALEDEV	10
#define	PPSCALEDEV	100

/* Value of the PP alpha gamma */
#define	PPALPHA		1.1
#define PPBETA		1

#define PPALPHASCLAE 0.1

/* Use uniform or gamma value priors*/
#define PPUNIFORM

#define MINRATE 1.0e-16
#define MAXRATE	1000
/*	#define MAXRATE	100 */
#define MINBL	0.0000001

#define MAX_NUMBER 50	/* maximum number of states */
#define MAX_NUM_PARAMS (((MAX_NUMBER) * (MAX_NUMBER)) - (MAX_NUMBER) + 1)
#define MAX_N (MAX_NUMBER * (MAX_NUMBER)) - 1

#define LOGFILEEXT		"log.txt"
#define UNKNOWNSTATE	'-'
#define SUMMARYFILEEXT	"sum.txt"
#define ESTDATAPOINT	"?"

#define ERRLH -999999

#define ZERORATENO		-1

//extern double LhPraxis(LhPraxisdouble *);

#define GAMMAMAX	10000
#define GAMMAMIN	0.05

#define	LOGFILEBUFFERSIZE	65536

#define MIN_DELTA	1E-07
#define MAX_DELTA	3

#define MIN_LAMBDA	1E-07
#define MAX_LAMBDA	1

#define MIN_KAPPA	1E-07
#define MAX_KAPPA	3

typedef enum
{
	CRUN,
	CRES,
	CUNRES,
	CRESALL,
	CUNRESALL,
	CPRIOR,
	CITTERS,
	CSAMPLE,
	CPRIORCAT,
	CMLTRIES,
	CINFO,
	CPRIORALL,
	CHELP,
	CNODE,
	CMRCA,
	CDELNODE,
	CADDTAXA,
	CDELTAXA,
	CEVENROOT,
	CLOGFILE,
	CMODEL,
	CRATEDEV,
	CPRESET,
	CSUMMARY,
	CBURNIN,
	CPIS,
	CKAPPA,
	CDELTA,
	CLAMBDA,
	CEXTTAXA,
	CTAXAINFO,
	CSAVETREES,
	CTESTCORREL,
	CSURFACE,
	CCOVARION,
	CREVJUMP,
	CEXIT,
	CFOSSIL,
	CNODEDATA,
	CALPHAZERO,
	CHYPERPRIOR,
	CHPRJ,
	CHPALL,
	CNODEBLDATA,
	CGAMMA,
	CCI,
	CDEPSITE,
	CHEADERS,
	CMODELFILE,
/*	CPREVAR, */
	CVARDATA,
	CRMODEL,
	CDATADEV,
	CCOMMENT,
	CNOSPERSITE,
	CSCHEDULE,
	CSTREEMOVE,
	CSETSEED,
	CMAKEUM,
	CPHYLOPLASTY,
	CUNKNOWN,
} COMMANDS;

static char    *COMMANDSTRINGS[] =
{
	"run",			"ru",
	"restrict",		"res",
	"unrestrict",	"unres",
	"restrictall",	"resall",
	"unrestrictall","unresall",
	"prior",		"pr",
	"iterations",	"it",
	"sample",		"sa",
	"priorcats",	"pcat",
	"mltries",		"mlt",
	"info",			"in",
	"priorall",		"pa",
	"help",			"he",
	"addnode",		"addn",
	"addmrca",		"mrca",
	"delnode",		"deln",
	"addtaxa",		"addtaxa",
	"deltaxa",		"deltaxa",
	"evenroot",		"er",
	"logfile",		"lf",
	"hiddenstate",	"hs",
	"ratedev",		"rd",
	"preset",		"ps",
	"summary",		"sum",
	"burnin",		"bi",
	"pis",			"pi",
	"kappa",		"ka",
	"delta",		"dl",
	"lambda",		"la",
	"excludetaxa",	"et",
	"taxainfo",		"ti",
	"savetrees",	"st",
	"testcorrel",	"tc",
	"surface",		"su",
	"covarion",		"cv",
	"revjump",		"rj",
	"exit",			"quit",
	"fossil",		"fo",
	"nodedata",		"nd",
	"alphazero",	"az",
	"hyperprior",	"hp",
	"revjumphp",	"rjhp",
	"hyperpriorall","hpall",
	"nodebldata",   "nbd",
	"gamma",		"ga",
	"confint",		"cf",
	"depsite",		"ds",
	"headers",		"hd",
	"modelfile",	"mf",
/*	"prevar",		"pv", */
	"vardata",		"vd",
	"rmodel",		"rm",
	"datadev",		"dd",
	"#",			"//",
	"fitnospersite","nps",
	"schedule",		"sch",
	"solotreemove",	"stree",
	"setseed",		"ss",
	"makeum",		"mum",
	"phyloplasty",	"pp", 
	""
};

static char		*MODELNAMES[] =
{
	"Multistates",
	"Discete Independant",
	"Discete Dependent",
	"Continuous Random Walk",
	"Continuous Directional",
	"Continuous Regression",
	"Independent Contrasts"
};

static char    *DEPPRAMS[] =
{
	"q12",
	"q13",
	"q21",
	"q24",
	"q31",
	"q34",
	"q42",
	"q43",
	""
};

static char    *INDEPPRAMS[] =
{
	"alpha1",
	"beta1",
	"alpha2",
	"beta2",
	""
};

static char    *DISTNAMES[] =
{
	"beta",
	"gamma",
	"uniform",
	"chi",
	"exp",
	""
};


static int	DISTPRAMS[] =
{
	2,
	2,
	2,
	2,
	1
};

typedef enum
{
	DISCRETE,
	CONTINUOUS
} DATATYPE;

typedef enum
{
	BETA=0,
	GAMMA=1,
	UNIFORM=2,
	CHI=3,
	EXP=4
} PRIORDIST;

typedef enum
{
	MULTISTATE,
	DESCINDEP,
	DESCDEP,
	CONTINUOUSRR,
	CONTINUOUSDIR,
	CONTINUOUSREG,
	CONTRASTM
} MODEL;

typedef enum
{
	ANALML,
	ANALMCMC,
} ANALSIS;

typedef enum
{
	RESNONE,
	RESCONST,
	RESRATE
} RESTYPES;

typedef enum
{
	PIUNI,
	PIEST,
	PIEMP,
	PINONE,
} PITYPES;

typedef enum
{
	PPNODE,
	PPBRANCH
} PLASTYTYPE;

typedef struct
{
	int	K;
	unsigned long *States;
} RANDSTATES;

typedef struct
{
	int		No;
	char*	Name;
	char**	DesDataChar;
	double*	ConData;
	int		Exclude;
	double	Dependant;

	int		EstData;
	char	*EstDataP;
	int		EstDepData;
	char	*RealData;
} TAXA;

typedef struct
{
	double	Data;
	double	Cont;
	double	Var;
	double	Err;
} CONTRAST;

struct INODE
{
	int		Tip;
	int		TipID;
	int		ID;

	double	Length;

	struct INODE	*Left;
	struct INODE	*Right;
	struct INODE	*Ans;

	char		Visited;

	double		**Partial;
	double		**GammaPartial;
	char		*Tag;

	TAXA		*Taxa;

	int			*Part;
	int			PSize;

	int			FossilState;

	CONTRAST*	Contrast;
};

typedef enum
{
	NODEREC,
	MRCA,
	FOSSIL,
} NODETYPE;

typedef struct INODE*	NODE;

struct RNODE
{
	NODETYPE	NodeType;
	char*		Name;

	int			NoOfTaxa;
	int			PresInTrees;
	int			FossilState;

	int			*TaxaID;
	TAXA		**Taxa;

	int			Hits;
	NODE*		TreeNodes;

	char**		ConData;
	
	struct		RNODE *Next;
};

typedef struct RNODE*	RECNODE;

typedef struct
{
	int		*NoOfBLVect;
	double	**TToTPath;
}  TAXADIST;

typedef struct
{
	MATRIX		*TrueV;
	MATRIX		*V;
	MATRIX		*InvV;
	MATRIX		*Sigma;
	MATRIX		*InvSigma;
	MATRIX		*KProd;
	MATRIX		*InvKProd;

	MATRIX		*TVT;
	MATRIX		*TVTTemp;

	MATRIX		*InvXVX;

	double		*Alpha;
	double		*Beta;
	double		*Z;
	double		*ZA;

	double		*ZATemp;

	double		*TVect1;
	double		*TVect2;
	double		*TVect3;
	double		*SVect;


	TAXADIST	*TaxaDist;

	double		LogDetOfV;
	double		LogDetOfSigma;

	double		*DepVect;

} CONVAR;

typedef	struct
{
	NODE		NodeList;
	NODE		Root;

	CONVAR*		ConVars;
} TREE;

typedef struct
{
	PRIORDIST	Dist;
	int			RateNo;
	double		*DistVals;
	double		*HP;
	double		OffSet;
	int			UseHP;
	int			NoOfCats;
	char		*RateName;
} PRIORS;

typedef	struct
{
	/* Matrix Inver info */
	MATRIX		*vec;
	MATRIX		*inv_vec;
	double		*val;
	MATRIX		*Q;
	MATRIX		*A;

	MATRIX		*TempA;
	double		*TempVect1;
	double		*TempVect2;
	double		*TempVect3;
	double		*TempVect4;

} INVINFO;

typedef struct
{
	/* Used in FindInvV */
	double	*T1;
	int		*T2;
	MATRIX*	TMat;

	/* FindMLRagVals */
	MATRIX	*X;
	MATRIX	*TranX;
	MATRIX	*NX;
	double	*Y;

	/* FindRegVar */
	MATRIX	*XT;
	MATRIX	*RVX;
	MATRIX	*TempV1;
	MATRIX	*TempV2;

} TEMPCONVAR;

typedef struct
{
	int			NoOfTrees;
	int			NoOfTaxa;
	int			NoOfNodes;
	int			NoOfSites;
	int			NoOfStates;

	INVINFO*	InvInfo;

	MATRIX		*PLeft;
	MATRIX		*PRight;

	TAXA		*Taxa;
	TREE		*Tree;

	char		*SymbolList;

	int			ValidCData;
	int			ValidDData;

	char		**RemovedTaxa;
	int			NoOfRemovedTaxa;

	int			UseCovarion;

	int			NOSPerSite;
	int			MaxNOS;
	int			*NOSList;
	char		**SiteSymbols;
	TEMPCONVAR	*TempConVars;

	int			JStop;
} TREES;

typedef struct
{
	double	**Data;
	char	**TaxaNames;
	int		Site;
	int		NoPoints;
} VARDATA;

typedef struct
{
	NODE		Node;
	double		Scale;
	PLASTYTYPE	Type;
	int			NodeID;
} PLASTYNODE;

typedef struct
{
	int			NoNodes;
	PLASTYNODE **NodeList;
	
	int			NoTrees;
	double		**TrueBL;
	double		*ScaleBL;

	int			NoValidNode;
	NODE		*ValidNode;

	NODE		*TempList;
	int			NoTempList;

	double		Alpha;
} PLASTY;

typedef struct
{
	double*	Alpha;
	double* Sigma;

	double*	EstAlpha;
	double*	EstSigma;

	double	AlphaErr;
} CONTRASTR;

typedef struct
{
	MODEL		Model;
	ANALSIS		Analsis;

	int			NoOfRates;
	char		**RateName;

	RESTYPES	*ResTypes;
	int			*ResNo;
	double		*ResConst;
	
	double		RateDev;
	double		*RateDevList;

	
	double		EstDataDev;
	double		PPScaleDev;
//	PPSCALEDEV

	PRIORS		**Priors;

	int			PriorCats;

	int			Itters;
	int			Sample;
	int			BurnIn;
	int			MLTries;

	int			NoOfRecNodes;

	RECNODE		RecNode;
	RECNODE		*RecNodeList;

	TREES		*Trees;

	char		*LogFN;
	char		*TreeFN;
	char		*DataFN;
	FILE		*LogFile;
	FILE		*LogFileRead;
	FILE		*PPTree;
	FILE		*PPLog;
	char		*LogFileBuffer;
	char		**PassedOut;

	int			Summary;

	PITYPES		PiTypes;

	DATATYPE	DataType;

	int			UseKappa;
	int			UseDelta;
	int			UseLambda;
	int			UseGamma;

	int			EstKappa;
	int			EstDelta;
	int			EstLambda;
	int			EstGamma;

	double		FixKappa;
	double		FixDelta;
	double		FixLambda;
	double		FixGamma;

	int			InvertV;

	PRIORS		*PriorKappa;
	PRIORS		*PriorDelta;
	PRIORS		*PriorLambda;
	PRIORS		*PriorGamma;

	int			GammaCats;

	int			TestCorrel;

	int			UseCovarion;

	int			SetPraxis;
	int			KTM;
	int			LMaxFun;

	int			UseRJMCMC;

	PRIORS		*RJPrior;

	int			NodeData;
	int			NodeBLData;
	int			AlphaZero;

	double		HPDev;

	int			FindCF;
	char*		CFRate;

	int			DependantSite;
	int			Headers;

	char*		ModelFile;
	int			UseModelFile;

	int			UseVarData;
	char		*VarDataFile;
	VARDATA*	VarData;

	int			UseRModel;
	double		RModelP;

	int			NoEstDataSite;
	int			*EstDataSites;
	int			NoEstChanges;

	int			AnalyticalP;

	int			NOSPerSite;

	int			UseSchedule;
	char		*ScheduleFile;

	int			SoloTreeMove;
	long		Seed;
	int			MakeUM;

	int			UsePhyloPlasty;

} OPTIONS;


typedef struct
{
	int		NoOfRates;
	int		NoOfFullRates;
	int		NoOfPriors;

	double	*Rates;
	double	*FullRates;
	double	*Root;
	double	*Pis;

	int		TreeNo;
	double	Lh;
	double	LhPrior;
	double	LnHastings;
	double	LogJacobion;

	PRIORS	**Prios;

	double	*Means;

	double	*Beta;

	double	Delta;
	double	Lambda;
	double	Kappa;

	double	OnToOff;
	double	OffToOn;
	double	CoVarPis[2];

	int		*MappingVect;


	double	*GammaMults;
	int		GammaCats;
	double	Gamma;
	double	LastGamma;
	PRIORS	*GammaPrior;

	double	Numer;
	double	Donom;

	double	HMeanSum;
	int		HMeanCount;

	int		UseEstData;
	int		*EstDescData;
	double	*EstData;
	int		NoEstData;

	double	**FixedModels;
	int		NoOfModels;
	int		ModelNo;

	int		VarDataSite;

	RANDSTATES		*RandStates;

	PLASTY			*Plasty;
	CONTRASTR		*Contrast;
} RATES;

typedef struct
{
	int		N;
	double	Sum;
	double	SumSqrs;
} SUMMARYNO;

typedef struct
{
	SUMMARYNO	Lh;
	SUMMARYNO	*Rates;
	SUMMARYNO	*Root;
} SUMMARY;

#define NOOFOPERATORS	14

static char    *SHEDOP[] =
{
	"Rate",
	"CV",
	"Kappa",
	"Delta",
	"Labda",
	"Jump",
	"Prior Change",
	"Est Data",
	"Var Data",
	"Solo Tree Move",
	"PP Add / Remove",
	"PP Move", 
	"PP Change Scale",
	"PP Hyper Prior"
};

typedef enum
{
	SRATES=0,
	SCV=1,
	SKAPPA=2,
	SDELTA=3,
	SLABDA=4,
	SJUMP=5,
	SPPROR=6,
	SESTDATA=7,
	SVARDATA=8,
	SSOLOTREEMOVE=9,
	SPPADDREMOVE=10,
	SPPMOVE=11,
	SPPCHANGESCALE=12,
	SPPHYPERPRIOR=13,
} OPERATORS;

typedef struct
{
	int		Op;
	int		NoOfOpts;

	double	OptFreq[NOOFOPERATORS];
	int		Tryed[NOOFOPERATORS];
	int		Accepted[NOOFOPERATORS];
} SCHEDULE;

typedef struct
{
	int	NoOfGroups;
	int	*GroupSize;
	int	*GroupID;
	int **GroupPos;

	int	NoInZero;
	int	*ZeroPos;
} MAPINFO;

#endif
