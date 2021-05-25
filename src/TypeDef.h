#ifndef TYPEDEFS_H
#define TYPEDEFS_H


//#pragma warning(disable : 4996)
#include <stdio.h>
#include <assert.h>
#include <gsl/gsl_randist.h>

#include "Matrix.h"
#include "RandLib.h"
#include "StableDist.h"
#include "AutoTune.h"

//#define JNIRUN
//#define OPENMP_THR
//#define BIG_LH

// #define BTOCL
//#define BTLAPACK

// use the intel MLK lib
//#define USE_MKL

#ifdef USE_MKL
	#define MKL_INT int
	#include <mkl.h>
#endif

#ifdef BTLAPACK 
#ifndef USE_MKL
	#include <cblas.h>
#endif
#endif


#ifdef BTOCL
	#include "btocl_runtime.h"
#endif

//#define	CLIK_P
//#define	PUBLIC_BUILD

// Information for the phomem cog est runs
// #define	PHONEIM_RUN

//#define QUAD_DOUBLE

// use nlopt libray for ML esitmates, much better than praxis if available 
#define NLOPT
#ifdef NLOPT
	#include <nlopt.h>
#else
	// maximum number of states 
	#define MAX_NUMBER 50	
	#define MAX_NUM_PARAMS (((MAX_NUMBER) * (MAX_NUMBER)) - (MAX_NUMBER) + 1)
	#define MAX_N (MAX_NUMBER * (MAX_NUMBER)) - 1
#endif

#ifdef BIG_LH
	#include <gmp.h>
	#include <mpfr.h>

	#define DEF_ROUND GMP_RNDZ
	#define	DEF_PRE	256
#endif

#ifdef QUAD_DOUBLE
	#include <quadmath.h>
	#define QDOUBLE	__float128
#endif

#ifdef	OPENMP_THR
	#include <omp.h>
#endif

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


/* Program Control */
#define	RATE_CHANGE_UNI
#define RATE_CHANGE_NORM

// #define NO_GSL

//#define LH_UNDER_FLOW 9.9296E-153
#define LH_UNDER_FLOW 1.3839E-87

// Change only one rate at a time
//#define RATE_CHANGE_ONE

// Used to test if tree transforms are normalised, so there is no cahnge in tree lenght. 
#define		NORMALISE_TREE_CON_SCALING FALSE

/* If defined Sigma for "Independent Contrast: Full" usning MCMC is restricted to a value. */
//#define RES_SIGMA	0.1
//#define RES_ALPHA	0


/* Use uniform or gamma value priors*/
//#define PPUNIFORM
#define PPGAMMA

/* Cost of Uniform prior */
#define PPUNICOST		2.302585093
#define PPJACOBIAN		1


/* Max uniform vlaue */
#define PPMAXSCALE	10
#define	PPSCALEDEV	100

/*	Modify the prior cost adding a VarRates scalar */
#define VAR_RATES_PRIOR_SCALE 1

/* Value of the VarRates alpha beta*/
#define	VARRATES_ALPHA		1.1
#define VARRATES_BETA		1

#define VARRATES_HP_ALPHA_SCLAE 0.1

// Only allow one type (Node, BL, Kappa ect) of variables rates operator per node 
//#define VARRATES_ONE_OP_PER_NODE

// use ML paramter for indpedent contrast MCMC / Var Rates
//#define CONTRAST_ML_PARAM



// Minimum number of taxa to transform a node (kappa, lamabed ect)
//#define MIN_TAXA_VR_TRANS	5

// Minimum number of taxa to Var Rate a ndoe
#define MIN_TAXA_VR_NODE	0

// def min and max rates, min rate is enforced with user overrides 
#define RATE_MIN 1.0e-8
#define RATE_MAX 100

/*	#define MAXRATE	100 */
// Minium barnch legnth
#define MIN_BL	0.0000001

// Output File Extensions 
#define OUTPUT_EXT_LOG			".Log.txt"
#define OUTPUT_EXT_SCHEDULE		".Schedule.txt"
#define OUTPUT_EXT_ANC			".AncStates.txt"
#define OUTPUT_EXT_DUMMY_CODE	".DummyCode.txt"
#define OUTPUT_EXT_SIM			".Sim.txt"
#define OUTPUT_EXT_VAR_RATES	".VarRates.txt"
#define OUTPUT_EXT_STONES		".Stones.txt"
#define OUTPUT_EXT_TREES		".Output.trees"


#define UNKNOWNSTATE	'-'
#define ESTDATAPOINT	"?"

#define ERRLH -999999

#define ZERO_RATE_NO		-1

/* Use to set V to idenity for testing */
//#define ID_MATRIX

//extern double LhPraxis(LhPraxisdouble *);

#define	LOGFILEBUFFERSIZE	65536

#define DISPLAY_INFO	printf("BayesTraits V3.0 (%s)\nMark Pagel and Andrew Meade\nwww.evolution.reading.ac.uk\n\n\n",__DATE__);fflush(stdout);

#define MIN_DELTA	1E-07
#define MAX_DELTA	3

#define MIN_LAMBDA	1E-07
#define MAX_LAMBDA	1

#define MIN_KAPPA	1E-07
#define MAX_KAPPA	3

#define MIN_OU		1E-07
#define MAX_OU		50

#define MIN_GAMMA	1E-07
#define MAX_GAMMA	100

#define MIN_LOCAL_RATE 1E-07
#define MAX_LOCAL_RATE 100

#define	NORM_MEAN_BL	0.1

#define	MAX_VR_KAPPA	25
#define MAX_VR_DELTA	25
#define MAX_VR_LAMBDA	1
#define MAX_VR_OU		100

// Value for a normal distribution for a fat tail model
#define FAT_TAIL_NORMAL_VAL 2.0

// Minimum and maximum acceptance rates, set for Auto tuning  
#define MIN_VALID_ACC 0.2
#define MAX_VALID_ACC 0.4

#define MIN_NO_TAXA_RJ_LOCAL_TRANS 10

// Number of local scalars 
#define	NO_RJ_LOCAL_SCALAR	6

static char    *RJ_LOCAL_SCALAR_NAMES[] =
{
	"kappa",
	"lambda",
	"delta", 
	"ou",
	"node", 
	"branch"
};

typedef enum 
{
	VR_KAPPA,
	VR_LAMBDA, 
	VR_DELTA,
	VR_OU,
	VR_NODE,
	VR_BL,
} TRANSFORM_TYPE;

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
	CMLTOL, 
	CMLEVAL,
	CMLALG,
	CINFO,
	CPRIORALL,
	CHELP,
	CNODE,
	CMRCA,
	CADDTAXA,
	CDELTAXA,
	CEVENROOT,
	CLOGFILE,
	CPRESET,
	CSUMMARY,
	CBURNIN,
	CPIS,
	CKAPPA,
	CDELTA,
	CLAMBDA,
	CEXTTAXA,
	CSAVEINITIALTREES,
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
	CRMODEL,
	CCOMMENT,
	CNOSPERSITE,
	CSETSEED,
	CMAKEUM,
	CEQUALTREES,
	CPRECISION,
	CCORES,
	CSYMMETRICAL,
	CMCMCMLSTART,
	CCAPRJRATES,
	CSAVEMODELS,
	CLOADMODELS,
	COU,
	CVARRATES,
	CSTONES,
	CADDERR,
	CSHEDULE, 
	CRJDUMMY,
	CSCALETREES,
	CRJLOCALTRANSFORM,
	CFATTAILNORMAL,
	CADDTAG,
	CLOCALTRANSFORM,
	CDISTDATA,
	CNOLH,
	CSAVETREES,
	CCSCHED, 
	CADDTIMESLICE,
	CADDPATTERN,
	CSETMINTAXATRANS,
	CSETMINMAXRATE, 
	CNORMQMAT, 
	CNOSLICESAMPLESTEPS,
	CPISANCSTATES,
	CRJZERO,
	CLANDSCAPE,
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
	"mltol",		"mlto",
	"mlmaxeval",	"mlme",
	"mlalg",		"mla",
	"info",			"in",
	"priorall",		"pa",
	"help",			"he",
	"addnode",		"addn",
	"addmrca",		"mrca",
	"addtaxa",		"addtaxa",
	"deltaxa",		"deltaxa",
	"evenroot",		"er",
	"logfile",		"lf",
	"preset",		"ps",
	"summary",		"sum",
	"burnin",		"bi",
	"pis",			"pi",
	"kappa",		"ka",
	"delta",		"dl",
	"lambda",		"la", 
	"excludetaxa",	"et",
	"saveinitialtrees", "sit",
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
	"rmodel",		"rm",
	"#",			"//",
	"fitnospersite","nps",
	"seed",			"se",
	"makeum",		"mum",
	"equaltrees",	"eqt",
	"precision",	"pre",
	"cores",		"cor",
	"symmetrical",	"sym", 
	"mcmcmlstart",	"mls",
	"CapRJRates",	"cap", 
	"SaveModels",	"sm",
	"LoadModels",	"lm", 
	"OU",			"ou",
	"VarRates",		"vr",
	"Stones",		"ss",
	"AddErr",		"er",
	"Schedule",		"sh", 
	"RJDummy",		"rjd",
	"ScaleTrees",	"sct", 
	"RJLocaltransform", "rjlt", 
	"FatTailNormal", "ftn",
	"AddTag",		"at",
	"LocalTransform",	"lt",
	"DistData",			"dd",
	"NoLh",			"nl",
	"SaveTrees",	"savetree",
	"CustomSchedule", "csched", 
	"AddTimeSlice", "ats",
	"AddPattern", "ap",
	"SetMinTransTaxaNo", "smttn",
	"SetMinMaxRate", "smmr",
	"NormaliseQMatrix", "nqm",
	"NoSliceSampleSteps", "nsss",
	"PisAncStates",		"pas",
	"RJZero",			"rjz",
	"Landscape",		"LS",
	""
};

static char		*MODELNAMES[] =
{
	"MultiState",
	"Discrete: Independent",
	"Discrete: Dependent",
	"Continuous: Random Walk",
	"Continuous: Directional",
	"Continuous: Regression",
	"Independent Contrasts",
	"Independent Contrasts: Correlation",
	"Independent Contrasts: Regression",
	"Discrete: Covarion", 
	"Discrete: Heterogeneous",
	"Fat Tail",
	"Geo"
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

static char    *DEPCVPRAMS[] =
{
	"alpha1",
	"beta1",
	"alpha2",
	"beta2",
	"q12",
	"q13",
	"q21",
	"q24",
	"q31",
	"q34",
	"q42",
	"q43",
	"qDI",
	"qID"
	""
};

static char    *DEPHETROPRAMS[] =
{
	"alpha1",
	"beta1",
	"alpha2",
	"beta2",
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

#define NO_PRIOR_DIST 7

static char    *DISTNAMES[] =
{
	"gamma",
	"uniform",
	"chi",
	"exp",
	"sgamma",
	"lognormal",
	"normal"
};

static int	DISTPRAMS[] =
{
	2,
	2,
	1,
	1,
	2,
	2,
	2
};

typedef enum
{
	PDIST_GAMMA,
	PDIST_UNIFORM,
	PDIST_CHI,
	PDIST_EXP,
	PDIST_SGAMMA,
	PDIST_LOGNORMAL,
	PDIST_NORMAL
} PRIORDIST;

typedef enum
{
	DISCRETE,
	CONTINUOUS
} DATATYPE;

typedef enum
{
	M_MULTISTATE,
	M_DISC_INDEP,
	M_DISC_DEP,
	M_CONTINUOUS_RR,
	M_CONTINUOUS_DIR,
	M_CONTINUOUS_REG,
	M_CONTRAST,
	M_CONTRAST_CORREL,
	M_CONTRAST_REG,
	M_DISC_CV,
	M_DISC_HET,
	M_FATTAIL,
	M_GEO
} MODEL;

typedef enum
{
	MT_DISCRETE,
	MT_CONTINUOUS,
	MT_CONTRAST,
	MT_FATTAIL
} MODEL_TYPE;

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
	PI_UNI,
	PI_EMP,
	PI_NONE,
} PITYPES;

typedef enum
{
	 RJDUMMY_INTER,
	 RJDUMMY_INTER_SLOPE
} RJDUMMY_TYPE;

typedef enum
{
	NODEREC,
	MRCA,
	FOSSIL,
} NODETYPE;


typedef struct
{
	char *Name;

	int FixedTime;
	int FixedScale;

	double Time;
	double Scale;
	
} TIME_SLICE;

typedef struct
{
	int NoTimeSlices;
	TIME_SLICE	**TimeSlices;

} TIME_SLICES;

typedef struct
{
	// User supplied taxa number, may not be continues  
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
	TAXA *Taxa;

	int		*NoSites;
	double	**Data;

	int		Linked;
} DIST_DATA_TAXA;

typedef struct
{
	char *FName;

	int		NoSites;
	int		NoTaxa;
	DIST_DATA_TAXA	**DistDataTaxa;

} DIST_DATA;

typedef struct
{
	int NoTaxa;
	int NoSites;

	int **SiteMap;
} DIST_DATE_RATES;

typedef struct
{
	double	*Data;
	double	*Cont;
//	double	*Var;
//	double	*Err;

	double	Var;
	double	Err;


	double	*v;
} CONTRAST;

typedef struct
{
	CONTRAST	**Contrast;
	int			NoContrast;

	double		*GVar;
	double		*SumLogVar;
} CONDATA;

typedef struct
{
	double	*Ans;

	double	Data, Cont, Err, Var, v;
	
} FATTAILNODE;

typedef struct
{
	int NoTaxa;
	int *Taxa;

	int		Freq;
	double	Prob;
	size_t	PartID;
} PART;

struct INODE
{
	int		Tip;
	int		TipID;
	int		ID;

	double		Length;
	double		DistToRoot;
	double		UserLength;
	double		ScaleFactor;

	int			VPosX, VPosY;

	struct INODE	*Ans;

	struct	INODE	**NodeList;
	int				NoNodes;
	char			Visited;

	double		**Partial;
	double		**GammaPartial;
	char		*Tag;
	int			NoUnderFlow;
	
#ifdef BIG_LH
	mpfr_t		**BigPartial;
	mpfr_t		t1, t2, t3;
#endif

#ifdef QUAD_DOUBLE
	QDOUBLE		**BigPartial;
#endif

	TAXA		*Taxa;
	
	PART		*Part;
	
	int			*FossilMask;

	int			PatternNo;

	CONDATA		*ConData;
	FATTAILNODE	*FatTailNode;

	double		LandscapeBeta;
};

typedef struct INODE*	NODE;


typedef struct
{
	char	*Name;

	int		NoTaxa;
	char	**Taxa;

	NODE	*NodeList;

	PART	*Part;
} TAG;

typedef struct
{
	NODETYPE	NodeType;
	char		*Name;

	int			*FossilStates;
	int			NoFossilStates;
		
	int			Hits;
	
	char		**ConData;

	TAG			*Tag;

} RECNODE;


typedef struct 
{
	double	*AnsVect;

} FATTAILTREE;

typedef struct
{
	MATRIX		*TrueV;
	MATRIX		*V;
	MATRIX		*InvV;
	MATRIX		*Sigma;
	MATRIX		*InvSigma;

#ifdef BTOCL   // continuous: V matrix and Z,ZA for Kronecker product
	cl_mem      buffer_invV;
	cl_mem		buffer_invSigma;
	cl_mem		buffer_ZA;
	cl_mem		buffer_ZATemp;	
#endif

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

	double		LogDetOfV;
	double		LogDetOfSigma;

	double		*DepVect;

	double		*MultiVarNormState;
	double		*MultiVarNormTemp;
	
} CONVAR;

typedef	struct
{
	int				NoNodes;
	NODE			*NodeList;
	NODE			Root;

	NODE			**FNodes;
	int				*NoFNodes;
	int				NoFGroups;

	int				NoPNodes;
	NODE			*PNodes;

	int				NoContrast;

	CONVAR*			ConVars;
	
	FATTAILTREE*	FatTailTree;

	double			AveBL;
	
// information needed to traverse the tree and accumulate partial results
#ifdef BTOCL
    int* groups;
	int* groupsIdx;  // start/end
	int* children;
	int* childrenIdx;
	int height;
	int max_nchildren;
	int* parentInfo;
	int* isTip;
#endif
	
} TREE;

typedef struct
{
	char		*Name;
	PRIORDIST	Dist;
	double		*DistVals;
	
	double		*HP;
	int			UseHP;
	
	int			Discretised;
	double		Width;

} PRIOR;

typedef	struct
{
	/* Matrix Inver info */
	MATRIX		*vec;
	MATRIX		*inv_vec;
	double		*val;
	MATRIX		*Q;
	MATRIX		*A;
#ifdef BTOCL    // discrete: Set PMatrix related
	// may want to separate InvInfo and PMatrix related 
	double*     vect_t;
	int*        vect_id;
	double*     vect_test;
	cl_mem      buffer_vec;  //  NOS*NOS  read-only
	cl_mem		buffer_inv_vec; // NOS*NOS  read-only
	cl_mem		buffer_val;  // NOS   read-only
	cl_mem		buffer_t;    // NoNodes  read-only
	cl_mem      buffer_id;    // NoNodes read only
	cl_mem		buffer_temp; // MaxNodes*NOS*NOS read-write
#endif

	MATRIX		*TempA;
	double		*TempVect1;
	double		*TempVect2;
	double		*TempVect3;
	double		*TempVect4;

	int			NoThreads;
	MATRIX		**As;
	double		**Ets;

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
	int			NoTrees;
	int			NoTaxa;
	int			NoSites;
	int			NoUserSites;
	int			NoStates;

	INVINFO**	InvInfo;

	double		*PMem;
	MATRIX		**PList;
	int			MaxNodes;

	TAXA		**Taxa;
	TREE		**Tree;

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

	double		*PMean;
	double		*PSD;

#ifdef BTOCL
	// discrete: SetPMatrix related 
	cl_mem 		buffer_pmatrix; // MaxNodes*NOS*NOS  write-only NO read/write!
	cl_mem      buffer_exp_eigen;  // MaxNodes*NOS Read/Write
	cl_mem      buffer_error;
	double* check_pmatrix;
	// discrete: compute partialLH related
	cl_mem      buffer_partialLh;
	cl_mem		buffer_groups;
	cl_mem		buffer_groupsIdx;
	cl_mem		buffer_children;
	cl_mem		buffer_childrenIdx;
	cl_mem		buffer_plhFactor; // MaxNodes*NoSites*NOS
	cl_mem      buffer_debug_plhFactor;
	int*        perror;
	double* previous_plh;
	double* plhFactor;
	double* temp_plh;
	double* debug_plhFactor;
	// new version
	cl_mem      buffer_parentInfo;
	cl_mem      buffer_isTip;
	int max_nchildren;
#endif
	
	PART	**PartList;
	size_t	NoParts;
} TREES;

typedef struct
{
//	NODE			Node;
	double			Scale;
	TRANSFORM_TYPE	Type;
	long long		NodeID;

	int				Fixed;
	int				UserSupplied;

	PART			*Part;
	NODE			*NodeList;

} VAR_RATES_NODE;

typedef struct
{
	int				NoNodes;
	VAR_RATES_NODE **NodeList;

	NODE		*TempList;
	int			NoTempList;

	double		Alpha;
} VARRATES;

typedef struct
{
	MATRIX *Uy, *Ux, *TUx, *InvUx;
	MATRIX *Prod1, *Prod2, *Prod3;
	double	*TempDVect;
	int		*TempIVect;
	
} REG_BETA_SPACE;

typedef struct
{
	double*	Alpha;
	double* Sigma;

	MATRIX_INVERT	*SigmaInvInfo;
	MATRIX			*SigmaMat;
	double			*SigmaInvVec;
	
	double	RegSigma;
	double	RegAlpha;
	double	*RegBeta;

	double	GlobalVar;

	int		NoSites;

	REG_BETA_SPACE	*RegSapce;
} CONTRASTR;

typedef struct
{
	double	*Power;
	double	*MLh;
	double	LastLh;
	double	Diff;
	double	Sum;
	double	Scalar, Length;

	int		NoStones;
	double	Alpha, Beta;
	
	int			ItPerStone;
	long long	ItStart;

	int		SampleFreq;
	int		Started;
	int		N;
} STONES;

typedef struct
{
	char			*Name;
	TRANSFORM_TYPE	Type;
	double			Scale;
	int				Est;
	TAG				**TagList;
	int				NoTags;
} LOCAL_TRANSFORM;


typedef struct
{
	long long	Iteration;
	double		*Frequencies;
	int			Default;
} CUSTOM_SCHEDULE;

typedef struct
{
	char	*Name;
	TAG		**TagList;
	int		NoTags;
} PATTERN;

typedef struct
{
	double		Beta;
	PART		*Part;
	NODE		*NodeList;
	long long	NodeID;
} LANDSCAPE_NODE;

typedef struct
{
	int				NoNodes;
	LANDSCAPE_NODE	**NodeList;
} LANDSCAPE;

typedef struct
{
	MODEL		Model;
	ANALSIS		Analsis;
	MODEL_TYPE	ModelType;

	int			NoOfRates;
	char		**RateName;

	int			DefNoRates;
	char		**DefRateNames;

	RESTYPES	*ResTypes;
	int			*ResNo;
	double		*ResConst;
	
	
	int			PriorCats;

	long long	Itters;
	int			Sample;
	long long	BurnIn;
	
	int			MLTries;
	int			MLMaxEVals;
	double		MLTol;
	char		*MLAlg;

	int			NoOfRecNodes;
	RECNODE		**RecNodeList;

	TREES		*Trees;

//	char		*LogFN;
	char		*BaseOutputFN;
	char		*TreeFN;
	char		*DataFN;
	FILE		*LogFile;
	FILE		*LogFileRead;
	FILE		*OutTrees;
	FILE		*VarRatesLog;
	FILE		*LogFatTail;
	char		*LogFileBuffer;
	char		**PassedOut;

	int			Summary;

	PITYPES		PiTypes;

	DATATYPE	DataType;

	int			UseKappa;
	int			UseDelta;
	int			UseLambda;
	int			UseGamma;
	int			UseOU;

	int			EstKappa;
	int			EstDelta;
	int			EstLambda;
	int			EstGamma;
	int			EstOU;

	double		FixKappa;
	double		FixDelta;
	double		FixLambda;
	double		FixGamma;
	double		FixOU;

	int			InvertV;

	PRIOR		**AllPriors;
	int			NoAllPriors;

	char		*SaveInitialTrees;

	int			GammaCats;

	int			TestCorrel;

	int			UseCovarion;

	int			SetPraxis;
	int			KTM;
	int			LMaxFun;

	int			UseRJMCMC;
	int			CapRJRatesNo;
	int			MCMCMLStart;
	
	int			NodeData;
	int			NodeBLData;
	int			AlphaZero;

	double		HPDev;

	int			FindCF;
	char*		CFRate;

	
	int			UseRModel;
	double		RModelP;

	int			NoEstDataSite;
	int			*EstDataSites;
	int			NoEstChanges;

	int			AnalyticalP;

	int			NOSPerSite;

	int			UseSchedule;

	long		Seed;
	int			MakeUM;

//	int			UseVarRates;

	int			UseEqualTrees;
	int			ETreeBI;

	int			Precision;
	int			Cores;

	int			SaveModels;
	char		*SaveModelsFN;

	int			LoadModels;
	char		*LoadModelsFN;

	STONES		*Stones;

	int			RJDummy;
	FILE		*RJDummyLog;

	double		RJDummyBetaDev;

	double		ScaleTrees;

	int			UseRJLocalScalar[NO_RJ_LOCAL_SCALAR];
	int			FatTailNormal;

	int			EstData;

	int			NoTags;
	TAG			**TagList;

	int				NoLocalTransforms;
	LOCAL_TRANSFORM	**LocalTransforms;

	DIST_DATA	*DistData;
	int			UseDistData;

	int			NoLh;

	int			SaveTrees;


	int				NoCShed;
	CUSTOM_SCHEDULE	**CShedList;

	TIME_SLICES		*TimeSlices;

	PATTERN		**PatternList;
	int			NoPatterns;

	int			MinTransTaxaNo;
	
	double		RateMin, RateMax;

	int			NormQMat;

	int			NoSliceSampleSteps;

	int			UsePisInAncStates;
		
	int			RJZero;

	int			UseLandscape;

} OPTIONS;

typedef struct
{
	int			NoModels;
	INVINFO**	ModelInv;
	
	int			MListSize;
	int			*MList;
} HETERO;

typedef struct
{
	int		NoModels;
	int		NoParam;
	char	*FName;
	double	**ModelP;
} MODELFILE;

typedef struct
{
	RJDUMMY_TYPE	Type;
	NODE			Node;
	double			*Beta;
	long long		Iteration;
} DUMMYCODE;

typedef struct
{
	DUMMYCODE	**DummyList;
	int			NoDummyCode;
	int			NoMaxDummy;
	double		*DummyBeta;
} RJDUMMY;

typedef struct
{
	// Total number of steps to use
	int			NoSteps;

	// Current number of horistoal slices
	int			NoSlices;

	// X,Y posititions of slice
	double		*SliceX;
	double		*SliceY;

	// start / end paires of 
	double		*SliceMin;
	double		*SliceMax;

} SLICESAMPLER;

typedef struct
{
	double		*Alpha;
	double		*Scale;

	double		*AnsVect;

	double		*SiteLh;

	double		*SiteMin;
	double		*SiteMax;
	double		*SiteSD;

	SLICESAMPLER	**SliceSamplers;
	
	int				NoSD;
	STABLEDIST**	SDList;

} FATTAILRATES;

typedef struct
{
	// Number of rate to esimate, or max num if RJ is used. 
	int		NoOfRates;

	// Total number of rates to use, inc, including resections and constatns. 
	int		NoOfFullRates;

	// RJ number of rates used currently
	int		NoOfRJRates;
	
	// Values for the rates being estimated. 
	double	*Rates;

	// The names of the paramters being estmated
	char	**RateNames;

	// Values for the all rates, inc constants and resections.  
	double	*FullRates;
	
	// Base frequencies , can all be set to 1 for back capability
	double	*Pis;

	// The current tree being evaluated
	int		TreeNo;
	
	// Lh of the rates
	double	Lh;

	double	LhPrior;
	double	LnHastings;
	double	LnJacobion;

	int			NoPatterns;

	PRIOR		**Priors;
	int			NoPriors;	

	double		*Means;

	double		*Beta;

	double	Delta;
	double	Lambda;
	double	Kappa;
	double	OU;

	double	OnToOff;
	double	OffToOn;
	double	CoVarPis[2];

	int		*MappingVect;


	double	*GammaMults;
	int		GammaCats;
	double	Gamma;
		
	int		HMeanCount;

#ifndef BIG_LH
	double	HMeanSum;
#else
	mpfr_t	HMeanSum;
#endif

	int		UseEstData;
	int		*EstDescData;
	double	*EstData;
	int		*EstDataSiteNo;
	int		NoEstData;

	MODELFILE		*ModelFile;
	int				ModelNo;

	RANDSTATES		*RS;
	RANDSTATES		**RSList;
	gsl_rng			*RNG;

	VARRATES		*VarRates;
	CONTRASTR		*Contrast;
	HETERO			*Hetero;
	RJDUMMY			*RJDummy;
	FATTAILRATES	*FatTailRates;

	int				UseLocalTransforms;
	int				EstLocalTransforms;
	
	int				NoLocalTransforms;
	LOCAL_TRANSFORM	**LocalTransforms;

	DIST_DATE_RATES	*DistDataRates;

	TIME_SLICES		*TimeSlices;

	int				AutoAccept;
	int				CalcLh;

	double			GlobablRate;
	double			NormConst;

	LANDSCAPE		*Landscape;
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

#define NO_SCHEDULE_OPT	28

static char    *SHEDOP[] =
{
	"Rate",
	"CV",
	"Kappa",
	"Delta",
	"Labda",
	"Rev Jump",
	"Prior Change",
	"Est Data",
	"Solo Tree Move",
	"RJ Local Transform Add / Remove",
	"Local Transform Move",
	"Local Transform Change Scale",
	"Local Transform Hyper Prior",
	"Change Hetero",
	"Tree Move",
	"OU",
	"Gamma",
	"RJ Dummy Add / Remove",
	"RJ Dummy Move Node",
	"RJ Dummy Change Beta",
	"Fat Tail Ans All",
	"Fat Tail Ans",
	"Geo Move All",
	"Local Rates",
	"Data Dist",
	"Time Slice - Time",
	"Time Slice - Scale",
	"Global Rate"
};

typedef enum
{
	S_RATES,
	S_CV,
	S_KAPPA,
	S_DELTA,
	S_LABDA,
	S_JUMP,
	S_PPROR,
	S_EST_DATA,
	S_SOLO_TREE_MOVE,
	S_VARRATES_ADD_REMOVE,
	S_VARRATES_MOVE,
	S_VARRATES_CHANGE_SCALE,
	S_VARRATES_HYPER_PRIOR,
	S_HETERO, 
	S_TREE_MOVE,
	S_OU,
	S_GAMMA_MOVE,
	S_RJ_DUMMY, 
	S_RJ_DUMMY_MOVE,
	S_RJ_DUMMY_CHANG_EBETA,
	S_FAT_TAIL_ANS_ALL,
	S_FAT_TAIL_ANS,
	S_GEO_MOVE_ALL,
	S_LOCAL_RATES,
	S_DATA_DIST,
	S_TIME_SLICE_TIME,
	S_TIME_SLICE_SCALE,
	S_GLOBAL_RATE
} OPERATORS;

typedef struct
{
	int		Op;
	int		NoOfOpts;

	int		GNoAcc, GNoTried;
	int		SNoAcc, SNoTried;
	
	double		*OptFreq;
	int			*Tryed;
	int			*Accepted;
	AUTOTUNE	*DataDevAT;
	AUTOTUNE	*VarRateAT;

	int			NoParm;
	AUTOTUNE	**RateDevATList;
	int			*PTried;
	int			*PAcc;
	int			PNo;

	int			NoVarRatesOp;
	double		*FreqVarRatesOp;
	TRANSFORM_TYPE	*VarRatesOp;


	AUTOTUNE	*KappaAT;
	AUTOTUNE	*DeltaAT;
	AUTOTUNE	*LambdaAT;
	AUTOTUNE	*OUAT;

	AUTOTUNE	*GammaAT;

	AUTOTUNE	*RJDummyBetaAT;

	AUTOTUNE	**FullATList;
	int			NoFullATList;

	AUTOTUNE	*LocalRatesAT;

	AUTOTUNE	*TimeSliceTimeAT;
	AUTOTUNE	*TimeSliceScaleAT;

	AUTOTUNE	*CurrentAT;

	AUTOTUNE	*GlobalRateAT;


	int				NoCShed;
	CUSTOM_SCHEDULE	**CShedList;
	double			*DefShed;

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
