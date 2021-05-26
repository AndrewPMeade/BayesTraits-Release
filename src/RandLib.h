#if !defined RANDLIB
#define RANDLIB

#pragma warning(disable : 4996)

#define	RSMallocErr() RSMallocErrFull(__FILE__, __LINE__)


#define		DEFSTATSIZE	25
#define		MAGIC		7

#define		PSIZE		1000

/* initial 25 seeds, change as you wish */
static long DEFRANDSTAT[DEFSTATSIZE] =
{ 
	0x95f24dab, 0x0b685215, 0xe76ccae7, 0xaf3ec239, 0x715fad23,
	0x24a590ad, 0x69e4b5ef, 0xbf456141, 0x96bc1b7b, 0xa7bdf825,
	0xc1de75b7, 0x8858a9c9, 0x2da87693, 0xb657f9dd, 0xffdc8a9f,
	0x8121da71, 0x8b823ecb, 0x885d05f5, 0x4e20cd47, 0x5a9ad5d9,
	0x512c0c03, 0xea857ccd, 0x4cc1d30f, 0x8891a8a1, 0xa6b7aadb
};

/* this is magic vector `a', don't change */
static unsigned long mag01[2]=
{ 
	0x0, 0x8ebfd028 
};

static const long A = 48271L;
static const long M = 2147483647L;
static const long Q = 2147483647L / 48271L;
static const long R = 2147483647L % 48271L;

typedef struct
{
	int				Size;
	long			*States;
	long			Seed;
	long			PID;
	int				Pos;
	unsigned long	TempSeed;
} RANDSTATES;

RANDSTATES*		CreateRandStates(void);
RANDSTATES*		CreateSeededRandStates(long Seed);

unsigned long	GetProcID(void);


unsigned long	GetSeed(void);
RANDSTATES*		CreateRandDefStates(void);
void			FreeRandStates(RANDSTATES* RS);

double			RandDouble(RANDSTATES* RS);
unsigned long	RandUSLong(RANDSTATES* RS);
long			RandLong(RANDSTATES* RS);
unsigned int	RandUSInt(RANDSTATES* RS);
int				RandInt(RANDSTATES*	RS);

double			RandUniDouble(RANDSTATES* RS, double Min, double Max);

int				RandIntBetween(RANDSTATES* RS, int l, int t);
int				RandIntRange(RANDSTATES* RS, int Range);

int				RandInProportion(RANDSTATES* RS, double *List, int No);

int				RandPoisson(RANDSTATES* RS, double ExpectedValue);
double			NegExp(RANDSTATES* RS, double ExpectedValue);
double			RandNormal(RANDSTATES* RS, double Mean, double Std);
void			DirichletRandomVariable (RANDSTATES* RS, double *alp, double *z, int n);
void			FindPriod(void);

void			SaveStates(RANDSTATES* RS, char *FName);
RANDSTATES*		LoadStates(char* FName);
void			DumpStates(FILE* Str, RANDSTATES *RS);

void			SaveDieHardTest(RANDSTATES* RS, char *FName, int NoLong);
#endif
