#if !defined (TRESSPASSER_HEADDER)
#define TRESSPASSER_HEADDER

static char    *TREETAGS[] =
{
	"utree",
	"rtree",
	"tree", 
	""
};

typedef struct
{
	char*	Name;
	int		No;
	int*	Data;
} NTAXA;

typedef struct
{
	int		No;
	int		*TaxaNo;
	int		Freq;
	double	Prob;
} PARTITION;

struct NINODE
{
	
	int				NoOfNodes;
	struct NINODE**	NodeList;
	struct NINODE*	Left;
	struct NINODE*	Right;
	struct NINODE*	Ans;
	
	double			Length;
	int				Tip;
	int				TaxaID;
	NTAXA*			Taxa;

	PARTITION		*Part;	

	int				Touched;
	int				State;
};

typedef struct NINODE*	NNODE;

typedef struct
{
	char*		Tag;
	int			NoOfNodes;

	int			Binary;

	NNODE		Root;
	NNODE*		NodeList;

	PARTITION**	PartList;
	int			NoOfParts;
} NTREE;

typedef struct
{
	int		NoOfTrees;
	int		NoOfTaxa;
	NTAXA*	Taxa;
	NTREE*	Trees;

	PARTITION**	PartList;
	int			NoOfParts;

	NTREE*	ConTrees;
	int		Sites;
} NTREES;

NTREES*	LoadNTrees(char* FileName, char** Err);
void	FreeNTrees(NTREES* Trees);

void	FindNoOfNodes(NNODE N, int *No);
void	FlattenNodes(NTREE *Tree);
void	MakeConsensus(NTREES *Trees);
void	ResolveTreesDet(NTREES *Trees, double Len);

void	PrintNTreeHeadder(FILE* Str, NTREES *Trees);
void	PrintNTree(FILE* Str, NNODE N);
void	PrintNTrees(FILE* Str, NTREES *Trees);
void	SaveNTrees(char* FileName, NTREES *Trees);

int		BinaryTree(NTREE *Tree);

double	GetTreeLength(NTREES* Trees, int TreeNo);
double	GetPathLength(NNODE Node);
double	AveBL(NTREE* Trees);

int		GetTaxaID(NTREES* Trees, char *TaxaName);
void	MakeUltrametric(NTREE* Tree, double Length);

void	NNI(NTREE* T);
void	SetBinaryTree(NTREE* Tree);

#endif
