#if !defined GENLIB
#define GENLIB

#pragma warning(disable : 4996)


#ifndef TRUE
	#define TRUE 1
#endif

#ifndef FALSE
	#define	FALSE 0
#endif

#define	BUFFERSIZE	1048576

#define	MallocErr() MallocErrFull(__FILE__, __LINE__)

typedef struct
{
	char	*FileName;
	int		NoOfLines;
	int		MaxLine;
	int		Size;
	char	**Data;
} TEXTFILE;

typedef struct
{
	int		NoOfSeq;
	char	**Seq;
	char	**Tags;
} FASTA;

typedef struct
{
	double	**Data;
	int		NoOfLines;
	int		*NoPerLine;
} NUMFILE;

TEXTFILE*	LoadTextFile(char* Name, char DelComments);
void		FreeTextFile(TEXTFILE* TextFile);


FASTA*		LoadFasta(char *FileName);
void		SaveFasta(char *FileName, FASTA* Fasta);
void		WriteFasta(FILE* Str, char* Tag, char *Seq, int CharPerLine);
void		FreeFasta(FASTA* Fasta);

NUMFILE*	LoadNumFile(char *FileName);
void		FreeNumFile(NUMFILE *NumFile);

FILE*		OpenWrite(char *FileName);
FILE*		OpenRead(char *FileName);
void		MallocErrFull(char* FileName, int LineNo);
char*		StrMake(char* Str);
void		MakeUpper(char* Str);
void		MakeLower(char* Str);
void		RemoveChar(char c, char* String);
int			MakeArgv(char*	string, char *argv[], int argvsize);


int			IsValidInt(char* Str);
int			IsValidDouble(char* Str);

double*		LoadDouble(char *FileName, int *No);

void		RemoveChar(char c, char* String);
void		ReplaceChar(char Rep, char With, char* String);

char*		FormatInt(int No, int Size);
void		revstr(char * String);
void		MakeComplement(char* DNA);

void 		Swap(void** a, void** b);
void		PrintFixSize(char *String, int Size, FILE* Str);

#endif
