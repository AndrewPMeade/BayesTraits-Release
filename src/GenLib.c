/*
*  BayesTriats 4.0
*
*  copyright 2022
*
*  Andrew Meade
*  School of Biological Sciences
*  University of Reading
*  Reading
*  Berkshire
*  RG6 6BX
*
* BayesTriats is free software: you can redistribute it and/or modify
* it under the terms of the GNU General Public License as published by
* the Free Software Foundation, either version 3 of the License, or
* (at your option) any later version.
* 
* This program is distributed in the hope that it will be useful,
* but WITHOUT ANY WARRANTY; without even the implied warranty of
* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
* GNU General Public License for more details.
* 
* You should have received a copy of the GNU General Public License
* along with this program.  If not, see <http://www.gnu.org/licenses/>
*
*/



#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <ctype.h>
#include <time.h>
#include <stdint.h>

#include "GenLib.h"


void		MallocErrFull(char* FileName, int LineNo)
{
	fprintf(stdout, "Memory allocation error in file %s line %d\n", FileName, LineNo);
	fprintf(stdout, "Two main causes of memory allocation error, 1) running out of usable memory, 2) a programming error.\n");

	exit(1);
}

void*	smalloc(size_t n, char* FName, unsigned long LineNo)
{
	void *Ret;

	Ret = malloc(n);

	if(Ret == NULL)
		MallocErrFull(FName, LineNo);
	
	return Ret;
}

void*	sbmalloc(size_t n, char* FName, unsigned long LineNo)
{
	void *Ret;
	
	Ret = smalloc(n, FName, LineNo);

	memset(Ret, 0, n);

	return Ret;
}


char*		StrMake(const char* Str)
{
	char*	Ret;

	Ret = (char*)SMalloc(sizeof(char) * (strlen(Str) + 1));
	
	strcpy(Ret, Str);

	return Ret;
}

FILE*		OpenRead(char *FileName)
{
	FILE*	Ret;

	Ret = fopen(FileName, "r");
	if(Ret == NULL)
	{
		fprintf(stderr, "Could not open file %s for reading\n", FileName);
		exit(1);
	}

	return Ret;
}

FILE*		OpenWrite(char *FileName)
{
	FILE*	Ret;

	Ret = fopen(FileName, "w");
	if(Ret == NULL)
	{
		fprintf(stderr, "Could not open file %s for writting\n", FileName);
		exit(1);
	}

	return Ret;
}

FILE*		OpenAppend(char *FileName)
{
	FILE*	Ret;

	Ret = fopen(FileName, "a");
	if(Ret == NULL)
	{
		fprintf(stderr, "Could not open file %s for appending\n", FileName);
		exit(1);
	}

	return Ret;
}

char*		GetFileNameWithExt(char *Base, char *Ext)
{
	char *Buffer;
	Buffer = (char*)SMalloc(strlen(Base) + strlen(Ext) + 2);

	sprintf(Buffer, "%s%s", Base, Ext);

	return Buffer;

}

FILE*		OpenWriteWithExt(char *Base, char *Ext)
{
	FILE *Ret;
	char *Buffer;

	Buffer = GetFileNameWithExt(Base, Ext);

	Ret = OpenWrite(Buffer);

	free(Buffer);

	return Ret;
}

FILE*		OpenAppendWithExt(char *Base, char *Ext)
{
	FILE *Ret;
	char *Buffer;

	Buffer = GetFileNameWithExt(Base, Ext);

	Ret = OpenAppend(Buffer);

	free(Buffer);

	return Ret;
}

/*
	PC	:	CRLF	013 010
	Mac	:	CR		013
	Unix:	LF		010
*/

char	IsNewLine(char *Char)
{
	char	C1;
	char	C2;

	if(*Char == '\0')
		return 'E';


	C1 = *Char;
	Char++;
	C2 = *Char;

	if((C1 == 13) && (C2 == 10))
		return 'P';

	if(C1 == 13)
		return 'M';

	if(C1 == 10)
		return 'U';

	return 'N';
}

int		FindNoOfNL(char* FileBuffer, int FileSize)
{
	int		BIndex;
	char	Line;
	int		Ret;

	Ret=0;
	for(BIndex=0;BIndex<FileSize;BIndex++)
	{
		Line = IsNewLine(&FileBuffer[BIndex]);

		if(Line != 'N')
			Ret++;

		if(Line == 'P')
			BIndex++;
	}

	Ret++;
	return Ret;
}

int		FindLineLen(char *Start)
{
	int	Ret;

	Ret =0;
	while(IsNewLine(Start) == 'N')
	{
		Start++;
		Ret++;
	}

	return Ret;
}

char**	ProcessBinaryFile(char* FileBuffer, int FileSize, int* NoOfLines)
{
	char	*P;
	int		LineIndex;
	int		LineLen;
	char	NL;
	char	**Ret;

	(*NoOfLines) = FindNoOfNL(FileBuffer, FileSize);

	Ret = (char**)malloc(sizeof(char*) * (*NoOfLines));
	if(Ret  == NULL)
		MallocErr();

	P = FileBuffer;
	for(LineIndex=0;LineIndex<*NoOfLines;LineIndex++)
	{
		LineLen = FindLineLen(P);
		Ret[LineIndex] = (char*)malloc(sizeof(char) * (LineLen + 1));
		if(Ret[LineIndex] == NULL)
			MallocErr();

		memcpy(Ret[LineIndex], P, sizeof(char) * LineLen);
		Ret[LineIndex][LineLen] = '\0';

		P += LineLen;

		NL = IsNewLine(P);
		if((NL == 'U') || (NL == 'M'))
			P++;

		if(NL == 'P')
			P+=2;
	}

	return Ret;
}

char*	RemoveComment(char* FileBuffer, int* FileSize, char *FileName)
{
	char*	Block;
	char*	Ret;
	int		BIndex;
	int		DIndex;
	int		InCom;

	Block = (char*)malloc(sizeof(char) * (*FileSize)+1);
	if(Block == NULL)
		MallocErr();

	InCom = FALSE;
	BIndex=0;
	for(DIndex=0;DIndex<(*FileSize)+1;DIndex++)
	{
		if(FileBuffer[DIndex] == '[')
		{
			if(InCom == TRUE)
			{
				printf("File %s has nested comments\n", FileName);
				exit(0);
			}
			InCom = TRUE;
		}

		if(InCom != TRUE)
		{
			Block[BIndex] = FileBuffer[DIndex];
			BIndex++;
		}

		if(FileBuffer[DIndex] == ']')
		{
			if(InCom == FALSE)
			{
				printf("Found close comment with out opening in %s\n", FileName);
				exit(0);
			}
			InCom = FALSE;
		}
	}

	BIndex--;
	*FileSize = BIndex;

	free(FileBuffer);

	Ret = (char*)malloc(sizeof(char) * (*FileSize) + 1);
	if(Ret == NULL)
		MallocErr();

	memcpy(Ret, Block, sizeof(char) * (*FileSize) + 1);
	free(Block);

	return Ret;

}

void		FindTextFileMaxLine(TEXTFILE* TextFile)
{
	int	Index;
	int	Len;

	if(TextFile->NoOfLines == 0)
	{
		TextFile->MaxLine = 0;
		return;
	}

	TextFile->MaxLine = (int)strlen(TextFile->Data[0]);
	for(Index=1;Index<TextFile->NoOfLines;Index++)
	{
		Len = (int)strlen(TextFile->Data[Index]);
		if(Len > TextFile->MaxLine)
			TextFile->MaxLine = Len;
	}
}

TEXTFILE*	LoadTextFile(char* Name, char DelComments)
{
	TEXTFILE*	Ret;
	FILE*		InFile;
	char*		Buffer;
	int			NoRead;


	Ret = (TEXTFILE*)malloc(sizeof(TEXTFILE));
	if(Ret == NULL)
		MallocErr();

	Ret->FileName	= StrMake(Name);
	Ret->Size		= 0;
	Ret->NoOfLines	= 0;
	Ret->MaxLine	= 0;
	Ret->Data		= NULL;

	InFile = fopen(Name, "rb");
	if(InFile == NULL)
	{
		fprintf(stdout, "Could not open file %s for reading\n", Name);
		fprintf(stderr, "Could not open file %s for reading\n", Name);
		exit(1);
	}

	Buffer = (char*)malloc(sizeof(char) * BUFFERSIZE);
	if(Buffer == NULL)
		MallocErr();

	do
	{
		NoRead = (int)fread(&Buffer[0], sizeof(char), BUFFERSIZE, InFile);
		Ret->Size += NoRead;
	} while(NoRead == BUFFERSIZE);

	free(Buffer);

	Buffer = (char*)malloc(sizeof(char) * (Ret->Size + 1));
	if(Buffer == NULL)
		MallocErr();

	rewind(InFile);
	fread(&Buffer[0], sizeof(char), Ret->Size, InFile);
	fclose(InFile);
	Buffer[Ret->Size] = '\0';

	if(DelComments == TRUE)
		Buffer = RemoveComment(Buffer, &Ret->Size, Name);

	Ret->Data = ProcessBinaryFile(Buffer, Ret->Size, &Ret->NoOfLines);

	FindTextFileMaxLine(Ret);

	free(Buffer);
	return Ret;
}

void	FreeTextFile(TEXTFILE* TextFile)
{
	int	Index;

	free(TextFile->FileName);

	for(Index=0;Index<TextFile->NoOfLines;Index++)
		free(TextFile->Data[Index]);
	free(TextFile->Data);

	free(TextFile);
}

int		IsValidFastaStart(char* Line)
{
	int	Index;

	Index = 0;
	while((Line[Index] == ' ' ) || (Line[Index] == '\t'))
		Index++;

	if(Line[Index] == '>')
		return TRUE;

	return FALSE;
}

char*		MergeStrings(char* Old, char* Add)
{
	char*	Ret;

	Ret = (char*)malloc(sizeof(char) * (strlen(Old) + strlen(Add) + 1));
	if(Ret == NULL)
		MallocErr();
	strcpy(Ret, Old);
	strcat(Ret, Add);

	free(Old);

	return Ret;
}

char*	RemoveFastaTag(char* Str)
{
	char*	DelChar;
	char*	Ret;

	DelChar = Str;

	while(*Str != '>')
		Str++;
	Str++;

	Ret = (char*)malloc(sizeof(char) * (strlen(Str) + 1));
	if(Ret==NULL)
		MallocErr();

	strcpy(Ret, Str);

	free(DelChar);

	return Ret;
}

FASTA*		LoadFasta(char *FileName)
{
	FASTA*		Ret;
	TEXTFILE*	File;
	int			Index;
	int			No;

	Ret = (FASTA*)malloc(sizeof(FASTA));
	if(Ret == NULL)
		MallocErr();

	File = LoadTextFile(FileName, FALSE);

	Ret->NoOfSeq = 0;
	for(Index=0;Index<File->NoOfLines;Index++)
		if(IsValidFastaStart(File->Data[Index]) == TRUE)
			Ret->NoOfSeq++;

	Ret->Seq = (char**)malloc(sizeof(char*) * Ret->NoOfSeq);
	Ret->Tags= (char**)malloc(sizeof(char*) * Ret->NoOfSeq);
	if((Ret->Seq == NULL) || (Ret->Tags == NULL))
		MallocErr();

	for(Index=0;Index<Ret->NoOfSeq;Index++)
	{
		Ret->Seq[Index]	= StrMake("\0");
		Ret->Tags[Index]= StrMake("\0");
	}

	Index = 0;
	while(IsValidFastaStart(File->Data[Index]) == FALSE)
		Index++;

	No= 0;
	Ret->Tags[0] = MergeStrings(Ret->Tags[0], File->Data[Index]);
	Index++;
	for(;Index<File->NoOfLines;Index++)
	{
		if(IsValidFastaStart(File->Data[Index]) == TRUE)
		{
			No++;
			Ret->Tags[No] = MergeStrings(Ret->Tags[No], File->Data[Index]);
		}
		else
		{
			if(strlen(File->Data[Index]) > 0)
				Ret->Seq[No] = MergeStrings(Ret->Seq[No], File->Data[Index]);
		}
	}

	for(Index=0;Index<Ret->NoOfSeq;Index++)
		Ret->Tags[Index] = RemoveFastaTag(Ret->Tags[Index]);


	FreeTextFile(File);

	return Ret;
}


void FreeFasta(FASTA* Fasta)
{
	int	Index;

	for(Index=0;Index<Fasta->NoOfSeq;Index++)
	{
		free(Fasta->Seq[Index]);
		free(Fasta->Tags[Index]);
	}

	free(Fasta->Seq);
	free(Fasta->Tags);
	free(Fasta);
}

void		MakeUpper(char* Str)
{
	int Size;
	int	Index;

	Size = (int)strlen(Str);
	for(Index=0;Index<Size;Index++)
		Str[Index] = toupper(Str[Index]);
}

void		MakeLower(char* Str)
{
	int Size;
	int	Index;

	Size = (int)strlen(Str);
	for(Index=0;Index<Size;Index++)
		Str[Index] = tolower(Str[Index]);
}

int	IsSapce(char C)
{
	if(C == ' ' || C == '\t')
		return 1;

	return 0;
}

int	MakeArgv(char*	string, char *argv[], int argvsize)
{
	char*	p = string;
	int		i;
	int		argc = 0;

	for(i=0; i<argvsize; i++)
	{
		/* Skip Leading whitespace */
		while(IsSapce(*p)) p++;

		if(*p != '\0')
			argv[argc++] = p;
		else
		{
			argv[argc] = 0;
			break;
		}

		/* Scan over arg */
		while(*p != '\0' && !IsSapce(*p))
			p++;

		/* Terminate argv */
		if(*p != '\0' && i < argvsize - 1)
			*p++ = '\0';
	}

	return argc;
}

int			IsChar(char Break, char C)
{
	if(C == Break)
		return FALSE;

	return TRUE;
}

int			MakeArgvChar(char*	string, char *argv[], int argvsize, char Break)
{
	char*	p = string;
	int		i;
	int		argc = 0;

	for(i=0; i<argvsize; i++)
	{
		/* Skip Leading whitespace */
		while(*p == Break) p++;

		if(*p != '\0')
			argv[argc++] = p;
		else
		{
			argv[argc] = 0;
			break;
		}

		/* Scan over arg */
		while(*p != '\0' && !(*p == Break))
			p++;

		/* Terminate argv */
		if(*p != '\0' && i < argvsize - 1)
			*p++ = '\0';

	}

	return argc;
}

int		IsValidDouble(char* Str)
{
	int	Point;

	if(atof(Str) != 0.0)
		return TRUE;

	/* Skip Leading whitespace */
	while(isspace(*Str)) Str++;

	Point = FALSE;

	if(*Str == '-')
		Str++;

	while(*Str)
	{
		if(*Str == '.')
		{
			if(Point == FALSE)
				Point = TRUE;
			else
				return FALSE;
		}
		else
		{
			if(!((*Str >= '0') && (*Str <= '9')))
				return FALSE;
		}

		Str++;
	}

	return TRUE;
}

int		IsValidInt(char* Str)
{
	if(atoi(Str) != 0)
		return TRUE;

	if(*Str == '-')
		Str++;

	if(*Str == '+')
		Str++;

	while(*Str)
	{
		if(!((*Str >= '0') && (*Str <= '9')))
			return FALSE;
		Str++;
	}

	return TRUE;
}


double*		LoadDouble(char *FileName, int *No)
{
	TEXTFILE*	Data;
	int			Index;
	double		*Ret;


	Ret	= NULL;
	Data = LoadTextFile(FileName, TRUE);

	*No = 0;
	for(Index=0;Index<Data->NoOfLines;Index++)
	{
		if(IsValidDouble(Data->Data[Index]) == TRUE)
			(*No)++;
	}

	if(*No == 0)
	{
		*No = -1;
		FreeTextFile(Data);
		return NULL;
	}

	Ret = (double*)malloc(sizeof(double) * (*No));
	if(No == NULL)
		MallocErr();

	*No = 0;
	for(Index=0;Index<Data->NoOfLines;Index++)
	{
		if(IsValidDouble(Data->Data[Index]) == TRUE)
		{
			Ret[*No] = atof(Data->Data[Index]);
			(*No)++;
		}

	}

	FreeTextFile(Data);

	return Ret;
}

void	ReplaceChar(char Rep, char With, char* String)
{
	int	Index;
	int	Size;

	Size = (int)strlen(String);

	for(Index=0;Index<Size;Index++)
	{
		if(String[Index] == Rep)
			String[Index] = With;
	}

}

void	RemoveChar(char c, char* String)
{
	int	Size;
	int	Index;
	int	SIndex;

	Size = (int)strlen(String);

	for(Index=0;Index<Size;Index++)
	{
		if(String[Index] == c)
		{
			for(SIndex=Index;SIndex<Size;SIndex++)
				String[SIndex] = String[SIndex+1];
			Index--;
		}
	}
}

char*	FormatInt(int No, int Size)
{
	char*	Buffer;
	char*	Ret;
	int		Index;
	int		Len;

	Buffer = (char*)malloc(sizeof(char) * BUFFERSIZE);
	if(Buffer == NULL)
		MallocErr();

	sprintf(Buffer, "%d", No);
	if((signed)strlen(Buffer) >= Size)
	{
		Ret = (char*)malloc(sizeof(char) * (strlen(Buffer) + 1));
		if(Ret == NULL)
			MallocErr();
		sprintf(Ret, "%d", No);
		free(Buffer);
		return Ret;
	}

	Ret = (char*)malloc(sizeof(char) * (Size + 1));
	if(Ret == NULL)
		MallocErr();

	Len = (int)strlen(Buffer);
	for(Index=0;Index<Size - Len;Index++)
		Ret[Index] = '0';
	Ret[Index] = '\0';
	strcat(Ret, Buffer);

	free(Buffer);

	return Ret;
}

void		WriteFasta(FILE* Str, char* Tag, char *Seq, int CharPerLine)
{
	int	C;

	C = 0;
	fprintf(Str, "> %s\n", Tag);

	while(*Seq!= '\0')
	{
		fprintf(Str, "%c", *Seq);
		C++;
		if(C == CharPerLine)
		{
			fprintf(Str, "\n");
			C=0;
		}
		Seq++;
	}

	fprintf(Str, "\n");
}

void		SaveFasta(char *FileName, FASTA* Fasta)
{
	FILE*	OutFile;
	int		Index;

	OutFile = OpenWrite(FileName);

	for(Index=0;Index<Fasta->NoOfSeq;Index++)
	{
		WriteFasta(OutFile, Fasta->Tags[Index], Fasta->Seq[Index], 80);
	}

	fclose(OutFile);
}

void    CharSwap(char * A, char * B)
{
    char    C;

    C   = *A;
    *A  = *B;
    *B  = C;
}

void    revstr(char * String)
{
    char    *Start;
    char    *End;
    int     Size;

    Size = (int)strlen(String);

    Start = String;
    End   = String + (Size - 1);

    while(Start < End)
    {
        CharSwap(Start, End);
        Start++;
        End--;
    }
}

void	MakeComplement(char*	DNA)
{
	revstr(DNA);

	while(*DNA != '\0')
	{
		*DNA = toupper(*DNA);
		switch(*DNA)
		{
			case 'A':	*DNA = 'T'; break;
			case 'T':	*DNA = 'A'; break;
			case 'C':	*DNA = 'G'; break;
			case 'G':	*DNA = 'C'; break;
			case 'N':	*DNA = 'N'; break;
		}
		DNA++;
	}
}

/*
typedef struct
{
	double	**Data;
	int		NoOfLines;
	int		*NoPerLine;
} NUMFILE;

NUMFILE*	LoadNumFile(char *FileName);
void		FreeNumFile(NUMFILE *NumFile);
*/

NUMFILE*	LoadNumFile(char *FileName)
{
	NUMFILE*	Ret;
	TEXTFILE*	InFile;
	int			Index;
	char		**Passed;
	char		*Buffer;
	int			Tokes;
	int			LIndex;
	int			No;

	InFile = LoadTextFile(FileName, FALSE);

	Ret = (NUMFILE*)malloc(sizeof(NUMFILE));
	if(Ret == NULL)
		MallocErr();

	Ret->NoOfLines = InFile->NoOfLines;

	Ret->Data = (double**)malloc(sizeof(double*) * Ret->NoOfLines);
	if(Ret->Data == NULL)
		MallocErr();

	Ret->NoPerLine = (int*)malloc(sizeof(int) * Ret->NoOfLines);
	if(Ret->NoPerLine == NULL)
		MallocErr();

	Buffer = (char*)malloc(sizeof(char) * (InFile->MaxLine + 1));
	if(Buffer == NULL)
		MallocErr();

	Passed = (char**)malloc(sizeof(char*) * InFile->MaxLine);
	if(Passed == NULL)
		MallocErr();

	No = 0;
	for(Index=0;Index<InFile->NoOfLines;Index++)
	{
		strcpy(Buffer, InFile->Data[Index]);

		Tokes = MakeArgv(Buffer, Passed, InFile->MaxLine);

		if(Tokes > 0)
		{
			Ret->NoPerLine[No] = Tokes;

			Ret->Data[No] = (double*)malloc(sizeof(double) * Tokes);
			if(Ret->Data[No] == NULL)
				MallocErr();

			for(LIndex=0;LIndex<Tokes;LIndex++)
			{
				if(IsValidDouble(Passed[LIndex]) == FALSE)
				{
					printf("%s::%d Error: Token %s is not a valid number, Line %d, Token number %d\n", __FILE__, __LINE__, Passed[LIndex], Index, LIndex);
					exit(0);
				}

				Ret->Data[No][LIndex] = atof(Passed[LIndex]);
			}

			No++;
		}
	}

	Ret->NoOfLines = No;

	free(Buffer);
	free(Passed);

	FreeTextFile(InFile);

	return Ret;
}

void		FreeNumFile(NUMFILE *NumFile)
{
	int	Index;

	for(Index=0;Index<NumFile->NoOfLines;Index++)
		free(NumFile->Data[Index]);

	free(NumFile->NoPerLine);

	free(NumFile);
}

void	PrintFixSize(char *String, int Size, FILE* Str)
{
	int		Index;
	int		StrSize;

	StrSize = (int)strlen(String);

	fprintf(Str, "%s", String);

	for(Index=0;Index<Size-StrSize;Index++)
		fprintf(Str, " ");
}

void Swap(void** a, void** b)
{
	void* t;

	t = *a;
	*a = *b;
	*b = t;
}


void	GotoFileEnd(FILE *File, char *Buffer, int Size)
{
	while(fgets(Buffer, Size, File) != NULL);
}

void	PrintTime(FILE* Str)
{
	time_t	*Now;
	struct tm *Time;

	Now = (time_t*)malloc(sizeof(time_t));

	time(Now);
	Time = localtime(Now);

	fprintf(Str, "%02d/%02d/%d ", Time->tm_mday, Time->tm_mon+1, Time->tm_year+1900);
	fprintf(Str, "%d:%02d:%02d", Time->tm_hour, Time->tm_min, Time->tm_sec);

	free(Now);
}

void**	AddToList(int *No, void** OldList, void* Item)
{
	void	**Ret;

	if(*No == 0)
	{
		Ret = (void**)malloc(sizeof(void*));
		if(Ret == NULL)
			MallocErr();
		Ret[0] = Item;
		(*No)++;
		return Ret;
	}

	Ret = (void**)malloc(sizeof(void*) * (*No + 1));
	if(Ret == NULL)
		MallocErr();

	memcpy(Ret, OldList, *No * sizeof(void*));
	free(OldList);


	Ret[*No] = Item;

	(*No)++;

	return Ret;
}


void double2HexString(double a) 
{ 
//   char *buf; // double is 8-byte long, so we have 2*8 + terminating \0 
   char *d2c; 
//	char *n; 
   int i; 

  // buf = (char*)malloc(sizeof(double)+sizeof(double));

 //  n = buf;

   d2c = (char *) &a; 

   for(i = 0; i < 8; i++) 
   { 
   //   sprintf(n, "%02X", *d2c++); 
	  printf("%x", *d2c++); 
 //     n += 2; 
   }  
  // *(n) = '\0'; 

//   return buf;
} 

char		IntToHex(int i)
{
	if(i<10)
		return '0'+i;
	return 'A'+(i-10);
}

void		PrintDoubleHex(FILE *Str, double D)
{
	unsigned char *h, in;
	int Index;
	unsigned char t,b;

	h = (unsigned char*)&D;
	
	for(Index=0;Index<sizeof(double);Index++)
	{

		in = h[Index];

		b = in & 0x0f;
		t = (in >> 4) & 0x0f;

		printf("%c%c", IntToHex(b), IntToHex(t));
	}
}

void		CalcRSqr(double *x, double *y, int Size, double *R2, double *Slope, double *Intercept)
{
	double	sumOfX, sumOfY, sumOfXSq, sumOfYSq, ssX, ssY, sumCoVar;
	double	MeanX, MeanY, RDenom, RNum, sCo;
	int Index;
	
	sumOfX = sumOfY = sumOfXSq = sumOfYSq = ssX = ssY = sumCoVar = 0;
	MeanX =  MeanY = 0;

	for(Index=0;Index<Size;Index++)
	{
		sumCoVar += (x[Index] * y[Index]);
		
		sumOfX += x[Index];
		sumOfY += y[Index];

		sumOfXSq += (x[Index] * x[Index]);
		sumOfYSq += (y[Index] * y[Index]);
	}

	ssX = sumOfXSq - ((sumOfX * sumOfX) / Size);
	ssY = sumOfYSq - ((sumOfY * sumOfY) / Size);

	RNum = (Size * sumCoVar) - (sumOfX * sumOfY);
	
	RDenom = (Size * sumOfXSq - (pow(sumOfX, 2))) * (Size * sumOfYSq - (pow(sumOfY, 2)));

	sCo = sumCoVar - ((sumOfX * sumOfY) / Size);
	
	MeanX = sumOfX / Size;
	MeanY = sumOfY / Size;

	*Slope = sCo / ssX;
	*Intercept = MeanY - (*Slope * MeanX);
	*R2 = RNum / sqrt(RDenom);
	*R2 = *R2 * *R2;
}

int			CountChar(char *Str, char C)
{
	int Ret;

	Ret = 0;
	while(*Str != '\0')
	{
		if(*Str == C)
			Ret++;
		Str++;
	}

	return Ret;
}

void*		CloneMem(size_t Size, void *Mem)
{
	void* Ret;

	Ret = SMalloc(Size);
	memcpy(Ret, Mem, Size);

	return Ret;
}

int StrICmp(char const *a, char const *b)
{
	int d;

	for (;; a++, b++) {
		d = tolower(*a) - tolower(*b);
		if (d != 0 || !*a)
			return d;
	}
}


void	NormaliseVector(double *Vect, int Size)
{
	double	SF;
	int		Index;

	SF = 0;
	for(Index=0;Index<Size;Index++)
		SF += Vect[Index];

	SF = 1 / SF;

	for(Index=0;Index<Size;Index++)
		Vect[Index] *= SF;
}


union DoubleStr 
{
    double Double;
    uint8_t i8[sizeof(double) + 1];
};


void	PrintDoubleHexStr(double Double)
{
	union DoubleStr DS;
	int i;

	DS.Double = Double;

	for(i = 0; i < (int) sizeof(DS.i8) - 1; ++i) 
	{
		printf("%02X", DS.i8[i]);

	}
}

char*	DoubleToHexStr(double Double)
{
	char *Ret, *P;
	union DoubleStr DS;
	size_t Size;
	int i;

	DS.Double = Double;

	Size = sizeof(DS.i8) * 2;

	Ret = (char*)SMalloc(Size);
	memset(Ret, 0, Size);
	P = Ret;

	for(i = 0; i < (int) sizeof(DS.i8) - 1; ++i) 
	{
		sprintf(P, "%02X", DS.i8[i]);
		P+=2;		
	}

	return Ret;
}

double	HexStrToDouble(char *Str)
{	
	union DoubleStr DS;
	int j, i, k;
	
	j = 0;


	for(i = 0; i < (int) strlen(Str); i += 2) 
	{
		uint8_t val = 0;
		for(k = 0; k < 2; ++k) 
		{
			char ch = Str[i + k];

			val *= 16;
			if (ch >= '0' && ch <= '9') 
			{
				val += ch - '0';
			}
			else if (ch >= 'A' && ch <= 'F') 
			{
				val += ch - 'A' + 10;
			}
		}
		DS.i8[j++] = val;
	}

	return DS.Double;
}

int		FileExists(char *FName)
{
	FILE *FIn;

	FIn = fopen(FName, "rb");
	if(FIn == NULL)
		return FALSE;

	fclose(FIn);

	return TRUE;

}

FILE*	OpenWithExt(int Append, char *Base, char *Ext)
{
	FILE *File;

	if(Append == FALSE)
		File = OpenWriteWithExt(Base, Ext);
	else
		File = OpenAppendWithExt(Base, Ext);

	return File;
}
