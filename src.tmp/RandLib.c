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


#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <time.h>

#ifdef _WIN32
	#include <windows.h>
#else
	#include <unistd.h>
#endif

#include "RandLib.h"

void		RSMallocErrFull(char* FileName, int LineNo)
{
	fprintf(stderr, "Memory allocation error in file %s line %d\n", FileName, LineNo);
	fprintf(stdout, "Memory allocation error in file %s line %d\n", FileName, LineNo);

	exit(1);
}

RANDSTATES*	AllocRandStates(int Size)
{
	RANDSTATES*	Ret;

	Ret = (RANDSTATES*)malloc(sizeof(RANDSTATES));
	if(Ret == NULL)
		RSMallocErr();

	Ret->States = (long*)malloc(sizeof(long) * Size);
	if(Ret->States == NULL)
		RSMallocErr();

	Ret->Seed	= 0;
	Ret->Pos	= 0;
	Ret->PID	= 0;
	Ret->Size	= Size;

	return Ret;
}

RANDSTATES*	CreateRandDefStates(void)
{
	RANDSTATES*	Ret;

	Ret = AllocRandStates(DEFSTATSIZE);

	memcpy(Ret->States, &DEFRANDSTAT[0], sizeof(unsigned long) * Ret->Size);

	return Ret;
}

unsigned long RandomLong(RANDSTATES* RS)
{
	long TmpSeed;

	TmpSeed = A * (RS->TempSeed % Q ) - R * ( RS->TempSeed / Q );

	if(TmpSeed >= 0)
		RS->TempSeed = TmpSeed;
	else
		RS->TempSeed = TmpSeed + M;

	return RS->TempSeed;
}

RANDSTATES*	CreateSeededRandStates(long Seed)
{
	RANDSTATES*	Ret;
	int	Index;

	Ret = AllocRandStates(DEFSTATSIZE);

	Ret->Seed = Seed;
	Ret->TempSeed = Seed;

	for(Index=0;Index<Ret->Size;Index++)
		Ret->States[Index] = RandomLong(Ret);

	return Ret;
}

RANDSTATES**	CreateRandStatesList(RANDSTATES* RS, int No)
{
	RANDSTATES** Ret;
	int Index;

	Ret = (RANDSTATES**)malloc(sizeof(RANDSTATES*) * No);
	if(Ret == NULL)
		RSMallocErr();
	
	for(Index=0;Index<No;Index++)
		Ret[Index] = CreateSeededRandStates(RandomLong(RS));
	
	return Ret;
}

unsigned long	GetProcID(void)
{
	#ifdef _WIN32
		return GetCurrentProcessId();
	#else
		return getpid();
	#endif
}


unsigned long ReverseUSLong(unsigned long x)
{
	unsigned long h;
	int i;

	h = 0;

	for(i = 0; i < 32; i++)
	{
		h = (h << 1) + (x & 1);
		x >>= 1;
	}

	return h;
}


unsigned long	GetSeed(void)
{
	unsigned long Seed;
	unsigned long Pid;

	Seed = (unsigned)time(NULL);
	Pid = GetProcID();
	Pid = ReverseUSLong(Pid);

	return Seed ^ Pid;
}

RANDSTATES*	CreateRandStates(void)
{
	RANDSTATES*	Ret;
	unsigned long Seed;

	Seed = GetSeed();

	Ret = CreateSeededRandStates(Seed);

	return Ret;
}

void		FreeRandStates(RANDSTATES* RS)
{
	free(RS->States);
	free(RS);
}

/* A C-program for TT800 : July 8th 1996 Version */
/* by M. Matsumoto, email: matumoto@math.keio.ac.jp */
/* genrand() generate one pseudorandom number with double precision */
/* which is uniformly distributed on [0,1]-interval */
/* for each call.  One may choose any initial 25 seeds */
/* except all zeros. */

/* See: ACM Transactions on Modelling and Computer Simulation, */
/* Vol. 4, No. 3, 1994, pages 254-266. */

unsigned long RandUSLong(RANDSTATES*	RS)
{
	unsigned long y;
	int kk;

	/* generate N words at one time */
	if (RS->Pos == RS->Size)
	{/*
		for (kk=0;kk<RS->Size-MAGIC;kk++)
			RS->States[kk] = RS->States[kk+MAGIC] ^ (RS->States[kk] >> 1) ^ mag01[(unsigned long)RS->States[kk] % 2];

		for (;kk<RS->Size;kk++)
			RS->States[kk] = RS->States[kk+(MAGIC-RS->Size)] ^ (RS->States[kk] >> 1) ^ mag01[(unsigned long)RS->States[kk] % 2];
		*/

		for (kk=0;kk<RS->Size-MAGIC;kk++)
			RS->States[kk] = (unsigned long)RS->States[kk+MAGIC] ^ ((unsigned long)RS->States[kk] >> 1) ^ mag01[(unsigned long)RS->States[kk] % 2];

		for (;kk<RS->Size;kk++)
			RS->States[kk] = (unsigned long)RS->States[kk+(MAGIC-RS->Size)] ^ ((unsigned long)RS->States[kk] >> 1) ^ mag01[(unsigned long)RS->States[kk] % 2];

		RS->Pos = 0;
    }

	y = RS->States[RS->Pos];
	y ^= (y << 7) & 0x2b5b2500;		/* s and b, magic vectors */
	y ^= (y << 15) & 0xdb8b0000;	/* t and c, magic vectors */
	y &= 0xffffffff;				/* you may delete this line if word size = 32 */

	y ^= (y >> 16);					/* added to the 1994 version */
	RS->Pos++;

/*	To make rand doubles between 0 and 1 */
/*	return( (double) y / (unsigned long) 0xffffffff); */
	return y;
}

long RandLong(RANDSTATES*	RS)
{
	return (long)RandUSLong(RS);
}

double RandDouble(RANDSTATES* RS)
{
	unsigned long y;

	y = RandUSLong(RS);

//	0 - 1 Exclusive
	return (((unsigned)y + 1.0)/0x100000002);

	//	0 - 1 Inclusive
//	return( (double) y / (unsigned long) 0xffffffff);
}


unsigned int RandUSInt(RANDSTATES* RS)
{
	return (unsigned int)RandUSLong(RS);
}

int RandInt(RANDSTATES*	RS)
{
	return (int)RandUSLong(RS);
}


int			RandIntBetween(RANDSTATES* RS, int l, int t)
{
	unsigned long y;

	y = RandUSLong(RS);

	return  (y%(t-l+1)+l);
}

double RandNormal(RANDSTATES* RS, double Mean, double Std)
{
 /* gives a distribution with mean 0 and std 1.

  To change the mean to M, simply add M to whatever
  is returned

  To change the std to S, simply multiply whatever is returned
  by S. Do the mult first.

  Eg: this returns Z (the thing in the return line)

  for mean M and std S, instead return M + S*Z
*/
	double a, b;
	double Norm;
	double pi = 3.14159265358979323846;

	a = RandDouble(RS);
	b = RandDouble(RS);

	Norm = sqrt(-2.0*log(a)) * cos(2*pi*b);
	Norm = Norm * Std;
	Norm += Mean;
	return Norm;
}

void	SaveDieHardTest(RANDSTATES* RS, char *FName, int NoLong)
{
	FILE*	Out;
	int		Index;
	unsigned long *OList;

	Out = fopen(FName, "wb");
	if(Out == NULL)
	{
		printf("Could not open file %s for writting.\n", FName);
		exit(0);
	}

	OList = (unsigned long*)malloc(sizeof(unsigned long) * NoLong);
	if(OList == NULL)
		RSMallocErr();


	for(Index=0;Index<NoLong;Index++)
		OList[Index]  = RandUSLong(RS);
/*		Do not use as doubles are set to be between 0 - 1 so are non-random, can be used to validate diehard*/
/*		OList[Index]  = (unsigned long)RandDouble(RS); */

	fwrite(OList, sizeof(unsigned long), NoLong, Out);

	free(OList);
	fclose(Out);
}

int	RandPoisson(RANDSTATES* RS, double ExpectedValue)
{
    double Limit;
    double Product;
    int Count;

	Limit = exp(-ExpectedValue);
	Product = RandDouble(RS);

    for( Count = 0; Product > Limit; Count++ )
        Product *= RandDouble(RS);

    return Count;
}


double NegExp(RANDSTATES* RS, double ExpectedValue)
{
    return - ExpectedValue * log(RandDouble(RS));
}

int		CheckRest(RANDSTATES* RS, unsigned long	*CList)
{
	int i;

	for(i=1;i<PSIZE;i++)
	{
		if(CList[i] != RandUSLong(RS))
			return 0;
	}

	return 1;
}


void	FindPriod(void)
{
	RANDSTATES*		RS;
	int				i;
	unsigned long	*CList;
	unsigned long	Count;


	RS = CreateRandStates();

	CList = (unsigned long*)malloc(sizeof(unsigned long) * PSIZE);
	if(CList == NULL)
		RSMallocErr();
	for(i=0;i<PSIZE;i++)
		CList[i] = RandUSLong(RS);

	Count = 0;
	for(;;)
	{
		if(Count + 1 < Count)
		{
			printf("%lu\n", Count);
			fflush(stdout);
		}
		Count++;

		if(RandUSLong(RS) == CList[0])
		{
			if(CheckRest(RS, CList) == 1)
			{
				printf("found repeat at %lu\n", Count);
				free(CList);
				return;
			}
		}
	}
}

void	SaveStates(RANDSTATES* RS, char *FName)
{
	FILE*	Out;

	Out = fopen(FName, "wb");

	if(Out == NULL)
	{
		printf("Could not open file %s to write to. %s::%d\n", FName, __FILE__, __LINE__);
		exit(0);
	}

	fwrite(&RS->Size, sizeof(int), 1, Out);
	fwrite(&RS->Seed, sizeof(long), 1, Out);
	fwrite(&RS->PID, sizeof(long), 1, Out);
	fwrite(&RS->Pos, sizeof(int), 1, Out);
	fwrite(&RS->States[0], sizeof(long), RS->Size, Out);

	fclose(Out);
}

RANDSTATES*	LoadStates(char* FName)
{
	RANDSTATES	*Ret;
	int			Size;
	FILE		*In;
	size_t		wr;

	In = fopen(FName, "rb");

	wr = fread(&Size, sizeof(int), 1, In);
	Ret = AllocRandStates(Size);

	wr = fread(&Ret->Seed, sizeof(long), 1, In);
	wr = fread(&Ret->PID, sizeof(long), 1, In);
	wr = fread(&Ret->Pos, sizeof(int), 1, In);
	wr = fread(&Ret->States[0], sizeof(long), Ret->Size, In);

	fclose(In);

	return Ret;
}

void	DumpStates(FILE* Str, RANDSTATES *RS)
{
	int Index;

	fprintf(Str, "Size:\t%d\n", RS->Size);
	fprintf(Str, "Seed:\t%lu\n", RS->Seed);
	fprintf(Str, "PID:\t%lu\n", RS->PID);
	fprintf(Str, "Tempseed:\t%lu\n", RS->TempSeed);
	fprintf(Str, "Pos\t%d\n", RS->Pos);

	for(Index=0;Index<RS->Size;Index++)
		fprintf(Str, "\t%d\t%lu\n", Index, RS->States[Index]);
}

double RndGamma1 (RANDSTATES* RS, double  s)
{
	double	r, x=0.0, lsmall, w;
	static double a, p, uf, ss=10.0, d;

	lsmall=1e-37 ;
	if (s != ss)
		{
		a  = 1.0 - s;
		p  = a / (a+s*exp(-a));
		uf = p * pow(lsmall/a,s);
		d  = a * log(a);
		ss = s;
		}
	for (;;)
		{
		r = RandDouble(RS);
		if (r > p)
			x = a - log((1.0 - r) / (1.0 - p)), w = a * log(x) - d;
		else if (r>uf)
			x = a * pow(r/p,1/s), w = x;
		else
			return (0.0);
		r = RandDouble(RS);
		if (1.0-r <= w && r > 0.0)
		if (r*(w+1.0) >= 1.0 || -log(r) <= w)
			continue;
		break;
		}
	return (x);
}

double	RndGamma2 (RANDSTATES* RS, double s)
{
	double	r ,d, f, g, x;
	static double b, h, ss=0;

	if (s!=ss)
		{
		b  = s-1.0;
		h  = sqrt(3.0*s-0.75);
		ss = s;
		}
	for (;;)
		{
		r = RandDouble(RS);
		g = r-r*r;
		f = (r-0.5)*h/sqrt(g);
		x = b+f;
		if (x <= 0.0)
			continue;
		r = RandDouble(RS);
		d = 64*r*r*g*g*g;
		if (d*x < x-2.0*f*f || log(d) < 2*(b*log(x/b)-f))
			break;
		}
	return (x);
}

double RndGamma (RANDSTATES* RS, double s)
{

	double r=0.0;

	if (s <= 0.0)
		puts ("jgl gamma..");
	else if (s < 1.0)
		r = RndGamma1 (RS, s);
	else if (s > 1.0)
		r = RndGamma2 (RS, s);
	else
		r = -log(RandDouble(RS));
	return (r);
}

void DirichletRandomVariable (RANDSTATES* RS, double *alp, double *z, int n)
{

	int		i;
	double	sum;

	sum = 0.0;
	for(i=0; i<n; i++)
		{
		z[i] = RndGamma(RS, alp[i]) / 1.0;
		sum += z[i];
		}
	for(i=0; i<n; i++)
		z[i] /= sum;
}

/* Generate an integer in proprtion to values in List */
int RandInProportion(RANDSTATES* RS, double *List, int No)
{
	double	Sum, SSF;
	int		Index;
	double	Point;

	Sum = 0;
	for(Index=0;Index<No;Index++)
		Sum += List[Index];

	Point = RandDouble(RS);
	SSF = 0;
	for(Index=0;Index<No;Index++)
	{
		if(Point <= (SSF + List[Index]) / Sum)
			return Index;

		SSF += List[Index];
	}

	return 0;
}

int	RandIntRange(RANDSTATES* RS, int Range)
{
	unsigned long USL;

	USL = RandUSLong(RS);

	USL = USL % Range;
	return (int)USL;
}

double			RandUniDouble(RANDSTATES* RS, double Min, double Max)
{
	double Ret;

	Ret = RandDouble(RS) * (Max - Min);
	Ret += Min;
	return Ret;
}

double			RandExp(RANDSTATES* RS, double Mean)
{
	double Ret;

	Mean = 1.0 / Mean;

	Ret = log(1.0-RandDouble(RS));
	Ret = Ret / (-Mean);

	return Ret;	
}

void		ShuffleIntList(RANDSTATES* RS, int *List, size_t Size)
{
	size_t Index, Pos;
	int		Temp;

	if(Size < 2)
		return;

	for(Index=0;Index<Size-1;Index++)
	{
		Pos = RandUSInt(RS) % (Size - Index);
		Pos += Index;

		Temp = List[Pos];
		List[Pos] = List[Index];
		List[Index] = Temp;
	}
}

void		ShuffleDoubleList(RANDSTATES* RS,double *List,size_t Size)
{
	size_t Index,Pos;
	double Temp;

	if(Size < 2)
		return;

	for(Index=0;Index<Size-1;Index++)
	{
		Pos = RandUSInt(RS) % (Size - Index);
		Pos += Index;

		Temp = List[Pos];
		List[Pos] = List[Index];
		List[Index] = Temp;
	}
}

