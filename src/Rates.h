#ifndef RATES_HEADDER
#define RATES_HEADDER

RATES*		CreatRates(OPTIONS *Opt, TREES *Trees, long Seed);
void		FreeRates(RATES *Rates, TREES *Trees);

void		MapMCMCConRates(RATES* Rates, OPTIONS *Opt, TREES *Trees);
void		MapRates(RATES* Rates, OPTIONS *Opt, TREES *Trees);

void		InitHMean(RATES* Rates, OPTIONS *Opt);
double		GetHMean(OPTIONS *Opt, RATES *Rates);

void		CopyRates(RATES *A, RATES *B, OPTIONS *Opt, TREES *Trees);


void		PrintRatesHeadder(FILE* Str, OPTIONS *Opt, TREES *Trees);
void		PrintRates(FILE* Str, RATES* Rates, OPTIONS *Opt, SCHEDULE* Shed, TREES *Trees);

double		ChangeRateExp(double Value, double Dev, gsl_rng *RNG, double *LnHastings);

int			FindNoOfRates(OPTIONS *Opt);
int			FindRatePos(int Rate, OPTIONS *Opt);
void		MutateRates(OPTIONS* Opt, RATES* Rates, TREES *Trees, SCHEDULE* Shed, size_t It);
void		MutateHetero(RATES *Rates);

SUMMARY*	CreatSummary(OPTIONS *Opt, TREES *Trees);
void		FreeSummary(SUMMARY*	Summary);
void	PrintSummaryHeadder(FILE* Str, SUMMARY	*Summary, OPTIONS *Opt, TREES *Trees);
void	PrintSummary(FILE* Str, SUMMARY	*Summary, OPTIONS *Opt, TREES *Trees);
void		UpDataSummary(SUMMARY *Summary, RATES* Rates, OPTIONS *Opt, TREES *Trees);

//SCHEDULE*	CreatSchedule(OPTIONS *Opt);
void		PrintShed(OPTIONS* Opt, SCHEDULE* Shed, FILE* Str);
void		PrintShedHeadder(OPTIONS* Opt, SCHEDULE* Shed, FILE* Str);
void		BlankSchedule(SCHEDULE*	Shed);

char		RJModelType(int *ModelStr);

void FindRSquared(RATES* Rates, OPTIONS *Opt, double *R2, double *SSE, double *SST, TREES *Trees);

int			FindNoEstDataPoint(TREES *Trees);

void		PrintAutoTune(FILE* Str, OPTIONS *Opt, SCHEDULE* Shed);
//void		PrintAutoTuneHeader(FILE* Str, OPTIONS *Opt);

int			FindNoConRates(OPTIONS *Opt, TREES *Trees);
void		SetEstDataFromPrior(RATES *Rates);

double*		GetEmpPis(OPTIONS *Opt, TREES *Trees);
int			ModelDep(MODEL Model);

#endif
