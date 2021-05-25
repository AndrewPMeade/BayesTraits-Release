#ifndef RATES_HEADDER
#define RATES_HEADDER

RATES*		CreatRates(OPTIONS *Opt);
void		FreeRates(RATES *Rates, TREES *Trees);

void		MapMCMCConRates(RATES* Rates, OPTIONS *Opt);
void		MapRates(RATES* Rates, OPTIONS *Opt);

void		InitHMean(RATES* Rates, OPTIONS *Opt);
double		GetHMean(OPTIONS *Opt, RATES *Rates);

void		CopyRates(RATES *A, RATES *B, OPTIONS *Opt);

void		PrintRatesHeadder(FILE* Str, OPTIONS *Opt);
void		PrintRates(FILE* Str, RATES* Rates, OPTIONS *Opt, SCHEDULE* Shed);


int			FindNoOfRates(OPTIONS *Opt);
int			FindRatePos(int Rate, OPTIONS *Opt);
void		MutateRates(OPTIONS* Opt, RATES* Rates, SCHEDULE* Shed, long long It);
void		MutateHetero(RATES *Rates);



SUMMARY*	CreatSummary(OPTIONS *Opt);
void		FreeSummary(SUMMARY*	Summary);
void		PrintSummaryHeadder(FILE* Str, SUMMARY	*Summary, OPTIONS *Opt);
void		PrintSummary(FILE* Str, SUMMARY	*Summary, OPTIONS *Opt);
void		UpDataSummary(SUMMARY *Summary, RATES* Rates, OPTIONS *Opt);

//SCHEDULE*	CreatSchedule(OPTIONS *Opt);
void		PrintShed(OPTIONS* Opt, SCHEDULE* Shed, FILE* Str);
void		PrintShedHeadder(OPTIONS* Opt, SCHEDULE* Shed, FILE* Str);
void		BlankSchedule(SCHEDULE*	Shed);

char		RJModelType(int *ModelStr);

void		FindRSquared(RATES* Rates, OPTIONS *Opt, double *R2, double *SSE, double *SST);

int			FindNoEstDataPoint(TREES *Trees);

void		PrintAutoTune(FILE* Str, OPTIONS *Opt, SCHEDULE* Shed);
//void		PrintAutoTuneHeader(FILE* Str, OPTIONS *Opt);

int			FindNoConRates(OPTIONS *Opt);

void		SetEstDataFromPrior(RATES *Rates);

double*		GetEmpPis(OPTIONS *Opt);
#endif
