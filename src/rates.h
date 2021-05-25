#if !defined (RATESHEADDER)
#define RATESHEADDER

RATES*		CreatRates(OPTIONS *Opt);
void		MapRates(RATES* Rates, OPTIONS *Opt);

void		PrintRatesHeadder(FILE* Str, OPTIONS *Opt);
void		PrintRates(FILE* Str, RATES* Rates, OPTIONS *Opt);
void		CopyRates(RATES *A, RATES *B, OPTIONS *Opt);

int			FindNoOfRates(OPTIONS *Opt);
void		MutateRates(OPTIONS* Opt, RATES* Rates, SCHEDULE*	Shed);
void		FreeRates(RATES *Rates);


SUMMARY*	CreatSummary(OPTIONS *Opt);
void		FreeSummary(SUMMARY*	Summary);
void		PrintSummaryHeadder(FILE* Str, SUMMARY	*Summary, OPTIONS *Opt);
void		PrintSummary(FILE* Str, SUMMARY	*Summary, OPTIONS *Opt);
void		UpDataSummary(SUMMARY *Summary, RATES* Rates, OPTIONS *Opt);

SCHEDULE*	CreatSchedule(OPTIONS *Opt);
void		PrintShed(OPTIONS* Opt, SCHEDULE* Shed, FILE* Str);
void		PrintShedHeadder(OPTIONS* Opt, SCHEDULE* Shed, FILE* Str);
void		BlankSchedule(SCHEDULE*	Shed);

char		RJModelType(int *ModelStr);

double	FindRSquared(RATES* Rates, OPTIONS *Opt);
#endif


