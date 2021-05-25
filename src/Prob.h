#ifndef PROB_H
#define PROB_H

double		PDFNorm(double X, double Mean, double Var);
double		PDFExp(double X, double Mean);
double		PDFGamma(double X, double Shape, double Scale);
double		PDFBeta(double X, double Alpha, double Beta);
double		PDFChi(double X, double Alpha, double Beta);
double		PDFInvGamma(double X, double Alpha, double Beta);
double		PDFSGamma(double x, double Alpha, double Beta);

double		CDFNorm(double X, double Mean, double Var);
double		CDFExp(double X, double Alpha);
double		CDFGamma(double X, double Shape, double Scale);
double		CDFBeta(double X, double Alpha, double Beta);
double		CDFChi(double X, double Alpha, double Beta);
double		CDFInvGamma(double X, double Alpha, double Beta);

//void		ProbTest(void);

#endif