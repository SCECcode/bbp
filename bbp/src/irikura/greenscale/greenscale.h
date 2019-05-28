#include <stdio.h>
#include <math.h>
#include <complex>
#include <vector>
#include <stdlib.h>
#include <string.h>

struct FLT{
		double LengthElement;
		double WidthElement;
		double Ds;
		double M0s;
		double myu;
		double Alpha_p;
		double AomegaCs;
		double Beta_s;
		double Sigma_s;
		double EqMoment;
		double lamda;
		double omegaCs;
		double Rs;
		double omegaSpq;
		double omegaDpq;
		double Rpq;
		double Dpq;
 	    double Ss;
		float  vr;                // ”j‰ó“`”d‘¬“x(m)
		double Qdepend1;          // •Ï‰»Q’l
		double Qdepend2;          // Q’l‚Ìü”g”•Ï‰»‚Ìæ”€
		double Qdepend3;          // Q’l‚Ì•Ï‰»‹«ŠEü”g”
		double Qdepend4;          // Q’l‚»‚Ì‘¼
        double Qp;
		double offsettime;         // ”j‰óŠÔ‚ÌƒIƒtƒZƒbƒgŠÔ(sec)
        double fmax;

};


std::complex<double> ScaleFactorS( double omega, double Qvalue, double EDpq, double EomegaDpq, double Rpq, double Rs, double Sigmapq,struct FLT mf ,struct FLT subfault );
std::complex<double> ScaleFactorP( double omega, double Qp, double EDpq, double EomegaDpq, double Rpq, double Rs, double Sigmapq,struct FLT mf ,struct FLT subfault );
void init( FLT *mainf, FLT *sbf );
void MakeFaultPara( FLT *subf );
void CFastFourierTransform( std::vector< std::complex<double> >& cx, int ind );
void GetGousei( std::vector<double> *sum, std::vector<double> *wv, int m, int vrshift);
int seismic_intensity(int nt, double dt, double xr[], double& eji);
void wvfast(int n, double* ar, double* ai, int ind);
double a0_func(double crit, double *x, int nt, double dt, double offset);
double rtbis(double (*func)(double, double *, int, double, double), double x1, double x2, double xacc, double *xx, int nt, double dt, double offset);
void Distance(double dLonA, double dLatA,double dHgtA, double dLonB,double dLatB,double dHgtB,double &dRet);
void SetXYZ(double dLon,double dLat,double dHgt,double& dX, double& dY, double& dZ );
void Ssekibun(std::vector<double> *indata, std::vector<double> *outdata, int n, double s_time, double term, double gain);
void tsekibun( std::complex<double> sekibunc, double dt, int nn, double v0 );
void GetGousei( std::vector<double> *sum, std::vector<double> *wv, int m , int vrshift);
