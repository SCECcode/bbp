/***********************************************************

wcc_resamp_arbdt:   Decimate or interpolate time series
                    to new (arbitrary) sampling rate

                    copyright (c) 1999
                     Robert W. Graves
                 Woodward-Clyde Consultants
                   566 El Dorado Street
                    Pasadena, CA 91101

                    tel (626) 449-7650
                    fax (626) 449-3536

***********************************************************/

#include "include.h"
#include "structure.h"
#include "function.h"
#include "getpar.h"

#define TAP_PERC 0.05

int main(int ac,char **av)
{
struct statdata head1;
float *s1;
float *p;
int nt1, i, j;

int it, ntmax, ntap;
char dt_str[128];
float df, fac, arg;
double dnt, double_dt;
float single_dt;

int order = 4;
float nyq_perc = 1.0;
float tap_perc = TAP_PERC;

int resamp = 0;
int ntpad, ntrsmp;
float fnt, *space;
int ntout = -1;

double dt_tol = 1.0001;
double tol = 1.0e-02;

char infile[1024];
char outfile[1024];

int inbin = 0;
int outbin = 0;

int use_fftw = 1;
int use_double = 0;

sprintf(infile,"stdin");
sprintf(outfile,"stdout");

setpar(ac,av);
getpar("infile","s",infile);
getpar("outfile","s",outfile);

mstpar("newdt","f",&single_dt);
mstpar("newdt","s",&dt_str);
double_dt = 1.000000*atof(dt_str);

getpar("use_fftw","d",&use_fftw);
getpar("use_double","d",&use_double);

getpar("nyq_perc","f",&nyq_perc);
getpar("tap_perc","f",&tap_perc);
getpar("order","d",&order);
getpar("ntout","d",&ntout);
getpar("inbin","d",&inbin);
getpar("outbin","d",&outbin);
endpar();

s1 = NULL;
s1 = read_wccseis(infile,&head1,s1,inbin);

wcc_resamp_arbdt(ac, av, &s1, &head1);

write_wccseis(outfile,&head1,s1,outbin);
}
