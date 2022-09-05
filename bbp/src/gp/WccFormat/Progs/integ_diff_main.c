#include "include.h"
#include "structure.h"
#include "function.h"
#include "getpar.h"

void t1t2(float *, int, float *, float *, float *, float *, float *, int);
void gelim_double(double *, int, double *);
void invmat3x3(float *, float *);
void invmat4x4(float *, float *);
void baseline(float *, int, float *, int);
void get_trend(float *, int, float *, float *);
void detrend(float *, int, float *, float *);
void remean(float *, int, float *, float *, float *);
void demean(float *, int, float *, float *);
void tapr(float *, int, int, int, int, int);
void retrend(float *, int, float *);
void differ(float *, int, float *, float *);
void integrate(float *, int, float *, float *);

int main (ac, av)
int ac;
char **av;
{
struct statdata shead;
float *seis;
char filein[256];
char fileout[256];
char title[128], header[128];

int inbin = 0;
int outbin = 0;

setpar(ac, av);
mstpar("filein","s",filein);
mstpar("fileout","s",fileout);
getpar("inbin","d",&inbin);
getpar("outbin","d",&outbin);
endpar();

seis = NULL;
seis = read_wccseis(filein,&shead,seis,inbin);

integ_diff(ac, av, seis, &shead);

write_wccseis(fileout,&shead,seis,outbin);
}

