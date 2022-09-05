#include "include.h"
#include "structure.h"
#include "function.h"
#include "getpar.h"

#define TAP_PERC 0.05

void borch_ampf(float *,float *,int,float *,float *,float *,float *,float *,float *,float *,float *,float *);
void cb2008_ampf(float *,float *,int,float *,float *,float *,float *,float *,float *,float *,float *,float *,float *,float *);
void cb2014_ampf(float *,float *,int,float *,float *,float *,float *,float *,float *,float *,float *,float *,float *,float *);
void bssa2014_ampf(float *,float *,int,float *,float *,float *,float *,float *,float *,float *,float *,float *,float *,float *);
void ampfac(struct complex *,float *,int);
void zero(float *, int);
void norm(float *, float *, int);
void taper_norm(float *, float *, int, float *);
int getnt_p2(int);
void getpeak(float *, int, float *);

int main(int ac,char **av)
{
struct statdata head1;
float *s1;

char infile[128];
char outfile[128];

int inbin = 0;
int outbin = 0;

setpar(ac,av);

mstpar("infile","s",infile);
mstpar("outfile","s",outfile);
getpar("inbin","d",&inbin);
getpar("outbin","d",&outbin);

endpar();

s1 = NULL;
s1 = read_wccseis(infile,&head1,s1,inbin);

wcc_siteamp14(ac, av, &s1, &head1);

write_wccseis(outfile,&head1,s1,outbin);
}

