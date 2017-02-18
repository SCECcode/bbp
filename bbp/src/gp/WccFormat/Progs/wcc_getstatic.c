#include "include.h"
#include "structure.h"
#include "function.h"

void integrate(float *,int,float *);

int size_float = sizeof(float);
int size_int = sizeof(int);

main(int ac, char **av)
{
struct statdata head1;
float *s1;
int i, itst, itend;
float sum;

int n_integ = 0;
float tstart = -99.0;
float tlen = -99.0;

char infile[128];
char outfile[128];

int inbin = 0;
int outbin = 0;

int print2file = 1;

setpar(ac,av);
mstpar("infile","s",infile);
getpar("print2file","d",&print2file);
if(print2file)
   mstpar("outfile","s",outfile);
getpar("n_integ","d",&n_integ);
getpar("tstart","f",&tstart);
getpar("tlen","f",&tlen);
getpar("inbin","d",&inbin);
getpar("outbin","d",&outbin);
endpar();

s1 = NULL;
s1 = read_wccseis(infile,&head1,s1,inbin);

if(n_integ) /*  integrate  */
   {
   for(i=0;i<n_integ;i++)
      integrate(s1,head1.nt,&head1.dt);
   }

itst = head1.nt - 1;
if(tstart > 0.0)
   itst = (int)(0.5 + tstart/head1.dt);

if(itst > head1.nt - 1)
   itst = head1.nt - 1;

itend = head1.nt;
if(tlen > 0.0)
   itend = itst + (int)(0.5 + tlen/head1.dt) + 1;

if(itend > head1.nt)
   itend = head1.nt;

sum = 0.0;
for(i=itst;i<itend;i++)
   sum = sum + s1[i];

sum = sum/(itend - itst);

for(i=0;i<head1.nt;i++)
   s1[i] = sum;

if(print2file)
   write_wccseis(outfile,&head1,s1,outbin);
else
   fprintf(stdout,"%13.5e\n",sum);
}

void integrate(float *s,int nt,float *dt)
{
float s0, s1;
int it;

/*
   trapezoid rule 0
*/

s0 = 0.5*(s[0])*(*dt);
for(it=1;it<nt;it++)
   {
   s1 = 0.5*(s[it-1] + s[it])*(*dt) + s0;
   s[it-1] = s0;
   s0 = s1;
   }
s[nt-1] = s1;
}
