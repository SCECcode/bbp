#include "include.h"
#include "structure.h"
#include "function.h"

int size_float = sizeof(float);
int size_int = sizeof(int);

main(int ac,char **av)
{
struct statdata head0;
float *s0;
float t0, t1, amax, time, tmax;
int it;

char infile[128];

float pstart = 0.05;
float pend = 0.95;
int inbin = 0;

sprintf(infile,"stdin");

setpar(ac, av);
getpar("infile","s",infile);
getpar("pstart","f",&pstart);
getpar("pend","f",&pend);
getpar("inbin","d",&inbin);
endpar();

s0 = NULL;
s0 = read_wccseis(infile,&head0,s0,inbin);

s0[0] = head0.dt*sqrt(s0[0]*s0[0]);
amax = 0.0;
tmax = 0.0;
for(it=1;it<head0.nt;it++)
   {
   if(s0[it] > amax)
      {
      amax = s0[it];
      tmax = it*head0.dt;
      }
   if(-s0[it] > amax)
      {
      amax = -s0[it];
      tmax = it*head0.dt;
      }

   s0[it] = s0[it-1] + head0.dt*sqrt(s0[it]*s0[it]);
   }

amax = s0[head0.nt-1];
for(it=0;it<head0.nt;it++)
   s0[it] = s0[it]/amax;

t0 = -99.9;
t1 = -99.9;
for(it=0;it<head0.nt;it++)
   {
   time = it*head0.dt;

   if(s0[it] > pstart && t0 < 0.0)
      t0 = time;

   if(s0[it] > pend && t1 < 0.0)
      t1 = time;
   }

printf("%10.2f %10.2f %10.2f %10.2f\n",t0,t1,(t1-t0)/(pend-pstart),tmax - t0 + pstart*(t1-t0)/(pend-pstart));
}
