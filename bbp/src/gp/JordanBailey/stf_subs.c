#include "include.h"
#include "structure.h"
#include "function.h"

void do_cnvlv(float *s,float *u,int nt,float *stf,int nstf)
{
float *sv, *sn, *se, *uv, *un, *ue;
int it, k, kend;

sv = s;
sn = s + nt;
se = s + 2*nt;

uv = u;
un = u + nt;
ue = u + 2*nt;

for(it=0;it<nt;it++)
   {
   kend = it + 1;
   if(kend > nstf)
      kend = nstf;

   for(k=0;k<kend;k++)
      {
      sv[it] = sv[it] + stf[k]*uv[it-k];
      //if (k==5 && it%100==0) {
	//printf("it=%d, k=%d, sv[it]: %f, stf[k]: %f, uv[it-k]: %f\n", it, k, sv[it], stf[k], uv[it-k]);
      //}
      sn[it] = sn[it] + stf[k]*un[it-k];
      se[it] = se[it] + stf[k]*ue[it-k];
      }
   }
}

void sum_nostf(float *s,float *u,float *slip,int nt)
{
float *sv, *sn, *se, *uv, *un, *ue;
int it;

sv = s;
sn = s + nt;
se = s + 2*nt;

uv = u;
un = u + nt;
ue = u + 2*nt;

for(it=0;it<nt;it++)
   {
   sv[it] = sv[it] + (*slip)*uv[it];
   sn[it] = sn[it] + (*slip)*un[it];
   se[it] = se[it] + (*slip)*ue[it];
   }
}
