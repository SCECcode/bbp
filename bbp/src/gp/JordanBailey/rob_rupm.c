#include "include.h"
#include "structure.h"
#include "function.h"

void get_rrmpars(struct rob *rrm,int i,int j,float *xp,float *zp,float *rt,float *vs,float *rk,float *tsf)
{
float xx, zz;
int ip;

ip = i + j*(rrm->nstk);

xx = *xp - rrm->shyp;
zz = *zp - rrm->dhyp;

*rt = -1.0;
if(rrm->vrup[ip] > 0.0)
   *rt = sqrt(xx*xx + zz*zz)/rrm->vrup[ip] + rrm->tsfac[ip];

*rk = rrm->rake[ip];
*vs = rrm->slip[ip];
*tsf = rrm->tsfac[ip];
}

void rob_stf(struct rob *rrm,int i,int j,float *s,float *u,float *stf,int nt,float *dt,float *z0)
{
FILE *fpw;
int it, nstf;
int ip, it0, it1, it2;
float tr, amp, a0;
float sum;

float alpha = 0.1;      /* 1st triangle has pulse width = 2*alpha*trise */
float betadeep = 0.2;       /* 2nd triangle has amplitude = beta*A (z0>dmax)*/
float betashal = 0.5;       /* 2nd triangle has amplitude = beta*A (z0<dmin)*/
float beta, dbdd;

float dmin = 4.0;
float dmax = 6.0;

dbdd = (betadeep - betashal)/(dmax-dmin);

if((*z0) >= dmax)
   beta = betadeep;
else if((*z0) < dmax && (*z0) > dmin)
   beta = betadeep - (dmax-(*z0))*dbdd;
else
   beta = betashal;

zapit(stf,nt);

ip = i + j*(rrm->nstk);

tr = rrm->trise[ip];

alpha = alpha*tr;

it0 = (int)((alpha)/(*dt) + 0.5);
if(it0 < 2)
   it0 = 2;
it1 = (int)((tr)/(*dt) + 0.5);
if(it1 < 4)
   it1 = 4;

it2 = (2 - beta)*it0;

a0 = 1.0;
amp = a0/(float)(it0);

for(it=0;it<it0;it++)
   stf[it] = it*amp;

for(it=it0;it<it2;it++)
   stf[it] = (2*it0-it)*amp;

amp = beta*a0/(float)(it1-it2);

for(it=it2;it<it1;it++)
   stf[it] = beta*a0 + (it2-it)*amp;

nstf = nt-1;
while(stf[nstf] == (float)(0.0) && nstf)
   nstf--;

if(nstf == 0)
   {
   sum_nostf(s,u,&(rrm->slip[ip]),nt);
   return;
   }

if(nstf < nt-1)
   nstf = nstf + 2;;

sum = 0.0;
for(it=0;it<nstf;it++)
   sum = sum + (*dt)*stf[it];

if(sum <= 0.0)
   return;

/* scale STF by slip and add factor of dt to prenormalize convolution */
sum = (*dt)*(rrm->slip[ip])/sum;
for(it=0;it<nstf;it++)
   stf[it] = stf[it]*sum;

/*
if(i==1 && j==0)
   {
   fpw = fopen("stf_file","w");
   fprintf(fpw,"stf stf\n");

   fprintf(fpw,"%d %13.5e\n",nstf,(*dt));
   for(it=0;it<nstf;it++)
      fprintf(fpw,"%13.5e\n",stf[it]);
   fclose(fpw);
   }
*/

do_cnvlv(s,u,nt,stf,nstf);
}

void read_rob(struct rob *rrm,char *rfile,float *tsf)
{
FILE *fopfile(), *fpr;
float elon, elat, stk, dip, dtop, shyp, dhyp;
float xmax, xavg, sf;
int i, j, k, it;
char string[256], *sptr;

fpr = fopfile(rfile,"r");

fgets(string,256,fpr);
sscanf(string,"%f %f %d %d %f %f",&rrm->elon,
				  &rrm->elat,
				  &rrm->nstk,
                                  &rrm->ndip,
                                  &rrm->flen,
                                  &rrm->fwid);

fgets(string,256,fpr);
sscanf(string,"%f %f %f %f %f",&rrm->stk,
				  &rrm->dip,
				  &rrm->dtop,
                                  &rrm->shyp,
                                  &rrm->dhyp);

rrm->dlen = (rrm->flen)/(rrm->nstk);
rrm->dwid = (rrm->fwid)/(rrm->ndip);

rrm->as    = (float *) check_malloc ((rrm->nstk)*(rrm->ndip)*sizeof(float));
rrm->dd    = (float *) check_malloc ((rrm->nstk)*(rrm->ndip)*sizeof(float));
rrm->slip  = (float *) check_malloc ((rrm->nstk)*(rrm->ndip)*sizeof(float));
rrm->rake  = (float *) check_malloc ((rrm->nstk)*(rrm->ndip)*sizeof(float));
rrm->trise = (float *) check_malloc ((rrm->nstk)*(rrm->ndip)*sizeof(float));
rrm->vrup  = (float *) check_malloc ((rrm->nstk)*(rrm->ndip)*sizeof(float));
rrm->tsfac = (float *) check_malloc ((rrm->nstk)*(rrm->ndip)*sizeof(float));

xmax = 0.0;
xavg = 0.0;
for(j=0;j<(rrm->ndip);j++)
   {
   for(i=0;i<(rrm->nstk);i++)
      {
      k = i + j*(rrm->nstk);

      fgets(string,256,fpr);
      sscanf(string,"%f %f %f %f %f %f",&rrm->as[k],
                                        &rrm->dd[k],
                                        &rrm->slip[k],
                                        &rrm->rake[k],
                                        &rrm->trise[k],
                                        &rrm->vrup[k]);

      xavg = xavg + rrm->slip[k];
      if(rrm->slip[k] > xmax)
	 xmax = rrm->slip[k];
      }
   }
fclose(fpr);

xavg = xavg/(float)((rrm->nstk)*(rrm->ndip));
if((xmax-xavg) != (float)(0.0))
   sf = 1.0/(xmax-xavg);
else
   *tsf = sf = 0.0;

for(j=0;j<(rrm->ndip);j++)
   {
   for(i=0;i<(rrm->nstk);i++)
      {
      k = i + j*(rrm->nstk);

      rrm->tsfac[k] = sf*(rrm->slip[k]-xavg)*(*tsf);
      }
   }
}

int gen_rob_stf(struct rob *rrm,int i,int j,float *stf,int nt,float *dt,float *z0)
{
int it, nstf;
int ip, it0, it1, it2;
float tr, amp, a0;
float sum;

float alpha = 0.1;      /* 1st triangle has pulse width = 2*alpha*trise */
float betadeep = 0.2;       /* 2nd triangle has amplitude = beta*A (z0>dmax)*/
float betashal = 0.5;       /* 2nd triangle has amplitude = beta*A (z0<dmin)*/
float beta, dbdd;

float dmin = 4.0;
float dmax = 6.0;

dbdd = (betadeep - betashal)/(dmax-dmin);

if((*z0) >= dmax)
   beta = betadeep;
else if((*z0) < dmax && (*z0) > dmin)
   beta = betadeep - (dmax-(*z0))*dbdd;
else
   beta = betashal;

zapit(stf,nt);

ip = i + j*(rrm->nstk);

tr = rrm->trise[ip];

alpha = alpha*tr;

it0 = (int)((alpha)/(*dt) + 0.5);
if(it0 < 2)
   it0 = 2;
it1 = (int)((tr)/(*dt) + 0.5);
if(it1 < 4)
   it1 = 4;

it2 = (2 - beta)*it0;

a0 = 1.0;
amp = a0/(float)(it0);

for(it=0;it<it0;it++)
   stf[it] = it*amp;

for(it=it0;it<it2;it++)
   stf[it] = (2*it0-it)*amp;

amp = beta*a0/(float)(it1-it2);

for(it=it2;it<it1;it++)
   stf[it] = beta*a0 + (it2-it)*amp;

nstf = nt-1;
while(stf[nstf] == (float)(0.0) && nstf)
   nstf--;

if(nstf == 0)
   return(0);

if(nstf < nt-1)
   nstf = nstf + 2;;

sum = 0.0;
for(it=0;it<nstf;it++)
   sum = sum + (*dt)*stf[it];

if(sum <= 0.0)
   return(0);

/* scale STF by slip */
sum = (rrm->slip[ip])/sum;
for(it=0;it<nstf;it++)
   stf[it] = stf[it]*sum;

return(nstf);
}
