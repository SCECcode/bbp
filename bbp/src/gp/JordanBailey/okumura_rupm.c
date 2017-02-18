#include "include.h"
#include "structure.h"
#include "function.h"

void get_ormpars(struct okumura *orm,int i,int j,float *xp,float *zp,float *rt,float *vs)
{
float dx, dz, dxdz;
float w0, w1, w2, w3;
int ip0, ip1, ip2, ip3;
float one = 1.0;

ip0 = i + j*(orm->nstk);

*rt = orm->rupt[ip0];
*vs = orm->slip[ip0];
}

void okumura_stf(struct okumura *orm,int i,int j,float *s,float *u,float *stf,int nt,float *dt)
{
FILE *fpw;
int it, nstf;
int ip0, it0, it1;
float td, tb, tr, ts, alpha, beta, gamma, ep, tt;
float sum, da, c0;

float quart = 0.25;
float half = 0.5;
float fone = 1.0;
float ftwo = 2.0;
float twop5 = 2.5;
float fthree = 3.0;
float p95 = 0.95;
float p125 = 1.25;
float p150 = 1.5;

float sp, vs;
float amax = 0.0;

zapit(stf,nt);

ip0 = i + j*(orm->nstk);

sp = orm->slip[ip0];
vs = orm->sv[ip0];

td = 0.1;
tb = p125*td;
tr = orm->rist[ip0];
ts = p150*tr;

if(tr <= (*dt))  /* no STF needed */
   {
   sum_nostf(s,u,&sp,nt);
   return;
   }

if(tb >= tr)
   {
   tb = p95*tr;
   td = tb/p125;
   }

alpha = ftwo/td;
gamma = half/td;
ep = (twop5*tb - fthree*td)/(fone - td/tb);
beta = alpha*tb*(fone - gamma*tb)*sqrt(tb - ep);
c0 = beta/sqrt(tr - ep);
da = -c0/(ts - tr);

it0 = (int)(tb/(*dt) + 1.0);
it1 = (int)(tr/(*dt) + 1.0);
nstf = (int)(ts/(*dt) + 1.0);

/*
sum = 0.0;
for(it=0;it<it0;it++)
   {
   tt = it*(*dt);
   stf[it] = alpha*tt*(fone - gamma*tt);

   sum = sum + (*dt)*stf[it];
   if(stf[it] > amax)
      amax = stf[it];
   }

for(it=it0;it<it1;it++)
   {
   tt = it*(*dt);
   stf[it] = beta/sqrt(tt - ep);

   sum = sum + (*dt)*stf[it];
   if(stf[it] > amax)
      amax = stf[it];
   }

for(it=it1;it<nstf;it++)
   {
   tt = it*(*dt);
   stf[it] = c0 + (tt - tr)*da;

   sum = sum + (*dt)*stf[it];
   if(stf[it] > amax)
      amax = stf[it];
   }
*/

sum = 0.0;
for(it=0;it<nstf;it++)
   {
   tt = it*(*dt);

   if(tt <= tb)
      stf[it] = alpha*tt*(fone - gamma*tt);
   else if(tt <= tr)
      stf[it] = beta/sqrt(tt - ep);
   else if(tt < ts)
      stf[it] = c0 + (tt - tr)*da;

   sum = sum + (*dt)*stf[it];
   if(stf[it] > amax)
      amax = stf[it];
   }

if(sum <= 0.0)
   return;

/*
*/
fprintf(stderr,"D= %13.5e Dvs= %13.5e r= %13.5f\n",sp,100*vs*sum/amax,(sp-100*vs*sum/amax)/sp);

/*
if(i==0 && j==0)
   {
   fpw = fopen("stf_file","w");
   fprintf(fpw,"okumura stf %d %d %13.5e %13.5e %13.5e %13.5e %13.5e\n",it0,it1,vs,td,tb,tr,ts);
   fprintf(fpw,"%d %13.5e\n",nstf,(*dt));
   for(it=0;it<nstf;it++)
      fprintf(fpw,"%13.5e\n",stf[it]);
   fclose(fpw);
   }
*/

/* scale STF by slip and add factor of dt to prenormalize convolution */
sum = (*dt)*sp/sum;
for(it=0;it<nstf;it++)
   stf[it] = stf[it]*sum;

do_cnvlv(s,u,nt,stf,nstf);
}

void read_okumura(struct okumura *orm,char *rfile,float *len2)
{
FILE *fopfile(), *fpr;
int i, j, k, kp;
char string[256];

fpr = fopfile(rfile,"r");

fgets(string,256,fpr);
sscanf(string,"%d %d %f %f %f %f %f",&orm->nstk,
                                     &orm->ndip,
                                     &orm->dlen,
                                     &orm->dwid,
                                     &orm->shypo,
                                     &orm->dhypo,
                                     &orm->vrup);

orm->flen = (orm->nstk)*(orm->dlen);
orm->fwid = (orm->ndip)*(orm->dwid);

*len2 = 0.5*(orm->flen);

orm->shypo = orm->shypo - (*len2);

orm->as   = (float *) check_malloc ((orm->nstk)*(orm->ndip)*sizeof(float));
orm->dd   = (float *) check_malloc ((orm->nstk)*(orm->ndip)*sizeof(float));
orm->slip = (float *) check_malloc ((orm->nstk)*(orm->ndip)*sizeof(float));
orm->sv   = (float *) check_malloc ((orm->nstk)*(orm->ndip)*sizeof(float));
orm->rist = (float *) check_malloc ((orm->nstk)*(orm->ndip)*sizeof(float));
orm->rupt = (float *) check_malloc ((orm->nstk)*(orm->ndip)*sizeof(float));

for(j=0;j<(orm->ndip);j++)
   {
   for(i=0;i<(orm->nstk);i++)
      {
      k = i + j*(orm->nstk);

      fgets(string,256,fpr);
      sscanf(string,"%f %f %f %f %f %f",&orm->as[k],
                                              &orm->dd[k],
                                              &orm->slip[k],
                                              &orm->sv[k],
                                              &orm->rist[k],
                                              &orm->rupt[k]);

      orm->as[k] = orm->as[k] - (*len2); /* move origin to top center */
      if(orm->dd[k] < 0.0)               /* make positive down-dip */
	 orm->dd[k] = -orm->dd[k];

      }
   }
fclose(fpr);
}
