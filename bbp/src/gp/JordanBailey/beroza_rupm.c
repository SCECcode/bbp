#include "include.h"
#include "structure.h"
#include "function.h"

void get_brmpars(struct beroza *brm,int i,int j,float *xp,float *zp,float *rt,float *vs)
{
float dx, dz, dxdz;
float w0, w1, w2, w3;
int ip0, ip1, ip2, ip3;
float one = 1.0;

ip0 = i + j*(brm->npstk);
ip1 = ip0 + 1;
ip2 = ip0 + (brm->npstk);
ip3 = ip2 + 1;;

dx = (*xp - brm->as[ip0])/(brm->as[ip1] - brm->as[ip0]);
dz = (*zp - brm->dd[ip0])/(brm->dd[ip2] - brm->dd[ip0]);
dxdz = dx*dz;

w0 = one - dx - dz + dxdz;
w1 = dx - dxdz;
w2 = dz - dxdz;
w3 = dxdz;

*rt = w0*brm->rupt[ip0] + w1*brm->rupt[ip1] + w2*brm->rupt[ip2] + w3*brm->rupt[ip3];
*vs = w0*brm->slip[ip0] + w1*brm->slip[ip1] + w2*brm->slip[ip2] + w3*brm->slip[ip3];
}

void beroza_stf(struct beroza *brm,int i,int j,float *s,float *u,float *stf,int nt,float *dt,float *z0)
{
FILE *fpw;
int it, k, kend, nstf;
int ip0, ip1, ip2, ip3;
float tp, tr, t2, invtp, tp2, tt;
float sum, *sv, *sn, *se, c0, c1;
float *uv, *un, *ue;

float quart = 0.25;
float half = 0.5;
float one = 1.0;
float two = 2.0;
float p95 = 0.95;
float p105 = 1.05;

float sp, vs;
float amax = 0.0;

float alpha = 0.1;      /* 1st triangle has pulse width = 2*alpha*trise */
float betadeep = 0.2;       /* 2nd triangle has amplitude = beta*A (z0>dmax)*/
float betashal = 0.5;       /* 2nd triangle has amplitude = beta*A (z0<dmin)*/
float beta, dbdd;

float dmin = 4.0;
float dmax = 6.0;

zapit(stf,nt);

ip0 = i + j*(brm->npstk);
ip1 = ip0 + 1;
ip2 = ip0 + (brm->npstk);
ip3 = ip2 + 1;

sp = quart*(brm->slip[ip0] + brm->slip[ip1] + brm->slip[ip2] + brm->slip[ip3]);
vs = quart*(brm->sv[ip0]   + brm->sv[ip1]   + brm->sv[ip2]   + brm->sv[ip3]  );

tp = quart*(brm->tdur[ip0] + brm->tdur[ip1] + brm->tdur[ip2] + brm->tdur[ip3]);
tr = quart*(brm->rist[ip0] + brm->rist[ip1] + brm->rist[ip2] + brm->rist[ip3]);
t2 = quart*(brm->t2[ip0]   + brm->t2[ip1]   + brm->t2[ip2]   + brm->t2[ip3]  );

if(brm->robstf == 1)
   {
   dbdd = (betadeep - betashal)/(dmax-dmin);

   if((*z0) >= dmax)
      beta = betadeep;
   else if((*z0) < dmax && (*z0) > dmin)
      beta = betadeep - (dmax-(*z0))*dbdd;
   else
      beta = betashal;

   tp = two*alpha*tr;
   t2 = (one - half*beta)*tp;
   }

if(tp <= (*dt) || tr <= (*dt) || t2 <= (*dt))  /* no STF needed */
   {
   sum_nostf(s,u,&sp,nt);
   return;
   }

if(t2 >= tp)
   t2 = p95*tp;
if(tr <= tp)
   tr = p105*tp;

tp2 = half*tp;
invtp = one/tp;
c0 = (one - invtp*t2);
c1 = c0/(tr - t2);
nstf = (int)(tr/(*dt) + 1.0);

sum = 0.0;
for(it=0;it<nstf;it++)
   {
   tt = it*(*dt);

   if(tt < tp2)
      stf[it] = invtp*tt;
   else if(tt < t2)
      stf[it] = one - invtp*tt;
   else if(tt < tr)
      stf[it] = c0 - c1*(tt - t2);

   sum = sum + (*dt)*stf[it];
   if(stf[it] > amax)
      amax = stf[it];
   }

if(sum <= 0.0)
   return;

/* scale STF by slip and add factor of dt to prenormalize convolution */
sum = (*dt)*sp/sum;
for(it=0;it<nstf;it++)
   stf[it] = stf[it]*sum;

/*
fprintf(stderr,"sp=%13.5e %13.5e vs=%13.5e %13.5e\n",sp,vs*(t2 - 0.5*tp + tr*(1.0 - t2/tp)),vs,sp*amax*sum);
*/

/*

if(i==0 && j==0)
   {
   fpw = fopen("stf_file","w");
   fprintf(fpw,"stf stf %13.5e %13.5e %13.5e\n",tr,tp,t2);
   fprintf(fpw,"%d %13.5e\n",nstf,(*dt));
   for(it=0;it<nstf;it++)
      fprintf(fpw,"%13.5e\n",stf[it]);
   fclose(fpw);
   }
*/

do_cnvlv(s,u,nt,stf,nstf);
}

void read_beroza(struct beroza *brm,char *rfile,float *len2)
{
FILE *fopfile(), *fpr;
int i, j, k, kp;
char string[256];

if(((brm->npstk) - 1)%(brm->inc_stk) != 0)
   {
   fprintf(stderr,"***input nodes along strike are not integer multiple of desired resampling-\n");
   fprintf(stderr,"     inc_stk   =%d\n",(brm->inc_stk));
   fprintf(stderr,"     (npstk-1) =%d\n",(brm->npstk)-1);
   fprintf(stderr,"     exiting...\n");
   exit(-1);
   }

if(((brm->npdip) - 1)%(brm->inc_dip) != 0)
   {
   fprintf(stderr,"***input nodes down dip are not integer multiple of desired resampling-\n");
   fprintf(stderr,"     inc_dip   =%d\n",(brm->inc_dip));
   fprintf(stderr,"     (npdip-1) =%d\n",(brm->npdip)-1);
   fprintf(stderr,"     exiting...\n");
   exit(-1);
   }

brm->npstk = (int)(((brm->npstk) - 1)/(brm->inc_stk)) + 1;
brm->npdip = (int)(((brm->npdip) - 1)/(brm->inc_dip)) + 1;

brm->as   = (float *) check_malloc ((brm->npstk)*(brm->npdip)*sizeof(float));
brm->dd   = (float *) check_malloc ((brm->npstk)*(brm->npdip)*sizeof(float));
brm->slip = (float *) check_malloc ((brm->npstk)*(brm->npdip)*sizeof(float));
brm->sv   = (float *) check_malloc ((brm->npstk)*(brm->npdip)*sizeof(float));
brm->rupt = (float *) check_malloc ((brm->npstk)*(brm->npdip)*sizeof(float));
brm->rist = (float *) check_malloc ((brm->npstk)*(brm->npdip)*sizeof(float));
brm->tdur = (float *) check_malloc ((brm->npstk)*(brm->npdip)*sizeof(float));
brm->t2   = (float *) check_malloc ((brm->npstk)*(brm->npdip)*sizeof(float));

fpr = fopfile(rfile,"r");
for(j=0;j<(brm->npdip);j++)
   {
   for(i=0;i<(brm->npstk);i++)
      {
      k = i + j*(brm->npstk);

      fgets(string,256,fpr);
      sscanf(string,"%f %f %f %f %f %f %f %f",&brm->dd[k],
                                                 &brm->as[k],
                                                 &brm->slip[k],
                                                 &brm->sv[k],
                                                 &brm->rupt[k],
                                                 &brm->rist[k],
                                                 &brm->tdur[k],
                                                 &brm->t2[k]);

      if((brm->generic_risetime) > 0.0)
         {
         brm->rist[k] = (brm->generic_risetime);
         brm->tdur[k] = (brm->generic_pulsedur);
         brm->t2[k] = (brm->generic_t2);
         }

      brm->as[k] = brm->as[k] - (*len2);

      if(i != (brm->npstk)-1) /* skip inc_stk-1 lines */
         {
         for(kp=0;kp<(brm->inc_stk)-1;kp++)
            fgets(string,256,fpr);
         }
      }

   if(j != (brm->npdip)-1) /* skip ((brm->inc_dip)-1)*(((brm->npstk)-1)*(brm->inc_stk) + 1) lines */
      {
      for(kp=0;kp<((brm->inc_dip)-1)*(((brm->npstk)-1)*(brm->inc_stk) + 1);kp++)
         fgets(string,256,fpr);
      }
   }
fclose(fpr);
}
