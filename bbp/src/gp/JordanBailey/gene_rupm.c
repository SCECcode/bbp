#include "include.h"
#include "structure.h"
#include "function.h"

void get_grmpars(struct gene *grm,int i,int j,float *xp,float *zp,float *rt,float *vs,float *rk)
{
float xx, zz;
int ip;

ip = i + j*(grm->nstk);

xx = *xp - grm->shypo;
zz = *zp - grm->dhypo;

*rt = sqrt(xx*xx + zz*zz)/grm->vrup[ip];
*rk = grm->rake[ip];
*vs = grm->slip[ip];
}

void gene_stf(struct gene *grm,int i,int j,float *s,float *u,float *stf,int nt,float *dt)
{
FILE *fpw;
int it, k, ntri, nstf;
int ip, itdel;
float t2, tt, amp;
float sum;

float half = 0.5;
float one = 1.0;

zapit(stf,nt);

ip = i + j*(grm->nstk);

t2 = 0.5*grm->trise;
ntri = (int)(grm->trise/(*dt) + one);

for(k=0;k<grm->nt[ip];k++)
   {
   itdel = (int)((k*grm->tdel)/(*dt) + half);
   amp = grm->swgt[k + grm->maxtw*ip];

   for(it=0;it<ntri;it++)
      {
      tt = it*(*dt);

      if(tt <= t2)
         stf[it+itdel] = stf[it+itdel] + amp*tt;
      else if(tt <= grm->trise)
         stf[it+itdel] = stf[it+itdel] + amp*(grm->trise - tt);
      }
   }

nstf = nt-1;
while(stf[nstf] == (float)(0.0) && nstf)
   nstf--;

if(nstf == 0)
   {
   sum_nostf(s,u,&(grm->slip[ip]),nt);
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
sum = (*dt)*(grm->slip[ip])/sum;
for(it=0;it<nstf;it++)
   stf[it] = stf[it]*sum;

/*
if(i==0 && j==0)
   {
   fpw = fopen("stf_file","w");
   fprintf(fpw,"stf stf");

   for(k=0;k<grm->nt[ip];k++)
      fprintf(fpw," %4.0f",grm->swgt[k + grm->maxtw*ip]);
   fprintf(fpw,"\n");

   fprintf(fpw,"%d %13.5e\n",nstf,(*dt));
   for(it=0;it<nstf;it++)
      fprintf(fpw,"%13.5e\n",stf[it]);
   fclose(fpw);
   }
*/

do_cnvlv(s,u,nt,stf,nstf);
}

void read_gene(struct gene *grm,char *rfile,float *len2)
{
FILE *fopfile(), *fpr;
int i, j, k, it;
char string[256], *sptr;

fpr = fopfile(rfile,"r");

fgets(string,256,fpr);
sscanf(string,"%d %d %f %f %f %f %f %f %d",&grm->nstk,
                                           &grm->ndip,
                                           &grm->flen,
                                           &grm->fwid,
                                           &grm->shypo,
                                           &grm->dhypo,
                                           &grm->trise,
                                           &grm->tdel,
                                           &grm->maxtw);

grm->dlen = (grm->flen)/(grm->nstk);
grm->dwid = (grm->fwid)/(grm->ndip);

*len2 = 0.5*(grm->flen);

grm->as   = (float *) check_malloc ((grm->nstk)*(grm->ndip)*sizeof(float));
grm->dd   = (float *) check_malloc ((grm->nstk)*(grm->ndip)*sizeof(float));
grm->slip = (float *) check_malloc ((grm->nstk)*(grm->ndip)*sizeof(float));
grm->rake = (float *) check_malloc ((grm->nstk)*(grm->ndip)*sizeof(float));
grm->vrup = (float *) check_malloc ((grm->nstk)*(grm->ndip)*sizeof(float));
grm->nt   = (int *) check_malloc ((grm->nstk)*(grm->ndip)*sizeof(int));
grm->swgt = (float *) check_malloc ((grm->maxtw)*(grm->nstk)*(grm->ndip)*sizeof(float));

for(j=0;j<(grm->maxtw)*(grm->nstk)*(grm->ndip);j++)
   grm->swgt[j] = 0.0;

for(j=0;j<(grm->ndip);j++)
   {
   for(i=0;i<(grm->nstk);i++)
      {
      k = i + j*(grm->nstk);

      fgets(string,256,fpr);
      sscanf(string,"%f %f %f %f %f %d",&grm->as[k],
                                        &grm->dd[k],
                                        &grm->slip[k],
                                        &grm->rake[k],
                                        &grm->vrup[k],
                                        &grm->nt[k]);

      sptr = string;
      sptr = skipval(6,sptr);

      for(it=0;it<grm->nt[k];it++)
	 {
         sscanf(sptr,"%f",&grm->swgt[it + (grm->maxtw)*k]);
         sptr = skipval(1,sptr);
	 }
      }
   }
fclose(fpr);
}
