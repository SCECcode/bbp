#include <errno.h>
#include <fcntl.h>
#include <math.h>
#include <stdio.h>
#include <unistd.h>
#include <stdlib.h>

#include <sys/file.h>
#ifndef __APPLE__
#include <sys/procfs.h>
#endif
#include <sys/resource.h>
#include <sys/signal.h>
#include <sys/stat.h>
#include <sys/syscall.h>
#include <sys/time.h>
#include <sys/types.h>
#include <string.h>

#include "getpar.h"

void *check_malloc(int);
void uncert(int,int *,float *,float *,int,float *,float *,float *,float *,float *,float *,float *);
FILE *fopfile(char*, char*);
char *skipval(int,char *);
char *getstr(char *,char *);
char *getflt(float *,char *);

float t95[] = {  6.3138, 2.9200, 2.3534, 2.1318, 2.0150, 1.9432, 1.8946,
                 1.8595, 1.8331, 1.8125, 1.7959, 1.7823, 1.7709, 1.7613,
                 1.7531, 1.7459, 1.7396, 1.7341, 1.7291, 1.7247, 1.7207,
                 1.7171, 1.7139, 1.7109, 1.7081, 1.7056, 1.7033, 1.7011,
                 1.6991, 1.6973, 1.6955, 1.6939, 1.6924, 1.6909, 1.6896,
                 1.6883, 1.6871, 1.6860, 1.6849, 1.6839, 1.6829, 1.6820,
                 1.6811, 1.6802, 1.6794, 1.6787, 1.6779, 1.6772, 1.6766,
                 1.6759, 1.6753, 1.6747, 1.6741, 1.6736, 1.6730 };

int main(int ac,char **av)
{
FILE *fpr, *fpw, *fopfile();
float *per, *resid, vs30, cdst, xcos, ycos, *tmin, *tmax;
float *bias, *sigma, *sigma0, *cl90m, *cl90p;
int nstat, nper, *nval, nstat_read, i;

float min_cdst = -1e+15;
float max_cdst =  1e+15;
float min_vs30 = -1e+15;
float max_vs30 =  1e+15;
float min_xcos = -1e+15;
float max_xcos =  1e+15;
float min_ycos = -1e+15;
float max_ycos =  1e+15;

char residfile[256], fileroot[256], comp[16], rdcomp[16];
char *sptr, string[2048];

setpar(ac,av);

mstpar("residfile","s",residfile);
mstpar("fileroot","s",fileroot);
mstpar("comp","s",comp);
mstpar("nstat","d",&nstat);
mstpar("nper","d",&nper);

getpar("min_cdst","f",&min_cdst);
getpar("max_cdst","f",&max_cdst);
getpar("min_vs30","f",&min_vs30);
getpar("max_vs30","f",&max_vs30);
getpar("min_xcos","f",&min_xcos);
getpar("max_xcos","f",&max_xcos);
getpar("min_ycos","f",&min_ycos);
getpar("max_ycos","f",&max_ycos);

endpar();

nval = (int *) check_malloc (nper*sizeof(float));
tmin = (float *) check_malloc (nstat*sizeof(float));
tmax = (float *) check_malloc (nstat*sizeof(float));

per = (float *) check_malloc (nper*sizeof(float));
bias = (float *) check_malloc (nper*sizeof(float));
sigma = (float *) check_malloc (nper*sizeof(float));
sigma0 = (float *) check_malloc (nper*sizeof(float));
cl90m = (float *) check_malloc (nper*sizeof(float));
cl90p = (float *) check_malloc (nper*sizeof(float));
resid = (float *) check_malloc (nstat*nper*sizeof(float));

fpr = fopfile(residfile,"r");

fgets(string,2048,fpr);
while(strncmp(string,"#",1) == 0)
   fgets(string,2048,fpr);

sptr = skipval(13,string);

for(i=0;i<nper;i++)
   {
   sptr = getflt(&per[i],sptr);
   nval[i] = 0;
   }

nstat_read = 0;
while(fgets(string,2048,fpr) != NULL)
   {
   sscanf(string,"%*s %*f %*s %*f %*f %*d %f %f %f %f %f %f %s",&vs30,&cdst,&xcos,&ycos,&tmin[nstat_read],&tmax[nstat_read],rdcomp);

   if((vs30 >= min_vs30 && vs30 <= max_vs30) &&
      (cdst >= min_cdst && cdst <= max_cdst) &&
      (xcos >= min_xcos && xcos <= max_xcos) &&
      (ycos >= min_ycos && ycos <= max_ycos) &&
      (strcmp(rdcomp,comp)==0) )
      {
      sptr = skipval(13,string);

      for(i=0;i<nper;i++)
         sptr = getflt(&resid[i+nstat_read*nper],sptr);

      nstat_read++;
      }

   if(nstat_read>nstat)
      {
      fprintf(stderr,"(nstat_read= %d) > (nstat= %d), exiting...\n",nstat_read,nstat);
      exit(-1);
      }
   }

fclose(fpr);

fprintf(stderr,"nstat_read= %d\n",nstat_read);

uncert(nstat_read,nval,tmin,tmax,nper,per,resid,bias,sigma,sigma0,cl90m,cl90p);

sprintf(string,"%s.bias",fileroot);
fpw = fopfile(string,"w");

for(i=0;i<nper;i++)
   fprintf(fpw,"%13.5e %13.5e\n",per[i],bias[i]);

fclose(fpw);

sprintf(string,"%s.sigma",fileroot);
fpw = fopfile(string,"w");

for(i=0;i<nper;i++)
   fprintf(fpw,"%13.5e %13.5e\n",per[i],sigma[i]);

fclose(fpw);

sprintf(string,"%s.sigma0",fileroot);
fpw = fopfile(string,"w");

for(i=0;i<nper;i++)
   fprintf(fpw,"%13.5e %13.5e\n",per[i],sigma0[i]);

fclose(fpw);

sprintf(string,"%s.m90",fileroot);
fpw = fopfile(string,"w");

for(i=0;i<nper;i++)
   fprintf(fpw,"%13.5e %13.5e\n",per[i],cl90m[i]);

fclose(fpw);

sprintf(string,"%s.p90",fileroot);
fpw = fopfile(string,"w");

for(i=0;i<nper;i++)
   fprintf(fpw,"%13.5e %13.5e\n",per[i],cl90p[i]);

fclose(fpw);
}

void *check_malloc(int len)
{
char *ptr;

ptr = (char *) malloc (len);

if(ptr == NULL)
   {
   fprintf(stderr,"*****  memory allocation error\n");
   exit(-1);
   }

return(ptr);
}

FILE *fopfile(char *name,char *mode)
{
FILE *fp;

if((fp = fopen(name,mode)) == NULL)
   {
   fprintf(stderr,"CAN'T FOPEN FILE = %s, MODE = %s\n", name, mode);
   exit(-1);
   }
return(fp);
}

char *skipval(int j,char *str)
{
while(j--)
   {
   while(str[0] != ' ' && str[0] != '\t' && str[0] != '\b' && str[0] != '\n')
      str++;

   while(str[0] == ' ' || str[0] == '\t' || str[0] == '\b' || str[0] == '\n')
      str++;
   }

return(str);
}

char *getstr(char *name,char *str)
{
sscanf(str,"%s",name);

while(str[0] != ' ' && str[0] != '\t' && str[0] != '\b' && str[0] != '\n')
   str++;

while(str[0] == ' ' || str[0] == '\t' || str[0] == '\b' || str[0] == '\n')
   str++;

return(str);
}

char *getflt(float *v,char *str)
{
sscanf(str,"%f",v);

while(str[0] != ' ' && str[0] != '\t' && str[0] != '\b' && str[0] != '\n')
   str++;

while(str[0] == ' ' || str[0] == '\t' || str[0] == '\b' || str[0] == '\n')
   str++;

return(str);
}

void uncert(int ns,int *nv,float *tmin,float *tmax,int np,float *per,float *res,float *b,float *sig,float *sig0,float *m90,float *p90)
{
int i, j;
float invn, invn1, ttfac;

/*

Bias is given by:
B = (1/n)*SUM(r[i])

Sigma is given by:
sigma = sqrt { 1/(n-1) SUM(r[i] - B)**2 }

NOTE:
SUM(r[i] - B)**2 = SUM(r[i]*r[i] - 2*r[i]*B + B*B)
	    = SUM(r[i]*r[i]) - 2*B*SUM(r[i]) + n*(B*B)
	    = SUM(r[i]*r[i]) - 2*B*(n*B) + n*(B*B)      where n*B = SUM(r[i])
            = SUM(r[i]*r[i]) - n*B*B

SO:
sigma = sqrt { 1/(n-1) SUM(r[i]*r[i]) - n*B*B }

*/

for(j=0;j<np;j++)
   {
   b[j] = 0.0;
   m90[j] = 0.0;
   nv[j] = 0;
   }

for(i=0;i<ns;i++)
   {
   for(j=0;j<np;j++)
      {
      if(per[j] >= tmin[i] && per[j] <= tmax[i])
         {
         b[j] = b[j] + res[j + i*np];
         m90[j] = m90[j] + res[j + i*np]*res[j + i*np];  /* temp for now */
	 nv[j] = nv[j] + 1;
	 }
      }
   }

for(j=0;j<np;j++)
   {
   if(nv[j] > 1)
      {
      invn = 1.0/(float)(nv[j]);
      invn1 = 1.0/(float)(nv[j]-1);

      if(nv[j] > 56)
         ttfac = 1.64*sqrt(invn);
      else
         ttfac = t95[nv[j]-2]*sqrt(invn);

      b[j] = b[j]*invn;
      sig[j] = sqrt(invn1*(m90[j] - nv[j]*b[j]*b[j]));  /* corrected for bias */
      sig0[j] = sqrt(invn*(m90[j]));
      m90[j] = b[j] - sig[j]*ttfac;
      p90[j] = b[j] + sig[j]*ttfac;
      }
   else
      {
      b[j] = 0.0;
      sig[j] = 0.0;
      sig0[j] = 0.0;
      m90[j] = 0.0;
      p90[j] = 0.0;
      }
   }
}
