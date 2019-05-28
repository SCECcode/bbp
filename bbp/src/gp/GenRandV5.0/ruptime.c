#include "include.h"
#include "structure.h"
#include "function.h"
#include "defs.h"

void get_range(vm,slay,sth,rlay,rth,p,r0,linc)
struct velmodel *vm;
double *sth, *rth, *p, *r0;
int linc, slay, rlay;
{
int i;
double denom, arg;
double invp2;

invp2 = 1.0/((*p)*(*p));

denom = sqrt(invp2*vm->invb2[slay] - 1.0);
*r0 = (*sth)/denom;

for(i=slay+linc;i!=rlay;i=i+linc)
   {
   denom = sqrt(invp2*vm->invb2[i] - 1.0);
   *r0 = *r0 + vm->th[i]/denom;
   }

denom = sqrt(invp2*vm->invb2[rlay] - 1.0);
*r0 = *r0 + (*rth)/denom;
}

void get_radtime(vm,slay,sth,rlay,rth,p,r0,tt,linc)
struct velmodel *vm;
double *sth, *rth, *p, *r0;
float *tt;
int linc, slay, rlay;
{
int i;
double r1, rad, arg;
double denom, invp2;

if(*p > 0.0)
   {
   arg = 1.0 - (*p)*(*p)*vm->vs[slay]*vm->vs[slay];
   denom = sqrt(arg);

   *r0 = (*sth)/denom;
   *tt = *r0/vm->vs[slay];

   for(i=slay+linc;i!=rlay;i=i+linc)
      {
      arg = 1.0 - (*p)*(*p)*vm->vs[i]*vm->vs[i];
      denom = sqrt(arg);

      rad = vm->th[i]/denom;
      *r0 = *r0 + rad;
      *tt = *tt + rad/vm->vs[i];
      }

   arg = 1.0 - (*p)*(*p)*vm->vs[rlay]*vm->vs[rlay];
   denom = sqrt(arg);

   rad = (*rth)/denom;
   *r0 = *r0 + rad;
   *tt = *tt + rad/vm->vs[rlay];
   }
else
   {
   *r0 = *sth;
   *tt = *r0/vm->vs[slay];

   for(i=slay+linc;i!=rlay;i=i+linc)
      {
      *r0 = *r0 + vm->th[i];
      *tt = *tt + vm->th[i]/vm->vs[i];
      }

   *r0 = *r0 + *rth;
   *tt = *tt + (*rth)/vm->vs[rlay];
   }
}

void get_headtime(mod,slay,sth,rlay,rth,rad,tt)
struct velmodel *mod;
double *sth, *rth, *rad;
float *tt;
int slay, rlay;
{
int vflag, j, jj, jst, jnd;
double inv2, rc, tinc, arg;

jst = rlay;
if(slay > rlay)
   jst = slay;

*tt = 1.0e+5;
for(jnd=jst+1;jnd<mod->nlay;jnd++)
   {
   jj = rlay;
   if(slay < rlay)
      jj = slay;

   vflag = 1;
   for(j=jj;j<jnd;j++)
      {
      if(mod->vs[j] > mod->vs[jnd])
	 vflag = -1;
      }

   if(vflag == 1)
      {
      tinc = (*rad)/mod->vs[jnd];
      inv2 = 1.0/(mod->vs[jnd]*mod->vs[jnd]);

      arg = 1.0/(mod->vs[slay]*mod->vs[slay]) - inv2;
      arg = sqrt(arg);
      tinc = tinc + (*sth)*arg;
      rc = (*sth)/(arg*mod->vs[jnd]);

      for(j=slay+1;j<jnd;j++)
         {
         arg = 1.0/(mod->vs[j]*mod->vs[j]) - inv2;
         arg = sqrt(arg);
         tinc = tinc + mod->th[j]*arg;
         rc = rc + mod->th[j]/(arg*mod->vs[jnd]);
         }

      for(j=rlay+1;j<jnd;j++)
         {
         arg = 1.0/(mod->vs[j]*mod->vs[j]) - inv2;
         arg = sqrt(arg);
         tinc = tinc + mod->th[j]*arg;
         rc = rc + mod->th[j]/(arg*mod->vs[jnd]);
         }

      arg = 1.0/(mod->vs[rlay]*mod->vs[rlay]) - inv2;
      arg = sqrt(arg);
      tinc = tinc + (*rth)*arg;
      rc = rc + (*rth)/(arg*mod->vs[jnd]);

      if(tinc < *tt && rc < (*rad))
         *tt = tinc;
      }
   }
}

void bisect_p(vm,slay,sth,rlay,rth,p,eps,tol,rng,linc)
struct velmodel *vm;
double *sth, *rth, *p, *rng, *eps;
float *tol;
int linc, slay, rlay;
{
double tp0, pp, pm, p0, r0, delr;
int i, ic;

int nc = 100;

p0 = 1.0/vm->vs[slay];
for(i=slay+linc;i!=rlay;i=i+linc)
   {
   tp0 = 1.0/vm->vs[i];
   if(tp0 < p0)
      p0 = tp0;
   }
tp0 = 1.0/vm->vs[rlay];
if(tp0 < p0)
   p0 = tp0;

*p = *eps;

get_range(vm,slay,sth,rlay,rth,p,&r0,linc);

if(r0 < *rng)  /* if not, then p=0 (vertical ray) */
   {
                   /* bracket range with ray parameter extremes */

   ic = 0;
   while(r0 < *rng && ic < nc)
      {
      ic++;
      *p = p0*(1.0 - (*eps)/(double)(ic*ic));
      get_range(vm,slay,sth,rlay,rth,p,&r0,linc);
      }

   pp = *p;
   pm = *eps;

   delr = r0 - *rng;

/*
   use bisection to converge to correct ray parameter
*/

   ic = 0;
   while(delr > *tol)
      {
      *p = 0.5*(pp + pm);

      if(*p == pp || *p == pm) /* beyond double precision accuracy */
         break;

      get_range(vm,slay,sth,rlay,rth,p,&r0,linc);
      if(r0 >= *rng)
         {
         delr = r0 - *rng;
         pp = *p;
         }
      else
         {
         delr = *rng - r0;
         pm = *p;
         }

      ic++;
      if(ic > nc)
         break;
      }
   }
else
   *p = 0.0;
}

void default_velmodel(struct velmodel *vm)
{
int i;

vm->nlay = 17;

vm->vp = (float *)check_malloc(vm->nlay*sizeof(float));
vm->vs = (double *)check_malloc(vm->nlay*sizeof(double));
vm->den = (float *)check_malloc(vm->nlay*sizeof(float));
vm->th = (float *)check_malloc(vm->nlay*sizeof(float));
vm->dep = (float *)check_malloc(vm->nlay*sizeof(float));
vm->mu = (float *)check_malloc(vm->nlay*sizeof(float));
vm->invb2 = (double *)check_malloc(vm->nlay*sizeof(double));

vm->th[ 0] = 0.002; vm->vp[ 0] = 1.70; vm->vs[ 0] = 0.35; vm->den[ 0] = 2.0;
vm->th[ 1] = 0.004; vm->vp[ 1] = 1.80; vm->vs[ 1] = 0.55; vm->den[ 1] = 2.1;
vm->th[ 2] = 0.006; vm->vp[ 2] = 1.80; vm->vs[ 2] = 0.80; vm->den[ 2] = 2.1;
vm->th[ 3] = 0.008; vm->vp[ 3] = 1.90; vm->vs[ 3] = 0.90; vm->den[ 3] = 2.1;
vm->th[ 4] = 0.010; vm->vp[ 4] = 2.00; vm->vs[ 4] = 1.00; vm->den[ 4] = 2.2;
vm->th[ 5] = 0.070; vm->vp[ 5] = 2.40; vm->vs[ 5] = 1.30; vm->den[ 5] = 2.2;
vm->th[ 6] = 0.200; vm->vp[ 6] = 3.30; vm->vs[ 6] = 1.80; vm->den[ 6] = 2.3;
vm->th[ 7] = 0.200; vm->vp[ 7] = 3.90; vm->vs[ 7] = 2.20; vm->den[ 7] = 2.3;
vm->th[ 8] = 0.200; vm->vp[ 8] = 4.15; vm->vs[ 8] = 2.40; vm->den[ 8] = 2.4;
vm->th[ 9] = 0.300; vm->vp[ 9] = 4.50; vm->vs[ 9] = 2.60; vm->den[ 9] = 2.4;
vm->th[10] = 2.000; vm->vp[10] = 5.00; vm->vs[10] = 2.90; vm->den[10] = 2.5;
vm->th[11] = 3.000; vm->vp[11] = 5.70; vm->vs[11] = 3.30; vm->den[11] = 2.5;
vm->th[12] = 5.000; vm->vp[12] = 6.10; vm->vs[12] = 3.50; vm->den[12] = 2.6;
vm->th[13] = 5.000; vm->vp[13] = 6.25; vm->vs[13] = 3.60; vm->den[13] = 2.7;
vm->th[14] = 5.000; vm->vp[14] = 6.45; vm->vs[14] = 3.70; vm->den[14] = 2.7;
vm->th[15] = 11.00; vm->vp[15] = 6.60; vm->vs[15] = 3.80; vm->den[15] = 2.8;
vm->th[16] = 999.0; vm->vp[16] = 7.80; vm->vs[16] = 4.50; vm->den[16] = 3.2;

for(i=0;i<vm->nlay;i++)
   {
   if(i==0)
      vm->dep[i] = vm->th[i];
   else
      vm->dep[i] = vm->dep[i-1] + vm->th[i];

   vm->mu[i] = vm->vs[i]*vm->vs[i]*vm->den[i]*1.0e+10;  /* in CMS units */
   }
}

void read_velmodel(char *vfile,struct velmodel *vm)
{
FILE *fpr, *fopfile();
int i;
char str[512];

fpr = fopfile(vfile,"r");

fgets(str,512,fpr);
sscanf(str,"%d",&vm->nlay);

vm->vp = (float *)check_malloc(vm->nlay*sizeof(float));
vm->vs = (double *)check_malloc(vm->nlay*sizeof(double));
vm->den = (float *)check_malloc(vm->nlay*sizeof(float));
vm->th = (float *)check_malloc(vm->nlay*sizeof(float));
vm->dep = (float *)check_malloc(vm->nlay*sizeof(float));
vm->mu = (float *)check_malloc(vm->nlay*sizeof(float));
vm->invb2 = (double *)check_malloc(vm->nlay*sizeof(double));

for(i=0;i<vm->nlay;i++)
   {
   fgets(str,512,fpr);
   sscanf(str,"%f %f %lf %f",&vm->th[i],&vm->vp[i],&vm->vs[i],&vm->den[i]);

   if((i+1) == vm->nlay)
      vm->th[i] = 10000.0;

   if(i==0)
      vm->dep[i] = vm->th[i];
   else
      vm->dep[i] = vm->dep[i-1] + vm->th[i];

   vm->mu[i] = vm->vs[i]*vm->vs[i]*vm->den[i]*1.0e+10;  /* in CMS units */
   }
fclose(fpr);
}

void read_Fvelmodel(char *vfile,struct velmodel *vm)
{
FILE *fpr, *fopfile();
int i, nr;
char str[512];

int nblock = 50;

if(strcmp(vfile,"NOT_PROVIDED") == 0)
   {
   vm->nlay = 1;
   vm->vp = (float *)check_malloc(vm->nlay*sizeof(float));
   vm->vs = (double *)check_malloc(vm->nlay*sizeof(double));
   vm->den = (float *)check_malloc(vm->nlay*sizeof(float));
   vm->th = (float *)check_malloc(vm->nlay*sizeof(float));
   vm->dep = (float *)check_malloc(vm->nlay*sizeof(float));
   vm->mu = (float *)check_malloc(vm->nlay*sizeof(float));
   vm->invb2 = (double *)check_malloc(vm->nlay*sizeof(double));

   vm->th[0] = 9999.9;
   vm->vp[0] = 6.23;
   vm->vs[0] = 3.60;
   vm->den[0] = 2.80;

   vm->dep[0] = vm->th[0];
   vm->mu[0] = vm->vs[0]*vm->vs[0]*vm->den[0]*1.0e+10;  /* in CMS units */
   }
else
   {
   fpr = fopfile(vfile,"r");

   fgets(str,512,fpr);

   if(str[0] == '#')
      {
      while(str[0] == '#')
         fgets(str,512,fpr);

      vm->nlay = nblock;
      vm->vp = (float *)check_malloc(vm->nlay*sizeof(float));
      vm->vs = (double *)check_malloc(vm->nlay*sizeof(double));
      vm->den = (float *)check_malloc(vm->nlay*sizeof(float));
      vm->th = (float *)check_malloc(vm->nlay*sizeof(float));
      vm->dep = (float *)check_malloc(vm->nlay*sizeof(float));
      vm->mu = (float *)check_malloc(vm->nlay*sizeof(float));
      vm->invb2 = (double *)check_malloc(vm->nlay*sizeof(double));

      i = 0;
      sscanf(str,"%f %f %lf %f",&vm->th[i],&vm->vp[i],&vm->vs[i],&vm->den[i]);

      vm->dep[i] = vm->th[i];
      vm->mu[i] = vm->vs[i]*vm->vs[i]*vm->den[i]*1.0e+10;  /* in CMS units */

      while(fgets(str,512,fpr) != NULL)
         {
         i++;

         if(i == vm->nlay)
            {
	    vm->nlay = vm->nlay + nblock;
	    vm->vp = (float *)check_realloc(vm->vp,vm->nlay*sizeof(float));
	    vm->vs = (double *)check_realloc(vm->vs,vm->nlay*sizeof(double));
	    vm->den = (float *)check_realloc(vm->den,vm->nlay*sizeof(float));
	    vm->th = (float *)check_realloc(vm->th,vm->nlay*sizeof(float));
	    vm->dep = (float *)check_realloc(vm->dep,vm->nlay*sizeof(float));
	    vm->mu = (float *)check_realloc(vm->mu,vm->nlay*sizeof(float));
	    vm->invb2 = (double *)check_realloc(vm->invb2,vm->nlay*sizeof(double));
	    }

         sscanf(str,"%f %f %lf %f",&vm->th[i],&vm->vp[i],&vm->vs[i],&vm->den[i]);

         vm->dep[i] = vm->dep[i-1] + vm->th[i];
         vm->mu[i] = vm->vs[i]*vm->vs[i]*vm->den[i]*1.0e+10;  /* in CMS units */
         }

      vm->nlay = i+1;
      vm->vp = (float *)check_realloc(vm->vp,vm->nlay*sizeof(float));
      vm->vs = (double *)check_realloc(vm->vs,vm->nlay*sizeof(double));
      vm->den = (float *)check_realloc(vm->den,vm->nlay*sizeof(float));
      vm->th = (float *)check_realloc(vm->th,vm->nlay*sizeof(float));
      vm->dep = (float *)check_realloc(vm->dep,vm->nlay*sizeof(float));
      vm->mu = (float *)check_realloc(vm->mu,vm->nlay*sizeof(float));
      vm->invb2 = (double *)check_realloc(vm->invb2,vm->nlay*sizeof(double));
      }
   else
      {
      sscanf(str,"%d",&vm->nlay);

      vm->vp = (float *)check_malloc(vm->nlay*sizeof(float));
      vm->vs = (double *)check_malloc(vm->nlay*sizeof(double));
      vm->den = (float *)check_malloc(vm->nlay*sizeof(float));
      vm->th = (float *)check_malloc(vm->nlay*sizeof(float));
      vm->dep = (float *)check_malloc(vm->nlay*sizeof(float));
      vm->mu = (float *)check_malloc(vm->nlay*sizeof(float));
      vm->invb2 = (double *)check_malloc(vm->nlay*sizeof(double));

      for(i=0;i<vm->nlay;i++)
         {
         fgets(str,512,fpr);
         sscanf(str,"%f %f %lf %f",&vm->th[i],&vm->vp[i],&vm->vs[i],&vm->den[i]);

         if(i==0)
            vm->dep[i] = vm->th[i];
         else
            vm->dep[i] = vm->dep[i-1] + vm->th[i];

         vm->mu[i] = vm->vs[i]*vm->vs[i]*vm->den[i]*1.0e+10;  /* in CMS units */
         }
      }

   fclose(fpr);
   }
}

void conv2vrup(struct velmodel *vm,struct velmodel *rvm,float *dip,float *ztop,float *wid,float *rvf,float *shal_vr)
{
int i, j, k;
float invsinA, dep, zbot;
char string[256];

float rperd = 0.017453293;
float dmin = 4.0;
float dmax = 6.0;
float rvfac;

rvm->nlay = vm->nlay;
rvm->vs = (double *)check_malloc(rvm->nlay*sizeof(double));
rvm->th = (float *)check_malloc(rvm->nlay*sizeof(float));
rvm->invb2 = (double *)check_malloc(rvm->nlay*sizeof(double));

i = 0;
dep = vm->th[0];
while(dep < (*ztop))
   {
   i++;
   dep = dep + vm->th[i];
   }

zbot = *ztop + (*wid)*sin((*dip)*rperd);
invsinA = 1.0/sin((*dip)*rperd);

if(dep >= dmax)
   rvfac = (*rvf);
else if(dep < dmax && dep > dmin)
   rvfac = (*rvf)*(1.0 - (1.0 - (*shal_vr))*(dmax-dep)/(dmax-dmin));
else
   rvfac = (*rvf)*(*shal_vr);

rvm->th[0] = invsinA*(dep - (*ztop));
rvm->vs[0] = rvfac*vm->vs[i];

j = i;
k = 0;
while(dep < zbot)
   {
   j++; k++;
   dep = dep + vm->th[j];

   if(dep >= dmax)
      rvfac = (*rvf);
   else if(dep < dmax && dep > dmin)
      rvfac = (*rvf)*(1.0 - (1.0 - (*shal_vr))*(dmax-dep)/(dmax-dmin));
   else
      rvfac = (*rvf)*(*shal_vr);

   rvm->th[k] = invsinA*vm->th[j];
   rvm->vs[k] = rvfac*vm->vs[j];
   }

rvm->nlay = k + 1;

for(i=0;i<rvm->nlay;i++)
   rvm->invb2[i] = 1.0/(rvm->vs[i]*rvm->vs[i]);
}

void conv2vrup_dd(struct velmodel *vm,struct velmodel *rvm,float *dip,float *ztop,float *wid,float *rvf,float *shal_vr,float *dmin,float *dmax)
{
int i, j, k;
float invsinA, dep, zbot;
char string[256];

float rperd = 0.017453293;
float rvfac;

rvm->nlay = vm->nlay;
rvm->vs = (double *)check_malloc(rvm->nlay*sizeof(double));
rvm->th = (float *)check_malloc(rvm->nlay*sizeof(float));
rvm->invb2 = (double *)check_malloc(rvm->nlay*sizeof(double));

i = 0;
dep = vm->th[0];
while(dep < (*ztop))
   {
   i++;
   dep = dep + vm->th[i];
   }

zbot = *ztop + (*wid)*sin((*dip)*rperd);
invsinA = 1.0/sin((*dip)*rperd);

if(dep >= (*dmax))
   rvfac = (*rvf);
else if(dep < (*dmax) && dep > (*dmin))
   rvfac = (*rvf)*(1.0 - (1.0 - (*shal_vr))*((*dmax)-dep)/((*dmax)-(*dmin)));
else
   rvfac = (*rvf)*(*shal_vr);

rvm->th[0] = invsinA*(dep - (*ztop));
rvm->vs[0] = rvfac*vm->vs[i];

j = i;
k = 0;
while(dep < zbot)
   {
   j++; k++;
   dep = dep + vm->th[j];

   if(dep >= (*dmax))
      rvfac = (*rvf);
   else if(dep < (*dmax) && dep > (*dmin))
      rvfac = (*rvf)*(1.0 - (1.0 - (*shal_vr))*((*dmax)-dep)/((*dmax)-(*dmin)));
   else
      rvfac = (*rvf)*(*shal_vr);

   rvm->th[k] = invsinA*vm->th[j];
   rvm->vs[k] = rvfac*vm->vs[j];
   }

rvm->nlay = k + 1;

for(i=0;i<rvm->nlay;i++)
   rvm->invb2[i] = 1.0/(rvm->vs[i]*rvm->vs[i]);
}

void conv2vrup_dd2(struct velmodel *vm,struct velmodel *rvm,float *dip,float *ztop,float *wid,float *rvf,float *shal_vr,float *dmin1,float *dmax1,float *dmin2,float *dmax2)
{
int i, j, k;
float invsinA, dep, zbot;
char string[256];

float rperd = 0.017453293;
float rvfac;

rvm->nlay = vm->nlay;
rvm->vs = (double *)check_malloc(rvm->nlay*sizeof(double));
rvm->th = (float *)check_malloc(rvm->nlay*sizeof(float));
rvm->invb2 = (double *)check_malloc(rvm->nlay*sizeof(double));

i = 0;
dep = vm->th[0];
while(dep < (*ztop))
   {
   i++;
   dep = dep + vm->th[i];
   }

zbot = *ztop + (*wid)*sin((*dip)*rperd);
invsinA = 1.0/sin((*dip)*rperd);

if(dep >= (*dmax1) && dep <= (*dmin2))
   rvfac = (*rvf);
else if(dep < (*dmax1) && dep > (*dmin1))
   rvfac = (*rvf)*(1.0 - (1.0 - (*shal_vr))*((*dmax1)-dep)/((*dmax1)-(*dmin1)));
else if(dep < (*dmax2) && dep > (*dmin2))
   rvfac = (*rvf)*(1.0 - (1.0 - (*shal_vr))*(dep-(*dmin2))/((*dmax2)-(*dmin2)));
else
   rvfac = (*rvf)*(*shal_vr);

rvm->th[0] = invsinA*(dep - (*ztop));
rvm->vs[0] = rvfac*vm->vs[i];

j = i;
k = 0;
while(dep < zbot)
   {
   j++; k++;
   dep = dep + vm->th[j];

   if(dep >= (*dmax1) && dep <= (*dmin2))
      rvfac = (*rvf);
   else if(dep < (*dmax1) && dep > (*dmin1))
      rvfac = (*rvf)*(1.0 - (1.0 - (*shal_vr))*((*dmax1)-dep)/((*dmax1)-(*dmin1)));
   else if(dep < (*dmax2) && dep > (*dmin2))
      rvfac = (*rvf)*(1.0 - (1.0 - (*shal_vr))*(dep-(*dmin2))/((*dmax2)-(*dmin2)));
   else
      rvfac = (*rvf)*(*shal_vr);

   rvm->th[k] = invsinA*vm->th[j];
   rvm->vs[k] = rvfac*vm->vs[j];
   }

rvm->nlay = k + 1;

for(i=0;i<rvm->nlay;i++)
   rvm->invb2[i] = 1.0/(rvm->vs[i]*rvm->vs[i]);
}

void get_rupt(vm,h,srcd,recd,srcr,recr,p,rad,tt)
struct velmodel *vm;
float *h, *srcr, *recr, *recd, *tt, *srcd;
double *p, *rad;
{
double sth, rth, rng;
float sdep, rdep;
float tol;
float tup, thead;
int k, slay, rlay, linc;

float tenth = 0.1;
double eps = 1.0e-12;

tol = tenth*(*h);

rng = *srcr - *recr;
if(rng < 0.0)
   rng = -rng;

k = 0;
sdep = vm->th[0];
while((*srcd) > sdep)
   {
   k++;
   sdep = sdep + vm->th[k];
   }
slay = k;

k = 0;
rdep = vm->th[0];
while((*recd) > rdep)
   {
   k++;
   rdep = rdep + vm->th[k];
   }
rlay = k;

sth = sdep - *srcd;
rth = rdep - *recd;
get_headtime(vm,slay,&sth,rlay,&rth,&rng,&thead);

if(slay != rlay)
   {
   if(sdep > rdep)
      {
      sth = vm->th[slay] - (sdep - *srcd);
      rth = rdep - *recd;
      linc = -1;
      }
   else
      {
      sth = sdep - *srcd;
      rth = vm->th[rlay] - (rdep - *recd);
      linc = 1;
      }

/*
   bisection method
*/
    bisect_p(vm,slay,&sth,rlay,&rth,p,&eps,&tol,&rng,linc);

         /* get path length and travel time for correct ray parameter */

   get_radtime(vm,slay,&sth,rlay,&rth,p,rad,&tup,linc);
   }
else
   {
   *rad = sqrt(rng*rng + ((*srcd)-(*recd))*((*srcd)-(*recd)));
   tup = (*rad)/vm->vs[slay];
   }

*tt = thead;
if(tup < thead)
   *tt = tup;
/*
else
   fprintf(stderr,"*** thead selected\n");

fprintf(stderr,"*** thd= %f tup= %f\n",thead,tup);
*/
}

void get_rspeed(struct velmodel *vm,float *rspd,struct pointsource *ps,int nx,int ny,float *smax,float *savg,float *rvfav,float *shal_vr,float *dmin1,float *dmax1,float *deep_vr,float *dmin2,float *dmax2,float *rvfmn,float *rvfmx,int scl_slip)
{
int ix, iy, k, ip;
float rfdep, rfslip, rfavg, dep;

float rfmin, rfmax;

rfmin = (*rvfmn);
rfmax = (*rvfmx);
rfavg = (*rvfav);

for(iy=0;iy<ny;iy++)
   {
   k = 0;
   dep = vm->th[0];
   while(dep < ps[iy*nx].dep && k < vm->nlay)
      {
      k++;
      dep = dep + vm->th[k];
      }

   if(ps[iy*nx].dep <= (*dmin1))
      rfdep = (*shal_vr);
   else if(ps[iy*nx].dep < (*dmax1) && ps[iy*nx].dep > (*dmin1))
      rfdep = (1.0 - (1.0 - (*shal_vr))*((*dmax1)-ps[iy*nx].dep)/((*dmax1)-(*dmin1)));
   else if(ps[iy*nx].dep >= (*dmax1) && ps[iy*nx].dep <= (*dmin2))
      rfdep = 1.0;
   else if(ps[iy*nx].dep < (*dmax2) && ps[iy*nx].dep > (*dmin2))
      rfdep = (1.0 - (1.0 - (*deep_vr))*(ps[iy*nx].dep-(*dmin2))/((*dmax2)-(*dmin2)));
   else
      rfdep = (*deep_vr);

   for(ix=0;ix<nx;ix++)
      {
      ip = ix + iy*nx;

      if(scl_slip)
         {
         if(ps[ip].slip < *savg)
            rfslip = rfmin + (rfavg-rfmin)*ps[ip].slip/(*savg);
         else
            rfslip = rfavg + (rfmax-rfavg)*(ps[ip].slip-(*savg))/((*smax)-(*savg));
	 }
      else
         rfslip = rfavg;

      rspd[ip] = rfslip*rfdep*vm->vs[k];
      }
   }
}

void get_rspeed_vsden(float *rspd,struct pointsource *ps,int nx,int ny,float *smax,float *savg,float *rvfav,float *shal_vr,float *dmin1,float *dmax1,float *deep_vr,float *dmin2,float *dmax2,float *rvfmn,float *rvfmx,int scl_slip)
{
int ix, iy, k, ip;
float rfdep, rfslip, rfavg, dep;

float rfmin, rfmax;

rfmin = (*rvfmn);
rfmax = (*rvfmx);
rfavg = (*rvfav);

for(iy=0;iy<ny;iy++)
   {
   if(ps[iy*nx].dep <= (*dmin1))
      rfdep = (*shal_vr);
   else if(ps[iy*nx].dep < (*dmax1) && ps[iy*nx].dep > (*dmin1))
      rfdep = (1.0 - (1.0 - (*shal_vr))*((*dmax1)-ps[iy*nx].dep)/((*dmax1)-(*dmin1)));
   else if(ps[iy*nx].dep >= (*dmax1) && ps[iy*nx].dep <= (*dmin2))
      rfdep = 1.0;
   else if(ps[iy*nx].dep < (*dmax2) && ps[iy*nx].dep > (*dmin2))
      rfdep = (1.0 - (1.0 - (*deep_vr))*(ps[iy*nx].dep-(*dmin2))/((*dmax2)-(*dmin2)));
   else
      rfdep = (*deep_vr);

   for(ix=0;ix<nx;ix++)
      {
      ip = ix + iy*nx;

      if(scl_slip)
         {
         if(ps[ip].slip < *savg)
            rfslip = rfmin + (rfavg-rfmin)*ps[ip].slip/(*savg);
         else
            rfslip = rfavg + (rfmax-rfavg)*(ps[ip].slip-(*savg))/((*smax)-(*savg));
	 }
      else
         rfslip = rfavg;

      rspd[ip] = rfslip*rfdep*ps[ip].vs;
      }
   }
}

void get_rspeed_rvfslip(struct velmodel *vm,float *rspd,struct pointsource *ps,float *slip_r,int nx,int ny,float *rvfav,float *shal_vr,float *dmin1,float *dmax1,float *deep_vr,float *dmin2,float *dmax2,float *rvfmn,float *rvfmx,float *rvf_sig)
{
int ix, iy, k, ip;
int isupsh, imin, imax;
float rfdep, rfslip, rfavg, dep;
float rfmin, rfmax;
float slp_sig, slp_avg, sigfac;

rfmin = (*rvfmn);
rfmax = (*rvfmx);
rfavg = (*rvfav);

slp_avg = 0.0;
for(iy=0;iy<ny;iy++)
   {
   for(ix=0;ix<nx;ix++)
      {
      ip = ix + iy*nx;

      slp_avg = slp_avg + slip_r[ip];
      }
   }
slp_avg = slp_avg/(nx*ny);

slp_sig = 0.0;
for(iy=0;iy<ny;iy++)
   {
   for(ix=0;ix<nx;ix++)
      {
      ip = ix + iy*nx;

      slp_sig = slp_sig + (slip_r[ip]-slp_avg)*(slip_r[ip]-slp_avg);
      }
   }
slp_sig = sqrt(slp_sig/(nx*ny));

sigfac = (*rvf_sig)/slp_sig;

isupsh = 0;
imin = 0;
imax = 0;
for(iy=0;iy<ny;iy++)
   {
   k = 0;
   dep = vm->th[0];
   while(dep < ps[iy*nx].dep && k < vm->nlay)
      {
      k++;
      dep = dep + vm->th[k];
      }

   if(ps[iy*nx].dep <= (*dmin1))
      rfdep = (*shal_vr);
   else if(ps[iy*nx].dep < (*dmax1) && ps[iy*nx].dep > (*dmin1))
      rfdep = (1.0 - (1.0 - (*shal_vr))*((*dmax1)-ps[iy*nx].dep)/((*dmax1)-(*dmin1)));
   else if(ps[iy*nx].dep >= (*dmax1) && ps[iy*nx].dep <= (*dmin2))
      rfdep = 1.0;
   else if(ps[iy*nx].dep < (*dmax2) && ps[iy*nx].dep > (*dmin2))
      rfdep = (1.0 - (1.0 - (*deep_vr))*(ps[iy*nx].dep-(*dmin2))/((*dmax2)-(*dmin2)));
   else
      rfdep = (*deep_vr);

   for(ix=0;ix<nx;ix++)
      {
      ip = ix + iy*nx;

      rfslip = rfavg + sigfac*(slip_r[ip]-slp_avg);
/*
fprintf(stderr,"%.5e %.5e\n",rfslip,slip_r[ip]);
*/

      if(rfslip > 1.0)
         isupsh++;

      if(rfslip > rfmax)
         {
         rfslip = rfmax;
	 imax++;
	 }
      if(rfslip < rfmin)
         {
         rfslip = rfmin;
	 imin++;
	 }

      rspd[ip] = rfslip*rfdep*vm->vs[k];
      }
   }

fprintf(stderr,"*** rvfsig stats: <(%.3f Vs)= %.5f\n",rfmin,(float)(imin)/(float)(nx*ny));
fprintf(stderr,"                  >(1.000 Vs)= %.5f\n",(float)(isupsh)/(float)(nx*ny));
fprintf(stderr,"                  >(%.3f Vs)= %.5f\n",rfmax,(float)(imax)/(float)(nx*ny));
}

void get_rslow(float *rspd,double *rslw,int nx,int ny,float *tsf,long *seed)
{
int ix, iy, ip;

float sf = 1.0;
float fzero = 0.0;

for(iy=0;iy<ny;iy++)
   {
   for(ix=0;ix<nx;ix++)
      {
      ip = ix + iy*nx;

      if(*tsf > 0.0)
         sf = exp(gaus_rand(tsf,&fzero,seed));

      rslw[ip] = sf/rspd[ip];
      }
   }
}

void get_rslow_stretch(float *rspd,int ns,int nd,double *rslw,int nsfd,int ndfd,int isoff,int idoff,float *tsf,long *seed)
{
int i, j, ip, iq;

float sf = 1.0;
float fzero = 0.0;

for(j=0;j<nd;j++)
   {
   for(i=0;i<ns;i++)
      {
      ip = i + j*ns;
      iq = i+isoff + (j+idoff)*nsfd;

      if(*tsf > 0.0)
         sf = exp(gaus_rand(tsf,&fzero,seed));

      rslw[iq] = sf/rspd[ip];
      }
   }

for(j=0;j<nd;j++)
   {
   for(i=0;i<isoff;i++)
      {
      ip = j*ns;
      iq = i + (j+idoff)*nsfd;

      if(*tsf > 0.0)
         sf = exp(gaus_rand(tsf,&fzero,seed));

      rslw[iq] = sf/rspd[ip];
      }
   }

for(j=0;j<nd;j++)
   {
   for(i=ns;i<nsfd;i++)
      {
      ip = (ns-1) + j*ns;
      iq = i + (j+idoff)*nsfd;

      if(*tsf > 0.0)
         sf = exp(gaus_rand(tsf,&fzero,seed));

      rslw[iq] = sf/rspd[ip];
      }
   }

for(j=0;j<idoff;j++)
   {
   for(i=0;i<nsfd;i++)
      {
      ip = i + idoff*nsfd;
      iq = i + j*nsfd;

      rslw[iq] = rslw[ip];
      }
   }

for(j=nd;j<ndfd;j++)
   {
   for(i=0;i<nsfd;i++)
      {
      ip = i + (nd-1)*nsfd;
      iq = i + j*nsfd;

      rslw[iq] = rslw[ip];
      }
   }
}

void get_rsegdelay(float *rspd,int nx,int ny,double *dh,int nb,float *gbnd,float *gwid,float *rvf)
{
int ib, ixs, ixe, ix, iy, ip;

for(ib=0;ib<nb;ib++)
   {
   ixs = (int)(1.0*(gbnd[ib]-0.5*gwid[ib])/(*dh) + 0.5);
   ixe = (int)(1.0*(gbnd[ib]+0.5*gwid[ib])/(*dh) + 0.5);

   if(ixs < 0)
      ixs = 0;
   if(ixe > nx-1)
      ixe = nx-1;

   if(ixs<ixe)
      {
      for(iy=0;iy<ny;iy++)
         {
         for(ix=ixs;ix<=ixe;ix++)
            {
            ip = ix + iy*nx;

	    rspd[ip] = rvf[ib]*rspd[ip];
            }
         }
      }
   }
}

void load_vsden(struct pointsource *ps,struct velmodel *vm,int nstk,int ndip)
{
int i, j, k, ip;

for(j=0;j<ndip;j++)
   {
   for(i=0;i<nstk;i++)
      {
      ip = i + j*nstk;

      k = 0;
      while(ps[ip].dep > vm->dep[k] && k < (vm->nlay)-1)
         k++;

      ps[ip].vs = vm->vs[k];
      ps[ip].den = vm->den[k];
      ps[ip].mu = vm->mu[k];
      }
   }
}
