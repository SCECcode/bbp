#include "include.h"
#include "structure.h"
#include "function.h"

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

   if(i==0)
      vm->dep[i] = vm->th[i];
   else
      vm->dep[i] = vm->dep[i-1] + vm->th[i];

   vm->mu[i] = vm->vs[i]*vm->vs[i]*vm->den[i]*1.0e+10;  /* in CMS units */
   vm->invb2[i] = 1.0/(vm->vs[i]*vm->vs[i]);
   }
fclose(fpr);
}

void conv2vrup(struct velmodel *vm,struct velmodel *rvm,float *dip,float *ztop,float *wid,float *rvf,float *shal_vr)
{
int i, j, k;
float invsinA, dep, zbot;

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
