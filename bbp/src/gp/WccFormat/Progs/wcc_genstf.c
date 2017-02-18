#include "include.h"
#include "structure.h"
#include "function.h"

#define AMP 3.0
#define PI 3.141592654

struct complex
   {
   float re;
   float im;
   };

void wfilt(struct complex *,int,float *,float *);
void randomize(float *,int,float *,float *,float *,long *,int);
double frand(void);
double sfrand(long *);

int size_float = sizeof(float);

main(int ac,char **av)
{
FILE *fpw, *fopfile();
struct statdata head1;
float *st;
int nt6, nt, it, i, j, ip;
float tsh, dt;
float *t0ptr, t0[60], ts[20], a0[20];

float tshift = 0.0;
int npeak = 1;

int bigN = 1;
int pos_only = 1;

int outbin = 0;

char title[128];
char name[16];
char comp[8];
char stype[16];

char outfile[256];

float rand = 0.0;
float fzero = 1.0;
long seed = 1;

float scale2slip = 1.0;

sprintf(stype,"tri");
sprintf(name,"stf");
sprintf(comp,"stf");
sprintf(title,"TITLE");

setpar(ac,av);
mstpar("dt","f",&dt);
mstpar("nt","d",&nt);
mstpar("outfile","s",outfile);
mstpar("npeak","d",&npeak);
mstpar("t0","vf",t0);
mstpar("ts","vf",ts);
mstpar("a0","vf",a0);

getpar("tshift","f",&tshift);
getpar("stype","s",stype);

if(strncmp("irikura",stype,7) == 0)
   mstpar("bigN","d",&bigN);

getpar("name","s",name);
getpar("comp","s",comp);
getpar("title","s",title);
getpar("rand","f",&rand);
getpar("pos_only","d",&pos_only);
getpar("fzero","f",&fzero);
getpar("seed","d",&seed);
getpar("outbin","d",&outbin);

getpar("scale2slip","f",&scale2slip);
endpar();

st = (float *) check_malloc (nt*size_float);

zap(st,nt);
t0ptr = t0;
for(ip=0;ip<npeak;ip++)
   {
   tsh = tshift + ts[ip];
   tsh = ts[ip];
   makesource(stype,st,&dt,nt,t0ptr,&a0[ip],&tsh,bigN);

   if(strncmp("trap",stype,4) == 0)
      t0ptr = t0ptr + 3;
   else if(strncmp("2tri",stype,4) == 0)
      t0ptr = t0ptr + 2;
   else
      t0ptr++;
   }

if(rand != (float)(0.0))
   randomize(st,nt,&rand,&dt,&fzero,&seed,pos_only);

if(tshift != (float)(0.0))
   apply_tsh(st,nt,&dt,&tshift);

normal(st,nt,&dt,&scale2slip);

strcpy(head1.stat,name);
strcpy(head1.comp,comp);
strcpy(head1.stitle,title);

head1.nt = nt;
head1.dt = dt;

head1.hr = 0;
head1.min = 0;
head1.sec = -tshift;

head1.edist = 0.0;
head1.az = 0.0;
head1.baz = 0.0;

write_wccseis(outfile,&head1,st,outbin);
}

makesource(type,s,dt,nt,x0,src_amp,tsh,bn)
float *s;
float *dt, *x0, *src_amp, *tsh;
int nt, bn;
char *type;
{
double exp(), exp_arg, cos_arg;
float t, t0, t2, tshift, alpha, beta, gamma;
float ep, em;
float td, tb, tr, ts, c0;
int it, itsh;
char str[128];
 
int ite, it2;
int it0, it1;
float da;

float tau, tau1, tau2, tau1x2, arg1, arg2;

float half = 0.5;
float one = 1.0;
float two = 2.0;

itsh = (*tsh)/(*dt);

if(strncmp("delta",type,5) == 0)
   s[itsh] = s[itsh] + *src_amp;
 
if(strncmp("gaus",type,4) == 0)
   {
   tshift = AMP*(*x0);
   alpha = -1.0/((*x0)*(*x0));
   for(it=0;it<nt;it++)
      {
      t = (it-1)*(*dt) - tshift;
      t2 = t*t;
      exp_arg = alpha*t2;
 
      s[it] = s[it] + (*src_amp)*exp(exp_arg);
      }
   }
 
if(strncmp("ucsb",type,4) == 0)
   {
   tau = (*x0);
   tau1 = 0.13*tau;
   tau2 = tau - tau1;
   tau1x2 = 2.0*tau1;

   it0 = (int)(tau/(*dt));
   for(it=0;it<it0;it++)
      {
      t = it*(*dt);

      alpha = 0.0;
      if(t < tau1)
         {
	 arg1 = PI*t/tau1;
	 arg2 = 0.5*arg1;
	 alpha = 0.7 - 0.7*cos(arg1) + 0.6*sin(arg2);
	 }
      else if(t < tau1x2)
         {
	 arg1 = PI*t/tau1;
	 arg2 = PI*(t - tau1)/tau2;
	 alpha = 1.0 - 0.7*cos(arg1) + 0.3*cos(arg2);
	 }
      else if(t < tau)
         {
	 arg1 = PI*(t - tau1)/tau2;
	 alpha = 0.3 + 0.3*cos(arg1);
	 }

      s[it+itsh] = s[it+itsh] + alpha;
      }
   }
 
if(strncmp("esg2006",type,7) == 0)
   {
   alpha = 4.0/(*x0);
   beta = 2.0*(*x0);
   gamma = alpha/sqrt(PI);

   for(it=0;it<nt-itsh;it++)
      {
      t = it*(*dt);

      arg1 = alpha*(t - beta);
      s[it+itsh] = s[it+itsh] + gamma*exp(-arg1*arg1);
      /*
      s[it+itsh] = s[it+itsh] + 0.5*(1.0 + erf(arg1));
      */
      }
   }

if(strncmp("tanh",type,4) == 0)
   {
   alpha = 2.0/(*x0);
   for(it=-itsh;it<nt-itsh;it++)
      {
      t = it*(*dt);
      exp_arg = t*alpha + 1.0;
      ep = exp(exp_arg);
      em = exp(-exp_arg);
 
      s[it+itsh] = s[it+itsh] + (*src_amp)*(alpha/((ep + em)*(ep + em)));
      }
   }

if(strncmp("tri",type,3) == 0)
   {
   ite = (int)((*x0)/(*dt) + 0.5);
   if(ite%2)
      ite++;

   it2 = ite/2;
   da = (float)((*src_amp)/it2);
 
   for(it=0;it<it2;it++)
      {
      s[it+itsh] = s[it+itsh] + it*da;
      s[(ite+itsh)-it] = s[(ite+itsh)-it] + it*da;
      }
   s[it2+itsh] = s[it2+itsh] + it2*da;
   }

if(strncmp("atri",type,4) == 0)
   {
   if(type[4] == '\0')
      alpha = 0.5;
   else
      alpha = 0.01*atof(type+4);

   ite = (int)((*x0)/(*dt) + 0.5);
   if(ite%2)
      ite++;

   it2 = ite*alpha;
   alpha = (*src_amp)*(float)(2.0/ite);

   da = alpha/(float)(it2);
   for(it=0;it<=it2;it++)
      s[it+itsh] = s[it+itsh] + it*da;

   da = alpha/(float)(ite-it2);
   for(it=it2+1;it<=ite;it++)
      s[it+itsh] = s[it+itsh] + (ite-it)*da;
   }

if(strncmp("gabor",type,5) == 0)
   {
   alpha = 11.0;
   beta = 2.0*PI/(*x0);
   tshift = 0.45*alpha*(*x0);
   for(it=0;it<nt;it++)
      {
      t = (it)*(*dt) - tshift;
      exp_arg = t*beta/alpha;
      exp_arg = -exp_arg*exp_arg;
      cos_arg = beta*t + PI/2.0;
 
      if(t <= tshift)
         s[it+itsh] = s[it+itsh] + (*src_amp)*exp(exp_arg)*cos(cos_arg);
      }
   }

if(strncmp("rick",type,4) == 0)
   {
   tshift = AMP*(*x0);
   alpha = -1.0/((*x0)*(*x0));
   beta = 2.0/(*x0);
   beta = PI/(*x0);
   for(it=0;it<nt;it++)
      {
      t = (it-1)*(*dt) - tshift;
      t2 = t*t;
      exp_arg = alpha*t2;
      cos_arg = beta*t;
 
      s[it+itsh] = s[it+itsh] + (*src_amp)*exp(exp_arg)*cos(cos_arg);
      }
   }

if(strncmp("zahr",type,4) == 0)
   {
   tshift = 0.5*AMP*(*x0);
   beta = PI/(*x0);
   alpha = beta*beta;
   for(it=0;it<nt;it++)
      {
      t = (it-1)*(*dt) - tshift;
      t2 = t*t;
      exp_arg = alpha*t2;
 
      s[it+itsh] = s[it+itsh] + (*src_amp)*(two*exp_arg - one)*exp(-exp_arg);
      }
   }

if(strncmp("cos",type,3) == 0)
   {
   alpha = 2.0*PI/(*x0);
   ite = (*x0)/(*dt);
   for(it=0;it<ite;it++)
      {
      cos_arg = it*(*dt)*alpha;

      s[it+itsh] = s[it+itsh] + (*src_amp)*(one - cos(cos_arg));
      }
   }

if(strncmp("expcos",type,6) == 0)
   {
   ite = (*x0)/(*dt);
   for(it=0;it<ite;it++)
      {
      cos_arg = it*(*dt)/(*x0);

      s[it+itsh] = s[it+itsh] + (*src_amp)*cos_arg*(one + cos(cos_arg*PI))*exp(-cos_arg);
      }
   }

if(strncmp("expcos2",type,7) == 0)
   {
   ite = (*x0)/(*dt);
   for(it=0;it<ite;it++)
      {
      cos_arg = it*(*dt)/(*x0);

      s[it+itsh] = s[it+itsh] + (*src_amp)*(one-cos(2.0*cos_arg*PI))*exp(-cos_arg);
      }
   }

if(strncmp("scec",type,4) == 0)
   {
   alpha = 1.0/((*x0)*(*x0));
   for(it=0;it<nt-itsh;it++)
      {
      exp_arg = it*(*dt);

      s[it+itsh] = s[it+itsh] + (*src_amp)*exp_arg*alpha*(exp(-exp_arg/(*x0)));
      }
   }

if(strncmp("brune",type,5) == 0)
   {
   t0 = 0.1600*(*x0);
   t0 = 0.2108*(*x0);  /* x0 is t95 ~ rise time */
   t0 = 0.1250*(*x0);
   for(it=0;it<nt-itsh;it++)
      {
      exp_arg = it*(*dt)/t0;

      s[it+itsh] = s[it+itsh] + (*src_amp)*(exp_arg/t0)*(exp(-exp_arg));
      }
   }

if(strncmp("sqrtT",type,5) == 0)
   {
   for(it=0;it<nt-itsh;it++)
      {
      if(it == 0)
	 t = 0.01*(*dt);
      else
         t = it*(*dt);

      s[it+itsh] = s[it+itsh] + (one - t)/(sqrt(t));
      }
   }

if(strncmp("miyatake",type,8) == 0)
   {
   td = 0.1;
   tb = 1.25*td;
   tr = (*x0);
   ts = 1.5*(*x0);

   alpha = two*(*src_amp)/td;
   gamma = half/td;
   ep = (2.5*tb - 3.0*td)/(one - td/tb);
   beta = alpha*tb*(one - gamma*tb)*sqrt(tb - ep);
   c0 = beta/sqrt(tr - ep);
   da = -c0/(ts - tr);

   it0 = (int)(tb/(*dt)) + 1;
   it1 = (int)(tr/(*dt)) + 1;
   it2 = (int)(ts/(*dt)) + 1;

   if(it0 > nt-itsh)
      it0 = nt-itsh;
   if(it1 > nt-itsh)
      it1 = nt-itsh;
   if(it2 > nt-itsh)
      it2 = nt-itsh;

   /*
   for(it=0;it<it0;it++)
      {
      t = it*(*dt);
      s[it+itsh] = s[it+itsh] + alpha*t*(one - gamma*t);
      }

   for(it=it0;it<it1;it++)
      {
      t = it*(*dt);
      s[it+itsh] = s[it+itsh] + beta/sqrt(t - ep);
      }

   for(it=it1;it<it2;it++)
      {
      t = it*(*dt);
      s[it+itsh] = s[it+itsh] + c0 + (t - tr)*da;
      }
*/

   /*
   */
   for(it=0;it<it2;it++)
      {
      t = it*(*dt);

      if(t <= tb)
         s[it+itsh] = s[it+itsh] + alpha*t*(one - gamma*t);
      else if(t <= tr)
         s[it+itsh] = s[it+itsh] + beta/sqrt(t - ep);
      else
         s[it+itsh] = s[it+itsh] + c0 + (t - tr)*da;
      }
   }

if(strncmp("2triGOOD",type,8) == 0)
   {
   if(type[4] == '\0')  /* set at 10% of total rise time */
      alpha = 0.1*(*x0);
   else if(strncmp("-p",type+4,2) == 0)  /* percentage is given after '-p' */
      alpha = 0.01*(*x0)*atof(type+6);
   else
      alpha = atof(type+4);

   it0 = (int)((alpha)/(*dt) + 0.5);
   if(it0 < 2)
      it0 = 2;
   it1 = (int)((*x0)/(*dt) + 0.5);
   if(it1 < 4)
      it1 = 2;

   it2 = it0 + it0/2;

   alpha = 1.0;
   da = alpha/(float)(it0);

   for(it=0;it<it0;it++)
      s[it+itsh] = s[it+itsh] + it*da;

   for(it=it0;it<it2;it++)
      s[it+itsh] = s[it+itsh] + (2*it0-it)*da;

   alpha = 0.5*alpha;
   da = alpha/(float)(it1-it2);

   for(it=it2;it<it1;it++)
      s[it+itsh] = s[it+itsh] + alpha + (it2-it)*da;
   }

if(strncmp("2tri",type,4) == 0)
   {
   /* DEFAULT: set t0 at 10% of total rise time, h at 0.2*A */

   alpha = 0.1*(*x0);
   beta = 0.2;

   if(strncmp("-p",type+4,2) == 0)       /* t0 percentage given after '-p' */
      {
      alpha = 0.01*(*x0)*atof(type+6);

      if(strncmp("-h",type+8,2) == 0)    /* ampl. percentage given after '-h' */
         beta = 0.01*atof(type+10);
      }

   if(strncmp("-h",type+4,2) == 0)       /* ampl. percentage given after '-h' */
      {
      beta = 0.01*atof(type+6);

      if(strncmp("-p",type+8,2) == 0)    /* t0 percentage given after '-p' */
         alpha = 0.01*(*x0)*atof(type+10);
      }

   it0 = (int)((alpha)/(*dt) + 0.5);
   if(it0 < 2)
      it0 = 2;
   it1 = (int)((*x0)/(*dt) + 0.5);
   if(it1 < 4)
      it1 = 2;

   it2 = (2 - beta)*it0;

   alpha = 1.0;
   da = alpha/(float)(it0);

   for(it=0;it<it0;it++)
      s[it+itsh] = s[it+itsh] + it*da;

   for(it=it0;it<it2;it++)
      s[it+itsh] = s[it+itsh] + (2*it0-it)*da;

   alpha = beta*alpha;
   da = alpha/(float)(it1-it2);

   for(it=it2;it<it1;it++)
      s[it+itsh] = s[it+itsh] + alpha + (it2-it)*da;
   }

if(strncmp("trap",type,4) == 0)
   {
   it0 = (int)((x0[0])/(*dt) + 0.5);
   if(it0 < 0)
      it0 = 0;
   it1 = (int)((x0[1])/(*dt) + 0.5);
   if(it1 < 0)
      it1 = 0;
   it2 = (int)((x0[2])/(*dt) + 0.5);
   if(it2 < 0)
      it2 = 0;

   da = (float)(it0)/(float)(it2);
   it1 = it0 + it1;
   it2 = it1 + it2;

   for(it=0;it<it0;it++)
      s[it+itsh] = s[it+itsh] + it;

   for(it=it0;it<it1;it++)
      s[it+itsh] = s[it+itsh] + it0;

   for(it=it1;it<it2;it++)
      s[it+itsh] = s[it+itsh] + (it2-it)*da;
   }

if(strncmp("irikura",type,7) == 0)
   {
   it0 = (int)((x0[0])/(*dt) + 0.5);
   beta = 1.0/(float)(it0);
   alpha = (x0[0])/((bn-1.0)*(*dt))*(1.0 - exp(-1.0));
   alpha = 1.0/alpha;

   for(it=0;it<it0;it++)
      s[it+itsh] = s[it+itsh] + alpha*exp(-it*beta);

   s[itsh] = s[itsh] + 1.0;
   /*
   */
   }
}

zap(s,n)
float *s;
int n;
{
int i;

for(i=0;i<n;i++)
   s[i] = 0.0;
}

apply_tsh(s,nt,dt,tsh)
float *s, *dt, *tsh;
int nt;
{
int it, itsh;

itsh = (*tsh)/(*dt);

if(itsh > 0)
   {
   for(it=nt-1;it>=itsh;it--)
      s[it] = s[it-itsh];

   for(it=0;it<itsh;it++)
      s[it] = 0.0;
   }
else if(itsh < 0)
   {
   for(it=0;it<nt+itsh;it++)
      s[it] = s[it-itsh];
   }

*tsh = (*dt)*itsh;
}

normal(s,n,dt,scl)
float *s, *dt, *scl;
int n;
{
int i;
float sum = 0.0;

for(i=0;i<n;i++)
   sum = sum + s[i]*(*dt);

sum = (*scl)/sum;
if(sum < 0.0)
   sum = -sum;

for(i=0;i<n;i++)
   s[i] = s[i]*sum;
}

void randomize(float *st,int nt,float *rand,float *dt,float *f0,long *seed,int pos_only)
{
int it, ntp2, tapl;
float rfac, rmax, *rt;
float arg;

float amax = -99.0;
int itzero = 0;
int nzflag = 0;

for(it=0;it<nt;it++)
   {
   if(st[it] > amax)
      amax = st[it];

   if(st[it] == (float)(0.0) && nzflag)
      itzero = it;

   if(st[it] != (float)(0.0))
      nzflag = 1;
   else
      nzflag = 0;
   }

ntp2 = 2;
while(itzero > ntp2)
   ntp2 = ntp2*2;

rt = (float *) check_malloc (ntp2*sizeof(float));

for(it=0;it<ntp2;it++)
   rt[it] = sfrand(seed);

forfft(rt,ntp2,-1);
wfilt((struct complex *) rt,ntp2,dt,f0);
invfft(rt,ntp2,1);

tapl = (int)(0.5/((*f0)*(*dt)));
for(it=0;it<tapl;it++)
   {
   arg = it*PI/tapl;

   rt[it] = 0.5*(1.0 - cos(arg))*rt[it];
   rt[itzero-it-1] = 0.5*(1.0 - cos(arg))*rt[itzero-it-1];
   }

rmax = 0.0;
for(it=0;it<itzero;it++)
   {
   if(rt[it] > rmax)
      rmax = rt[it];
   if(-rt[it] > rmax)
      rmax = -rt[it];
   }

rfac = (*rand)*amax/rmax;

for(it=0;it<itzero;it++)
   {
   st[it] = st[it] + rfac*rt[it];
   if(st[it] < 0.0 && pos_only)
      st[it] = 0.0;

   /*
   st[it] = rfac*rt[it];
   */
   }
}

void wfilt(struct complex *rc,int n,float *dt,float *f0)
{
int i;
float t2, df, f, amp, af;
float norm;

float pi = 3.14159;

norm = 4.0*pi/(float)(n);

t2 = 1.0/((*f0)*(*f0));
df = 1.0/(n*(*dt));
for(i=1;i<n/2;i++)
   {
   amp = sqrt(rc[i].re*rc[i].re + rc[i].im*rc[i].im);

   f = i*df;
   af = 1.0/(1.0 + f*f*t2);

   if(amp == 0.0)
      af = 0.0;
   else
      af = af/amp;

   rc[i].re = norm*af*rc[i].re;
   rc[i].im = norm*af*rc[i].im;
   }
rc[0].re = norm;
rc[0].re = 0.0;
rc[0].im = 0.0;
}

static  long    frandx = 1;

/* frand() returns a uniform distribution of random numbers
 * in the range -1.0 -> 1.0.
 */
double frand(void)
{
frandx = (frandx * 1103515245 + 12345) & 0x7fffffff;
return((double)(frandx)/1073741824.0 - 1.0);
}

/* sfrand() returns a uniform distribution of random numbers
 * in the range -1.0 -> 1.0.
 */
double sfrand(long *seed)
{
*seed = ((*seed) * 1103515245 + 12345) & 0x7fffffff;
return((double)(*seed)/1073741824.0 - 1.0);
}
