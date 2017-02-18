#include "include.h"
#include "structure.h"
#include "function.h"
#include "defs.h"

int gen_2tri_stfOLD(float *slip,float *trise,float *stf,int nt,float *dt,float *z0)
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
float dmax = 8.0;
float rtfac;
float rtfac0 = 1.0;

dbdd = (betadeep - betashal)/(dmax-dmin);

if((*z0) >= dmax)
   beta = betadeep;
else if((*z0) < dmax && (*z0) > dmin)
   beta = betadeep - (dmax-(*z0))*dbdd;
else
   beta = betashal;

zapit(stf,nt);

if((*z0) >= dmax)
   rtfac = 1.0;
else if((*z0) < dmax && (*z0) > dmin)
   rtfac = 1.0 + rtfac0*(dmax-(*z0))/(dmax-dmin);
else
   rtfac = 1.0 + rtfac0;

tr = (*trise)*rtfac;

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
sum = (*slip)/sum;
for(it=0;it<nstf;it++)
   stf[it] = stf[it]*sum;

return(nstf);
}

int gen_brune_stf(float *slip,float *t0,float *stf,int nt,float *dt,float *z0)
{
int it, nstf;
float t95, tend, sfac, tfac;
float sum;

zapit(stf,nt);

t95 = 1.745*exp(1.0)*(*t0);
tend = 3.0*t95;

nstf = (int)((tend)/(*dt) + 0.5);
if(nstf > nt)
   nstf = nt;

if(nstf == 0)
   return(0);

sfac = (*slip)/(*t0);
tfac = (*dt)/(*t0);
for(it=0;it<nstf;it++)
   stf[it] = sfac*(it*tfac)*exp(-it*tfac);

sum = 0.0;
for(it=0;it<nstf;it++)
   sum = sum + (*dt)*stf[it];

if(sum <= 0.0)
   return(0);

return(nstf);
}

int gen_2tri_stf(float *slip,float *trise,float *stf,int nt,float *dt,float *z0)
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

tr = (*trise);
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
sum = (*slip)/sum;
for(it=0;it<nstf;it++)
   stf[it] = stf[it]*sum;

return(nstf);
}

int gen_ucsb_stf(float *slip,float *t0,float *stf,int nt,float *dt,float *z0)
{
int it, nstf;
float tau, tau1, tau2, tau1x2, arg1, arg2;
float sum, t, alpha;
float pi = 3.141592654;

zapit(stf,nt);

tau = (*t0);
tau1 = 0.13*tau;
tau2 = tau - tau1;
tau1x2 = 2.0*tau1;

nstf = (int)((tau)/(*dt) + 0.5);
if(nstf > nt)
   nstf = nt;

if(nstf == 0)
   return(0);

for(it=0;it<nstf;it++)
   {
   t = it*(*dt);

   alpha = 0.0;
   if(t < tau1)
      {
      arg1 = pi*t/tau1;
      arg2 = 0.5*arg1;
      alpha = 0.7 - 0.7*cos(arg1) + 0.6*sin(arg2);
      }
   else if(t < tau1x2)
      {
      arg1 = pi*t/tau1;
      arg2 = pi*(t - tau1)/tau2;
      alpha = 1.0 - 0.7*cos(arg1) + 0.3*cos(arg2);
      }
   else if(t < tau)
      {
      arg1 = pi*(t - tau1)/tau2;
      alpha = 0.3 + 0.3*cos(arg1);
      }

   stf[it] = alpha;
   }

sum = 0.0;
for(it=0;it<nstf;it++)
   sum = sum + (*dt)*stf[it];

if(sum <= 0.0)
   return(0);

/* scale STF by slip */
sum = (*slip)/sum;
for(it=0;it<nstf;it++)
   stf[it] = stf[it]*sum;

return(nstf);
}

int gen_ucsbvT_stf(float *slip,float *t0,float *stf,int nt,float *dt,float *trat)
{
int it, nstf;
float tau, tau1, tau2, tau1x2, arg1, arg2;
float sum, t, alpha;
float pi = 3.141592654;

zapit(stf,nt);

tau = (*t0);
tau1 = (*trat)*tau;
tau2 = tau - tau1;
tau1x2 = 2.0*tau1;

nstf = (int)((tau)/(*dt) + 0.5);
if(nstf > nt)
   nstf = nt;

if(nstf == 0)
   return(0);

for(it=0;it<nstf;it++)
   {
   t = it*(*dt);

   alpha = 0.0;
   if(t < tau1)
      {
      arg1 = pi*t/tau1;
      arg2 = 0.5*arg1;
      alpha = 0.7 - 0.7*cos(arg1) + 0.6*sin(arg2);
      }
   else if(t < tau1x2)
      {
      arg1 = pi*t/tau1;
      arg2 = pi*(t - tau1)/tau2;
      alpha = 1.0 - 0.7*cos(arg1) + 0.3*cos(arg2);
      }
   else if(t < tau)
      {
      arg1 = pi*(t - tau1)/tau2;
      alpha = 0.3 + 0.3*cos(arg1);
      }

   stf[it] = alpha;
   }

sum = 0.0;
for(it=0;it<nstf;it++)
   sum = sum + (*dt)*stf[it];

if(sum <= 0.0)
   return(0);

/* scale STF by slip */
sum = (*slip)/sum;
for(it=0;it<nstf;it++)
   stf[it] = stf[it]*sum;

return(nstf);
}

int gen_ucsb2_stf(float *slip,float *t0,float *stf,int nt,float *dt,float *z0)
{
int it, nstf;
float tau, tau1, tau2, tau1x2, arg1, arg2;
float sum, t, alpha;
float pi = 3.141592654;

zapit(stf,nt);

tau = 2.0*(*t0);
tau1 = 0.5*0.13*tau;   /* keep peak at same place */
tau2 = tau - tau1;
tau1x2 = 2.0*tau1;

nstf = (int)((tau)/(*dt) + 0.5);
if(nstf > nt)
   nstf = nt;

if(nstf == 0)
   return(0);

for(it=0;it<nstf;it++)
   {
   t = it*(*dt);

   alpha = 0.0;
   if(t < tau1)
      {
      arg1 = pi*t/tau1;
      arg2 = 0.5*arg1;
      alpha = 0.7 - 0.7*cos(arg1) + 0.6*sin(arg2);
      }
   else if(t < tau1x2)
      {
      arg1 = pi*t/tau1;
      arg2 = pi*(t - tau1)/tau2;
      alpha = 1.0 - 0.7*cos(arg1) + 0.3*cos(arg2);
      }
   else if(t < tau)
      {
      arg1 = pi*(t - tau1)/tau2;
      alpha = 0.3 + 0.3*cos(arg1);
      }

   stf[it] = alpha;
   }

sum = 0.0;
for(it=0;it<nstf;it++)
   sum = sum + (*dt)*stf[it];

if(sum <= 0.0)
   return(0);

/* scale STF by slip */
sum = (*slip)/sum;
for(it=0;it<nstf;it++)
   stf[it] = stf[it]*sum;

return(nstf);
}

int gen_ucsbT_stf(float *slip,float *t0,float *stf,int nt,float *dt,float *beta)
{
int it, nstf;
float tau, tau1, tau2, tau1x2, arg1, arg2;
float sum, t, alpha;
float pi = 3.141592654;

zapit(stf,nt);

tau = (*beta)*(*t0);
tau1 = 0.13*tau/(*beta);   /* keep peak at same place */
tau2 = tau - tau1;
tau1x2 = 2.0*tau1;

nstf = (int)((tau)/(*dt) + 0.5);
if(nstf > nt)
   nstf = nt;

if(nstf == 0)
   return(0);

for(it=0;it<nstf;it++)
   {
   t = it*(*dt);

   alpha = 0.0;
   if(t < tau1)
      {
      arg1 = pi*t/tau1;
      arg2 = 0.5*arg1;
      alpha = 0.7 - 0.7*cos(arg1) + 0.6*sin(arg2);
      }
   else if(t < tau1x2)
      {
      arg1 = pi*t/tau1;
      arg2 = pi*(t - tau1)/tau2;
      alpha = 1.0 - 0.7*cos(arg1) + 0.3*cos(arg2);
      }
   else if(t < tau)
      {
      arg1 = pi*(t - tau1)/tau2;
      alpha = 0.3 + 0.3*cos(arg1);
      }

   stf[it] = alpha;
   }

sum = 0.0;
for(it=0;it<nstf;it++)
   sum = sum + (*dt)*stf[it];

if(sum <= 0.0)
   return(0);

/* scale STF by slip */
sum = (*slip)/sum;
for(it=0;it<nstf;it++)
   stf[it] = stf[it]*sum;

return(nstf);
}

int gen_esg2006_stf(float *slip,float *trise,float *stf,int nt,float *dt,float *z0)
{
int it, nstf;
int ip, it0, it1, it2;
float tr, amp, a0, tt;
float sum, arg1;
float alpha, beta, gamma;
float pi = 3.14159265;

zapit(stf,nt);

tr = (*trise);
alpha = 4.0/tr;
beta= 2.0*tr;
gamma = alpha/sqrt(pi);

nstf = (int)((2.0*beta)/(*dt) + 0.5);
if(nstf > nt)
   nstf = nt;

if(nstf == 0)
   return(0);

for(it=0;it<nstf;it++)
   {
   tt = it*(*dt);

   arg1 = alpha*(tt - beta);
   stf[it] = exp(-arg1*arg1);
   }

for(it=0;it<nstf;it++)
   sum = sum + (*dt)*stf[it];

if(sum <= 0.0)
   return(0);

/* scale STF by slip */
sum = (*slip)/sum;
for(it=0;it<nstf;it++)
   stf[it] = stf[it]*sum;

return(nstf);
}

int gen_cos_stf(float *slip,float *t0,float *stf,int nt,float *dt,float *z0)
{
int it, nstf;
float tfac, sum;
float pi = 3.14159265;
float one = 1.0;

zapit(stf,nt);

nstf = (int)((*t0)/(*dt));
if(nstf > nt)
   nstf = nt;

if(nstf == 0)
   return(0);

tfac = 2.0*pi*(*dt)/(*t0);
for(it=0;it<nstf;it++)
   stf[it] = (one - cos(it*tfac));

sum = 0.0;
for(it=0;it<nstf;it++)
   sum = sum + (*dt)*stf[it];

if(sum <= 0.0)
   return(0);

/* scale STF by slip */
sum = (*slip)/sum;
for(it=0;it<nstf;it++)
   stf[it] = stf[it]*sum;

return(nstf);
}

int gen_seki_stf(float *slip,float *t0,float *stf,int nt,float *dt,float *z0)
{
int it, nstf;
float tfac, sum, t, htan, ep, em, arg;

zapit(stf,nt);

nstf = (int)(1.5*(*t0)/(*dt));
if(nstf > nt)
   nstf = nt;

if(nstf == 0)
   return(0);

tfac = 2.0/(*t0);
for(it=0;it<nstf;it++)
   {
   t = it*(*dt);
   arg = 2.0*(t*tfac - 1.0);
   ep = exp(arg);
   em = exp(-arg);
   htan = (ep - em)/(ep + em);

   stf[it] = (*slip)*tfac*(1.0 - htan*htan);
   }

sum = 0.0;
for(it=0;it<nstf;it++)
   sum = sum + (*dt)*stf[it];

if(sum <= 0.0)
   return(0);

/* scale STF by slip */
sum = (*slip)/sum;
for(it=0;it<nstf;it++)
   stf[it] = stf[it]*sum;

return(nstf);
}
