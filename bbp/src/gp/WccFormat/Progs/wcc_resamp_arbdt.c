/***********************************************************

wcc_resamp_arbdt:   Decimate or interpolate time series
                    to new (arbitrary) sampling rate

                    copyright (c) 1999
                     Robert W. Graves
                 Woodward-Clyde Consultants
                   566 El Dorado Street
                    Pasadena, CA 91101

                    tel (626) 449-7650
                    fax (626) 449-3536

***********************************************************/

#include "include.h"
#include "structure.h"
#include "function.h"
#include "fftw3.h"
#include "getpar.h"

void resample_fftw(float *,int,float *,int,int,int,double *,int,float *,float *);
void resample_fftwf(float *,int,float *,int,int,int,float *,int,float *,float *);
void resample(float *,int,float *,int,int,int,float *,float *,int,float *,float *);
void zapit(float *,int);
void taper_norm(float *,float *,int,float *);
void norm(float *,float *,int);
double nt_tol(float,int);
double nt_tol_d(double,int);

extern void fourg_(float *, int *, int *, float *);

#define TAP_PERC 0.05

int main(int ac,char **av)
{
struct statdata head1;
float *s1;
float *p;
int nt1, i, j;

int it, ntmax, ntap;
char dt_str[128];
float df, fac, arg;
double dnt, double_dt;
float single_dt;

int order = 4;
float nyq_perc = 1.0;
float tap_perc = TAP_PERC;

int resamp = 0;
int ntpad, ntrsmp;
float fnt, *space;
int gnt;
int ntout = -1;

double dt_tol = 1.00001;
double tol = 1.0e-02;

char infile[1024];
char outfile[1024];

int inbin = 0;
int outbin = 0;

int use_fftw = 1;
int use_double = 0;

sprintf(infile,"stdin");
sprintf(outfile,"stdout");

setpar(ac,av);
getpar("infile","s",infile);
getpar("outfile","s",outfile);

mstpar("newdt","f",&single_dt);
mstpar("newdt","s",&dt_str);
double_dt = 1.000000*atof(dt_str);

getpar("use_fftw","d",&use_fftw);
getpar("use_double","d",&use_double);

getpar("nyq_perc","f",&nyq_perc);
getpar("tap_perc","f",&tap_perc);
getpar("order","d",&order);
getpar("ntout","d",&ntout);
getpar("inbin","d",&inbin);
getpar("outbin","d",&outbin);
endpar();

s1 = NULL;
s1 = read_wccseis(infile,&head1,s1,inbin);

fprintf(stderr,"***nt=%d dt=%f\n",head1.nt,head1.dt);

if(double_dt <= 0.0 || (double_dt >= head1.dt/dt_tol && double_dt <= head1.dt*dt_tol))
   {
   resamp = 0;

   if(ntout < 0)
      ntout = head1.nt;

   if(ntout > head1.nt)
      {
      s1 = (float *) check_realloc(s1,ntout*sizeof(float));

      ntap = (head1.nt)*(tap_perc);

      df = 3.14159/(float)(ntap);
      for(it=(head1.nt)-ntap;it<head1.nt;it++)
         {
         arg = (it-((head1.nt)-(ntap+1)))*df;
         fac = 0.5*(1.0 + cos(arg));
         s1[it] = s1[it]*fac;
         }

      for(it=head1.nt;it<ntout;it++)
         s1[it] = 0.0;

      head1.nt = ntout;
      }
   }
else                /* need to resample time history */
   {
   ntpad = 2*head1.nt;
   dnt = (double)(1.000000*ntpad)*(double)(head1.dt)/double_dt;
   gnt = (int)(dnt + 0.5);
   ntrsmp = (int)(dnt);

   while(nt_tol_d(dnt,gnt) > tol || (ntrsmp%2) == 1)
      {
      ntpad = ntpad + 2;
      dnt = (double)(1.000000*ntpad)*(double)(head1.dt)/double_dt;
      gnt = (int)(dnt + 0.5);
      ntrsmp = (int)(dnt);
      }

   fprintf(stderr,"*** ntpad=%d ntrsmp=%d\n",ntpad,ntrsmp);

   if(single_dt < head1.dt)
      {
      resamp = 1;
      ntmax = 2*ntrsmp;
      }
   else
      {
      resamp = -1;
      ntmax = 2*ntpad;
      }
   s1 = (float *) check_realloc(s1,ntmax*sizeof(float));

   if(use_fftw == 0)
      {
      space = (float *) check_malloc (ntmax*sizeof(float));
      resample(s1,head1.nt,&head1.dt,resamp,ntpad,ntrsmp,&single_dt,space,order,&nyq_perc,&tap_perc);
      }
   else
      {
      if(use_double == 0)
         resample_fftwf(s1,head1.nt,&head1.dt,resamp,ntpad,ntrsmp,&single_dt,order,&nyq_perc,&tap_perc);
      else
         resample_fftw(s1,head1.nt,&head1.dt,resamp,ntpad,ntrsmp,&double_dt,order,&nyq_perc,&tap_perc);
      }

   if(ntout < 0)
      ntout = (int)(head1.nt*head1.dt/double_dt);

   if(ntout > ntmax/2)
      {
      s1 = (float *) check_realloc(s1,ntout*sizeof(float));

      ntap = (ntmax/2)*(tap_perc);

      df = 3.14159/(float)(ntap);
      for(it=(ntmax/2)-ntap;it<ntmax/2;it++)
         {
         arg = (it-((ntmax/2)-(ntap+1)))*df;
         fac = 0.5*(1.0 + cos(arg));
         s1[it] = s1[it]*fac;
         }

      for(it=ntmax/2;it<ntout;it++)
         s1[it] = 0.0;
      }

   head1.nt = ntout;
   head1.dt = double_dt;
   }

fprintf(stderr,"***nt=%d dt=%f\n",head1.nt,head1.dt);

write_wccseis(outfile,&head1,s1,outbin);
}

void resample(float *s,int nt,float *dt,int isamp,int ntpad,int ntrsmp,float *newdt,float *p,int ord,float *perc, float *tp)
{
float df, f, f0, fl, fl2, fac;
int i, j;

float one = 1.0;

int minus = -1;
int plus = 1;

taper_norm(s,dt,nt,tp);
zapit(s+nt,(ntpad)-(nt));

for(i=ntpad-1;i>=0;i--)
   {
   s[2*i] = s[i];
   s[2*i + 1] = 0.0;
   }

fourg_(s,&ntpad,&minus,p);

if(isamp > 0)
   zapit(s+ntpad,2*ntrsmp-ntpad);
else if(isamp < 0)
   {
   if(ord)  /* lowpass at 100*(*perc) % of new Nyquist */
      {
      f0 = (*perc)/(2.0*(*newdt));
      df = 1.0/(ntrsmp*(*newdt));
      for(i=1;i<ntrsmp/2;i++)
         {
         f = i*df;

         fl = f/f0;
         fl2 = fl*fl;
         fl = fl2;
         for(j=1;j<ord;j++)
            fl = fl*fl2;

         fac = one/(one + fl);

         s[2*i] = fac*s[2*i];
         s[2*i + 1] = fac*s[2*i + 1];
         }
      }

   s[ntrsmp] = s[ntrsmp+1] = 0.0; /* zero nyquist */
   }

for(i=1;i<ntrsmp/2;i++)  /* replicate with complex-conjugate */
   {
   s[2*(ntrsmp-i)] = s[2*i];
   s[2*(ntrsmp-i) + 1] = -s[2*i + 1];
   }

fourg_(s,&ntrsmp,&plus,p);

for(i=0;i<ntrsmp;i++)
   s[i] = s[2*i];

norm(s,newdt,ntrsmp);
}

void zapit(float *s,int n)
{
while(n--)
   {
   s[0] = 0.0;
   s++;
   }
}

void taper_norm(float *g,float *dt,int nt,float *tp)
{
float fac, df, arg;
int i;
int ntap;

ntap = nt*(*tp);

for(i=0;i<nt-ntap;i++)
   g[i] = g[i]*(*dt);

df = 3.14159/(float)(ntap);
for(i=nt-ntap;i<nt;i++)
   {
   arg = (i-(nt-(ntap+1)))*df;
   fac = (*dt)*0.5*(1.0 + cos(arg));
   g[i] = g[i]*fac;
   }
}

void norm(float *g,float *dt,int nt)
{
float fac;

fac = 1.0/((*dt)*nt);
while(nt--)
   {
   g[0] = g[0]*fac;
   g++;
   }
}

double nt_tol(float fnt,int gnt)
{
double diff;

diff = ((double)(fnt) - (double)(gnt));
if(diff < 0.0)
   diff = -diff;

/*
fprintf(stderr,"diff= %15.10e\n",diff);
*/

return(diff);
}

double nt_tol_d(double dnt,int gnt)
{
double diff;

diff = (dnt - (double)(gnt));
if(diff < 0.0)
   diff = -diff;

return(diff);
}

void resample_fftw(float *s,int nt,float *dt,int isamp,int ntpad,int ntrsmp,double *newdt,int ord,float *perc, float *tp)
{
double df, f, f0, fl, fl2, fac;
int i, j;

int minus = -1;
int plus = 1;

fftw_plan plan;
fftw_complex *dbl_cs;

taper_norm(s,dt,nt,tp);
for(i=nt;i<ntpad;i++)
   s[i] = 0.0;

dbl_cs = (fftw_complex *) check_malloc(sizeof(fftw_complex)*ntpad);
plan = fftw_plan_dft_1d(ntpad,dbl_cs,dbl_cs,FFTW_FORWARD,FFTW_ESTIMATE);

for(i=0;i<ntpad;i++)
   {
   dbl_cs[i][0] = s[i];
   dbl_cs[i][1] = 0.0;
   }

fftw_execute(plan);

dbl_cs = (fftw_complex *) check_realloc(dbl_cs,sizeof(fftw_complex)*ntrsmp);

if(isamp > 0)
   {
   for(i=ntpad/2;i<ntrsmp;i++)
      {
      dbl_cs[i][0] = 0.0;
      dbl_cs[i][1] = 0.0;
      }
   }
else if(isamp < 0)
   {
   if(ord)  /* lowpass at 100*(*perc) % of new Nyquist */
      {
      f0 = (*perc)/(2.0*(*newdt));
      df = 1.0/(ntrsmp*(*newdt));
      for(i=1;i<ntrsmp/2;i++)
         {
         f = i*df;

         fl = f/f0;
         fl2 = fl*fl;
         fl = fl2;
         for(j=1;j<ord;j++)
            fl = fl*fl2;

         fac = 1.0/(1.0 + fl);

         dbl_cs[i][0] = fac*dbl_cs[i][0];
         dbl_cs[i][1] = fac*dbl_cs[i][1];
         }
      }

   /* zero nyquist */
   dbl_cs[ntrsmp/2][0] = 0.0;
   dbl_cs[ntrsmp/2][1] = 0.0;
   }

for(i=1;i<ntrsmp/2;i++)  /* replicate with complex-conjugate */
   {
   dbl_cs[(ntrsmp-i)][0] = dbl_cs[i][0];
   dbl_cs[(ntrsmp-i)][1] = -dbl_cs[i][1];
   }

plan = fftw_plan_dft_1d(ntrsmp,dbl_cs,dbl_cs,FFTW_BACKWARD,FFTW_ESTIMATE);
fftw_execute(plan);

fac = 1.0/((*newdt)*ntrsmp);
for(i=0;i<ntrsmp;i++)
   s[i] = fac*dbl_cs[i][0];

fftw_destroy_plan(plan);
fftw_free(dbl_cs);
}

void resample_fftwf(float *s,int nt,float *dt,int isamp,int ntpad,int ntrsmp,float *newdt,int ord,float *perc, float *tp)
{
float df, f, f0, fl, fl2, fac;
int i, j;

int minus = -1;
int plus = 1;

fftwf_plan plan;
fftwf_complex *cs;

taper_norm(s,dt,nt,tp);
for(i=nt;i<ntpad;i++)
   s[i] = 0.0;

cs = (fftwf_complex *) check_malloc(sizeof(fftwf_complex)*ntpad);
plan = fftwf_plan_dft_1d(ntpad,cs,cs,FFTW_FORWARD,FFTW_ESTIMATE);

for(i=0;i<ntpad;i++)
   {
   cs[i][0] = s[i];
   cs[i][1] = 0.0;
   }

fftwf_execute(plan);

cs = (fftwf_complex *) check_realloc(cs,sizeof(fftwf_complex)*ntrsmp);

if(isamp > 0)
   {
   for(i=ntpad/2;i<ntrsmp;i++)
      {
      cs[i][0] = 0.0;
      cs[i][1] = 0.0;
      }
   }
else if(isamp < 0)
   {
   if(ord)  /* lowpass at 100*(*perc) % of new Nyquist */
      {
      f0 = (*perc)/(2.0*(*newdt));
      df = 1.0/(ntrsmp*(*newdt));
      for(i=1;i<ntrsmp/2;i++)
         {
         f = i*df;

         fl = f/f0;
         fl2 = fl*fl;
         fl = fl2;
         for(j=1;j<ord;j++)
            fl = fl*fl2;

         fac = 1.0/(1.0 + fl);

         cs[i][0] = fac*cs[i][0];
         cs[i][1] = fac*cs[i][1];
         }
      }

   /* zero nyquist */
   cs[ntrsmp/2][0] = 0.0;
   cs[ntrsmp/2][1] = 0.0;
   }

for(i=1;i<ntrsmp/2;i++)  /* replicate with complex-conjugate */
   {
   cs[(ntrsmp-i)][0] = cs[i][0];
   cs[(ntrsmp-i)][1] = -cs[i][1];
   }

plan = fftwf_plan_dft_1d(ntrsmp,cs,cs,FFTW_BACKWARD,FFTW_ESTIMATE);
fftwf_execute(plan);

fac = 1.0/((*newdt)*ntrsmp);
for(i=0;i<ntrsmp;i++)
   s[i] = fac*cs[i][0];

fftwf_destroy_plan(plan);
fftwf_free(cs);
}
