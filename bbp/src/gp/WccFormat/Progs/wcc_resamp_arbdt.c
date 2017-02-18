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
#include "getpar.h"

void resample(float *,int,float *,int,int,int,float *,float *,int,float *,float *);
void zapit(float *,int);
void taper_norm(float *,float *,int,float *);
void norm(float *,float *,int);
double nt_tol(float,int);

#define TAP_PERC 0.20

int size_cx = sizeof(struct complex);
int size_float = sizeof(float);
int size_int = sizeof(int);

extern void fourg_(float *, int *, int *, float *);

int main(int ac,char **av)
{
struct statdata head1;
float *s1;
float *p;
int nt1, i, j;

int order = 4;
float nyq_perc = 1.0;
float tap_perc = TAP_PERC;

int resamp = 0;
int ntpad, ntrsmp;
float fnt, dtout, *space;
int gnt;
float newdt = -1.0;
int ntout = -1;

float tol = 1.0e-02;

char infile[256];
char outfile[256];

int inbin = 0;
int outbin = 0;

int pow2 = 0;

sprintf(infile,"stdin");
sprintf(outfile,"stdout");

setpar(ac,av);
getpar("infile","s",infile);
getpar("outfile","s",outfile);
mstpar("newdt","f",&newdt);
getpar("nyq_perc","f",&nyq_perc);
getpar("tap_perc","f",&tap_perc);
getpar("order","d",&order);
getpar("ntout","d",&ntout);
getpar("inbin","d",&inbin);
getpar("outbin","d",&outbin);
getpar("pow2","d",&pow2);
endpar();

s1 = NULL;
s1 = read_wccseis(infile,&head1,s1,inbin);

if(newdt < 0.0 || newdt == head1.dt)
   {
   resamp = 0;

   if(ntout < 0)
      ntout = head1.nt;
   if(ntout > head1.nt)
      ntout = head1.nt;
   dtout = head1.dt;
   }
else                /* need to resample time history */
   {
   ntpad = 2*head1.nt;
   fnt = ntpad*head1.dt/newdt;
   gnt = (int)(fnt + 0.5);

   //force power of 2
   /*if (pow2==1) {
	printf("Forcing to power of 2.\n");
	int tot_t = head1.dt*head1.nt;
	//continue...
	ntrsmp = (int)fnt;
	while(nt_tol(fnt,gnt)>tol || ( (ntrsmp & (ntrsmp-1))!=0) ){
		ntpad++;
		fnt = ntpad*head1.dt/newdt;
		gnt = (int)(fnt + 0.5);
		ntrsmp = (int)fnt;
	}
   } else {*/
   	while(nt_tol(fnt,gnt) > tol)
   	   {
   	   ntpad++;
   	   fnt = ntpad*head1.dt/newdt;
   	   gnt = (int)(fnt + 0.5);
   	   }
   //}

   ntrsmp = (int)(fnt);

   if(ntout < 0)
      ntout = (int)(head1.nt*head1.dt/newdt);
   dtout = newdt;

   fprintf(stderr,"*** ntpad=%d ntrsmp=%d\n",ntpad,ntrsmp);

   if(newdt < head1.dt)
      {
      resamp = 1;

      if(ntout > ntrsmp)
         ntout = ntrsmp;

      space = (float *) check_malloc (2*ntrsmp*sizeof(float));
      s1 = (float *) check_realloc(s1,2*ntrsmp*sizeof(float));
      }
   else
      {
      resamp = -1;

      if(ntout > ntpad)
         ntout = ntpad;

      space = (float *) check_malloc (2*ntpad*sizeof(float));
      s1 = (float *) check_realloc(s1,2*ntpad*sizeof(float));
      }
   }

fprintf(stderr,"***nt=%d dt=%f\n",head1.nt,head1.dt);

if(resamp != 0)
   resample(s1,head1.nt,&head1.dt,resamp,ntpad,ntrsmp,&newdt,space,order,&nyq_perc,&tap_perc);

head1.nt = ntout;
head1.dt = dtout;

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
