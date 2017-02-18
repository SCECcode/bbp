#include "include.h"
#include "structure.h"
#include "function.h"

#define          SMALL   1.0e-15

int size_float = sizeof(float);
int size_int = sizeof(int);

main(ac,av)
int ac;
char **av;
{
struct statdata shead1, shead2, shead3;
float *s1, *s2;
char infile1[512], infile2[512];

float dt, time1, time2, *p1, *p2;
float cvel, *cc, tdel;
int i, nt, itcc;
int istr, iend;

int inbin1 = 0;
int inbin2 = 0;

int itwin;
float twin = -99.9;

int printmax = 1;
int printall = 0;
float tval = -1.0e+15;

setpar(ac,av);
mstpar("infile1","s",infile1);
mstpar("infile2","s",infile2);
getpar("inbin1","d",&inbin1);
getpar("inbin2","d",&inbin2);
getpar("twin","f",&twin);
getpar("tval","f",&tval);
getpar("printmax","d",&printmax);
getpar("printall","d",&printall);
endpar();

s1 = NULL;
s1 = read_wccseis(infile1,&shead1,s1,inbin1);
s2 = NULL;
s2 = read_wccseis(infile2,&shead2,s2,inbin2);

time1 = shead1.hr*3600 + shead1.min*60 + shead1.sec;
time2 = shead2.hr*3600 + shead2.min*60 + shead2.sec;

if(shead1.dt != shead2.dt)
   {
   fprintf(stderr,"dt1 not equal to dt2, exiting...\n");
   exit(-1);
   }
else
   dt = shead1.dt;

if(shead1.nt < shead2.nt)
   nt = shead2.nt;
else
   nt = shead1.nt;

p1 = (float *) check_malloc (3*nt*size_float);
p2 = (float *) check_malloc (nt*size_float);

for(i=0;i<nt;i++)
   p1[i] = p1[i+2*nt] = 0.0;
for(i=0;i<shead1.nt;i++)
   p1[i+nt] = s1[i];
for(i=shead1.nt;i<nt;i++)
   p1[i+nt] = 0.0;

for(i=0;i<shead2.nt;i++)
   p2[i] = s2[i];
for(i=shead2.nt;i<nt;i++)
   p2[i] = 0.0;

if(twin < 0.0)
   twin = 2*nt*dt;

cc = (float *) check_malloc (2*nt*size_float);

xcor(p1,p2,nt,cc);

if(printmax)
   {
   itcc = maxcctime(cc,nt,&time1,&time2,&dt,&twin);
   tdel = time2 - time1 - (itcc-nt)*dt;

   printf("cc= %13.5e t= %13.8f\n",cc[itcc],tdel);
   }

if(tval > -1.0e+10)
   {
   itcc = nt + (time2 - time1 - tval)/dt;

   printf("cc= %13.5e t= %13.8f\n",cc[itcc],tval);
   }

if(printall)
   {
   istr = 0;
   iend = 2*nt;
   if(tval > -1.0e+10 && twin > 0)
      {
      istr = nt + (time2 - time1 - tval - twin)/dt;
      iend = nt + (time2 - time1 - tval + twin)/dt;
      }
   for(i=istr;i<iend;i++)
      printf("t= %f cc= %13.5e\n",time2 - time1 - (i-nt)*dt,cc[i]);
   }
}

xcor(s1,s2,nt,cc)
float *s1, *s2, *cc;
int nt;
{
float *s1ptr;
float zap = 0.0;
int it;

for(it=0;it<2*nt;it++)
   cc[it] = zap;

for(it=0;it<2*nt;it++)
   {
   s1ptr = s1 + it;
   cor(s1ptr,s2,nt,&cc[it],1);
   }
}

cor(s1,s2,span,p,it)
float *s1, *s2, *p;
int span, it;
{
int i;
int i2 = 0;
int i3 = 0;
float fac = 1.0;
float sum1 = SMALL;
float sum2 = SMALL;
float sum3 = SMALL;
float smin = 1.0e-10;

for(i=0;i<span;i=i+it)
   {
   sum1 = sum1 + s1[i]*s2[i];
   sum2 = sum2 + s1[i]*s1[i];
   sum3 = sum3 + s2[i]*s2[i];
   }

if(sum2 == 0.0 || sum3 == 0.0)
   *p = -1.0;
else
   {
   while(sum2 <= smin)  /* for very low sum2,sum3 adjust to prevent underflow */      {  
      sum2 = sum2/smin;
      i2++;
      }    
   while(sum3 <= smin)
      {  
      sum3 = sum3/smin;
      i3++;
      }    

   if(i2 || i3)
      fac = 1.0/sqrt((i2+i3)*smin);

   *p = (fac*sum1)/sqrt(sum2*sum3);
   }
}

maxcctime(cc,nt,t1,t2,dt,tw)
float *cc, *t1, *t2, *dt, *tw;
int nt;
{
int it, itx, its, ite;
float max = -10.0;

its = nt + ((*t2) - (*t1) - (*tw))/(*dt);
if(its < 0)
   its = 0;

ite = nt + ((*t2) - (*t1) + (*tw))/(*dt);
if(ite > 2*nt)
   ite = 2*nt;

for(it=its;it<ite;it++)
   {
   if(cc[it] > max)
      {  
      max = cc[it];
      itx = it;
      }  
   }   
return(itx);
}
