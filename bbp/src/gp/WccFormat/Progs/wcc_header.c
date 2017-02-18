#include "include.h"
#include "structure.h"
#include "function.h"

int size_float = sizeof(float);
int size_int = sizeof(int);

main(int ac, char **av)
{
struct statdata head1;
float *s1;
int i, j, nt1;
int nt = -99;
float epi = -1.0e+10;
float baz = -999.9;
float tst = -1.0e+10;
float max = -99.0;
float scale = 0.0;
float dt = -99.0;

char infile[128];
char outfile[128];

char stat[16];
char comp[4];
char title[128];

int inbin = 0;
int outbin = 0;

stat[0] = '\0';
comp[0] = '\0';
title[0] = '\0';

setpar(ac,av);
mstpar("infile","s",infile);
mstpar("outfile","s",outfile);
getpar("stat","s",stat);
getpar("comp","s",comp);
getpar("title","s",title);
getpar("nt","d",&nt);
getpar("dt","f",&dt);
getpar("epi","f",&epi);
getpar("tst","f",&tst);
getpar("baz","f",&baz);
getpar("max","f",&max);
getpar("inbin","d",&inbin);
getpar("outbin","d",&outbin);
endpar();

s1 = NULL;
s1 = read_wccseis(infile,&head1,s1,inbin);

for(i=0;i<head1.nt;i++)
   {
   if(s1[i] > scale)
      scale = s1[i];
   else if(-s1[i] > scale)
      scale = -s1[i];
   }

if(stat[0] != '\0')
   strcpy(head1.stat,stat);

if(comp[0] != '\0')
   strcpy(head1.comp,comp);

if(title[0] != '\0')
   strcpy(head1.stitle,title);

if(max > -99.0)
   max = max/scale;
else
   max = 1.0;

if(epi > -1.0e+10)
   head1.edist = epi;

if(baz > -999.0)
   head1.baz = baz;

if(tst > -1.0e+10)
   {
   head1.hr = (int)(tst/3600.0);
   head1.min = (int)((tst - 3600.0*head1.hr)/60.0);
   head1.sec = (float)(tst - 60.0*head1.min - 3600.0*head1.hr);
   }

if(nt > 0)
   {
   if(nt > head1.nt)
      {
      s1 = (float *) check_realloc(s1,nt*sizeof(float));
      for(i=head1.nt;i<nt;i++)
         s1[i] = 0.0;
      }
   head1.nt = nt;
   }
if(dt > 0.0)
   head1.dt = dt;

for(i=0;i<head1.nt;i++)
   s1[i] = s1[i]*max;

write_wccseis(outfile,&head1,s1,outbin);
}
