#include "include.h"
#include "structure.h"
#include "function.h"

main(ac,av)
int ac;
char **av;
{
FILE *fpw, *fopfile();
struct statdata head1;
float *s1, tst, tt;
int i;
char infile[128], outfile[128];

int inbin = 0;
float tshift = 0.0;
int zerotst = 0;
int itskip = 1;
float scale = 1.0;
float shift = 0.0;
float tstart = -1.0e+15;
float tend = 1.0e+15;

sprintf(infile,"stdin");
sprintf(outfile,"stdout");

setpar(ac,av);
getpar("infile","s",infile);
getpar("outfile","s",outfile);
getpar("inbin","d",&inbin);
getpar("scale","f",&scale);
getpar("shift","f",&shift);
getpar("tshift","f",&tshift);
getpar("zerotst","d",&zerotst);
getpar("itskip","d",&itskip);
getpar("tstart","f",&tstart);
getpar("tend","f",&tend);
endpar();
 
s1 = NULL;
s1 = read_wccseis(infile,&head1,s1,inbin);

if(strcmp(outfile,"stdout") == 0)
   fpw = stdout;
else
   fpw = fopfile(outfile,"w");

tst = 3600.0*head1.hr + 60.0*head1.min + head1.sec;
if(zerotst)
   tst = 0.0;

for(i=0;i<head1.nt;i=i+itskip)
   {
   tt = tshift+tst+i*head1.dt;
   if(tt >= tstart && tt <= tend)
      fprintf(fpw,"%13.5e %13.5e\n",tt,scale*(s1[i]+shift));
   }

fclose(fpw);
}
