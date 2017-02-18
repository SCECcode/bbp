#include "include.h"
#include "structure.h"
#include "function.h"

main(ac,av)
int ac;
char **av;
{
FILE *fpw, *fopfile();
struct statdata head1;
float *s1;
int i, nt6, j;
char infile[128], outfile[128], frmt[32];

int inbin = 0;
int nperline = 8;

sprintf(infile,"stdin");
sprintf(outfile,"stdout");

sprintf(frmt,"%%9.3f");

setpar(ac,av);
getpar("infile","s",infile);
getpar("outfile","s",outfile);
getpar("inbin","d",&inbin);
getpar("nperline","d",&nperline);
getpar("frmt","s",frmt);
endpar();
 
s1 = NULL;
s1 = read_wccseis(infile,&head1,s1,inbin);

if(strcmp(outfile,"stdout") == 0)
   fpw = stdout;
else
   fpw = fopfile(outfile,"w");

fprintf(fpw,"%-10s %3s %s\n",head1.stat,head1.comp,head1.stitle);
fprintf(fpw,"%d %12.5e %d %d %12.5e %12.5e %12.5e %12.5e\n",head1.nt,
                                           head1.dt,
                                           head1.hr,
                                           head1.min,
                                           head1.sec,
                                           head1.edist,
                                           head1.az,
                                           head1.baz);

nt6 = head1.nt/nperline;
for(i=0;i<nt6;i++)
   {
   for(j=0;j<nperline;j++)
      fprintf(fpw,frmt,s1[nperline*i + j]);

   fprintf(fpw,"\n");
   }

if(nperline*nt6 != head1.nt)
   {
   for(i=nperline*nt6;i<head1.nt;i++)
      fprintf(fpw,frmt,s1[i]);

   fprintf(fpw,"\n");
   }
fclose(fpw);
}
