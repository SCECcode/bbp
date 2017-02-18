#include "include.h"
#include "structure.h"
#include "function.h"
#include "getpar.h"

int main(int ac,char **av)
{
FILE *fpr;
struct statdata head1, head2, head3;
float *s1, *s2, *s3;
int it;
double tst, *tbuf;

char nsfile[256];
char ewfile[256];
char udfile[256];
char bbfile[256];
char str[1024];

char nsname[128];
char ewname[128];
char udname[128];
char units[128];

char tformat[16];
char dformat[16];
char frmt[256];

char stat[16];
char comp[4];
char title[128];

int wcc2bbp = 1;

int inbin = 0;
int outbin = 0;

int bufinc = 10000;
int bufcnt = 1;

stat[0] = '\0';
title[0] = '\0';

nsname[0] = '\0';
ewname[0] = '\0';
udname[0] = '\0';
units[0] = '\0';

sprintf(tformat,"%%.6e");
sprintf(dformat,"\t%%.6e");

setpar(ac,av);

mstpar("nsfile","s",nsfile);
mstpar("ewfile","s",ewfile);
mstpar("udfile","s",udfile);

getpar("wcc2bbp","d",&wcc2bbp);

if(wcc2bbp == 1)
   {
   getpar("nsname","s",nsname);
   getpar("ewname","s",ewname);
   getpar("udname","s",udname);

   getpar("units","s",units);

   getpar("tformat","s",tformat);
   getpar("dformat","s",dformat);
   }

if(wcc2bbp == 1)
   sprintf(bbfile,"stdout");
else
   sprintf(bbfile,"stdin");

getpar("bbfile","s",bbfile);

getpar("stat","s",stat);
getpar("title","s",title);
getpar("inbin","d",&inbin);
getpar("outbin","d",&outbin);

endpar();

if(wcc2bbp == 1)
   {
   s1 = NULL;
   s1 = read_wccseis(nsfile,&head1,s1,inbin);

   s2 = NULL;
   s2 = read_wccseis(ewfile,&head2,s2,inbin);

   s3 = NULL;
   s3 = read_wccseis(udfile,&head3,s3,inbin);

   if(strcmp(bbfile,"stdout") == 0)
      fpr = stdout;
   else
      fpr = fopfile(bbfile,"w");

   if(title[0] != '\0')
      fprintf(fpr,"# %s\n",title);

   if(units[0] == '\0')
      sprintf(units,"cm/s");
   if(nsname[0] == '\0')
      sprintf(nsname,"N-S(%s)",units);
   if(ewname[0] == '\0')
      sprintf(ewname,"E-W(%s)",units);
   if(udname[0] == '\0')
      sprintf(udname,"U-D(%s)",units);

   if(strcmp(units,"-1") != 0)
      fprintf(fpr,"#    time(sec)      %s      %s      %s\n",nsname,ewname,udname);

   sprintf(frmt,"%s%s%s%s\n",tformat,dformat,dformat,dformat);

   tst = head1.sec;
   for(it=0;it<head1.nt;it++)
      {
      /*
      fprintf(fpr,"%14.6f %14.6e %14.6e %14.6e\n",tst,s1[it],s2[it],s3[it]);
      */
      fprintf(fpr,frmt,tst,s1[it],s2[it],s3[it]);
      tst = tst + head1.dt;
      }
   fclose(fpr);
   }
else
   {
   if(strcmp(bbfile,"stdin") == 0)
      fpr = stdin;
   else
      fpr = fopfile(bbfile,"r");

   fgets(str,1024,fpr);
   while(str[0] == '#' || str[0] == '%')
      fgets(str,1024,fpr);

   tbuf = (double *)check_malloc(bufinc*sizeof(double));
   s1 = (float *)check_malloc(bufinc*sizeof(float));
   s2 = (float *)check_malloc(bufinc*sizeof(float));
   s3 = (float *)check_malloc(bufinc*sizeof(float));

   sscanf(str,"%lf %f %f %f",&tbuf[0],&s1[0],&s2[0],&s3[0]);
   it = 1;



   while(fgets(str,1024,fpr) != NULL)
      {
      if(it >= bufcnt*bufinc)
	 {
         bufcnt++;
         tbuf = (double *)check_realloc(tbuf,bufcnt*bufinc*sizeof(double));
         s1 = (float *)check_realloc(s1,bufcnt*bufinc*sizeof(float));
         s2 = (float *)check_realloc(s2,bufcnt*bufinc*sizeof(float));
         s3 = (float *)check_realloc(s3,bufcnt*bufinc*sizeof(float));
	 }

      sscanf(str,"%lf %f %f %f",&tbuf[it],&s1[it],&s2[it],&s3[it]);
      it++;
      }
   fclose(fpr);

   if(stat[0] != '\0')
      {
      strcpy(head1.stat,stat);
      strcpy(head2.stat,stat);
      strcpy(head3.stat,stat);
      }
   else
      {
      sprintf(head1.stat,"UNKNOWN");
      sprintf(head2.stat,"UNKNOWN");
      sprintf(head3.stat,"UNKNOWN");
      }

   sprintf(head1.comp,"000");
   sprintf(head2.comp,"090");
   sprintf(head3.comp,"ver");

   strcpy(head1.stitle,title);
   strcpy(head2.stitle,title);
   strcpy(head3.stitle,title);

   head1.nt = it;
   head2.nt = it;
   head3.nt = it;

   head1.dt = (tbuf[it-1]-tbuf[0])/(double)(it-1.0);
   head2.dt = (tbuf[it-1]-tbuf[0])/(double)(it-1.0);
   head3.dt = (tbuf[it-1]-tbuf[0])/(double)(it-1.0);

   tst = tbuf[0];

   head1.hr = (int)(tst/3600.0);
   head2.hr = (int)(tst/3600.0);
   head3.hr = (int)(tst/3600.0);

   head1.min = (int)((tst - 3600.0*head1.hr)/60.0);
   head2.min = (int)((tst - 3600.0*head1.hr)/60.0);
   head3.min = (int)((tst - 3600.0*head1.hr)/60.0);

   head1.sec = (float)(tst - 60.0*head1.min - 3600.0*head1.hr);
   head2.sec = (float)(tst - 60.0*head1.min - 3600.0*head1.hr);
   head3.sec = (float)(tst - 60.0*head1.min - 3600.0*head1.hr);

   head1.edist = 0.0;
   head2.edist = 0.0;
   head3.edist = 0.0;

   head1.az = 0.0;
   head2.az = 0.0;
   head3.az = 0.0;

   head1.baz = 0.0;
   head2.baz = 0.0;
   head3.baz = 0.0;

   write_wccseis(nsfile,&head1,s1,outbin);
   write_wccseis(ewfile,&head2,s2,outbin);
   write_wccseis(udfile,&head3,s3,outbin);
   }
   return(0);
}
