#include <errno.h>
#include <fcntl.h>
#include <math.h>
#include <stdio.h>
#include <unistd.h>
#include <stdlib.h>

#include <sys/file.h>
#ifndef __APPLE__
#include <sys/procfs.h>
#endif
#include <sys/resource.h>
#include <sys/signal.h>
#include <sys/stat.h>
#include <sys/syscall.h>
#include <sys/time.h>
#include <sys/types.h>
#include <string.h>

#include "function.h"
#include "getpar.h"

#define SLEN 1024

int main(int ac,char **av)
{
FILE *fpr, *fopfile();
float *per, *sa1, *sa2, *res1, *res2, *avgr, v[8];
float rpga1, rpga2, avgpga, rpgv1, rpgv2, avgpgv;
int i, np2, np;
float flo = 0.0;
float fhi = 0.0;

char statinfo[512], str[SLEN];
char eq[128], mag[16], stat[16], lon[16], lat[16];
char seqno[16], vs30[16], cd[16], xc[16], yc[16], tmin[16], tmax[16];

int respect_format = 1;
int sa_field = 7;        /* 6=AA, 7=pseudo AA for respect format */

int bbp_format = 0;

int print_header = 0;

char datafile1[256], datafile2[256];
char simfile1[256], simfile2[256];
char comp1[128], comp2[128];

float pga1 = -999.9;
float pga2 = -999.9;
float pgv1 = -999.9;
float pgv2 = -999.9;

float sim_pga1 = -999.9;
float sim_pga2 = -999.9;
float sim_pgv1 = -999.9;
float sim_pgv2 = -999.9;

datafile1[0] = '\0';
datafile2[0] = '\0';
simfile1[0] = '\0';
simfile2[0] = '\0';

sprintf(eq,"-999");
sprintf(mag,"-999");
sprintf(stat,"-999");
sprintf(lon,"-999");
sprintf(lat,"-999");
sprintf(seqno,"-999");
sprintf(vs30,"-999");
sprintf(cd,"-999");
sprintf(xc,"-999");
sprintf(yc,"-999");

setpar(ac,av);

getpar("datafile1","s",datafile1);
getpar("datafile2","s",datafile2);
getpar("simfile1","s",simfile1);
getpar("simfile2","s",simfile2);

getpar("comp1","s",comp1);
getpar("comp2","s",comp2);

getpar("eqname","s",eq);
getpar("mag","s",mag);
getpar("stat","s",stat);
getpar("lon","s",lon);
getpar("lat","s",lat);
getpar("seqno","s",seqno);
getpar("vs30","s",vs30);
getpar("cd","s",cd);
getpar("xcos","s",xc);
getpar("ycos","s",yc);
getpar("flo","f",&flo);
getpar("fhi","f",&fhi);

getpar("print_header","d",&print_header);

getpar("pga1","f",&pga1);
getpar("pga2","f",&pga2);
getpar("pgv1","f",&pgv1);
getpar("pgv2","f",&pgv2);

getpar("sim_pga1","f",&sim_pga1);
getpar("sim_pga2","f",&sim_pga2);
getpar("sim_pgv1","f",&sim_pgv1);
getpar("sim_pgv2","f",&sim_pgv2);

getpar("bbp_format","d",&bbp_format);
if(bbp_format == 1)
   respect_format = 0;

getpar("respect_format","d",&respect_format);

if(respect_format == 0 && bbp_format == 0)
   getpar("np","d",&np);

getpar("sa_field","d",&sa_field);

endpar();

if(fhi > 0.0)
   sprintf(tmin,"%.3f",1.0/fhi);
else
   sprintf(tmin,"-99999.999");

if(flo > 0.0)
   sprintf(tmax,"%.3f",1.0/flo);
else
   sprintf(tmax,"99999.999");

sprintf(statinfo,"%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s",eq,mag,stat,lon,lat,seqno,vs30,cd,xc,yc,tmin,tmax);

if(respect_format)
   {
   fpr = fopfile(datafile1,"r");

   fgets(str,SLEN,fpr);
   fgets(str,SLEN,fpr);
   fgets(str,SLEN,fpr);
   fgets(str,SLEN,fpr);

   sscanf(str,"%d",&np);

   fclose(fpr);
   }
else if(bbp_format == 1)
   {
   fpr = fopfile(datafile1,"r");

   fgets(str,SLEN,fpr);
   while(strncmp(str,"#",1) == 0)
      fgets(str,SLEN,fpr);

   np = 1;
   while(fgets(str,SLEN,fpr) != NULL)
      np++;

   fclose(fpr);
   }

per = (float *)check_malloc(np*sizeof(float));
sa1 = (float *)check_malloc(np*sizeof(float));
sa2 = (float *)check_malloc(np*sizeof(float));
res1 = (float *)check_malloc(np*sizeof(float));
res2 = (float *)check_malloc(np*sizeof(float));
avgr = (float *)check_malloc(np*sizeof(float));

if(bbp_format == 1)
   read_bbpfile(datafile1,per,sa1,sa2,np);
else
   {
   read_specfile(datafile1,per,sa1,np,sa_field,respect_format);
   read_specfile(datafile2,per,sa2,np,sa_field,respect_format);
   }

if(simfile1[0] == '\0' || (bbp_format == 0 && simfile2[0] == '\0'))
   {
   for(i=0;i<np;i++)
      {
      res1[i] = sa1[i];
      res2[i] = sa2[i];
      avgr[i] = sqrt(res1[i]*res2[i]);
      }

   rpga1 = pga1;
   rpga2 = pga2;
   avgpga = -999.9;
   if(rpga1 >= 0.0 && rpga2 >= 0.0)
      avgpga = sqrt(rpga1*rpga2);

   rpgv1 = pgv1;
   rpgv2 = pgv2;
   avgpgv = -999.9;
   if(rpgv1 >= 0.0 && rpgv2 >= 0.0)
      avgpgv = sqrt(rpgv1*rpgv2);
   }
else
   {
   if(bbp_format == 1)
      read_bbpfile(simfile1,per,res1,res2,np);
   else
      {
      read_specfile(simfile1,per,res1,np,sa_field,respect_format);
      read_specfile(simfile2,per,res2,np,sa_field,respect_format);
      }

   for(i=0;i<np;i++)
      {
      if(res1[i] != 0.0)
         res1[i] = log(sa1[i]/res1[i]);
      else
         res1[i] = -99;

      if(res2[i] != 0.0)
         res2[i] = log(sa2[i]/res2[i]);
      else
         res2[i] = -99;

      avgr[i] = 0.5*(res1[i]+res2[i]);
      }

   if(pga1 > 0.0 && sim_pga1 > 0.0)
      rpga1 = log(pga1/sim_pga1);
   else
      rpga1 = -999.9;

   if(pga2 > 0.0 && sim_pga2 > 0.0)
      rpga2 = log(pga2/sim_pga2);
   else
      rpga2 = -999.9;

   avgpga = 0.5*(rpga1+rpga2);

   if(pgv1 > 0.0 && sim_pgv1 > 0.0)
      rpgv1 = log(pgv1/sim_pgv1);
   else
      rpgv1 = -999.9;

   if(pgv2 > 0.0 && sim_pgv2 > 0.0)
      rpgv2 = log(pgv2/sim_pgv2);
   else
      rpgv2 = -999.9;

   avgpgv = 0.5*(rpgv1+rpgv2);
   }

if(print_header)
   {
   fprintf(stdout,"%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t\%s","EQ",
                                                               "Mag",
                                                               "stat",
                                                               "lon",
                                                               "lat",
                                                               "stat_seq_no",
                                                               "Vs30",
                                                               "close_dist",
                                                               "Xcos",
                                                               "Ycos",
                                                               "T_min",
                                                               "T_max",
                                                               "comp");
   for(i=0;i<np;i++)
      fprintf(stdout,"\t%.5e",per[i]);

   if(pga1 >= 0.0)
      fprintf(stdout,"\t%s","pga");

   if(pgv1 >= 0.0)
      fprintf(stdout,"\t%s","pgv");

   fprintf(stdout,"\n");
   }

fprintf(stdout,"%s",statinfo);
fprintf(stdout,"\t%s",comp1);
for(i=0;i<np;i++)
   fprintf(stdout,"\t%.5e",res1[i]);

if(pga1 >= 0.0)
   fprintf(stdout,"\t%.5e",rpga1);

if(pgv1 >= 0.0)
   fprintf(stdout,"\t%.5e",rpgv1);

fprintf(stdout,"\n");

fprintf(stdout,"%s",statinfo);
fprintf(stdout,"\t%s",comp2);
for(i=0;i<np;i++)
   fprintf(stdout,"\t%.5e",res2[i]);

if(pga2 >= 0.0)
   fprintf(stdout,"\t%.5e",rpga2);

if(pgv2 >= 0.0)
   fprintf(stdout,"\t%.5e",rpgv2);

fprintf(stdout,"\n");

fprintf(stdout,"%s",statinfo);
fprintf(stdout,"\t%s","avgh");
for(i=0;i<np;i++)
   fprintf(stdout,"\t%.5e",avgr[i]);

if(pga1 >= 0.0 && pga2 >= 0.0)
   fprintf(stdout,"\t%.5e",avgpga);
if(pgv1 >= 0.0 && pgv2 >= 0.0)
   fprintf(stdout,"\t%.5e",avgpgv);

/* this prints avgh of two comps
if(pga1 >= 0.0 && pga2 >= 0.0)
   fprintf(stdout,"\t%.5e",sqrt(pga1*pga2));

if(pgv1 >= 0.0 && pgv2 >= 0.0)
   fprintf(stdout,"\t%.5e",sqrt(pgv1*pgv2));
*/

/* this prints max of both components
if(pga1 >= 0.0 && pga1 >= pga2)
   fprintf(stdout,"\t%.5e",pga1);
else if(pga2 >= 0.0 && pga2 >= pga1)
   fprintf(stdout,"\t%.5e",pga2);

if(pgv1 >= 0.0 && pgv1 >= pgv2)
   fprintf(stdout,"\t%.5e",pgv1);
else if(pgv2 >= 0.0 && pgv2 >= pgv1)
   fprintf(stdout,"\t%.5e",pgv2);
*/

fprintf(stdout,"\n");
 return(0);
}

void *check_malloc(int len)
{
char *ptr;

ptr = (char *) malloc (len);

if(ptr == NULL)
   {
   fprintf(stderr,"*****  memory allocation error\n");
   exit(-1);
   }

return(ptr);
}

FILE *fopfile(char *name,char *mode)
{
FILE *fp;

if((fp = fopen(name,mode)) == NULL)
   {
   fprintf(stderr,"CAN'T FOPEN FILE = %s, MODE = %s\n", name, mode);
   exit(-1);
   }
return(fp);
}

void read_specfile(char *file,float *per,float *sa,int np,int sa_field,int flag)
{
FILE *fpr;
float v[8];
int i, np2;
char str[SLEN];

fpr = fopfile(file,"r");

if(flag)
   {
   fgets(str,SLEN,fpr);
   fgets(str,SLEN,fpr);
   fgets(str,SLEN,fpr);
   fgets(str,SLEN,fpr);

   sscanf(str,"%d",&np2);

   if(np2 != np)
      {
      fprintf(stderr,"No. periods not equal, exiting...\n");
      exit(-1);
      }

   for(i=0;i<np;i++)
      {
      fgets(str,SLEN,fpr);
      sscanf(str,"%f %f %f %f %f %f %f %f",&v[0],&v[1],&v[2],&v[3],&v[4],&v[5],&v[6],&v[7]);
      per[i] = v[1];
      sa[i] = v[sa_field-1];
      }
   }
else
   {
   for(i=0;i<np;i++)
      {
      fgets(str,SLEN,fpr);
      sscanf(str,"%f %f",&per[i],&sa[i]);
      }
   }

fclose(fpr);
}

void read_bbpfile(char *file,float *per,float *sa1,float *sa2,int np)
{
FILE *fpr;
int i, nr;
char str[SLEN];

fpr = fopfile(file,"r");

fgets(str,SLEN,fpr);
while(strncmp(str,"#",1) == 0)
   fgets(str,SLEN,fpr);

for(i=0;i<np;i++)
   {
   nr = sscanf(str,"%f %f %f",&per[i],&sa1[i],&sa2[i]);
   if(nr < 3)
      {
      fprintf(stderr,"Error in file= %s\n",file);
      fprintf(stderr,"found %d columns, expecting at least 3, exiting...\n",nr);
      exit(-1);
      }

   if(i < np-1)
      fgets(str,SLEN,fpr);
   }
fclose(fpr);
}
