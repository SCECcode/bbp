#include "include.h"
#include "structure.h"
#include "function.h"
#include "getpar.h"

#define SLEN 1024

int main(int ac,char **av)
{
FILE *fpr;
float *per1, *sa1, *per2, *sa2, *per3, *sa3;
float v[8];
int ip, np1, np2, np3;

char nsfile[SLEN];
char ewfile[SLEN];
char udfile[SLEN];
char bbfile[SLEN];
char str[SLEN];

char nsname[128];
char ewname[128];
char udname[128];
char units[128];

char tformat[16];
char dformat[16];
char frmt[256];

char title[128];

int sa_column = 7; /* 3=rel disp, 4=rel vel, 5=pseudo rel vel, 6=absolute acc, 7=pseudo absolute acc; NGA=7 */

title[0] = '\0';

nsname[0] = '\0';
ewname[0] = '\0';
udname[0] = '\0';
units[0] = '\0';

sprintf(tformat,"%%.6e");
sprintf(dformat,"\t%%.6e");
sprintf(bbfile,"stdout");

setpar(ac,av);

mstpar("nsfile","s",nsfile);
mstpar("ewfile","s",ewfile);
mstpar("udfile","s",udfile);

getpar("nsname","s",nsname);
getpar("ewname","s",ewname);
getpar("udname","s",udname);

getpar("units","s",units);

getpar("tformat","s",tformat);
getpar("dformat","s",dformat);

getpar("bbfile","s",bbfile);
getpar("title","s",title);
getpar("sa_column","d",&sa_column);

endpar();

fpr = fopfile(nsfile,"r");

fgets(str,SLEN,fpr);
fgets(str,SLEN,fpr);
fgets(str,SLEN,fpr);
fgets(str,SLEN,fpr);

sscanf(str,"%d",&np1);

per1 = (float *)check_malloc(np1*sizeof(float));
sa1 = (float *)check_malloc(np1*sizeof(float));

for(ip=0;ip<np1;ip++)
   {
   fgets(str,SLEN,fpr);
   sscanf(str,"%f %f %f %f %f %f %f %f",&v[0],&v[1],&v[2],&v[3],&v[4],&v[5],&v[6],&v[7]);
   per1[ip] = v[1];
   sa1[ip] = v[sa_column-1];
   }
fclose(fpr);

fpr = fopfile(ewfile,"r");

fgets(str,SLEN,fpr);
fgets(str,SLEN,fpr);
fgets(str,SLEN,fpr);
fgets(str,SLEN,fpr);

sscanf(str,"%d",&np2);

per2 = (float *)check_malloc(np2*sizeof(float));
sa2 = (float *)check_malloc(np2*sizeof(float));

for(ip=0;ip<np2;ip++)
   {
   fgets(str,SLEN,fpr);
   sscanf(str,"%f %f %f %f %f %f %f %f",&v[0],&v[1],&v[2],&v[3],&v[4],&v[5],&v[6],&v[7]);
   per2[ip] = v[1];
   sa2[ip] = v[sa_column-1];
   }
fclose(fpr);

fpr = fopfile(udfile,"r");

fgets(str,SLEN,fpr);
fgets(str,SLEN,fpr);
fgets(str,SLEN,fpr);
fgets(str,SLEN,fpr);

sscanf(str,"%d",&np3);

per3 = (float *)check_malloc(np3*sizeof(float));
sa3 = (float *)check_malloc(np3*sizeof(float));

for(ip=0;ip<np3;ip++)
   {
   fgets(str,SLEN,fpr);
   sscanf(str,"%f %f %f %f %f %f %f %f",&v[0],&v[1],&v[2],&v[3],&v[4],&v[5],&v[6],&v[7]);
   per3[ip] = v[1];
   sa3[ip] = v[sa_column-1];
   }
fclose(fpr);

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
   fprintf(fpr,"#  period(sec)      %s      %s      %s\n",nsname,ewname,udname);

sprintf(frmt,"%s%s%s%s\n",tformat,dformat,dformat,dformat);

for(ip=0;ip<np1;ip++)
   fprintf(fpr,frmt,per1[ip],sa1[ip],sa2[ip],sa3[ip]);

fclose(fpr);
}
