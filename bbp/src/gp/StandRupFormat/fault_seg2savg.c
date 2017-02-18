#include "include.h"
#include "structure.h"
#include "function.h"
#include "defs.h"

void moment_sum(struct velmodel *vm,float *dtop,float *dip,float *flen,float *dwid,int nwid,float *msum);

main(int ac,char **av)
{
FILE *fopfile(), *fpr, *fpw;
float dip, dtop, flen, fwid, msum, dwid;
int i, nseg, nlen, nwid;

char *rchar, infile[512], str[1024], velfile[512];
struct velmodel vmod;

float moment = -1.0;
float mag = -1.0;

sprintf(infile,"stdin");
sprintf(velfile,"NOT_PROVIDED");

setpar(ac, av);
getpar("infile","s",infile);
getpar("velfile","s",velfile);
getpar("moment","f",&moment);
if(moment < 0.0)
   mstpar("mag","f",&mag);
endpar();

if(mag > 0.0)
   moment = exp(1.5*(mag+10.7)*log(10.0));

read_Fvelmodel(velfile,&vmod);

if(strcmp(infile,"stdin") == 0)
   fpr = stdin;
else
   fpr = fopfile(infile,"r");

fgets(str,1024,fpr);
while(strncmp(str,"#",1) == 0)
   {
   rchar = fgets(str,1024,fpr);
   if(rchar == NULL)
      {
      fprintf(stderr,"Unexpected EOF in %s, exiting...\n",infile);
      exit(-99);
      }
   }

sscanf(str,"%d",&nseg);

msum = 0.0;
for(i=0;i<nseg;i++)
   {
   fgets(str,1024,fpr);
   sscanf(str,"%f %f %f %f %d %d",&dtop,&dip,&flen,&fwid,&nlen,&nwid);

   dwid = fwid/nwid;
   moment_sum(&vmod,&dtop,&dip,&flen,&dwid,nwid,&msum);
   }
fclose(fpr);

fprintf(stdout,"average_slip= %.2f moment= %13.5e msum= %13.5e\n",moment/msum,moment,msum);
}

void moment_sum(struct velmodel *vm,float *dtop,float *dip,float *flen,float *dwid,int nwid,float *msum)
{
int j, iw;
float sinD, dep, tmom;

float pi = 3.141592654;

sinD = sin((*dip)*pi/180.0);

tmom = 0.0;
for(iw=0;iw<nwid;iw++)
   {
   dep = *dtop + (iw+0.5)*(*dwid)*sinD;

   j = 0;
   while(vm->dep[j] < dep && j < vm->nlay)
      j++;

   tmom = tmom + (*dwid)*vm->mu[j];
   }

*msum = *msum + tmom*(*flen)*1.0e+10;
}
