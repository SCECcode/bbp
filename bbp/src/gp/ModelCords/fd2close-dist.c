#include <errno.h>
#include <fcntl.h>
#include <math.h>
#include <stdio.h>
#include <unistd.h>

#include <sys/file.h>
#include <sys/procfs.h>
#include <sys/resource.h>
#include <sys/signal.h>
#include <sys/stat.h>
#include <sys/syscall.h>
#include <sys/time.h>
#include <sys/types.h>

char *check_malloc();

struct fc
   {
   float x;
   float y;
   float z;
   };

struct sc
   {
   float x;
   float y;
   float z;
   char name[16];
   };

main(ac,av)
int ac;
char **av;
{
FILE *fpw, *fpr, *fopfile();
struct sc *sc;
struct fc *fc;
float xx, yy, zz, r, cdmin;
int i, j, ns, nf;
char outfile[128], statfile[128], faultfile[128], str[512];

setpar(ac,av);
mstpar("statfile","s",statfile);
mstpar("faultfile","s",faultfile);
mstpar("outfile","s",outfile);
endpar();

fpr = fopfile(statfile,"r");
fscanf(fpr,"%d",&ns);

sc = (struct sc *) check_malloc (ns*sizeof(struct sc));

for(i=0;i<ns;i++)
   fscanf(fpr,"%f %f %f %s",&sc[i].x,&sc[i].y,&sc[i].z,sc[i].name);

fclose(fpr);

fpr = fopfile(faultfile,"r");
fscanf(fpr,"%d",&nf);

fc = (struct fc *) check_malloc (nf*sizeof(struct fc));

for(j=0;j<nf;j++)
   {
   fgets(str,512,fpr);
   sscanf(str,"%f %f %f",&fc[j].x,&fc[j].y,&fc[j].z);
   }

fclose(fpr);

fpw = fopfile(outfile,"w");

for(i=0;i<ns;i++)
   {
   cdmin = 1.0e+10;
   for(j=0;j<nf;j++)
      {
      xx = fc[j].x - sc[i].x;
      yy = fc[j].y - sc[i].y;
      zz = fc[j].z - sc[i].z;

      r = sqrt(xx*xx + yy*yy + zz*zz);
      if(r < cdmin)
	 cdmin = r;
      }

   fprintf(fpw,"%13.5e %s\n",cdmin,sc[i].name);
   }

fclose(fpw);
}

FILE *fopfile(name,mode)
char *name, *mode;
{
FILE *fp, *fopen();

if((fp = fopen(name,mode)) == NULL)
   {
   fprintf(stderr,"CAN'T FOPEN FILE = %s, MODE = %s\n", name, mode);
   exit(-1);
   }
return(fp);
}

char *check_malloc(len)
int len;
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
