#include <errno.h>
#include <fcntl.h>
#include <math.h>
#include <stdio.h>
#include <unistd.h>
#include <stdlib.h>

#include <sys/file.h>
#include <sys/procfs.h>
#include <sys/resource.h>
#include <sys/signal.h>
#include <sys/stat.h>
#include <sys/syscall.h>
#include <sys/time.h>
#include <sys/types.h>

#define         ERAD            6378.139
#define         RPERD           0.017453292

main(int ac,char **av)
{
FILE *fpw, *fpr, *fopfile();
float xlen, ylen, alpha, ref_lon, ref_lat, rlon, rlat;
float dx, dy, xf, yf;

double g0, b0;
double amat[9], ainv[9];

float ref_rad = ERAD;
int input_xy = 1;

setpar(ac,av);
mstpar("ref_lon","f",&ref_lon);
mstpar("ref_lat","f",&ref_lat);
mstpar("alpha","f",&alpha);
mstpar("xlen","f",&xlen);
mstpar("ylen","f",&ylen);
getpar("ref_rad","f",&ref_rad);
getpar("input_xy","d",&input_xy);
endpar();

gen_matrices(amat,ainv,&alpha,&ref_lon,&ref_lat);

g0 = (double)(0.5*ylen)/(double)(ref_rad);
b0 = (double)(0.5*xlen)/(double)(ref_rad);

/*
fprintf(stderr,"erad= %f g0= %f b0= %f\n",ref_rad,g0,b0);
fprintf(stderr,"a0= %f a1= %f a2= %f\n",amat[0],amat[1],amat[2]);
fprintf(stderr,"a0= %f a1= %f a2= %f\n",amat[3],amat[4],amat[5]);
fprintf(stderr,"a0= %f a1= %f a2= %f\n",amat[6],amat[7],amat[8]);
fprintf(stderr,"a0= %f a1= %f a2= %f\n",ainv[0],ainv[1],ainv[2]);
fprintf(stderr,"a0= %f a1= %f a2= %f\n",ainv[3],ainv[4],ainv[5]);
fprintf(stderr,"a0= %f a1= %f a2= %f\n",ainv[6],ainv[7],ainv[8]);
*/

if(input_xy)
   {
   while(scanf("%f %f",&xf,&yf) == 2)
      {
      gcproj(&xf,&yf,&rlon,&rlat,&ref_rad,&g0,&b0,amat,ainv,0);
      printf("%.5f %.5f\n",rlon,rlat);

      gcproj(&xf,&yf,&rlon,&rlat,&ref_rad,&g0,&b0,amat,ainv,1);
      printf("%.5f %.5f\n",xf,yf);
      }
   }
else
   {
   while(scanf("%f %f",&rlon,&rlat) == 2)
      {
      gcproj(&xf,&yf,&rlon,&rlat,&ref_rad,&g0,&b0,amat,ainv,1);
      printf("%.5f %.5f\n",xf,yf);

      gcproj(&xf,&yf,&rlon,&rlat,&ref_rad,&g0,&b0,amat,ainv,0);
      printf("%.5f %.5f\n",rlon,rlat);
      }
   }
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

void *check_realloc(void *ptr,size_t len)
{
ptr = (char *) realloc (ptr,len);

if(ptr == NULL)
   {
   fprintf(stderr,"*****  memory reallocation error\n");
   exit(-1);
   }

return(ptr);
}

void *check_malloc(size_t len)
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
