#include <errno.h>
#include <fcntl.h>
#include <math.h>
#include <stdio.h>
#include <unistd.h>

#include <sys/fault.h>
#include <sys/file.h>
#include <sys/procfs.h>
#include <sys/resource.h>
#include <sys/signal.h>
#include <sys/stat.h>
#include <sys/syscall.h>
#include <sys/time.h>
#include <sys/types.h>

char *check_malloc();

struct coords
   {
   int xp;
   int yp;
   };

main(ac,av)
int ac;
char **av;
{
struct coords *points;
int ix, iy;
int x0, y0, x1, y1;
int ip, np, i;
float x, y, r, cosA, sinA;
char str[64], name[32];

sprintf(name,"s");

setpar(ac, av);
mstpar("x0","d",&x0);
mstpar("y0","d",&y0);
mstpar("x1","d",&x1);
mstpar("y1","d",&y1);
getpar("name","s",name);
endpar();

x = (float)(x1 - x0);
y = (float)(y1 - y0);
r = sqrt(x*x + y*y);
np = 1 + (int)(r + 0.5);

cosA = x/r;
sinA = y/r;

points = (struct coords *) check_malloc (np*sizeof(struct coords));

points[0].xp = x0;
points[0].yp = y0;
i = 0;
for(ip=1;ip<np;ip++)
   {
   r = (float)(ip);
   ix = x0 + (int)(r*cosA + 0.5);
   iy = y0 + (int)(r*sinA + 0.5);

   if(ix != points[i].xp && iy != points[i].yp)
      {
      i++;
      points[i].xp = ix;
      points[i].yp = iy;
      }
   }

np = i + 1;
for(ip=0;ip<np;ip++)
   {
   if(ip < 10)
      sprintf(str,"%s00%1d",name,ip);
   else if(ip < 100)
      sprintf(str,"%s0%2d",name,ip);
   else
      sprintf(str,"%s%3d",name,ip);

   printf("%5d %5d %s\n",points[ip].xp,points[ip].yp,str);
   }
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
