#include <errno.h>
#include <fcntl.h>
#include <math.h>
#include <stdio.h>
#include <stdarg.h>
#include <unistd.h>
#include <stdlib.h>

#include <sys/file.h>
#include <sys/resource.h>
#include <sys/signal.h>
#include <sys/stat.h>
#include <sys/syscall.h>
#include <sys/time.h>
#include <sys/types.h>

void subr(float *);

main(int ac,char **av)
{
float xlon, lon;
int ffault = 0;

setpar(ac,av);

getpar("ffault","d",&ffault);

if(ffault == 2)
   subr(&xlon);

fprintf(stderr,"xlon= %13.5f lon= %13.5f\n",xlon,lon);

getpar("lon","f",&lon);

endpar();

fprintf(stderr,"xlon= %13.5f lon= %13.5f\n",xlon,lon);
}

void subr(float *x)
{
float lon;

mstpar("lon","f",&lon);

fprintf(stderr,"*x= %13.5f\n",lon);

*x = lon;
}
