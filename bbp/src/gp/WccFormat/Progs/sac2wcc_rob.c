#include "include.h"
#include "structure.h"
#include "function.h"
#include "sac.h"

int size_float = sizeof(float);
int size_int = sizeof(int);

int main(int ac, char **av)
{
struct statdata head1;
float *s1, time_sec;
SACHEAD sachd;
int i, j, nt1;

char infile[128];
char outfile[128];

char stat[16];
char comp[4];
char title[128];

int event_yr = -999;
int event_jday = -999;
int event_hr = -999;
int event_min = -999;
double event_sec = -1.0e+15;

int use_header_cmpaz = 0;
int outbin = 0;

int swap_bytes = 0;

stat[0] = '\0';
comp[0] = '\0';
title[0] = '\0';

setpar(ac,av);
mstpar("infile","s",infile);
mstpar("outfile","s",outfile);
getpar("swap_bytes","d",&swap_bytes);
getpar("stat","s",stat);
getpar("comp","s",comp);
getpar("title","s",title);
getpar("outbin","d",&outbin);
getpar("event_yr","d",&event_yr);
getpar("event_jday","d",&event_jday);
getpar("event_hr","d",&event_hr);
getpar("event_min","d",&event_min);
getpar("event_sec","F",&event_sec);
getpar("use_header_cmpaz","d",&use_header_cmpaz);
endpar();

if((s1 = read_sac(infile,&sachd,swap_bytes)) == NULL)
   {
   exit(-1);
   }

if(stat[0] != '\0')
   strcpy(head1.stat,stat);

if(comp[0] != '\0')
   strcpy(head1.comp,comp);

if(use_header_cmpaz != 0)
   sprintf(head1.comp,"%3.0f",sachd.cmpaz);

if(title[0] != '\0')
   strcpy(head1.stitle,title);

if(event_yr < 0)
   event_yr = sachd.nzyear;
if(event_jday < 0)
   event_jday = sachd.nzjday;
if(event_hr < 0)
   event_hr = sachd.nzhour;
if(event_min < 0)
   event_min = sachd.nzmin;
if(event_sec < -1.0e+10)
   event_sec = (double)(sachd.nzsec) + 0.001*((double)(sachd.nzmsec));

head1.nt = sachd.npts;
head1.dt = sachd.delta;

time_sec = 365*24*3600.0*(sachd.nzyear - event_yr)
	      + 24*3600.0*(sachd.nzjday - event_jday)
	      + 3600.0*(sachd.nzhour - event_hr)
	      + 60.0*(sachd.nzmin - event_min)
	      + (sachd.nzsec - event_sec)
	      + 0.001*sachd.nzmsec
	      + sachd.b;

/*
fprintf(stderr,"b=%g nzsec=%d nzmsec=%d esec=%.20f t=%.20f\n",sachd.b,sachd.nzsec,sachd.nzmsec,event_sec,time_sec);
*/

head1.hr = (int)(time_sec/3600.0);
head1.min = (int)(time_sec/60.0) - 60*head1.hr;
head1.sec = time_sec - 3600.0*head1.hr - 60.0*head1.min;

head1.edist = sachd.dist;
head1.az = 0.0;
head1.baz = 0.0;

write_wccseis(outfile,&head1,s1,outbin);
}
