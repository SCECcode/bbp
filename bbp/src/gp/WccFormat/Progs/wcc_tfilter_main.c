/**********************************************************************/
/*                                                                    */
/*           wcc_tfilter - written by RWG 03/92                       */
/*                                                                    */
/*           N-th order high-pass, low-pass and band-pass             */
/*           Butterworth filters for time series.                     */
/*                                                                    */
/**********************************************************************/

#include "include.h"
#include "structure.h"
#include "function.h"
#include "getpar.h"

#define 	PI 		3.14159265
#define         MAXFILES        500

void hp_filter(struct complex* q,struct complex* p,int n,float* alpha,float* beta,int sgn);
void lp_filter(struct complex* q,struct complex* p,int n,float* alpha,float* beta,int sgn);
void czero(struct complex* p,int n);

int main (ac, av)
int ac;
char **av;
{
char infilebuf[MAXFILES*256];
char outfilebuf[MAXFILES*256];

struct statdata shead1;
char *infile[MAXFILES], *outfile[MAXFILES], *string, *readline(), filelist[256];
char str[512];
struct complex *q, *p;
float *s1;
int j, it, i;
int nshft, inchar, outchar, nstat, nt6;

int inbin = 0;
int outbin = 0;

FILE *fopfile(), *fmake_or_open(), *fpr, *fpw;
char outpath[256];

if(ac == 1)
   {
   printf("**** %s:\n",av[0]);
   printf("     'flo' is low pass frequency cutoff (pass for f < flo)\n");
   printf("     'fhi' is high pass frequency cutoff (pass for f > fhi)\n\n");
   printf("     This means that the passband will be:  'fhi' < f < 'flo'.\n\n");
   printf("     For a high-pass only filter, set fhi=f0, flo=1.0e+10;\n");
   printf("     for a low-pass only filter, set fhi=0.0, flo=f0;\n");
   printf("     where 'f0' is the cutoff frequency in either case.\n");
   exit(1);
   }

setpar(ac, av);
mstpar("filelist","s",filelist);
mstpar("outpath","s",outpath);
getpar("inbin","d",&inbin);
getpar("outbin","d",&outbin);
endpar();

/*  read in input data filenames */
 
fpr = fopfile(filelist,"r");

i = 0;                       
infile[0] = infilebuf;
outfile[0] = outfilebuf;
string = readline(fpr);
while(string != NULL)
   {
   nshft = 0;
   inchar = getname(string,infile[i],&nshft);
   outchar = getname(&string[nshft],outfile[i],&nshft);

   if(outfile[i][0] == '\0')
      {
      j = 0;                                              
      while(infile[i][j] != '\0')
         j++;
 
      outchar = 0;
      while(j >= 0 && infile[i][j] != '/')
	 {
         j--;
	 outchar++;
	 }
      j++;  
 
      strcpy(outfile[i],infile[i]+j);
      }

   i++;
   if(i == MAXFILES)
      break;
   else
      {
      infile[i] = infile[i-1] + inchar;
      outfile[i] = outfile[i-1] + outchar;
      }

   string = readline(fpr);
   }
nstat = i;
fclose(fpr);
fflush(stdout);
if(!nstat) {
  printf("No entries in file %s.\n");
  fflush(stdout);
  exit(0);
}

s1 = NULL;
p = NULL;
q = NULL;
for(i=0;i<nstat;i++)
   {
   s1 = read_wccseis(infile[i],&shead1,s1,inbin);

   wcc_tfilter(ac, av, s1, &shead1);

   /* make sure output directory exists */
 
   makedir(outpath);

   set_fullpath(str,outpath,outfile[i]);
   printf("Writing to %s.\n", str);
   write_wccseis(str,&shead1,s1,outbin);
   }
   return(0);
}

