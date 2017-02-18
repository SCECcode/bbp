#include "include.h"
#include "structure.h"
#include "function.h"
#include "getpar.h"

int size_float = sizeof(float);
int size_int = sizeof(int);

int main(ac,av)
int ac;
char **av;
{
struct statdata head1;
float *s1, amax;
int i;
char infile[128];

float max = -1.0e+20;
float min = 1.0e+20;
int inbin = 0;
int outbin = 0;
int keepsign = 0;
float scale = 1.0;

sprintf(infile,"stdin");

setpar(ac,av);
getpar("infile","s",infile);
getpar("inbin","d",&inbin);
getpar("keepsign","d",&keepsign);
getpar("scale","f",&scale);
endpar();
 
s1 = NULL;
s1 = read_wccseis(infile,&head1,s1,inbin);

for(i=0;i<head1.nt;i++)
   {
   if(s1[i] > max)
      max = s1[i];
   if(s1[i] < min)
      min = s1[i];
   }
 
if(max >= 0.0 && min < 0.0)
   {
   if(max > -min)
      amax = max;
   else
      {
      amax = -min;
      if(keepsign)
	 amax = -amax;
      }
   }

else if(max >= 0.0 && min >= 0.0)
   amax = max;

else if(max < 0.0 && min < 0.0)
   {
   amax = -min;
   if(keepsign)
      amax = -amax;
   }

printf("%10.2f %13.5e %s\n",head1.edist,scale*amax,head1.stat);
}
