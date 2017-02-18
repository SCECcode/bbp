#include        "include.h"
#include        "structure.h"
#include        "function.h"

int size_float = sizeof(float);
float float_swap(char *);

main(int ac, char **av)
{
FILE *fpr, *fopfile();
struct seisheader seishead;
float *s;
int i, ix;

int swap_bytes = 0;
char cbuf[512];

int outbin = 0;
int ngf;
int fd;
char infile[128];
char outfile[128];

setpar(ac, av);
mstpar("infile","s",infile);
mstpar("outfile","s",outfile);
getpar("swap_bytes","d",&swap_bytes);
endpar();

fd = opfile_ro(infile);
fpr = fopfile(outfile,"w");

reed(fd,&ngf,sizeof(int));
if(swap_bytes)
   swap_in_place(1,(char *)(&ngf));

fprintf(fpr,"%d\n",ngf);

i = 0;
while(i<ngf)
   {
   reed(fd,&seishead,sizeof(struct seisheader));

   if(swap_bytes)
      {
      swap_in_place(1,(char *)(&seishead.indx));
      swap_in_place(1,(char *)(&seishead.ix));
      swap_in_place(1,(char *)(&seishead.iy));
      swap_in_place(1,(char *)(&seishead.iz));
      swap_in_place(1,(char *)(&seishead.nt));
      swap_in_place(1,(char *)(&seishead.dt));
      swap_in_place(1,(char *)(&seishead.h));
      swap_in_place(1,(char *)(&seishead.modelrot));
      swap_in_place(1,(char *)(&seishead.modellat));
      swap_in_place(1,(char *)(&seishead.modellon));
      }

   fprintf(fpr,"%d\n",seishead.indx);
   fprintf(fpr,"%d\n",seishead.ix);
   fprintf(fpr,"%d\n",seishead.iy);
   fprintf(fpr,"%d\n",seishead.iz);
   fprintf(fpr,"%d\n",seishead.nt);
   fprintf(fpr,"%f\n",seishead.dt);
   fprintf(fpr,"%f\n",seishead.h);
   fprintf(fpr,"%f\n",seishead.modelrot);
   fprintf(fpr,"%f\n",seishead.modellat);
   fprintf(fpr,"%f\n",seishead.modellon);
   fprintf(fpr,"%s\n",seishead.name);

   i++;
   }

s = (float *)check_malloc(3*ngf*seishead.nt*sizeof(float));

reed(fd,s,3*ngf*seishead.nt*sizeof(float));

if(swap_bytes)
   swap_in_place(3*ngf*seishead.nt,(char *)(s));

for(i=1;i<3*ngf*seishead.nt;i=i+3)
   fprintf(fpr,"%13.5e\n",s[i]);

close(fd);
}

long long_swap(char *cbuf)
{
union
   {
   char cval[4];
   long lval;
   } l_union;

l_union.cval[3] = cbuf[0];
l_union.cval[2] = cbuf[1];
l_union.cval[1] = cbuf[2];
l_union.cval[0] = cbuf[3];

return(l_union.lval);
}

float float_swap(char *cbuf)
{
union
   {
   char cval[4];
   float fval;
   } f_union;

f_union.cval[3] = cbuf[0];
f_union.cval[2] = cbuf[1];
f_union.cval[1] = cbuf[2];
f_union.cval[0] = cbuf[3];

return(f_union.fval);
}

void swap_in_place(int n,char *cbuf)
{
char cv;

while(n--)
   {
   cv = cbuf[0];
   cbuf[0] = cbuf[3];
   cbuf[3] = cv;

   cv = cbuf[1];
   cbuf[1] = cbuf[2];
   cbuf[2] = cv;

   cbuf = cbuf + 4;
   }
}
