#include        "include.h"
#include        "structure.h"
#include        "function.h"

float float_swap(char *);
long  long_swap(char *);

main(int ac, char **av)
{
int fdr, fdw, ntswap, i, ndat;
float *fbuf;
char *cbuf;
char infile[128];
char outfile[128];

int nt = -1;

setpar(ac, av);
mstpar("infile","s",infile);
mstpar("outfile","s",outfile);
getpar("nt","d",&nt);
endpar();

fdr = opfile_ro(infile);

if(nt < 0)
   {
   reed(fdr,&nt,sizeof(int));

   ndat = 8*nt + 6;
   cbuf = (char *) check_malloc(ndat*sizeof(float));
   fbuf = (float *) check_malloc(ndat*sizeof(float));
   }
else
   {
   ndat = 8*nt + 6;
   cbuf = (char *) check_malloc(ndat*sizeof(float));
   fbuf = (float *) check_malloc(ndat*sizeof(float));

   reed(fdr,&nt,sizeof(int));
   }

reed(fdr,cbuf,ndat*sizeof(float));

close(fdr);

ntswap = long_swap((char *)(&nt));

for(i=0;i<ndat;i++)
   fbuf[i] = float_swap(cbuf+(i*4));

fdw = croptrfile(outfile);

rite(fdw,&ntswap,sizeof(int));
rite(fdw,fbuf,ndat*sizeof(float));

close(fdw);
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
