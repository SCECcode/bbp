#include        "include.h"
#include        "structure.h"
#include        "function.h"

main(int ac, char **av)
{
struct statdata gfh;
float *gf;
int k, isub, j, nc, nt;
int nbyte;
float dt, rng, tst;

int fdr;
char infile[128];
char name[128];
char outfile[128];

int outbin = 0;

static char *comp[] = {"000","090","ver"};

gfh.stitle[0] = '\0';
gfh.hr = 0;
gfh.min = 0;
gfh.sec = 0.0;
gfh.edist = 0.0;
gfh.az = 0.0;
gfh.baz = 0.0;

setpar(ac, av);
mstpar("infile","s",infile);
mstpar("isub","d",&isub);
getpar("outbin","d",&outbin);
endpar();

sprintf(gfh.stat,"%d",isub);

fdr = opfile_ro(infile);

for(k=0;k<isub;k++)
   {
   reed(fdr,&nbyte,sizeof(int));
   reed(fdr,&nc,sizeof(int));
   reed(fdr,&nbyte,sizeof(int));

   for(j=0;j<nc;j++)
      {
      reed(fdr,&nbyte,sizeof(int));
      reed(fdr,&rng,sizeof(float));
      reed(fdr,&tst,sizeof(float));
      reed(fdr,&nbyte,sizeof(int));

      reed(fdr,&nbyte,sizeof(int));
      reed(fdr,&nt,sizeof(int));
      reed(fdr,&dt,sizeof(float));
      reed(fdr,&nbyte,sizeof(int));

      gf = (float *) check_realloc(gf,nt*sizeof(float));

      reed(fdr,&nbyte,sizeof(int));
      reed(fdr,gf,nt*sizeof(float));
      reed(fdr,&nbyte,sizeof(int));
      }
   }

reed(fdr,&nbyte,sizeof(int));
reed(fdr,&nc,sizeof(int));
reed(fdr,&nbyte,sizeof(int));

for(j=0;j<nc;j++)
   {
   reed(fdr,&nbyte,sizeof(int));
   reed(fdr,&rng,sizeof(float));
   reed(fdr,&tst,sizeof(float));
   reed(fdr,&nbyte,sizeof(int));

   reed(fdr,&nbyte,sizeof(int));
   reed(fdr,&nt,sizeof(int));
   reed(fdr,&dt,sizeof(float));
   reed(fdr,&nbyte,sizeof(int));

   gf = (float *) check_realloc(gf,nt*sizeof(float));

   reed(fdr,&nbyte,sizeof(int));
   reed(fdr,gf,nt*sizeof(float));
   reed(fdr,&nbyte,sizeof(int));

   gfh.nt = nt;
   gfh.dt = dt;
   gfh.edist = rng;

   gfh.hr = (int)(tst/3600.0);
   gfh.min = (int)((tst - 3600.0*gfh.hr)/60.0);
   gfh.sec = (float)(tst - 60.0*gfh.min - 3600.0*gfh.hr);

   strcpy(gfh.comp,comp[j]);
   sprintf(outfile,"%s.%s",gfh.stat,gfh.comp);

   write_wccseis(outfile,&gfh,gf,outbin);
   }
close(fdr);
}
