#include        "include.h"
#include        "structure.h"
#include        "function.h"

main(int ac,char **av)
{
struct tsheader tshead;
float *val;
int fdr, fdw, len, i;
char infile[512], outfile[512];

int nxts = -1;
int nyts = -1;
int nzts = -1;
int ntts = -1;
float dtts = -1.0;

int swap_bytes = 0;

setpar(ac, av);
mstpar("infile","s",infile);
getpar("swap_bytes","d",&swap_bytes);
endpar();

fdr = opfile_ro(infile);

reed(fdr,&tshead,sizeof(struct tsheader));

if(swap_bytes)
   {
   swap_in_place(1,(char *)(&tshead.ix0));
   swap_in_place(1,(char *)(&tshead.iy0));
   swap_in_place(1,(char *)(&tshead.iz0));
   swap_in_place(1,(char *)(&tshead.it0));
   swap_in_place(1,(char *)(&tshead.nx));
   swap_in_place(1,(char *)(&tshead.ny));
   swap_in_place(1,(char *)(&tshead.nz));
   swap_in_place(1,(char *)(&tshead.nt));
   swap_in_place(1,(char *)(&tshead.dx));
   swap_in_place(1,(char *)(&tshead.dy));
   swap_in_place(1,(char *)(&tshead.dz));
   swap_in_place(1,(char *)(&tshead.dt));
   swap_in_place(1,(char *)(&tshead.modelrot));
   swap_in_place(1,(char *)(&tshead.modellat));
   swap_in_place(1,(char *)(&tshead.modellon));
   }

fprintf(stderr,"ix0= %d\n",tshead.ix0);
fprintf(stderr,"iy0= %d\n",tshead.iy0);
fprintf(stderr,"iz0= %d\n",tshead.iz0);
fprintf(stderr,"it0= %d\n",tshead.it0);
fprintf(stderr,"nx= %d\n",tshead.nx);
fprintf(stderr,"ny= %d\n",tshead.ny);
fprintf(stderr,"nz= %d\n",tshead.nz);
fprintf(stderr,"nt= %d\n",tshead.nt);
fprintf(stderr,"dx= %13.5e\n",tshead.dx);
fprintf(stderr,"dy= %13.5e\n",tshead.dy);
fprintf(stderr,"dz= %13.5e\n",tshead.dz);
fprintf(stderr,"dt= %13.5e\n",tshead.dt);
fprintf(stderr,"modelrot= %13.5f\n",tshead.modelrot);
fprintf(stderr,"modellat= %13.5f\n",tshead.modellat);
fprintf(stderr,"modellon= %13.5f\n",tshead.modellon);

close(fdr);
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
