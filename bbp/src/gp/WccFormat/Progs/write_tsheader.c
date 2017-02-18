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
mstpar("outfile","s",outfile);
getpar("nxts","d",&nxts);
getpar("nyts","d",&nyts);
getpar("nzts","d",&nzts);
getpar("ntts","d",&ntts);
getpar("dtts","f",&dtts);
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

if(nxts > 0)
   tshead.nx = nxts;
if(nyts > 0)
   tshead.ny = nyts;
if(nzts > 0)
   tshead.nz = nzts;

if(ntts > 0)
   tshead.nt = ntts;
else
   ntts = tshead.nt;

if(dtts > 0)
   tshead.dt = dtts;

fprintf(stderr,"nxts= %d nyts= %d nzts= %d ntts= %d dtts= %13.6e\n",nxts,nyts,nzts,ntts,dtts);
fprintf(stderr,"hdnx= %d hdny= %d hdnz= %d hdnt= %d hddt= %13.5e\n",tshead.nx,tshead.ny,tshead.nz,tshead.nt,tshead.dt);

len = 3*tshead.nx*tshead.ny*tshead.nz*sizeof(float);
val = (float *) check_malloc (len);

fdw = croptrfile(outfile);

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

rite(fdw,&tshead,sizeof(struct tsheader));

for(i=0;i<ntts;i++)
   {
   fprintf(stderr,"%5d of %5d\n",i+1,ntts);

   reed(fdr,val,len);
   rite(fdw,val,len);
   }

close(fdr);
close(fdw);
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
