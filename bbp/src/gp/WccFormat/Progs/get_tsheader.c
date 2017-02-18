#include        "include.h"
#include        "structure.h"
#include        "function.h"

int size_float = sizeof(float);
float float_swap(char *);
void split_bytes(char *cbuf,char *c1,char *c2);

main(int ac, char **av)
{
struct tsheader tsh;
float s1[10];
int i, it, ixp, iyp;
int np1_byte, off1, off2, off3;
int fdr, fdw1, fdw2, fdw3;
char in_tsfile[128], out_tsfile[128];
off_t off;

int swap_bytes = 0;
int inbin = 0;
int zero_tsfile = 0;
int nt = -1;

it = -1;

setpar(ac, av);
mstpar("in_tsfile","s",in_tsfile);
getpar("it","d",&it);
endpar();

fdr = opfile_ro(in_tsfile);
reed(fdr,&tsh,sizeof(struct tsheader));

if(swap_bytes)
   {
   swap_in_place(1,(char *)(&tsh.ix0));
   swap_in_place(1,(char *)(&tsh.iy0));
   swap_in_place(1,(char *)(&tsh.iz0));
   swap_in_place(1,(char *)(&tsh.it0));
   swap_in_place(1,(char *)(&tsh.nx));
   swap_in_place(1,(char *)(&tsh.ny));
   swap_in_place(1,(char *)(&tsh.nz));
   swap_in_place(1,(char *)(&tsh.nt));
   swap_in_place(1,(char *)(&tsh.dx));
   swap_in_place(1,(char *)(&tsh.dy));
   swap_in_place(1,(char *)(&tsh.dz));
   swap_in_place(1,(char *)(&tsh.dt));
   swap_in_place(1,(char *)(&tsh.modelrot));
   swap_in_place(1,(char *)(&tsh.modellat));
   swap_in_place(1,(char *)(&tsh.modellon));
   }

fprintf(stderr,"ix= %d iy= %d iz= %d it= %d\n",tsh.ix0,tsh.iy0,tsh.iz0,tsh.it0);
fprintf(stderr,"nx= %d ny= %d nz= %d nt= %d\n",tsh.nx,tsh.ny,tsh.nz,tsh.nt);
fprintf(stderr,"dx= %12.4e dy= %12.4e dz= %12.4e dt= %12.4e\n",tsh.dx,tsh.dy,tsh.dz,tsh.dt);
fprintf(stderr,"lon= %10.4f lat= %10.4f rot= %10.4f\n",tsh.modellon,tsh.modellat,tsh.modelrot);

if(it<0)
   off = (3*tsh.nx*tsh.ny*tsh.nt-5)*sizeof(float);
else
   off = it*sizeof(float);

fprintf(stderr,"off= %20.0f\n",1.0*off);
fprintf(stderr,"off= %lld\n",off);
/*
split_bytes((char *)(&off),(char *)(&off1),(char *)(&off2));
*/
fprintf(stderr,"off2= %d off1= %d\n",off2,off1);
fprintf(stderr,"size= %d\n",sizeof(off_t));
fprintf(stderr,"size= %d\n",sizeof(off));

lseek(fdr,off,SEEK_CUR);

reed(fdr,s1,10*sizeof(float));
close(fdr);

for(i=0;i<10;i++)
   fprintf(stderr,"s1[%d]= %13.5e\n",i,s1[i]);
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

void split_bytes(char *cbuf,char *c1,char *c2)
{
c1[0] = cbuf[0];
c1[1] = cbuf[1];
c1[2] = cbuf[2];
c1[3] = cbuf[3];

c2[0] = cbuf[4];
c2[1] = cbuf[5];
c2[2] = cbuf[6];
c2[3] = cbuf[7];
}
