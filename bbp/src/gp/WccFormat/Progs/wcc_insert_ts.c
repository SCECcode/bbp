#include        "include.h"
#include        "structure.h"
#include        "function.h"

int size_float = sizeof(float);
float float_swap(char *);

main(int ac, char **av)
{
struct statdata head1, head2, head3;
struct tsheader tshead;
float *s1, *s2, *s3;
int i, it, ixp, iyp;
off_t np1_byte, off1, off2, off3;
int fdr, fdw1, fdw2, fdw3;
char seisfile1[128], seisfile2[128], seisfile3[128];
char in_tsfile[128], out_tsfile[128];

int swap_bytes = 0;
int inbin = 0;
int zero_tsfile = 0;
float dt = -1.0;
int nt = -1;

setpar(ac, av);

mstpar("out_tsfile","s",out_tsfile);
getpar("inbin","d",&inbin);
getpar("zero_tsfile","d",&zero_tsfile);
getpar("swap_bytes","d",&swap_bytes);
getpar("nt","d",&nt);
getpar("dt","f",&dt);

if(zero_tsfile == 0)
   {
   mstpar("seisfile1","s",seisfile1);
   mstpar("seisfile2","s",seisfile2);
   mstpar("seisfile3","s",seisfile3);

   mstpar("ixp","d",&ixp);
   mstpar("iyp","d",&iyp);
   }
else
   mstpar("in_tsfile","s",in_tsfile);

endpar();

if(zero_tsfile == 1)
   {
   fdr = opfile_ro(in_tsfile);
   reed(fdr,&tshead,sizeof(struct tsheader));
   close(fdr);

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

   if(nt > 0)
      tshead.nt = nt;
   if(dt > 0)
      tshead.dt = dt;

   s1 = (float *) check_malloc(3*tshead.nt*size_float);
   for(it=0;it<tshead.nt;it++)
      s1[it] = 0.0;

   fdw1 = croptrfile(out_tsfile);

   fprintf(stderr,"nx= %d ny= %d nt= %d\n",tshead.nx,tshead.ny,tshead.nt);

   rite(fdw1,&tshead,sizeof(struct tsheader));
   for(i=0;i<tshead.nx*tshead.ny*tshead.nz;i++)
      rite(fdw1,s1,3*tshead.nt*size_float);

   close(fdw1);
   }
else
   {
   fdw1 = opfile(out_tsfile);
   reed(fdw1,&tshead,sizeof(struct tsheader));

   fdw2 = opfile(out_tsfile);
   reed(fdw2,&tshead,sizeof(struct tsheader));

   fdw3 = opfile(out_tsfile);
   reed(fdw3,&tshead,sizeof(struct tsheader));

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

   fprintf(stderr,"nx= %d ny= %d nt= %d\n",tshead.nx,tshead.ny,tshead.nt);

   s1 = NULL;
   s1 = read_wccseis(seisfile1,&head1,s1,inbin);

   s2 = NULL;
   s2 = read_wccseis(seisfile2,&head2,s2,inbin);

   s3 = NULL;
   s3 = read_wccseis(seisfile3,&head3,s3,inbin);

   np1_byte = (3*tshead.nx*tshead.ny - 1)*size_float;

   off1 = (ixp + iyp*tshead.nx)*size_float;
   lseek(fdw1,off1,SEEK_CUR);

   off2 = (tshead.nx*tshead.ny + ixp + iyp*tshead.nx)*size_float;
   lseek(fdw2,off2,SEEK_CUR);

   off3 = (2*tshead.nx*tshead.ny + ixp + iyp*tshead.nx)*size_float;
   lseek(fdw3,off3,SEEK_CUR);

   for(it=0;it<tshead.nt;it++)
      {
      /*
      fprintf(stderr,"it= %d, 1 ... ",it);
      fprintf(stderr," 2 ... ");
      fprintf(stderr," 3 ... ");
      fprintf(stderr," DONE.\n");
      */
      rite(fdw1,&s1[it],size_float);
      lseek(fdw1,np1_byte,SEEK_CUR);

      rite(fdw2,&s2[it],size_float);
      lseek(fdw2,np1_byte,SEEK_CUR);

      rite(fdw3,&s3[it],size_float);
      lseek(fdw3,np1_byte,SEEK_CUR);
      }

   close(fdw1);
   close(fdw2);
   close(fdw3);
   }
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
