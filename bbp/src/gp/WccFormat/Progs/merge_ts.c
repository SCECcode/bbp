#include        "include.h"
#include        "structure.h"
#include        "function.h"

float float_swap(char *);

main(int ac, char **av)
{
FILE *fpr, *fopfile();
struct tsheader tshead;
struct tsheader_proc *tshead_p;
int i, j, ix, iy, it;
int dyts, nproc;
float h, *val;

int swap_bytes = 0;
char cbuf[512];

off_t off, cur_off, head_off, nrite, nxblen;

int fdr, fdw;
char str[512];
char *filebuf;
char **infile;
char filelist[128];
char outfile[128];

setpar(ac, av);
mstpar("filelist","s",filelist);
mstpar("outfile","s",outfile);
mstpar("nproc","d",&nproc);
mstpar("h","f",&h);
endpar();

tshead_p = (struct tsheader_proc *) check_malloc (nproc*sizeof(struct tsheader_proc));
filebuf = (char *) check_malloc (256*nproc*sizeof(char));
infile = (char **) check_malloc (nproc*sizeof(char *));

fpr = fopfile(filelist,"r");

i = 0;
while(fscanf(fpr,"%s",str) != EOF && i < nproc)
   {
   infile[i] = filebuf + i*256;
   strcpy(infile[i],str);

   fprintf(stderr,"%d %s\n",i,infile[i]);

   fdr = opfile_ro(infile[i]);
   reed(fdr,&tshead_p[i],sizeof(struct tsheader_proc));
   close(fdr);

   i++;
   }

if(i != nproc)
   {
   fprintf(stderr,"(%d) entries in filelist != nproc(%d), exiting ...\n",i,nproc);
   exit(-1);
   }

tshead.ix0 = tshead_p[0].ix0;
tshead.iy0 = tshead_p[0].iy0;
tshead.iz0 = tshead_p[0].iz0;
tshead.it0 = tshead_p[0].it0;
tshead.nx = tshead_p[0].nx;
tshead.ny = tshead_p[0].ny;
tshead.nz = tshead_p[0].nz;
tshead.nt = tshead_p[0].nt;
tshead.dx = tshead_p[0].dx;
tshead.dy = tshead_p[0].dy;
tshead.dz = tshead_p[0].dz;
tshead.dt = tshead_p[0].dt;
tshead.modelrot = tshead_p[0].modelrot;
tshead.modellat = tshead_p[0].modellat;
tshead.modellon = tshead_p[0].modellon;

val = (float *) check_malloc (3*tshead.nx*tshead.ny*sizeof(float));
for(iy=0;iy<3*tshead.nx*tshead.ny;iy++)
   val[iy] = 0;

fdw = croptrfile(outfile);

rite(fdw,&tshead,sizeof(struct tsheader));;

for(it=0;it<tshead.nt;it++)
   rite(fdw,val,3*tshead.nx*tshead.ny*sizeof(float));

lseek(fdw,sizeof(struct tsheader),SEEK_SET);
head_off = sizeof(struct tsheader);
cur_off = sizeof(struct tsheader);

dyts = (int)(tshead.dy/h + 0.5);

fprintf(stderr,"nx= %d\n",tshead.nx);
fprintf(stderr,"ny= %d\n",tshead.ny);
fprintf(stderr,"nt= %d\n",tshead.nt);
fprintf(stderr,"dy= %lg\n",tshead.dy);
fprintf(stderr,"h= %lg\n",h);
fprintf(stderr,"dyts= %d\n",dyts);

nxblen = (off_t)(tshead.nx)*sizeof(float);

for(i=0;i<nproc;i++)
   {
   fprintf(stderr,"%d of %d\n",i+1,nproc);

   iy = tshead_p[i].iyleft;
   while(iy%dyts != 0)
      iy++;

   iy = iy/dyts;

   /*
   fprintf(stderr,"%5d: %5d %5d %5d, %5d\n",i,tshead_p[i].iyleft,tshead_p[i].iyright,tshead_p[i].localny,iy);
   */

   if(tshead_p[i].localny > 0)
      {
      fdr = opfile_ro(infile[i]);
      lseek(fdr,sizeof(struct tsheader_proc),SEEK_SET);

      for(it=0;it<tshead.nt;it++)
         {
         for(j=0;j<3;j++)
            {
	    reed(fdr,val,tshead.nx*tshead_p[i].localny*sizeof(float));

	    off = head_off + (off_t)((3*it+j)*tshead.ny + iy)*nxblen - cur_off;

	    lseek(fdw,off,SEEK_CUR);
            nrite = rite(fdw,val,tshead.nx*tshead_p[i].localny*sizeof(float));
	    cur_off = nrite + off + cur_off;
            }
         }

      close(fdr);
      }
   }
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
