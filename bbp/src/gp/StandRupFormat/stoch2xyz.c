#include "include.h"
#include "structure.h"
#include "function.h"
#include "defs.h"

main(int ac,char **av)
{
FILE *fopfile(), *fpr, *fpw;
int nseg;
float xp, yp, dx, dy, val;
int nx, ny, ix, iy, i, nb, ne;
char infile[256], type[64], str[512], outfile[256];

float xoff = 0.0;
float yoff = 0.0;

sprintf(infile,"stdin");
sprintf(outfile,"stdout");
sprintf(type,"slip");
nseg = 0;

setpar(ac,av);
getpar("infile","s",infile);
getpar("outfile","s",outfile);
getpar("type","s",type);
getpar("nseg","d",&nseg);
getpar("xoff","f",&xoff);
getpar("yoff","f",&yoff);
endpar();

if(strcmp(infile,"stdin") == 0)
   fpr = stdin;
else
   fpr = fopfile(infile,"r");

if(strcmp(outfile,"stdout") == 0)
   fpw = stdout;
else
   fpw = fopfile(outfile,"w");

fgets(str,512,fpr);
sscanf(str,"%d",&nseg);

for(i=0;i<nseg;i++)
   {
   fgets(str,512,fpr);
   sscanf(str,"%*f %*f %d %d %f %f",&nx,&ny,&dx,&dy);
   fgets(str,512,fpr);

   if(strcmp(type,"slip") == 0)
      {
      nb = 0;
      ne = 2*ny;
      }
   if(strcmp(type,"tinit") == 0)
      {
      nb = 2*ny;
      ne = 0;
      }

   for(iy=0;iy<nb;iy++)
      fgets(str,512,fpr);

   for(iy=0;iy<ny;iy++)
      {
      yp = (iy + 0.5)*dy;
      for(ix=0;ix<nx;ix++)
         {
         xp = (ix + 0.5)*dx + xoff;
         fscanf(fpr,"%f",&val);

	 fprintf(fpw,"%13.5e %13.5e %13.5e\n",xp,yp,val);
	 }
      }

   xoff = xoff + nx*dx;
   fgets(str,512,fpr); /* get rouge newline character */

   for(iy=0;iy<ne;iy++)
      fgets(str,512,fpr);
   }

fclose(fpr);
fclose(fpw);
}
