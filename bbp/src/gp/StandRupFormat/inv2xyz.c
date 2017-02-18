#include "include.h"
#include "structure.h"
#include "function.h"
#include "defs.h"

main(int ac,char **av)
{
FILE *fpr, *fpw, *fopfile();
int nx, ny, nwin, ix, iy, i;
float flen, fwid, dx, dy, dt;
float xp, yp, *tp, slip, rake;
char infile[256], outfile[256], str[512];

sprintf(infile,"stdin");
sprintf(outfile,"stdout");

setpar(ac,av);
getpar("infile","s",infile);
getpar("outfile","s",outfile);
mstpar("flen","f",&flen);
mstpar("fwid","f",&fwid);
mstpar("dx","f",&dx);
mstpar("dy","f",&dy);
endpar();

nx = (int)(flen/dx + 0.5);
ny = (int)(fwid/dy + 0.5);

tp = (float *)check_malloc(nx*ny*sizeof(float));

if(strcmp(infile,"stdin") == 0)
   fpr = stdin;
else
   fpr = fopfile(infile,"r");

fgets(str,512,fpr);
sscanf(str,"%d %*d %*f %*f %f",&nwin,&dt);

for(ix=0;ix<nx;ix++)
   {
   for(iy=0;iy<ny;iy++)
      {
      tp[ix+iy*nx] = -99;
      fgets(str,512,fpr);
      for(i=0;i<nwin;i++)
         {
	 fscanf(fpr,"%f %f",&slip,&rake);
	 if(slip > 0.01 && tp[ix+iy*nx] < 0.0)
	    tp[ix+iy*nx] = (i-1)*dt;
	 }
      fgets(str,512,fpr);  /* get rouge newline */
      }
   }
fclose(fpr);

if(strcmp(outfile,"stdout") == 0)
   fpw = stdout;
else
   fpw = fopfile(outfile,"w");

for(iy=0;iy<ny;iy++)
   {
   yp = (iy+0.5)*dy;
   for(ix=0;ix<nx;ix++)
      {
      xp = (ix+0.5)*dx;
      fprintf(fpw,"%13.5e %13.5e %13.5e\n",xp,yp,tp[ix+iy*nx]);
      }
   }
fclose(fpw);
}
