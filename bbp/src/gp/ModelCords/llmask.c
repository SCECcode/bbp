#include "include.h"

main(ac,av)
int ac;
char **av;
{
FILE *fpw, *fpr, *fopfile();
float lonp, latp;
int nx, ny, npts, i, j, ix, iy;
int jp, ip, nedge;
char string[512], infile[512], outfile[512];
char edgefile[512];
float elat[5000], elon[5000];

sprintf(infile,"stdin");
sprintf(outfile,"stdout");

setpar(ac,av);
getpar("infile","s",infile);
getpar("outfile","s",outfile);
mstpar("edgefile","s",edgefile);
endpar();

fpr = fopfile(edgefile,"r");

i = 0;
while(fscanf(fpr,"%f %f",&elon[i],&elat[i]) != EOF)
   i++;

fclose(fpr);
nedge = i;

if(strcmp(infile,"stdin") == 0)
   fpr = stdin;
else
   fpr = fopfile(infile,"r");

if(strcmp(outfile,"stdout") == 0)
   fpw = stdout;
else
   fpw = fopfile(outfile,"w");

while(fgets(string,512,fpr) != NULL)
   {
   sscanf(string,"%f %f",&lonp,&latp);

   if(inside(&lonp,&latp,elon,elat,nedge) == 1)
      fprintf(fpw,"%s",string);
   }
fclose(fpr);
fclose(fpw);
}

inside(xp,yp,xv,yv,nv)
float *xp, *yp, *xv, *yv;
int nv;
{
float x1, y1, x2, y2, xx;
int i, flag;
int nleft = 0;

x2 = xv[nv-1];
y2 = yv[nv-1];
for(i=0;i<nv;i++)
   {
   x1 = x2;
   y1 = y2;
   x2 = xv[i];
   y2 = yv[i];

   if(y1 == (*yp) && y2 == (*yp))
      {
      if((x1 < (*xp) && x2 >= (*xp)) || (x2 < (*xp) && x1 >= (*yp)))
         return(1);
      }
   else if((y1 < (*yp) && y2 >= (*yp)) || (y2 < (*yp) && y1 >= (*yp)))
      {
      xx = x1 + ((*yp)-y1)*(x2-x1)/(y2-y1);
      if(xx == (*xp))
         return(1);
      else if(xx > (*xp))
         nleft++;
      }
   }

if(nleft == 0)
   flag = 0;
else
   flag = nleft%2;

return(flag);
}
 
FILE *fopfile(name,mode)
char *name, *mode;
{
FILE *fp, *fopen();

if((fp = fopen(name,mode)) == NULL)
   {
   fprintf(stderr,"CAN'T FOPEN FILE = %s, MODE = %s\n", name, mode);
   exit(-1);
   }
return(fp);
}
