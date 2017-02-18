#include "include.h"
#include "structure.h"
#include "function.h"
#include "defs.h"

main(int ac,char **av)
{
int nseg;
char infile[256], type[64], outfile[256];

struct standrupformat srf;

int vmax2slip = 0;
int inbin = 0;
float maxslip = 0.0;
int kp = -1;

sprintf(infile,"stdin");
sprintf(outfile,"stdout");
sprintf(type,"slip");
nseg = 0;

setpar(ac,av);
getpar("infile","s",infile);
getpar("outfile","s",outfile);
getpar("type","s",type);
getpar("nseg","d",&nseg);
getpar("inbin","d",&inbin);
getpar("vmax2slip","d",&vmax2slip);
getpar("maxslip","f",&maxslip);
getpar("kp","d",&kp);
endpar();

read_srf(&srf,infile,inbin);
if(vmax2slip == 1)
   get_vmax2slip(outfile,&srf,type,nseg);
else
   {
   write_maxsvf(outfile,&srf,type,nseg,&maxslip,kp);
   fprintf(stderr,"maxslip= %f\n",maxslip);
   }
}
