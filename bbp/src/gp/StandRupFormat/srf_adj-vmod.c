#include "include.h"
#include "structure.h"
#include "function.h"
#include "defs.h"

main(int ac,char **av)
{
FILE *fpr;
char infile[2048], outfile[2048], velfile[2048], sbuf[2048];
int i, it;

struct standrupformat srf1;
struct srf_apointvalues *apval_ptr1;

double tmom, xmom, slip;
float sfac, *stfp;

int inbin = 0;
int outbin = 0;

int print_command = 1;

sprintf(infile,"stdin");
sprintf(outfile,"stdout");

setpar(ac,av);

getpar("infile","s",infile);
getpar("inbin","d",&inbin);
getpar("outfile","s",outfile);
getpar("outbin","d",&outbin);

mstpar("velfile","s",velfile);

getpar("print_command","d",&print_command);

endpar();

read_srf(&srf1,infile,inbin);

if(print_command && atof(srf1.version) >= 2.0)
   load_command_srf(&srf1,ac,av);

apval_ptr1 = srf1.srf_apnts.apntvals;

fpr = fopfile(velfile,"r");

xmom = 0.0;
tmom = 0.0;
for(i=0;i<srf1.srf_apnts.np;i++)
   {
   slip = sqrt(apval_ptr1[i].slip1*apval_ptr1[i].slip1 +
                     apval_ptr1[i].slip2*apval_ptr1[i].slip2 +
                     apval_ptr1[i].slip3*apval_ptr1[i].slip3);

   tmom = tmom + slip*apval_ptr1[i].area*apval_ptr1[i].vs*apval_ptr1[i].vs*apval_ptr1[i].den;

   fgets(sbuf,2048,fpr);
   sscanf(sbuf,"%f %f %f",&(apval_ptr1[i].vp),&(apval_ptr1[i].vs),&(apval_ptr1[i].den));

   xmom = xmom + slip*apval_ptr1[i].area*apval_ptr1[i].vs*apval_ptr1[i].vs*apval_ptr1[i].den;
   }
fclose(fpr);

fprintf(stderr,"original moment= %.5e dyne-cm (Mw= %.2f)\n",tmom,(2.0/3.0)*log10(tmom)-10.73);

sfac = tmom/xmom;
fprintf(stderr,"slip scale factor= %.5e\n",sfac);

for(i=0;i<srf1.srf_apnts.np;i++)
   {
   apval_ptr1[i].slip1 = sfac*apval_ptr1[i].slip1;
   if(apval_ptr1[i].nt1 > 0)
      {
      stfp = apval_ptr1[i].stf1;
      for(it=0;it<apval_ptr1[i].nt1;it++)
         stfp[it] = sfac*stfp[it];
      }

   apval_ptr1[i].slip2 = sfac*apval_ptr1[i].slip2;
   if(apval_ptr1[i].nt2 > 0)
      {
      stfp = apval_ptr1[i].stf2;
      for(it=0;it<apval_ptr1[i].nt2;it++)
         stfp[it] = sfac*stfp[it];
      }

   apval_ptr1[i].slip3 = sfac*apval_ptr1[i].slip3;
   if(apval_ptr1[i].nt3 > 0)
      {
      stfp = apval_ptr1[i].stf3;
      for(it=0;it<apval_ptr1[i].nt3;it++)
         stfp[it] = sfac*stfp[it];
      }
   }

write_srf(&srf1,outfile,outbin);
/*
*/
}
