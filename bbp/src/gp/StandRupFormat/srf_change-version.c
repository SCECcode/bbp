#include "include.h"
#include "structure.h"
#include "function.h"
#include "defs.h"

main(int ac,char **av)
{
char infile[256], outfile[256], new_version[256];

struct standrupformat srf1;
struct srf_prectsegments *prseg_ptr1;
struct srf_apointvalues *apval_ptr1;

int inbin = 0;
int outbin = 0;

int i, ioff;

sprintf(infile,"stdin");
sprintf(outfile,"stdout");
sprintf(new_version,"1.0");

setpar(ac,av);
getpar("infile","s",infile);
getpar("inbin","d",&inbin);
getpar("outfile","s",outfile);
getpar("outbin","d",&outbin);
getpar("new_version","s",new_version);
endpar();

read_srf(&srf1,infile,inbin);

strcpy(srf1.version,new_version);

/*
sprintf(srf1.srf_hcmnt.cbuf,"#\n");
sprintf((srf1.srf_hcmnt.cbuf)+MAXLINE,"# ");
sprintf((srf1.srf_hcmnt.cbuf)+2*MAXLINE,"#\n");

ioff = 2;
for(i=0;i<ac;i++)
   {
   sprintf((srf1.srf_hcmnt.cbuf)+ioff+MAXLINE,"%s ",av[i]);
   while(srf1.srf_hcmnt.cbuf[ioff+MAXLINE] != '\0' && ioff < MAXLINE-2)
      ioff++;
   }
sprintf((srf1.srf_hcmnt.cbuf)+ioff+MAXLINE,"\n");
*/

write_srf(&srf1,outfile,outbin);
}
