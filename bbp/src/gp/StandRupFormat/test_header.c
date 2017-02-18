#include "include.h"
#include "structure.h"
#include "function.h"
#include "defs.h"

#define         MAXLINE            10000

main(int ac,char **av)
{
struct standrupformat srf;
FILE *fpr, *fopfile();
char str[MAXLINE], infile[256], pword[32];
int iseg, nseg, np_seg, np_tot;

sprintf(infile,"stdin");

setpar(ac,av);
getpar("infile","s",&infile);
endpar();

if(strcmp(infile,"stdin") == 0)
   fpr = stdin;
else
   fpr = fopfile(infile,"r");

fgets(str,MAXLINE,fpr);
sscanf(str,"%s",srf.version);

fgets(str,MAXLINE,fpr);
while(str[0] == '#')
   {
   fprintf(stderr,"COMMENT LINE: %s",str);
   if(fgets(str,MAXLINE,fpr) == NULL)
      break;
   }

sscanf(str,"%s",pword);

if(strncmp(pword,"PLANE",5) == 0)
   {
   fprintf(stderr,"pword= %s, I'm going to get the PLANE information...\n",pword);
   }

if(strncmp(pword,"POINTS",6) != 0)
   {
   fgets(str,MAXLINE,fpr);
   while(str[0] == '#')
      {
      fprintf(stderr,"COMMENT LINE: %s",str);
      if(fgets(str,MAXLINE,fpr) == NULL)
         break;
      }

   sscanf(str,"%s",pword);
   }

if(strcmp(pword,"POINTS") == 0)
   {
   fprintf(stderr,"pword= %s, I'm going to get the POINT SOURCE information...\n",pword);
   }
else if(strcmp(pword,"POINTS_SEG") == 0)
   {
   nseg = 0;
   np_tot = 0;
   while(strcmp(pword,"POINTS_SEG") == 0)
      {
      nseg ++;
      sscanf(str,"%*s %d %d",&np_seg,&iseg);
      np_tot = np_tot + np_seg;
      fprintf(stderr,"SEGMENT %d: np_seg= %d np_tot= %d\n",nseg,np_seg,np_tot);
      /*  
      {
          }
	  */

      if(fgets(str,MAXLINE,fpr) == NULL)
         break;
      else
         sscanf(str,"%s",pword);
      }
   }
else
   {
   fprintf(stderr,"*** No valid data found.  Last input line read is:\n");
   fprintf(stderr,"%s",str);
   }
}
